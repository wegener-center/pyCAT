# (C) Wegener Center for Climate and Global Change, University of Graz, 2015
#
# This file is part of pyCAT.
#
# pyCAT is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License version 3 as published by the
# Free Software Foundation.
#
# pyCAT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with pyCAT. If not, see <http://www.gnu.org/licenses/>.

import abc

class BiasCorrection(object):
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, observation, model, scenarios, reference_period,
                 correction_period, window):
        """
        Args:

        * observation (pycat.io.dataset.Dataset):
            the observation dataset

        * model (pycat.io.dataset.Dataset):
            the model dataset

        * scenarios (pycat.io.dataset.Dataset or list of those)
            the scenarios that shall be bias corrected

        * reference_period (tuple of datetime.datetime):
            the reference period that observations and model share

        * correction_period (tuple of datetime.datetime):
            the period for which the correction shall be done

        * window (int):
            the window size in days that is considered for generating
            large enough samples of the observation and model distribution
        """
        self.obs = observation
        obs_phenomenon = {
            'units': self.obs._orig_units,
            'standard_name': self.obs._orig_standard_name,
            'long_name': self.obs._orig_long_name,
            'var_name': self.obs._orig_var_name
        }
        self.mod = model
        self.sce = scenarios
        self.obs.period = reference_period
        self.mod.period = reference_period

        try:
            obs_extent = self.obs.extent
        except AttributeError:
            obs_extent = self.obs._orig_extent

        self.mod.extent = obs_extent
        self.mod.adjustments = obs_phenomenon

        # make the scenarios list if they are not already
        if hasattr(scenarios, '__iter__'):
            self.sce = scenarios
        else:
            self.sce = [ scenarios ]

        for sce in self.sce:
            sce.extent = obs_extent
            sce.adjustments = obs_phenomenon
            if correction_period:
                sce.period = correction_period
                
        self.window = window

        # initialize the regridder
        # obs_cube = self.obs.get_cube()

    @abc.abstractmethod
    def correct(self):
        """
        function that implements the correction method
        """
        return

class QuantileMapping(BiasCorrection):
    def __init__(self, observation, model, scenarios, reference_period,
                 correction_period=None, window=15):
        super(self.__class__, self).__init__(
            observation, model, scenarios, reference_period, correction_period, window)
        

    def correct(self):
        """
        implementation of the correction method

        at the moment, the quantile mapped data is written to the
        temporary directory of the scenario dataset with a simple
        filename
        this will be changed in the future
        """
        import os
        from utils import generate_day_constraint_with_window
        import iris
        from iris.analysis import Linear
        import numpy as np
        from numpy import ma
        from statsmodels.tools.tools import ECDF

        regridder = None
        for day in xrange(self.mod.days_in_year[self.mod.calendar]):
            # check if obs and mod have calendars with same number of days
            obs_cube = None
            if self.obs.days_in_year[self.obs.calendar] != \
               self.mod.days_in_year[self.mod.calendar]:
                obs_day = self.obs.days_in_year[self.obs.calendar] * day /\
                          self.mod.days_in_year[self.mod.calendar]
                day_constraint, window_constraint = \
                    generate_day_constraint_with_window(obs_day, self.window, self.obs.calendar)
                obs_cube = self.obs.get_cube(window_constraint)
                
            day_constraint, window_constraint = \
                generate_day_constraint_with_window(day, self.window, self.mod.calendar)
            mod_cube = self.mod.get_cube(window_constraint)
            obs_cube = obs_cube or self.obs.get_cube(window_constraint)

            regridder = regridder or Linear().regridder(mod_cube, obs_cube)
            mod_cube = regridder(mod_cube)

            obs_first_time_slice = obs_cube.slices_over(obs_cube.coords(axis='T', dim_coords=True)).next()
            try:
                cells_with_value = np.where(~obs_first_time_slice.data.mask)
            except AttributeError:
                cells_with_value = np.where(np.isfinite(obs_first_time_slice.data))

            sce_cubes = []
            for sce in self.sce:
                sce_cube = regridder(sce.get_cube(day_constraint))
                try:
                    # mask the scenario data if the observation has a mask
                    sce_cube.data = ma.masked_array(
                        sce_cube.data, mask=[ obs_first_time_slice.data.mask ]*sce_cube.coord('time').shape[0])
                except AttributeError:
                    pass
                sce_cubes.append(sce_cube)
                
            # loop over the cells
            for cell in xrange(cells_with_value[0].shape[0]):
                i, j = cells_with_value[0][cell], cells_with_value[1][cell]
                obs_data = obs_cube.data[:,i,j]
                mod_data = mod_cube.data[:,i,j]
                mod_ecdf = ECDF(mod_data)

                for sce_cube in sce_cubes:
                    sce_data = sce_cube[:,i,j].data
                    p = mod_ecdf(sce_data)*100
                    corr = np.percentile(obs_data, p) - \
                           np.percentile(mod_data, p)
                    sce_cube.data[:,i,j] += corr
                

            # save the cubes temporarily
            for i, sce_cube in enumerate(sce_cubes):
                tmp_dir = self.sce[i].tmp_directory
                filename = 'qm_{}_scenario-{:d}_day-{:03d}.nc'.format(sce_cube.var_name, i, day)
                iris.save(sce_cube, os.path.join(tmp_dir, filename))
