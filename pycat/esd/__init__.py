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
# You should have received a copy of the GNU General Public License
# along with pyCAT. If not, see <http://www.gnu.org/licenses/>.

from tempfile import gettempdir

import numpy as np
from iris.analysis import Linear

from .methods import quantile_mapping, scaled_distribution_mapping


class BiasCorrector(object):

    """
    Base class for all bias correction classes
    """

    def __init__(self, call_func, observation, model, scenarios,
                 reference_period=None, correction_period=None,
                 time_unit='day', work_dir=gettempdir(),
                 interpolator=Linear(), save_regridded=False):
        """
        Args:

        * call_func (callable):
            | *call signature*: (obs_cube, ref_cube, sce_cubes,
                                 \*args, \**kwargs)

        * observation (:class:`.io.Dataset`):
            the observation dataset

        * model (:class:`.io.Dataset`):
            the model dataset

        * scenarios (:class:`.io.Dataset` or list of those)
            the scenarios that shall be bias corrected

        Kwargs:

        * reference_period (tuple of :class:`datetime.datetime`):
            the reference period that observations and model share;
            if not given, take it from the observation and model
            reference dataset, respectively

        * correction_period (tuple of :class:`datetime.datetime`):
            the period for which the correction shall be done;
            if not given, take it from the scenario dataset

        * time_unit (str):
            correction will be performed on daily (day) or
            monthly (month) basis

        * interpolator:
            an available interpolation scheme for regridding to the
            observational dataset. Currently available schemes are:

            * :class:`iris.analysis.Linear` (default)

            * :class:`iris.analysis.Nearest`

            * :class:`iris.analysis.AreaWeighted`

        * work_dir (path):
            directory where intermediate files will be written

        * save_regridded (boolean):
            wheter regridded data shall be stored to disk (default: False)
        """
        self.call_func = call_func
        self.obs = observation
        obs_phenomenon = {
            'units': self.obs._orig_units,
            'standard_name': self.obs._orig_standard_name,
            'long_name': self.obs._orig_long_name,
            'var_name': self.obs._orig_var_name
        }
        self.mod = model
        self.sce = scenarios

        # set the reference period
        if reference_period:
            self.obs.period = reference_period
            self.mod.period = reference_period

        # set the spatial extent of the observation to the model
        self.mod.extent = self.obs.extent
        self.mod.adjustments = obs_phenomenon

        # make the scenarios list if they are not already
        if hasattr(scenarios, '__iter__'):
            self.sce = scenarios
        else:
            self.sce = [scenarios]

        for sce in self.sce:
            sce.extent = self.obs.extent
            sce.adjustments = obs_phenomenon
            if correction_period:
                sce.period = correction_period

        self.time_unit = time_unit
        self.interpolator = interpolator
        self.work_dir = work_dir
        self.save_regridded = save_regridded

    def correct(self, unit_list=None, *args, **kwargs):
        """
        assemble data that is given to the actual correction function

        kwargs are passed to :meth:`call_func`

        Args:

        * unit_list (None, int, iterable):
            depending on self.time_unit this is interpreted as
            all days/months of year (None), single day/month (int) or
            list of days/months (iterable)
        """
        import os
        from .utils import (
            generate_day_constraint_with_window,
            generate_month_constraint)
        import iris
        from numpy import ma
        import collections

        regridder = None
        if isinstance(unit_list, int):
            unit_list = [unit_list]
        elif isinstance(unit_list, collections.Iterable):
            pass
        else:
            unit_list = self.time_unit == 'day' and \
                range(self.mod.days_in_year[self.mod.calendar]) or \
                range(1, 13)

        # padding for the filename day/month number
        padding = self.time_unit == 'day' and 3 or 2

        for unit in unit_list:
            obs_cube = None
            if self.time_unit == 'day':
                # check if obs and mod have calendars with same number of days
                if self.obs.days_in_year[self.obs.calendar] != \
                   self.mod.days_in_year[self.mod.calendar]:
                    obs_day = self.obs.days_in_year[self.obs.calendar] \
                        * unit / self.mod.days_in_year[self.mod.calendar]
                    single_constraint, window_constraint = \
                        generate_day_constraint_with_window(
                            obs_day, self.window, self.obs.calendar)
                    obs_cube = self.obs.get_cube(window_constraint)

                single_constraint, window_constraint = \
                    generate_day_constraint_with_window(
                        unit, self.window, self.mod.calendar)

            else:
                # a monthly method: single- and window-constraint are the same
                single_constraint, window_constraint = \
                    [generate_month_constraint(unit)] * 2

            # ok, got all constraints now it's the same for daily and monthly
            # extract data from obs, mod and sce
            obs_cube = obs_cube or self.obs.get_cube(window_constraint)
            mod_cube = self.mod.get_cube(window_constraint)

            regridder = regridder or self.interpolator.regridder(
                mod_cube, obs_cube)
            mod_cube = regridder(mod_cube)

            obs_first_time_slice = obs_cube.slices_over(
                obs_cube.coords(axis='T', dim_coords=True)).next()

            sce_cubes = iris.cube.CubeList()
            for sce in self.sce:
                sce_cube = regridder(sce.get_cube(single_constraint))
                try:
                    # mask the scenario data if the observation has a mask
                    sce_cube.data = ma.masked_array(
                        sce_cube.data,
                        mask=[obs_first_time_slice.data.mask]
                        * sce_cube.coord('time').shape[0])
                except AttributeError:
                    pass
                sce_cubes.append(sce_cube)

            fn_template = '{method}_{variable}_scenario-{scenario:d}' \
                          '_{startyear}-{endyear}_' \
                          '{time_unit}-{unit:0{padding}d}.nc'
            if self.save_regridded:
                for sce_number, sce_cube in enumerate(sce_cubes):
                    time_coord = sce_cube.coord('time')
                    startyear, endyear = [
                        date.year for date in
                        time_coord.units.num2date(
                            time_coord.points[np.array([0, -1])])]
                    filename = fn_template.format(
                        method='regridded', variable=sce_cube.var_name,
                        startyear=startyear, endyear=endyear,
                        scenario=sce_number, time_unit=self.time_unit,
                        unit=unit, padding=padding)
                    iris.save(sce_cube, os.path.join(self.work_dir, filename))

            # call the correction function
            self.call_func(obs_cube, mod_cube, sce_cubes, *args, **kwargs)

            # save the cubes into the temporary directory with a simple
            # filename
            for sce_number, sce_cube in enumerate(sce_cubes):
                time_coord = sce_cube.coord('time')
                startyear, endyear = [
                    date.year for date in
                    time_coord.units.num2date(
                        time_coord.points[np.array([0, -1])])]
                filename = fn_template.format(
                    method=self.call_func.__name__.strip('_'),
                    variable=sce_cube.var_name, scenario=sce_number,
                    startyear=startyear, endyear=endyear,
                    time_unit=self.time_unit, unit=unit, padding=padding)
                iris.save(sce_cube, os.path.join(self.work_dir, filename))


class QuantileMapping(BiasCorrector):

    """
    convenience class for quantile mapping
    """

    def __init__(self, observation, model, scenarios, window=15,
                 *args, **kwargs):
        super(QuantileMapping, self).__init__(
            quantile_mapping, observation, model, scenarios,
            time_unit='day', *args, **kwargs)
        self.window = window


class ScaledDistributionMapping(BiasCorrector):

    """
    convenience class for scaled distribution mapping
    """

    def __init__(self, observation, model, scenarios, *args, **kwargs):
        super(ScaledDistributionMapping, self).__init__(
            scaled_distribution_mapping, observation, model, scenarios,
            time_unit='month', *args, **kwargs)
