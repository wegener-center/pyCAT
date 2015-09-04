# -*- coding: utf-8 -*-

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
#
"""
A class that holds information about a climate dataset
"""

import iris
from iris.experimental.equalise_cubes import equalise_attributes
from cartopy.crs import Geodetic

import os
import glob
import numpy as np

class Dataset(object):
    """
    A dataset is a holding meta data to an iris cube
    """
    def __init__(self, directory, filename, constraints=None,
                 callback=None, tmp_directory=None):
        self.directory = directory
        self.filename = filename
        self.tmp_directory = tmp_directory or os.path.join(directory, 'tmp')
        # try to create the temporary directory
        # if there are other problems (e.g. permissions) except that it already
        # exists an exception will be raised
        try:
            os.makedirs(self.tmp_directory)
        except OSError:
            if not os.path.isdir(self.tmp_directory):
                raise

        self.cube_list = iris.load(os.path.join(directory, filename),
                                   constraints=constraints, callback=callback)
        self._orig_standard_name = self.cube_list[0].standard_name
        self._orig_var_name = self.cube_list[0].var_name
        self._orig_long_name = self.cube_list[0].long_name
        self.calendar = self.cube_list[0].coord('time').units.calendar
        self._orig_units = self.cube_list[0].units

        for ndim, coord in enumerate(self.cube_list[0].dim_coords):
            if coord.units.is_time_reference():
                self._timeaxis = ndim
                break

        # save the period of the entire dataset
        self._orig_period = (
            self.cube_list[0].dim_coords[self._timeaxis].units.num2date(
                self.cube_list[0].dim_coords[self._timeaxis].points[0]),
            self.cube_list[-1].dim_coords[self._timeaxis].units.num2date(
                self.cube_list[-1].dim_coords[self._timeaxis].points[-1])
        )

        # save the geographical extent
        x = self.cube_list[0].coord(axis="X", dim_coords=True)
        y = self.cube_list[0].coord(axis="Y", dim_coords=True)

        # save the coordinate system
        self._coord_system = x.coord_system.as_cartopy_crs()
        
        if not x.has_bounds():
            x.guess_bounds()
            y.guess_bounds()
            remove_bounds = True

        x_edges = np.concatenate((x.bounds[:,0], np.array([x.bounds[-1,-1]]*y.shape[0]),
                                  x.bounds[:,-1], np.array([x.bounds[0,0]]*y.shape[0])))
        y_edges = np.concatenate((np.array([y.bounds[0,0]]*x.shape[0]), y.bounds[:,0],
                                  np.array([y.bounds[-1,-1]]*x.shape[0]), y.bounds[:,-1]))

        ll_crs = Geodetic()
        lon_lat_edges = ll_crs.transform_points(
            x.coord_system.as_cartopy_crs(), x_edges, y_edges)[:,:-1]
        east, west = max(lon_lat_edges[:,0]), min(lon_lat_edges[:,0])
        north, south = max(lon_lat_edges[:,1]), min(lon_lat_edges[:,1])
            
        self._orig_extent = (north, east, south, west)
        if remove_bounds:
            x.bounds = None
            y.bounds = None


    def get_cube(self, extra_constraints=None):
        """
        return the cube of the dataset constrainted by period and extent
        (if they have been set) and extra constraints given in the call

        also all adjustments, i.e. units, standard_name, ... are applied

        Kwargs:
        
        * extra_constraints (iris.Constraint):
            will be applied to the cube

        Returns:

            the concatenated constrained cube of the dataset
        """
        constraints = extra_constraints
        try:
            start, end = self.period
            import datetime
            constraints &= iris.Constraint(
                time=lambda cell: start <= cell.point < end+datetime.timedelta(1)
            )
        except AttributeError:
            pass

        try:
            constraints &= self._extent_constraint()
        except AttributeError:
            pass
            
        with iris.FUTURE.context(cell_datetime_objects=True):
            cl = self.cube_list.extract(constraints)
        equalise_attributes(cl)
        
        merged_cube = self._merge_by_time(cl)
        try:
            for k, v in self.adjustments.iteritems():
                setattr(merged_cube, k, v)
        except AttributeError:
            pass
        
        return merged_cube


    def _extent_constraint(self):
        """
        if Dataset has an extent set return the geographical constraint

        Returns:
        
            iris.Constraint over the x and y axis
        """
        north, east, south, west = self.extent
        from matplotlib.path import Path
        poly = Path([[west, south], [east, south],
                     [east, north], [west, north]], closed=True)
        ll_crs = Geodetic()
        
        x = self.cube_list[0].coord(axis='X', dim_coords=True)
        y = self.cube_list[0].coord(axis='Y', dim_coords=True)

        dx = x.points[1] - x.points[0]
        dy = y.points[1] - y.points[0]

        # get the order of the spatial dimensions
        xdim = 1
        ydim = 0
        if self.cube_list[0].coord_dims(x) < self.cube_list[0].coord_dims(y):
            xdim = 0
            ydim = 1
            
        xgrid, ygrid = np.meshgrid(x.points, y.points)
        ll_field = ll_crs.transform_points(self._coord_system, xgrid, ygrid)[:,:,:2]
        ll_flat = ll_field.reshape((-1,2))
        mask = poly.contains_points(ll_flat).reshape(xgrid.shape)

        inside_indices = np.where(mask)
        minx, maxx = inside_indices[xdim].min(), inside_indices[xdim].max()
        miny, maxy = inside_indices[ydim].min(), inside_indices[ydim].max()

        west_bound, east_bound = x.points[minx]-dx, x.points[maxx]+dx
        south_bound, north_bound = y.points[miny]-dy, y.points[maxy]+dy

        return constraints & iris.Constraint(coord_values={
            x.standard_name: lambda cell: west_bound <= cell.point <= east_bound,
            y.standard_name: lambda cell: south_bound <= cell.point <= north_bound
        })

    def _merge_by_time(self, cl):
        """
        try to merge a cube_list by the time-axis and adjust units

        Args:

        * cl (iris.cube.CubeList):
            the already equalised cubelist

        Returns:

            a single concatenated cube with the proper units
        """
        try:
            units = self.units
        except AttributeError:
            units = cl[0].units
            
        min_dim = cl[0].coord(axis='T').shape[0]
        for cube in cl[1:]:
            if min_dim == 1:
                break
            min_dim = min(min_dim, cube.coord(axis='T').shape[0])

        if min_dim > 1:
            for cube in cl:
                cube.convert_units(units)
            ret = cl.concatenate_cube()
        else:
            cl_new = iris.cube.CubeList()
            for cube in cl:
                cube.convert_units(units)
                if cube.coord(axis='T').shape[0] == 1:
                    cl_new.append(cube)
                else:
                    for time in xrange(cube.coord(axis='T').shape[0]):
                        tmp = cube[time]
                        cl_new.append(tmp)

            ret = cl_new.merge_cube()
        return ret
            
    # define the time period
    @property
    def period(self):
        return self._period

    @period.setter
    def period(self, value):
        self._period = value

    @period.deleter
    def period(self):
        del self._period

    #define extent
    @property
    def extent(self):
        return self._extent

    @extent.setter
    def extent(self, value):
        self._extent = value

    @extent.deleter
    def extent(self):
        del self._extent

    # define units
    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, value):
        self._units = value

    @units.deleter
    def units(self):
        del self._units

    # dictionary for some adjustments of cube attributes
    @property
    def adjustments(self):
        return self._adjustments

    @adjustments.setter
    def adjustments(self, dictionary):
        self._adjustments = dictionary

    @adjustments.deleter
    def adjustments(self):
        del self._adjustments
        
