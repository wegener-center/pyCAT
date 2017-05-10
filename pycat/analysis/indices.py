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
"""
Calculation of climate indices
"""
import iris
import iris.coord_categorisation as ccat
import numpy as np
from iris.analysis import Aggregator
from iris.exceptions import CoordinateNotFoundError
from pycat.analysis.utils import _create_cube, _make_time_dimension


def consecutive_dry_days(cube, period='year', length=6, threshold=1.):
    """
    calculate consecutive dry days within an iris.cube.Cube

    Args:

    * cube (iris.cube.Cube):
        An iris.cube.Cube holding precipiation amount in mm/day
    * period (string):
        Period over that the CDD will be calculated. Can be 'year', 'season'
        or 'month'. If period is 'season' or 'month' the CDD will be averaged
        over the years

    Kwargs:

    * length (int):
        The number of days without rainfall that define a dry period

    * threshold (float):
        The upper limit of daily rainfall in mm that indicates
        'no precipitation'

    Returns:

        An iris.cube.CubeList that holds two iris.cube.Cubes with the longest
        period of dry days in the given period and the mean of the number of
        dry periods with respect to the given length
    """
    def _cdd_index(array, axis, threshold):
        """
        Calculate the consecutive dry days index.

        This function is used as an iris.analysis.Aggregator

        Args:

        * array (numpy.array or numpy.ma.array):
            array that holds the precipitation data

        * axis (int):
            the number of the time-axis

        * threshold (float):
            the threshold that indicates a precipiation-less day

        Returns:
            the aggregation result, collapsing the 'axis' dimension of
            the 'data' argument
        """
        from pycat.analysis.utils import (
            _get_max_true_block_length, _get_true_block_lengths)

        up_down = _get_true_block_lengths(array < threshold, axis)
        return _get_max_true_block_length(up_down)

    def _cdd_periods(array, axis, threshold, length):
        """
        Calculate the number of consecutive dry days periods.

        This function is used as an iris.analysis.Aggregator

        Args:

        * array (numpy.array or numpy.ma.array):
            array that holds the precipitation data

        * axis (int):
            the number of the time-axis

        * threshold (float):
            the threshold that indicates a precipiation-less day

        * length (int):
            number of days that a dry period must last

        Returns:
            the aggregation result, collapsing the 'axis' dimension
            of the 'data' argument
        """
        from pycat.analysis.utils import (
            _get_len_true_block_length, _get_true_block_lengths)

        up_down = _get_true_block_lengths(array < threshold, axis)
        return _get_len_true_block_length(up_down, length)

    # build the iris.analysis.Aggregators
    cdd_index = Aggregator('cdd_index', _cdd_index)
    cdd_periods = Aggregator('cdd_periods', _cdd_periods)

    # check if the cube already has the needed auxiliary coordinates
    if period == 'season':
        # add the season_year auxiliary coordinate
        try:
            years = np.unique(cube.coord('season_year').points)
        except CoordinateNotFoundError:
            ccat.add_season_year(cube, 'time')
            years = np.unique(cube.coord('season_year').points)
        constraint_year_key = 'season_year'
    else:
        # add calendar years
        try:
            years = np.unique(cube.coord('year').points)
        except CoordinateNotFoundError:
            ccat.add_year(cube, 'time')
            years = np.unique(cube.coord('year').points)
        constraint_year_key = 'year'

    if period in ['season', 'month']:
        try:
            index_period = np.unique(cube.coord('%s_number' % period).points)
        except CoordinateNotFoundError:
            cat = getattr(ccat, 'add_%s_number' % period)
            cat(cube, 'time')
            index_period = np.unique(cube.coord('%s_number' % period).points)

    # create time-axis of resulting cubes
    time_dimension = _make_time_dimension(
        cube.coord('time').units.num2date(cube.coord('time').points[0]),
        cube.coord('time').units.num2date(cube.coord('time').points[-1]),
        period=period)
    # create the empty resulting cubes
    dim_coords_and_dims = []
    slices = []
    for coord in cube.dim_coords:
        if coord.units.is_time_reference():
            dim_coords_and_dims.append(
                (time_dimension, cube.coord_dims(coord)))
            slices.append(0)
            time_axis = cube.coord_dims(coord)[0]
        else:
            dim_coords_and_dims.append((coord, cube.coord_dims(coord)))
            slices.append(slice(None, None, None))

    cdd_index_cube = _create_cube(
        long_name='Consecutive dry days is the greatest number of '
                  'consecutive days per time period with daily '
                  'precipitation amount below %s mm.' % threshold,
        var_name='consecutive_dry_days_index_per_time_period',
        units=iris.unit.Unit('1'),
        dim_coords_and_dims=dim_coords_and_dims)

    cdd_periods_cube = _create_cube(
        long_name='Number of cdd periods in given time period '
        'with more than %d days.' % length,
        var_name='number_of_cdd_periods_with_more_than_'
                 '%ddays_per_time_period' % length,
        units=iris.unit.Unit('1'),
        dim_coords_and_dims=dim_coords_and_dims)

    # differentiate between the considered period
    if period == 'year':
        # just run the aggregation over all given years resulting in
        # the maximum cdd length and the number of cdd periods for each year
        for year in years:
            tmp_cube = cube.extract(iris.Constraint(year=year))
            slices[time_axis] = year - years[0]
            cdd_index_data = tmp_cube.collapsed(
                'time', cdd_index, threshold=threshold).data
            cdd_periods_data = tmp_cube.collapsed(
                'time', cdd_periods, threshold=threshold, length=length).data

            cdd_index_cube.data[slices] = cdd_index_data
            cdd_periods_cube.data[slices] = cdd_periods_data

        return iris.cube.CubeList(
            (cdd_index_cube, cdd_periods_cube)
        )

    else:
        # run the aggregation over all seasons/months of all years
        # afterwards aggregate the seasons/month by MAX Aggregator
        # for the cdd_index and the MEAN Aggregator for cdd_periods
        for year in years:
            for p in index_period:
                constraint_dict = {'%s_number' % period: p,
                                   constraint_year_key: year}
                tmp_cube = cube.extract(iris.Constraint(**constraint_dict))
                if tmp_cube:
                    # the extraction can lead to empty cubes for seasons
                    # in the last year
                    time_index = (year - years[0]) * len(index_period) + p
                    # months numbers start at 1
                    if period == 'month':
                        time_index -= 1
                    slices[time_axis] = time_index
                    cdd_index_data = tmp_cube.collapsed(
                        'time', cdd_index, threshold=threshold).data
                    cdd_periods_data = tmp_cube.collapsed(
                        'time', cdd_periods, threshold=threshold,
                        length=length).data

                    cdd_index_cube.data[slices] = cdd_index_data
                    cdd_periods_cube.data[slices] = cdd_periods_data

        # aggregate over seasons/months
        cat = getattr(ccat, 'add_%s' % period)
        cat(cdd_index_cube, 'time')
        cat(cdd_periods_cube, 'time')

        cdd_index_mean = cdd_index_cube.aggregated_by(
            period, iris.analysis.MEAN)
        cdd_periods_mean = cdd_periods_cube.aggregated_by(
            period, iris.analysis.MEAN)

        cdd_index_mean.remove_coord('time')
        cdd_periods_mean.remove_coord('time')
        return iris.cube.CubeList(
            (cdd_index_mean, cdd_periods_mean)
        )
