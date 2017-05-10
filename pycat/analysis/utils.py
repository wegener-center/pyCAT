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

import datetime

import iris
import numpy as np
import numpy.ma as ma
from dateutil import parser
from dateutil.relativedelta import relativedelta
from iris.coords import DimCoord

try:
    from cf_units import Unit
except:
    # for iris<=1.9
    from iris.unit import Unit


def _get_max_true_block_length(up_down):
    """
    Calculate the maximum length of True blocks in an array.

    The True blocks in the array are defined by up/down changes from False
    to True and vice-versa.

    Args:

    * up_down (tuple of numpy.arrays):
        information about the up/down changes

    Returns:

        numpy.array of same dimension of up_down[0] holding the number
        of maximum of the True block lengths in the array
    """
    out_shape = up_down[0].shape[:-1]
    ret = ma.zeros(out_shape) - 1
    ret.fill_value = -1
    ret.mask = True
    for index in np.ndindex(out_shape):
        if not up_down[0][index].mask.any():
            start_idxs = ma.where(up_down[0][index])
            stop_idxs = ma.where(up_down[1][index])
            try:
                ret[index] = np.max(stop_idxs[-1] - start_idxs[-1] + 1)
            except ValueError:
                ret[index] = 0
    return ret


def _get_len_true_block_length(up_down, length):
    """
    Calculate the len of True blocks in an array that succeed the given length.

    The True blocks in the array are defined by up/down changes from False
    to True and vice-versa.

    Args:

    * up_down (tuple of numpy.arrays):
        information about the up/down changes

    * length (int or float):
        threshold for the length of blocks to be accounted for

    Returns:

        numpy.array of same dimension of up_down[0] holding the number
        of block lengths succeeding the given length
    """
    out_shape = up_down[0].shape[:-1]
    ret = ma.zeros(out_shape) - 1
    ret.fill_value = -1
    ret.mask = True
    for index in np.ndindex(out_shape):
        if not up_down[0][index].mask.any():
            start_idxs = ma.where(up_down[0][index])
            stop_idxs = ma.where(up_down[1][index])
            try:
                dry_blocks = stop_idxs[-1] - start_idxs[-1] + 1
                ret[index] = np.where(dry_blocks > length)[0].shape[0]
            except ValueError:
                ret[index] = 0
    return ret


def _get_true_block_lengths(array, axis=-1):
    """
    calculate the lengths of True blocks in an array over the given axis

    Args:

    * array (numpy.array):
        a boolean numpy.array in any dimension

    * axis (int):
        the axis over which the True blocks are calculated

    Returns:
        tuple of numpy.arrays holding the indices of the array from
        False to True and True to False, respectively
    """
    # roll the considered axis to the end
    a = np.rollaxis(array, axis, array.ndim)
    up = np.concatenate(
        (np.resize(a[..., 0], a.shape[:-1] + (1,)),
         np.logical_and(np.logical_not(a[..., :-1]), a[..., 1:])
         ), axis=a.ndim - 1)
    down = np.concatenate(
        (np.logical_and(a[..., :-1], np.logical_not(a[..., 1:])),
         np.resize(a[..., -1], a.shape[:-1] + (1,))
         ), axis=a.ndim - 1)

    if isinstance(a, ma.core.MaskedArray):
        up.mask = a.mask
    else:
        up = ma.masked_array(up, False)
    return up, down


def _make_time_dimension(start_date, end_date, period='year', align='center'):
    """
    create a temporal iris.coords.DimCoord

    Args:

    * start_date (string or datetime):
        the start date of the vector

    * end_date (string or datetime):
        the end date of the vector

    Kwargs:

    * period (string):
        create annual ('year'), seasonal ('season') or monthly ('month')
        intervals in days

    * align (string):
        put the datetime point at the first, center (default) or
        last of the period

    Returns:

    iris.coords.DimCoord with standard_name 'time' using a gregorian calendar
    """

    if not isinstance(start_date, datetime.datetime):
        start_date = parser.parse(start_date)
    if not isinstance(end_date, datetime.datetime):
        end_date = parser.parse(end_date)

    first_day_of_year = datetime.datetime(start_date.year, 1, 1)
    units = Unit('days since {:%F}'.format(first_day_of_year),
                 calendar='gregorian'
                 )

    if period == 'year':
        increment = relativedelta(years=1)
        start = datetime.datetime(start_date.year, 1, 1)
    elif period == 'season':
        increment = relativedelta(months=3)
        year = start_date.year
        month = start_date.month / 3 * 3
        if month in (0, 12):
            month = 12
            year -= 1
        start = datetime.datetime(year, month, 1)
        # add one month to the end_date in order to catch the december
        # of the last djf within the period
        end_date += relativedelta(months=1)
    elif period == 'month':
        increment = relativedelta(months=1)
        start = datetime.datetime(start_date.year, start_date.month, 1)
    else:
        raise ValueError("period must be one of year, season, month or day")

    if align == 'center':
        start = start + \
            ((start + increment - relativedelta(days=1)) - start) / 2
    elif align == 'last':
        start = start + increment - relativedelta(days=1)

    data = []
    while start < end_date:
        data.append((start - first_day_of_year).days)
        start += increment

    return DimCoord(np.array(data),
                    var_name='time',
                    standard_name='time',
                    units=units
                    )


def _create_cube(long_name='', var_name='', units='1',
                 dim_coords_and_dims=[], fill_value=-1):
    """
    Create an iris.cube.Cube given by its dimensions

    Kwargs:

    * long_name (string):
        Long description of the variable

    * var_name (string):
        Variable name

    * units (iris.unit.Unit or string):
        the unit of the variable

    * dims_coords_and_dims (list of iris.coords.DimCoord):
        the dimension of the variable

    Returns:
        An 'empty' iris.cube.Cube
    """
    shape = [x[0].shape[0] for x in dim_coords_and_dims]
    array = ma.ones(shape) * fill_value
    array.mask = True
    array.fill_value = fill_value

    if isinstance(units, str):
        units = Unit(units)
    return iris.cube.Cube(array, long_name=long_name, var_name=var_name,
                          units=units, dim_coords_and_dims=dim_coords_and_dims)
