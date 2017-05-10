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
# You should have received a copy of the GNU General Public License
# along with pyCAT. If not, see <http://www.gnu.org/licenses/>.
"""
Utility functions needed for generating temporal constraints on
:class:`iris.cube.Cubes <iris.cube.Cube>`
"""

import datetime as dt

from iris import Constraint
from iris.time import PartialDateTime


def generate_day_constraint_with_window(
        day_of_year, window, calendar='standard'):
    """
    generate two :class:`iris.Constraints <iris.Constraint>` for the time axis:

      1. for the exact day of the year over all years

      2. including all days over all years that lie within day_of_year ± window

    Args:

    * day_of_year (int):
        day of the year in the given calendar

    * window (int):
        the size of the temporal window around the given day (in days)

    * calendar (str):
        a supported calendar: standard (default), gregorian,
        proleptic_gregorian, noleap, 365_day, all_leap, 366_day, 360_day

    Returns:
        a 2-tuple of :class:`iris.Constraints <iris.Constraint>`
        on the time axis
    """
    if calendar in ['standard', 'gregorian', 'proleptic_gregorian',
                    'all_leap', '366_day']:
        # take a leap year to generate bounds
        start = dt.datetime(2000, 1, 1)
        year_start = PartialDateTime(month=1, day=1)
        year_end = PartialDateTime(month=12, day=31)
        begin = start + dt.timedelta(days=day_of_year - window)
        mid = start + dt.timedelta(days=day_of_year)
        end = start + dt.timedelta(days=day_of_year + window)
        begin = PartialDateTime(month=begin.month, day=begin.day)
        mid = PartialDateTime(month=mid.month, day=mid.day)
        end = PartialDateTime(month=end.month, day=end.day)
    elif calendar in ['noleap', '365_day']:
        # take a non-leap year to generate bounds
        start = dt.datetime(1999, 1, 1)
        year_start = PartialDateTime(month=1, day=1)
        year_end = PartialDateTime(month=12, day=31)
        begin = start + dt.timedelta(days=day_of_year - window)
        mid = start + dt.timedelta(days=day_of_year)
        end = start + dt.timedelta(days=day_of_year + window)
        begin = PartialDateTime(month=begin.month, day=begin.day)
        mid = PartialDateTime(month=mid.month, day=mid.day)
        end = PartialDateTime(month=end.month, day=end.day)
    elif calendar in ['360_day']:
        # construct the bounds manually
        year_start = PartialDateTime(month=1, day=1)
        year_end = PartialDateTime(month=12, day=30)
        month, day_of_month = divmod(day_of_year, 30)
        # add one as month and day counting starts at 1
        month += 1
        day_of_month += 1
        mid = PartialDateTime(month=month, day=day_of_month)
        start_day = day_of_month - window
        end_day = day_of_month + window
        begin = PartialDateTime(month=month, day=start_day) \
            if start_day >= 1 else \
            PartialDateTime(month=(month - 2) % 12 + 1,
                            day=start_day + 30)
        end = PartialDateTime(month=month, day=end_day) if end_day <= 30 else \
            PartialDateTime(month=month % 12 + 1, day=end_day - 30)
    else:
        raise ValueError("calendar '{}' not supported".format(calendar))

    day_constraint = Constraint(time=lambda cell: cell.point == mid)
    if begin.month <= end.month:
        window_constraint = Constraint(
            time=lambda cell: begin <= cell.point <= end)
    else:
        window_constraint = Constraint(
            time=lambda cell: year_start <= cell.point <= end or
            begin <= cell.point <= year_end)

    return day_constraint, window_constraint


def generate_year_constraint_with_window(year, window):
    """
    generate a :class:`iris.Constraint` on the time axis
    for specified year ± window

    Args:

    * year (int):
        centered year for the constraint

    * window (int):
        number of years around the given year

    Returns:

        an :class:`iris.Constraint` on the time-axis
    """
    first_year = PartialDateTime(year=year - window)
    last_year = PartialDateTime(year=year + window)
    return Constraint(time=lambda cell: first_year <= cell.point <= last_year)


def generate_month_constraint(month):
    """
    generate an :class:`iris.Constraint` on the time-axis for a specified month

    Args:

    * month (int):
       the desired month (1..jan, 12..dec)

    Returns:

       an :class:`iris.Constraint` on the time-axis
    """
    return Constraint(
        time=lambda cell: cell.point == PartialDateTime(month=month))
