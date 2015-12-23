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

def generate_day_constraint_with_window(day_of_year, window, calendar):
    """
    generate two `iris.Constraint`s for the time axis:
    1. for the exact day of the year over all years
    2. including all days over all years that lie within day_of_year ± window

    Args:

    * day_of_year (int):
        day of the year in the given calendar

    * window (int):
        the size of the temporal window around the given day (in days)

    * calendar (string):
        a supported calendar (standard, gregorian, proleptic_gregorian,
        noleap, 365_day, all_leap, 366_day, 360_day)

    Returns:
        a 2-tuple of iris.Constraint
    """
    from iris.time import PartialDateTime
    from iris import Constraint
    import datetime

    if calendar in ['standard', 'gregorian', 'proleptic_gregorian', 'all_leap', '366_day']:
        # take a leap year to generate bounds
        start = datetime.datetime(2000,1,1)
        year_start = PartialDateTime(month=1, day=1)
        year_end = PartialDateTime(month=12, day=31)
        begin = start + datetime.timedelta(days=day_of_year-window)
        mid = start + datetime.timedelta(days=day_of_year)
        end = start + datetime.timedelta(days=day_of_year+window)
        begin = PartialDateTime(month=begin.month, day=begin.day)
        mid = PartialDateTime(month=mid.month, day=mid.day)
        end = PartialDateTime(month=end.month, day=end.day)
    elif calendar in ['no_leap', '365_day']:
        # take a non-leap year to generate bounds
        start = datetime.datetime(1999,1,1)
        year_start = PartialDateTime(month=1, day=1)
        year_end = PartialDateTime(month=12, day=31)
        begin = start + datetime.timedelta(days=day_of_year-window)
        mid = start + datetime.timedelta(days=day_of_year)
        end = start + datetime.timedelta(days=day_of_year+window)
        begin = PartialDateTime(month=begin.month, day=begin.day)
        mid = PartialDateTime(month=mid.month, day=mid.day)
        end = PartialDateTime(month=end.month, day=end.day)
    elif calender in ['360_days']:
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
        begin = PartialDateTime(month=month, day=start_day) if start_day>=1 else \
                PartialDateTime(month=(month-2)%12 + 1, day=start_day+30)
        end = PartialDateTime(month=month, day=end_day) if end_day<=30 else \
              PartialDateTime(month=month%12 + 1, day=end_day-30)
    else:
        raise ValueError("calendar '{}' not supported".format(calendar))

    day_constraint = Constraint(time=lambda cell: cell.point == mid)
    if begin.month <= end.month:
        window_constraint = Constraint(time=lambda cell: begin <= cell.point <= end)
    else:
        window_constraint = Constraint(time=lambda cell: year_start <= cell.point <= end or
                                       begin <= cell.point <= year_end)

    return day_constraint, window_constraint

