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

import os
import warnings

import iris
import numpy as np
import pycat.esd.utils

try:
    from cf_units import Unit
except:
    from iris.unit import Unit

warnings.filterwarnings("ignore")


def _create_cube(calendar):
    """
    create a cube with time-vector from 1949-01-01 to 2050-12-31
    when asuming a gregorian calendar
    """
    length = 36525
    data = np.arange(length)
    cube = iris.cube.Cube(data)
    time = iris.coords.DimCoord(
        np.arange(length),
        standard_name="time",
        units=Unit("days since 1951-01-01", calendar=calendar),
        long_name="time",
        var_name="time",
    )
    cube.add_dim_coord(time, 0)
    return cube


def test_leap_year():
    calendar = "standard"
    leap_day = 59  # this is 29th of february in leap years
    cube = _create_cube(calendar)
    day_contraint, window_constraint = pycat.esd.utils.generate_day_constraint_with_window(
        leap_day, 1, calendar
    )
    c = cube.extract(day_contraint)

    # there must be 25 leap years
    assert c.shape[0] == 25


def test_qm():
    from pycat.io import Dataset
    from pycat.esd import QuantileMapping

    obs = Dataset("sample-data", "observation.nc")
    mod = Dataset("sample-data", "model*.nc")
    sce = Dataset("sample-data", "scenario*.nc")

    qm = QuantileMapping(obs, mod, sce)

    outfile = os.path.join(
        qm.work_dir, "quantile_mapping_tas_scenario-0_2021-2030_day-001.nc"
    )
    try:
        os.remove(outfile)
    except FileNotFoundError:
        pass

    qm.correct(1)

    assert os.path.exists(outfile)


def test_sdm():
    from pycat.io import Dataset
    from pycat.esd import ScaledDistributionMapping

    obs = Dataset("sample-data", "observation.nc")
    mod = Dataset("sample-data", "model*.nc")
    sce = Dataset("sample-data", "scenario*.nc")

    sdm = ScaledDistributionMapping(obs, mod, sce)

    outfile = os.path.join(
        sdm.work_dir, "scaled_distribution_mapping_tas_scenario-0_2021-2030_month-01.nc"
    )
    try:
        os.remove(outfile)
    except FileNotFoundError:
        pass

    sdm.correct(1)

    assert os.path.exists(outfile)
