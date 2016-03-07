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
this file provides the actual correction functions with the following singature:
obs_cube (Cube), mod_cube(Cube), sce_cubes(list of Cubes), *args, **kwargs
"""

import numpy as np

def _quantile_mapping(obs_cube, mod_cube, sce_cubes, *args, **kwargs):
    """
    Quantile Mapping
    ----------------

    apply quantile mapping to all scenario cubes using the distributions
    of obs_cube and mod_cube

    Args:

    * obs_cube (iris.cube.Cube):
        the observational data

    * mod_cube (iris.cube.Cube):
        the model data at the reference period

    * sce_cubes (iris.cube.CubeList):
        the scenario data that shall be corrected
    """
    from statsmodels.tools.tools import ECDF
    
    cell_iterator = np.nditer(obs_cube.data[0], flags=['multi_index'])
    while not cell_iterator.finished:
        index_list = list(cell_iterator.multi_index)
        cell_iterator.iternext()

        index_list.insert(0,0)
        index = tuple(index_list)
        if obs_cube.data.mask[index]:
            continue

        index_list[0] = slice(0,None,1)
        index = tuple(index_list)
        obs_data = obs_cube.data[index]
        mod_data = mod_cube.data[index]
        mod_ecdf = ECDF(mod_data)

        for sce_cube in sce_cubes:
            sce_data = sce_cube[index].data
            p = mod_ecdf(sce_data)*100
            corr = np.percentile(obs_data, p) - \
                   np.percentile(mod_data, p)
            sce_cube.data[index] += corr
