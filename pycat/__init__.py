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
"""
pyCAT
=====

A package for analyzing climate data

"""

import logging
import pycat.config

# pyCAT version
__version__ = '0.0.1-DEV'

# restrict imports when using "from pycat import *"
__all__ = ['site_configuration']

# initialize the site configuration
site_configuration = {}
