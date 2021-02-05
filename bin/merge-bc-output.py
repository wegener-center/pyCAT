#!/usr/bin/env python

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
The output of the bias correction is not continuous in time and thus the
netCDF files have to be merged in order to get proper time-series

- Quantile mapping writes out files for each day of the year for
  the entire correction period
- Scaled distribution mapping writes out file for each month for 
  each correction period
"""
import datetime as dt
import logging
import os
import sys

import iris
from iris.experimental.equalise_cubes import equalise_attributes


def getargs(test=None):
    from argparse import ArgumentParser

    parser = ArgumentParser(description='merge bc output')
    parser.add_argument('--infile', type=str, required=True,
                        help='glob name of input files')
    parser.add_argument('--outfile-base', type=str, required=True,
                        help='name of output file, _YYYY.nc will be appended')
    parser.add_argument('--start-year', type=int, required=False,
                        help='start year')
    parser.add_argument('--end-year', type=int, required=False,
                        help='end year')
    parser.add_argument('-v', '--verbose', dest="log_level", const=logging.INFO,
                        action='store_const', default=logging.WARNING,
                        help='be verbose')
    parser.add_argument('-d', '--debug', dest="log_level", const=logging.DEBUG,
                        action='store_const', help='print debug messages')

    args = parser.parse_args(test)

    # end constraint will be excluded
    if args.end_year:
        args.end_year += 1

    return args


if __name__ == '__main__':
    args = getargs()

    # set the logger
    logging.basicConfig(
        stream=sys.stderr,
        format="%(asctime)s %(levelname)s: %(message)s",
        level=args.log_level, datefmt='%F %T'
    )

    iris.config.warnings.filterwarnings('ignore')

    logging.info('Running {}'.format(' '.join(sys.argv)))
    # if start/end year are given read the input file constrained
    constraint = None
    try:
        start = dt.datetime(args.start_year, 1, 1)
        constraint &= iris.Constraint(time=lambda cell: start <= cell.point)
    except:
        logging.debug('No start year given. Take it from input')
    try:
        end = dt.datetime(args.end_year, 1, 1)
        constraint &= iris.Constraint(time=lambda cell: cell.point < end)
    except:
        logging.debug('No end year given. Take it from input')

    logging.debug('Reading input')
    cl_orig = iris.load(args.infile, constraints=constraint)

    if not cl_orig:
        logging.debug('No data for the specified date range - exit')
        sys.exit(1)

    # rearrange cubes
    logging.debug('Rearrange data by days')
    cl = iris.cube.CubeList()
    for cube in cl_orig:
        for dc in cube.slices_over('time'):
            cl.append(dc)

    logging.debug('Equalize cube attributes')
    equalise_attributes(cl)
    cube = cl.merge_cube()

    time = cube.coord('time')
    start_date, end_date = \
        time.units.num2date([time.points[0], time.points[-1]])

    try:
        if start.year < start_date.year:
            logging.warning('Wanted start year not in data: {} < {}'.format(
                start.year, start_date.year))
    except:
        pass

    try:
        if end_date.year < end.year - 1:
            logging.warning('Wanted end year not in data {} < {}'.format(
                end_date.year, end.year - 1))
    except:
        pass

    logging.info('Writing output files to {}'.format(
        os.path.dirname(args.outfile_base)
    ))
    for year in range(start_date.year, end_date.year + 1):
        fn = '{}_{:4d}.nc'.format(args.outfile_base, year)
        logging.debug(' {}'.format(os.path.basename(fn)))
        constraint = iris.Constraint(time=lambda cell: cell.point.year == year)
        iris.save(cube.extract(constraint), fn, zlib=True, complevel=9)
