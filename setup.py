#!/usr/bin/env python

from distutils.core import setup

setup(name='pyCAT',
      version='0.0.1.dev1',
      description='Climate Analysis Tool',
      license='GPLv3',
      author='Armin Leuprecht',
      author_email='armin.leuprecht@uni-graz.at',
      packages=['pycat', 'pycat.analysis', 'pycat.io', 'pycat.downscaling'],
      classifiers=[
          'Development Status :: 1 - Planning',
          'Topic :: Scientific/Engineering :: Physics',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.6',
          'Programming Language :: Python :: 2.7',
      ],
      keywords = ['climate', ],
      install_requires = [
          'numpy', 'scipy', 'NetCDF4', 'biggus', 'cython', 'pyshp', 'shapely',
          'statsmodels', 'cartopy', 'pyke', 'pillow', 'matplotlib'
      ]
      
     )
