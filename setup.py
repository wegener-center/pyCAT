#!/usr/bin/env python

from distutils.core import setup

from pycat import __version__ as version

setup(name='pyCAT',
      version=version,
      description='Climate Analysis Tool',
      license='GPLv3',
      author='Armin Leuprecht',
      author_email='armin.leuprecht@uni-graz.at',
      scripts=['bin/merge-bc-output.py'],
      packages=['pycat', 'pycat.analysis', 'pycat.io', 'pycat.esd'],
      classifiers=[
          'Development Status :: 1 - Planning',
          'Topic :: Scientific/Engineering :: Physics',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.7',
      ],
      keywords=['climate', ],
      tests_require=['nose'],
      )
