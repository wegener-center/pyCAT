Installing pyCAT from source
----------------------------

The latest release of pyCAT is available from
https://github.com/wegener-center/pyCAT

Building and running pyCAT requires a range of other libraries and
python modules. Once you have resolved all these dependencies (see
details below) change to the installation path and enter:

::

    python setup.py install

Requirements
~~~~~~~~~~~~

At the moment pyCAT requires Python 2.7 or Python 3.4.

Here you will find a list of external packages you will need to install
before building and running pyCAT.

Many of these packages are available via a Linux package manager such as
aptitude or yum. But we strongly encourage people to work within a
virtual environment and install the latest stable releases of the
packages using pip. In fact the easiest way uses conda as virtual
environment.

The proposed installation (using conda)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Download and install conda from http://conda.pydata.org/miniconda.html.

Once conda is running you just need to create an environment, add the
scitools-channel and install the dependencies

::

    ENV_NAME="env"
    conda create --name $ENV_NAME python=2.7
    conda config --add channels conda-forge
    source activate $ENV_NAME
    conda install --file conda-requirements.txt

Afterwards you can install pyCAT

::
   
    python setup.py install


Installation using python-virtualenv (the long way)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a debian-based Linux distribution enter:

::

    sudo aptitude install python-virtualenv

to install the python virtual environment. For RedHat and derived
distros you will use:

::

    sudo yum install python-virtualenv

Change to an empty directory and create the virtual environment with the
following command:

::

    virtualenv env

where *env* is the freely chooseable name of your environment. As the
installation of the python modules for the graphical environment is
rather tricky we suggest to use the system-wide installation of
python-gtk via aptitude or yum (which is probably already available) by
enabling these packages with the given switch (–system-site-packages).

Afterwards activate your freshly installed environment by

::

    source ./env/bin/activate

and install the following software using pip (iris dependencies):

::

    pip install distribute
    pip install python-dateutil
    pip install numpy
    pip install scipy
    pip install NetCDF4 (probably you have to set the CFLAGS=-I/usr/include/hdf5/serial)
    pip install biggus
    pip install cython
    pip install pyshp
    pip install shapely
    pip install statsmodels
    pip install cartopy
    pip install cf_units
    pip install pyke --allow-external pyke  --allow-unverified pyke
    pip install pillow

In order to use the matplotlib with the GTK backend you can install the
python matplotlib in your virtual environment by install all the GTK
stuff (which is rather laborius) or linking the required python packages
(e.g. on a debian jessie using python 2.7):

::

    ln -sf /usr/lib/python2.7/dist-packages/{glib,gobject,gtk-2.0,pygtk.py,pygtk.pth} $VIRTUAL_ENV/lib/python2.7/site-packages
    ln -sf /usr/lib/pymodules/python2.7/cairo $VIRTUAL_ENV/lib/python2.7/site-packages
    pip install matplotlib

This will install these packages along with some dependencies into
*env*.

Unfortunately *iris* is not available on pypi.python.org. Thus you will
have to download it from github:

::

    git clone https://github.com/SciTools/iris

and install it by changing into the iris directory and running:

::

    python setup.py install

Optional software
^^^^^^^^^^^^^^^^^

To use the powerful interactive shell and debugger you can install:

::

    pip install ipython
    pip install ipdb

To run the pyCAT test suite you will have to install:

::

    pip install nose

For building the full documentation you need:

::

    pip install Sphinx

The udunits2 library will help at the conversion of physical units. You
can install it via your package manager, e.g. aptitude:

::

    sudo aptitude install udunits-bin

    
