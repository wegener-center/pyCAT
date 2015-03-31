# Installing pyCAT from source

The latest release of pyCAT is available from
https://github.com/wegener-center/pyCAT

Building and running pyCAT requires a range of other libraries and python modules. Once you have resolved all these dependencies (see details below) change to the installation path and enter:

    python setup.py install


## Requirements

At the moment pyCAT requires Python 2 and is not compatible with Python 3.

Here you will find a list of external packages you will need to install before building and running pyCAT.

Many of these packages are available via a Linux package manager such as aptitude or yum. But we strongly encourage people to work within a virtual environment and install the latest stable releases of the packages using pip.

### The proposed installation

For a debian-based Linux distribution enter:

    sudo aptitude install python-virtualenv

to install the python virtual environment. For RedHat and derived distros
you will use:

    sudo yum install python-virtualenv

Change to an empty directory and create the virtual environment with the following command:

    virtualenv --system-site-packages env

where *env* is the freely chooseable name of your environment. As the installation of the python modules for the graphical environment is rather tricky we suggest to use the system-wide installation of python-gtk via aptitude or yum (which is probably already available) by enabling these packages with the given switch (--system-site-packages).

Afterwards activate your freshly installed environment by

    source ./env/bin/activate

and install the following software using pip with the --upgrade switch to allow the local installation parallel to the system installations:

    pip install --upgrade distribute
    pip install --upgrade numpy
    pip install --upgrade scipy
    pip install --upgrade biggus
    pip install --upgrade cython
    pip install --upgrade pyshp
    pip install --upgrade shapely
    pip install cartopy==0.11.0
    pip install --upgrade matplotlib
    pip install --upgrade pyke

This will install these packages along with some dependencies into *env*.

Unfortunately *iris* is not available on pypi.python.org. Thus you will have to download it from github:

    git clone https://github.com/SciTools/iris

and install it by changing into the iris directory and running:

    python setup.py install


### Optional software

To use the powerful interactive shell and debugger you can install:

    pip install --upgrade ipython
    pip install --upgrade ipdb
    
To run the pyCAT test suite you will have to install:

    pip install --upgrade nose

For building the full documentation you need:

    pip install --upgrade Sphinx

The udunits2 library will help at the conversion of physical units. You can install it via your package manager, e.g. aptitude:

    sudo aptitude install udunits-bin
