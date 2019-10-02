Installing pyCAT from source
----------------------------

The latest release of pyCAT is available from
https://github.com/wegener-center/pyCAT

Building and running pyCAT requires a range of other libraries and
python modules. Once you have resolved all these dependencies (see
details below) change to the download path and enter:

::

    python setup.py install

Requirements
~~~~~~~~~~~~

pyCAT requires Python 3.7.

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
conda-forge channel and install the dependencies

::

    ENV_NAME="pycat"
    conda create --name $ENV_NAME python=3.7
    conda config --add channels conda-forge
    conda activate $ENV_NAME
    conda install --file conda-requirements.txt

Afterwards you can install pyCAT

::
   
    python setup.py install


Installation using python-virtualenv
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

    virtualenv pycat

where *pycat* is the freely chooseable name of your environment.
Afterwards activate your freshly installed environment by

::

    source ./pycat/bin/activate

and install the following software using pip (iris dependencies):

::

    pip install -r requirements.txt

This will install these packages along with some dependencies into
*pyenv*.

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

    
