# BoltzTraP2

[![BoltzTraP Logo](http://www.icams.de/content/wp-content/uploads/2014/09/boltztrap_200x58.png)](https://www.imc.tuwien.ac.at//forschungsbereich_theoretische_chemie/forschungsgruppen/prof_dr_gkh_madsen_theoretical_materials_chemistry/boltztrap/)
[![TU Wien Logo](https://www.imc.tuwien.ac.at/fileadmin/tuw/main/images/TU-Logo.gif)](https://www.tuwien.ac.at/)

BoltzTraP2 is a modern implementation of the [smoothed Fourier interpolation algorithm](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.38.2721) for electronic bands that formed the base of the original and widely used [BoltzTraP](http://www.sciencedirect.com/science/article/pii/S0010465506001305) code. One of the most typical uses of BoltzTraP is the calculation of thermoelectric transport coefficients as functions of temperature and chemical potential in the rigid-band picture. However, many other features are available, including 3D plots of Fermi surfaces based on the reconstructed bands. For more information, check out the [BoltzTraP2 tutorial](https://gitlab.com/sousaw/BoltzTraP2/-/wikis/tutorial).

## Prerequisites

BoltzTraP2 is a Python module, with a small performance-critical portion written in C++ and [Cython](http://cython.org/). BoltzTraP2's runtime requirements are Python version 3.5 or higher, and the Python libraries [NumPy](http://www.numpy.org/), [SciPy](https://www.scipy.org/), [matplotlib](https://matplotlib.org/), [spglib](https://atztogo.github.io/spglib/), [NetCDF4](https://github.com/Unidata/netcdf4-python) and [ASE](https://wiki.fysik.dtu.dk/ase/). All of them can be easily obtained from the [Python Package Index](https://pypi.python.org/pypi) (PyPI), using tools such as pip. They may also be bundled with Python distributions aimed at scientists, like [Anaconda](https://anaconda.org/), and with a number of Linux distributions. If pip is used to install BoltzTraP2, dependencies should be resolved automatically.

If available, BoltzTraP2 will also make use of [pyFFTW](http://hgomersall.github.io/pyFFTW/) (for faster Fourier transforms), [colorama](https://github.com/tartley/colorama) (to colorize some console output) and [VTK](https://www.vtk.org/) (to generate 3D representations). Those packages are not required, but they are recommended to be able to access the full functionality of BoltzTraP2.

Furthermore, compiling BoltzTraP2 from its sources requires a C++ compiler, and the development headers and libraries for Python. Cython is **not** required for a regular compilation.

## Compiling and install BoltzTraP2

The easiest way to get BoltzTraP2 is to run:

    $ pip install BoltzTraP2

This should take care of downloading and installing the dependencies as well.

Users installing from source must install the dependencies first and then run:

    $ python setup.py install

from the source directory. For finer-grained control, please see the output of these commands:

    $ python setup.py --help
    $ python setup.py --help-commands
    $ python setup.py install --help

The BoltzTraP2 installer supports

    $ python setup.py develop

which install the module through a set of symbolic links to the source directory, allowing users to immediately tests the effects of their changes to the code.

## Running the tests

BoltzTraP2 comes with a comprehensive set of unit and integration tests of its core functionality. To run those, install [pytest](https://docs.pytest.org) (also available through pip), change to the source directory and use the command

    $ pytest -v tests