============
Installation
============

It is recommended to install the package under the `anaconda <https://docs.anaconda.com/anaconda/install/>`_ environment. Under the anaconda prompt, one can create a preferred directory and then run
 
  .. code-block:: bash

    git clone https://gitlab.com/yiwang62/BoltzTraP2.git
    cd BoltzTraP2/
    setenv CC g++

Edit the ``setup.py`` file, depending on your CC environment, use

    | required_compile_flags = ["-std=c++11"], or
    | required_compile_flags = ["-std=c++0x"]

Last, run

  .. code-block:: bash

    pip install -e . #for develop version (recomended)
    #if failed, try
    python setup.py install #or
    python setup.py develop #for develop version (recomended)
