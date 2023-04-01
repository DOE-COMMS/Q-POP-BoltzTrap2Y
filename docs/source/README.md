# README

## Quick run instruction

1. perform DFT calculations for, including
   a. 0-K electron energetics. These data contain the lattice structure, 0-K total energy, electronic density of states, electron energies in each *k*-mesh.
   b. interatomic force constants through supercell or linear response approach. These data will be used to calculate the lattice contribution to the thermodynamic properties.
2. invoke the BoltzTrap2Y package to calculate the kinetic/thermodynamic properties based on 0-K DFT results
3. when needed, invoke DFTTK to postprocess the interatomic force constants to get lattice contribution to the thermodynamic properties.
4. convert the data calculated by BoltzTrap2Y/DFTTK into the input to phase field simulation

## Run using Jupyter Notebook

To make ease of the package, we have prepared an extensive jupyter notebook script (which assumes that the step of DFT calculations have been completed) which can be excuted by the following steps

1. clicking the link [BoltzTrap2Y.ipynb](https://gitlab.com/yiwang62/BoltzTraP2/-/blob/20210126/BoltzTrap2Y.ipynb) followed by click the <img src="_static/download.png" /> icon in the right hand side of the web page. By default, a file named like ``BoltzTrap2Y.ipynb`` will be saved in your ``Downloands`` folder in the case of Windows computer. After that, go back to this page and run the codes using the free google notebook server by

2. clicking the link [jupyter notebook google](https://colab.research.google.com/notebooks/intro.ipynb) followed by uploading the downloaded code through clicking ``file->upload`` in the jupyter notebook google page.

3. clicking ``Runtime->Run all`` in the jupyter notebook google page


For more details, see the secction **Example by Jupyter Notebook**.

## Change logs in coding

As a numerical demonstration, we have implimented the computational procedure by python code through forking the [BoltzTraP2](https://gitlab.com/sousaw/BoltzTraP2) by @georg.madsen et al. The major extensions, which at present are only made for constant doping calculations, are briefed below:

   1. revised `interface.py` 
   2. extended `bandlib.py` into `bandlibEXT.py`.
   3. extended `io.py` into `ioEXT.py` 
   4. added ``remesh.py`` for calculations in low temperature (meshed increased to ~100,000 so that the calculated results can be good until a few of 10 K).
   5. improved the computational efficiency on the calculations of the chemical potential of electrons (Fermi levels)
   6. added `plotpng.py` for figure plot.


