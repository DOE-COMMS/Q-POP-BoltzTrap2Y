=========
Utilities
=========

Mixing two phases
-----------------

The mechanism to mix the DOSs between two structures and then calculate the chemical potential is given below

1. shift the DOS using Fermi energy as zero

  .. math::

    D_{AFM - b}^{'}\left( \varepsilon \right) = D_{AFM - b}\left( \varepsilon + \varepsilon_{F}^{AFM - b} \right)

    D_{AFM - a}^{'}\left( \varepsilon \right) = D_{AFM - a}\left( \varepsilon + \varepsilon_{F}^{AFM - a} \right)

2. mix the DOSs

  .. math::

    D_{\text{mix}}^{'}\left( \varepsilon \right) = (1 - x)*D_{AFM - b}^{'}\left( \varepsilon \right)+{x*D}_{AFM - a}^{'}\left( \varepsilon \right)

3. calculate the chemical potential of electrons

  .. math::

    \int_{- \infty}^{\infty}{\text{fD}_{\text{mix}}^{'}\left( \varepsilon \right)\text{dε}} = \int_{- \infty}^{0}{D_{\text{mix}}^{'}\left( \varepsilon \right)\text{dε}}

under Fermi distribution

  .. math::

    \mathbf{f} = \frac{1}{e^{\frac{\varepsilon - \mu}{k_{B}T}} + 1}

Usage:

  .. code-block:: bash

    python utilities/dosmixAPI.py -d0 dir0/ -d1 dir1/ -nC 11 -nT 101

The output to the text file ``thermo.out`` contains data as  functions of phase compostion ``x`` and ``T``. These data can be plotted following the section :ref:`Example by Jupyter Notebook`
 

