=================
Penn State's work
=================

Thermodynamic respect of certain thermoelectric quantities
----------------------------------------------------------

This branch is partially inspired from our recent understanding on the thermoelectric effects. The first is `on the thermodynamic understing of the Seebeck coefficent <https://doi.org/10.1103/PhysRevB.98.224101>`_ which shows that:

Thermoelectric effects, measured by the Seebeck coefficients, refer to the phenomena in which a temperature difference or gradient imposed across a thermoelectric material induces an electrical potential difference or gradient, and vice versa, enabling the direct conversion of thermal and electric energies. All existing first-principles calculations of Seebeck coefficients have been based on the Boltzmann kinetic transport theory. In this work, we present a fundamentally different method for the first-principles calculations of Seebeck coefficients without using any assumptions on the electron scattering mechanism, in contrast to the traditional theory by Cutler and Mott that shows the dependence of the Seebeck coefficient on the scattering mechanisms. It is shown that the Seebeck coefficient is a well-defined thermodynamic quantity that can be determined from the change in the chemical potential of electrons induced by the temperature change and thus can be computed solely based on the electronic density of states through first-principles calculations at different temperatures. The proposed approach is demonstrated using the prototype PbTe and SnSe thermoelectric materials.

The second is `Lorenz Number and Electronic Thermoelectric Figure of Merit: Thermodynamics and Direct DFT Calculations <https://arxiv.org/abs/2010.00664>`_ which shows that:

The Lorenz number (L) contained in the Wiedemann-Franz law represents the ratio of two kinetic parameters of electronic charge carriers: the electronic contribution to the thermal conductivity (K_el) and the electrical conductivity (sigma), and can be expressed as LT=K_el/sigma where T is temperature. We demonstrate that the Lorenz number simply equals to the ratio of two thermodynamic quantities: the electronic heat capacity (c_el) and the electrochemical capacitance (c_N) through LT=c_el/c_N , a purely thermodynamic quantity, and thus it can be calculated solely based on the electron density of states of a material. It is shown that our thermodynamic formulation for the Lorenz number leads to: i) the well-known Sommerfeld value L=pi^2/3(k_B/e)^2 at the low temperature limit, ii) the Drude value L=3/2(k_B/e)^2 at the high temperature limit with the free electron gas model, and iii) possible higher values than the Sommerfeld limit for semiconductors. It is also demonstrated that the purely electronic contribution to the thermoelectric figure-of-merit can be directly computed using high-throughput DFT calculations without resorting to the computationally more expensive Boltzmann transport theory to the electronic thermal conductivity and electrical conductivity.


Revised output for the dope module
----------------------------------

For the ``interpolation.dope.trace`` file, the collumns are made of

    | collum 1, :math:`\mu-E_f(eV)` - electron chemical potential
    | collum 2, :math:`T(K)` - temperature
    | collum 3, :math:`N(e/uc)` - number of charge carries due to doping
    | collum 4, :math:`DOS(ef)[1/(Ha*uc)]` -
    | collum 5, :math:`S(V/K)` - Seebeck coefficients by BTE theory
    | collum 6, :math:`\sigma/tau0[1/(ohm*m*s)]` - trace of electrical conductivity
    | collum 7, :math:`RH[m**3/C]` -
    | collum 8, :math:`kappae/tau0[W/(m*K*s)]` -
    | collum 9, :math:`C_{\mu}[J/(mole-atom*K)]` - constant voltage heat capacity
    | collum 10, :math:`chi[m**3/mol]` -
    | collum 11, :math:`C_{el}[J/(mole-atom*K)]` - heat capacity with constant number of eletrons
    | collum 12, :math:`S_e(V/K)` - Seebeck coefficients by thermodynamic understanding, Phys. Rev. B, 98 (2018) 224101.
    | collum 13, :math:`n_{eff}(e/cm^3)` - effective carrier concentration
    | collum 14, :math:`L(W*ohm/K**2)` - Lorenz number by thermodynamic understanding
    | collum 15, :math:`\sigma_h` - hole electrical conductivity
    | collum 16, :math:`\sigma_e` - electon electrical conductivity concentration
    | collum 17, :math:`N_h` - testing
    | collum 18, :math:`N_e` - testing
    | collum 19, :math:`n_h(e/cm^3)` - hole carrier concentration
    | collum 20, :math:`n_e(e/cm^3)` - electron carrier

Thermodynamic formulations
--------------------------

- Fermi distribution

  .. math::

    f = \frac{1}{e^{\frac{\varepsilon - \mu}{k_{B}T}} + 1}

where :math:`\mu` is chemical potential, also called Fermi level. This quantity is temperature dependent, to be determined by the
number of electrons in the system as

  .. math::

    N_{el} = \int_{- \infty}^{+ \infty}{fD\left( \varepsilon \right)d\varepsilon} = \int_{- \infty}^{\varepsilon_{F}}{D\left( \varepsilon \right)d\varepsilon}

where :math:`D\left( \varepsilon \right)` is the electron density of states defined as

  .. math::

    D\left( \varepsilon \right) = \frac{1}{V}\int_{}^{}{\sum_{i}^{}{\delta(\varepsilon - \varepsilon_{i}\mathbf{(}\mathbf{k}))}\frac{d\mathbf{k}}{8\pi^{3}}}

- Formulations on the thermodynamic theory on the electronic contribution under constant doping conditions

 - internal energy

  .. math::

    E_{el}\left( T \right) = \int_{}^{}{fD\left( \varepsilon \right)d\varepsilon} - \int_{- \infty}^{\varepsilon_{F}}{D\left( \varepsilon \right)d\varepsilon}

 - entropy

  .. math::

    S_{el} = {- k}_{B}\int_{}^{}{\left\lbrack flnf + \left( 1 - f \right)\ln\left( 1 - f \right) \right\rbrack D\left( \varepsilon \right)d\varepsilon}

 - free energy

  .. math::

    F_{el}\left( T \right) = E_{el}\left( T \right) - TS_{el}

 - heat capacity

  .. math::

    C_{el} = \frac{1}{k_{B}T^{2}}\left\{ \int_{}^{}{f\left( 1 - f \right)\left\lbrack \varepsilon - \mu\left( T \right) \right\rbrack^{2}}D\left( \varepsilon \right)d\varepsilon - \frac{\left\lbrack \int_{}^{}{f\left( 1 - f \right)\left\lbrack \varepsilon - \mu\left( T \right) \right\rbrack D(\varepsilon)d\varepsilon} \right\rbrack^{2}}{\int_{}^{}{f\left( 1 - f \right)D(\varepsilon)d\varepsilon}} \right\}

- Heat capacity under constant voltage condition


  .. math::

    C_{\mu} = T\left( \frac{\partial S_{el}}{\partial T} \right)_{\mu}\mathrm{=}\frac{1}{k_{B}T^{2}}\left\{ \int_{}^{}{f\left( 1 - f \right)\left\lbrack \varepsilon - \mu\left( T \right) \right\rbrack^{2}}D\left( \varepsilon \right)d\varepsilon \right\}

which is related to its constant doping countpart by

  .. math::

    C_{\mu} = C_{el} +  \frac{n_{eff}e^{2}}{k_{B}}S_{e}^{2}

where the effective carrier density is defined as

  .. math::

    n_{eff} = \int_{- \infty}^{\infty}{f(1 - f)\ D\left( \varepsilon \right)d\varepsilon}

while the Seebeck coefficient is thermodynamically determined by

  .. math::

    S_{e} = - \frac{1}{en_{eff}T}\int_{- \infty}^{\infty}{\left( \varepsilon - \mu \right)\left( 1 - f \right)fD(\varepsilon)d\varepsilon}

Finally, the Lorenz Number can thermodynamicallu be calculated by

  .. math::

    L = \frac{k_{B}}{e^{2}}\frac{C_{el}}{n_{eff}}

Electrical conductivity
-----------------------

  .. math::

    \mathbf{\sigma} = \frac{1}{k_{B}T}\int_{- \infty}^{\infty}{f(1 - f)\ \mathbf{X}\left( \varepsilon \right)d\varepsilon}

where :math:`\mathbf{X}` is called sransport distribution function, see Madsen, CPC 231}, 140 (2018) and Scheidemantel, PRB 68, 125210(2003)

  .. math::

    X^{\alpha\beta}\left( \varepsilon \right) = \frac{e^{2}}{V}\int_{}^{}{\sum_{i}^{}{v_{i}^{\alpha}\mathbf{(}\mathbf{k})v_{i}^{\beta}\mathbf{(}\mathbf{k}\mathbf{)}}\tau_{i,\mathbf{k}}\delta(\varepsilon - \varepsilon_{i}\mathbf{(}\mathbf{k}\mathbf{)})\frac{d\mathbf{k}}{8\pi^{3}}}

where   :math:`\tau_{i,\mathbf{k}}` is relaxation time and the electron group velocity is

  .. math::

    v_{i}^{\alpha}\mathbf{(}\mathbf{k}\mathbf{) =}\frac{\mathbf{\partial}\varepsilon_{i}\mathbf{(}\mathbf{k}\mathbf{)}}{\mathbf{\partial}k^{\mathbf{\alpha}}}


- effective carrier mobility

  .. math::

    \mathbf{M}\mathbf{=}\frac{\mathbf{\sigma}}{n_{eff}e}

Separation of electron and hole contributions
---------------------------------------------

- electron carrier concentration

  .. math::

    n = \int_{\varepsilon_{F}}^{+ \infty}{fDd\varepsilon}

- hole carrier concentration

  .. math::

    p = \int_{- \infty}^{\varepsilon_{F}}{fDd\varepsilon}

- electrical conductivity

  .. math::

    \mathbf{\sigma}_{h} = \frac{1}{k_{B}T}\int_{- \infty}^{\varepsilon_{F}}{f(1 - f)\ \mathbf{X}\left( \varepsilon \right)d\varepsilon}

    \mathbf{\sigma}_{e} = \frac{1}{k_{B}T}\int_{\varepsilon_{F}}^{\infty}{f(1 - f)\ \mathbf{X}\left( \varepsilon \right)d\varepsilon}

- mobility

  .. math::

    \mathbf{M}_{h} = \frac{\mathbf{\sigma}_{h}}{ep}

    \mathbf{M}_{e} = \frac{\mathbf{\sigma}_{e}}{en}
