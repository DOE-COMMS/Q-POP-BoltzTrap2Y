=================
Theoretical Brief
=================

In this section, we summarize the formualtions for the Seebeck coefficient and Lorenz Number between the Boltzmann transport equation (BTE) and the proposed thermodynamic equation (TE). For more details, see the section of :ref:`Penn State's work`

- Seebeck coefficient

  .. math::
    
    S_{BTE} = - \frac{1}{eT}\int_{}^{}\frac{\int_{}^{}{\left\lbrack \varepsilon - \mu \right\rbrack f(1 - f)\mathbf{X}\left( \varepsilon \right)d\varepsilon}}{\int_{}^{}{f(1 - f)\mathbf{X}\left( \varepsilon \right)d\varepsilon}}

vs

  .. math::

    S_{TE} = - \frac{1}{eT}\int_{}^{}\frac{\int_{}^{}{\left\lbrack \varepsilon - \mu \right\rbrack f(1 - f)D\left( \varepsilon \right)d\varepsilon}}{\int_{}^{}{f(1 - f)D\left( \varepsilon \right)d\varepsilon}}

- Lorenz number

  .. math::

    L_{BTE} = \frac{1}{e^{2}T^{2}}\frac{\left\{ \int_{}^{}{f\left( 1 - f \right)\left\lbrack \varepsilon - \mu \right\rbrack^{2}}\mathbf{X}\left( \varepsilon \right)d\varepsilon - \frac{\left\lbrack \int_{}^{}{f\left( 1 - f \right)\left\lbrack \varepsilon - \mu \right\rbrack\mathbf{X}\left( \varepsilon \right)d\varepsilon} \right\rbrack^{2}}{\int_{}^{}{f\left( 1 - f \right)\mathbf{X}\left( \varepsilon \right)d\varepsilon}} \right\}}{\int_{}^{}{f\left( 1 - f \right)\mathbf{X}\left( \varepsilon \right)d\varepsilon}}

vs

  .. math::

    L_{TE} = \frac{1}{e^{2}T^{2}}\frac{\left\{ \int_{}^{}{f\left( 1 - f \right)\left\lbrack \varepsilon - \mu \right\rbrack^{2}}D\left( \varepsilon \right)d\varepsilon - \frac{\left\lbrack \int_{}^{}{f\left( 1 - f \right)\left\lbrack \varepsilon - \mu \right\rbrack D(\varepsilon)d\varepsilon} \right\rbrack^{2}}{\int_{}^{}{f\left( 1 - f \right)D(\varepsilon)d\varepsilon}} \right\}}{\int_{}^{}{f\left( 1 - f \right)D(\varepsilon)d\varepsilon}};\ or

  .. math::

    L_{TE} = \frac{k_{B}}{e^{2}}\frac{C_{el}}{n_{eff}}


- Relation between :math:`\mathbf{X}\left( \varepsilon \right)` and :math:`D\left( \varepsilon \right)`

:math:`\mathbf{X}\left( \varepsilon \right)` in the equation for :math:`S_{BTE}` or :math:`L_{BTE}` for the Boltzmann transport theory is called the transport distribution function [see Madsen, CPC 231, 140 (2018);
Scheidemantel, PRB 68, 125210(2003)]

  .. math::

    X^{\alpha\beta}\left( \varepsilon \right) = \frac{e^{2}}{V}\int_{}^{}{\sum_{i}^{}{v_{i}^{\alpha}\mathbf{(}\mathbf{k})v_{i}^{\beta}\mathbf{(}\mathbf{k}\mathbf{)}}\tau_{i,\mathbf{k}}\delta(\varepsilon - \varepsilon_{i}\mathbf{(}\mathbf{k}\mathbf{)})\frac{d\mathbf{k}}{8\pi^{3}}}

where the electron group velocity is

  .. math::

    v_{i}^{\alpha}\mathbf{(}\mathbf{k}\mathbf{) =}\frac{\mathbf{\partial}\varepsilon_{i}\mathbf{(}\mathbf{k}\mathbf{)}}{\mathbf{\partial}k^{\mathbf{\alpha}}} 

and :math:`\tau_{i,\mathbf{k}}` is called the relaxation time.

  .. math::

    D\left( \varepsilon \right) = \frac{1}{V}\int_{}^{}{\sum_{i}^{}{\delta(\varepsilon - \varepsilon_{i}\mathbf{(}\mathbf{k}))}\frac{d\mathbf{k}}{8\pi^{3}}}

:math:`D\left( \varepsilon \right)` in the equation for :math:`S_{TE}` or :math:`L_{TE}` for the `thermodynamic theory <https://doi.org/10.1103/PhysRevB.98.224101>`_ represents the electron density of states

* It is observed that when the transport distribution function, :math:`\mathbf{X}\left( \varepsilon \right)`, is replaced by the the electron density of states, :math:`D\left( \varepsilon \right)`, the expressions for the Seebeck coefficient and Lorenz number due to Boltzmann transport equation completely reduced to those due to thermodynamic theory.