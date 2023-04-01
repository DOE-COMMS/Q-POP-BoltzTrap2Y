#    BoltzTraP2, a program for interpolating band structures and calculating
#                semi-classical transport coefficients.
#    Copyright (C) 2017-2020 Georg K. H. Madsen <georg.madsen@tuwien.ac.at>
#    Copyright (C) 2017-2020 Jesús Carrete <jesus.carrete.montana@tuwien.ac.at>
#    Copyright (C) 2017-2020 Matthieu J. Verstraete <matthieu.verstraete@ulg.ac.be>
#    Copyright (C) 2018-2019 Genadi Naydenov <gan503@york.ac.uk>
#    Copyright (C) 2020 Gavin Woolman <gwoolma2@staffmail.ed.ac.uk>
#    Copyright (C) 2020 Roman Kempt <roman.kempt@tu-dresden.de>
#
#    This file is part of BoltzTraP2.
#
#    BoltzTraP2 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    BoltzTraP2 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with BoltzTraP2.  If not, see <http://www.gnu.org/licenses/>.

# Test the Fermi-Dirac module

import numpy as np
import pytest

import BoltzTraP2
from BoltzTraP2.fd import *
from BoltzTraP2.units import BOLTZMANN
from BoltzTraP2.units import eV


def test_FD_0K():
    """The Fermi-Dirac occupancies at 0 K should be a step function."""
    # Test that the result is correct for mu = 0
    e = np.linspace(-5., 5., num=1001)
    r = FD(e, 0., 0.)
    assert np.allclose(r[e < 0.], 1.)
    assert np.allclose(r[e == 0.], 0.5)
    assert np.allclose(r[e > 0.], 0.)
    # Test that the result is correct for mu != 0
    mu = 3.
    e += mu
    r = FD(e, mu, 0.)
    assert np.allclose(r[e < mu], 1.)
    assert np.allclose(r[e == mu], 0.5)
    assert np.allclose(r[e > mu], 0.)


def test_FD_T():
    """The Fermi-Dirac occupancies should reproduce tabulated values."""
    kBT = BOLTZMANN * 300.
    r = FD(0., 0., kBT)
    assert np.allclose(r, 0.5)
    e = 30. * eV * 1e-3  # 30 meV
    r = FD(e, 0., kBT)
    assert np.allclose(r, 0.23858517713617)
    r = FD(-e, 0., kBT)
    assert np.allclose(r, 0.76141482286383)
    r = FD(e, 15. * eV * 1e-3, kBT)
    assert np.allclose(r, 0.358880600979885)
    r = FD(e, -15. * eV * 1e-3, kBT)
    assert np.allclose(r, 0.149226850069459)
    r = FD(np.inf, 0., kBT)
    assert r == 0.
    r = FD(-np.inf, 0., kBT)
    assert r == 1.


def test_dFDdx():
    """The Fermi-Dirac derivative should agree with tabulated values."""
    x = np.linspace(-20., 20., num=1001)
    r = dFDdx(x)
    assert np.allclose(r[x == 0.], -.25)
    assert np.all(r <= 0.)
    assert np.all(r[x != 0.] > -.25)
    assert r[0] == 0.
    assert r[-1] == 0.
    assert dFDdx(-1.)
    assert np.allclose(dFDdx(1.), -0.19661193)
    assert np.allclose(dFDdx(-1.), -0.19661193)


def test_dFDde():
    """The Fermi-Dirac derivative should scale with temperature correctly."""
    x = np.linspace(-20., 20., num=1001)
    assert np.allclose(dFDdx(x), dFDde(x * 33., 0., 33.) * 33.)
    assert np.allclose(dFDdx(x), dFDde(x * 33. + 10., 10., 33.) * 33.)
