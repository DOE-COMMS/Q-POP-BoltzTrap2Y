# -*- coding: utf-8 -*
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

# Test the code dealing with densities of states and Onsager coefficients

import os
import os.path
from collections import namedtuple

import ase
import numpy as np
import scipy as sp
import scipy.linalg as la
import pytest

import BoltzTraP2
import BoltzTraP2.io
import BoltzTraP2.dft
import BoltzTraP2.fite
import BoltzTraP2.bandlib
from BoltzTraP2.units import BOLTZMANN
from BoltzTraP2.units import BOLTZMANN_SI

mydir = os.path.abspath(os.path.dirname(__file__))
datadir = os.path.join(os.path.dirname(mydir), "data")


@pytest.fixture
def seed_pgenerator():
    """Return the pseudo-random number generator to a predictable state."""
    np.random.seed(1337)


@pytest.fixture
def parabolic_band():
    """Generate a set of electron energies from a single parabolic band."""
    # Take samples from a uniform grid
    L = 5.
    x = np.linspace(-L, L, num=201)
    y = np.linspace(-L, L, num=201)
    z = np.linspace(-L, L, num=201)
    x, y, z = points = np.meshgrid(x, y, z)
    # And discard all points outside of the inscribed sphere to avoid
    # biasing the DOS.
    k2 = x**2 + y**2 + z**2
    k2 = k2[k2 < L * L]
    # Return the squared modulus of each point as the energy, with a
    # prefactor that brings the values into a more reasonable range.
    return k2.reshape((-1, 1)) / k2.max() * 1e-1


@pytest.fixture()
def si_bands():
    """Load saved data for Si and reconstruct the bands."""
    data_noder = BoltzTraP2.dft.DFTData(os.path.join(datadir, "Si"),
                                        derivatives=False)
    data_withder = BoltzTraP2.dft.DFTData(os.path.join(datadir, "Si"),
                                          derivatives=True)
    equivalences = BoltzTraP2.sphere.get_equivalences(
        data_noder.atoms, data_noder.magmom, 2 * len(data_noder.kpoints))
    DataAndEquivalences = namedtuple(
        "DataAndEquivalences", ["data_noder", "data_withder", "equivalences"])
    si_data = DataAndEquivalences(data_noder, data_withder, equivalences)
    si_coefficients = np.load(os.path.join(mydir, "Si_fitde3D.npz"))
    return BoltzTraP2.fite.getBTPbands(si_data.equivalences,
                                       si_coefficients["noder"],
                                       si_data.data_noder.get_lattvec(),
                                       curvature=True)


def test_suggest_nbins(seed_pgenerator):
    """The histogram code should suggest the same number of bins as NumPy."""
    for i in range(1, 10):
        a = np.random.rand(i * 1000) * 10
        domain = (a.min(), a.max())
        h = np.histogram(a, range=domain, bins="auto")
        assert len(h[0]) == BoltzTraP2.bandlib._suggest_nbins(a, domain)


def test_DOS_parabolic(parabolic_band):
    """The DOS() function should work for a parabolic band."""
    npts = 60
    # Numerical estimate of the DOS
    dos = BoltzTraP2.bandlib.DOS(parabolic_band, npts=npts)
    # Test that the result has the right shape
    assert dos[0].size == npts
    assert dos[1].size == npts
    # Test that the result is normalized to 1
    de = dos[0][1] - dos[0][0]
    assert np.allclose(dos[1].sum() * de, 1.0)
    # Analytic result, unnormalized
    ref = np.sqrt(dos[0])
    # Renormalize the reference
    ref /= ref.sum() * de
    # Test that the difference in norm is small enough
    assert la.norm(dos[1] - ref) * de < 1e-3


def test_DOS_Si(si_bands):
    """The DOS() function should work for real Si data."""
    eband, vvband, cband = si_bands
    dos = BoltzTraP2.bandlib.DOS(eband.T)
    de = dos[0][1] - dos[0][0]
    assert np.allclose(dos[1].sum() * de, 6.0)
    reference = np.load(os.path.join(mydir, "Si_BTPdos.npz"))
    assert np.allclose(reference["dose"], dos[0])
    assert np.allclose(reference["dos"], dos[1])


def test_BTPDOS_Si(si_bands):
    """The BTPDOS() function should work for real Si data and several different
    scattering models and combinations of parameters.
    """
    reference = np.load(os.path.join(mydir, "Si_BTPdos.npz"))
    eband, vvband, cband = si_bands
    btpdos = BoltzTraP2.bandlib.BTPDOS(eband, vvband)
    assert np.allclose(btpdos[0], reference["dose"])
    assert np.allclose(btpdos[1], reference["dos"])
    assert np.allclose(btpdos[2], reference["vvdos"])
    assert btpdos[3] is None
    btpdos = BoltzTraP2.bandlib.BTPDOS(eband, vvband, cband)
    assert np.allclose(btpdos[0], reference["dose"])
    assert np.allclose(btpdos[1], reference["dos"])
    assert np.allclose(btpdos[2], reference["vvdos"])
    assert np.allclose(btpdos[3], reference["cdos"])
    btpdos = BoltzTraP2.bandlib.BTPDOS(eband,
                                       vvband,
                                       cband,
                                       scattering_model="uniform_lambda")
    reference = np.load(os.path.join(mydir, "Si_BTPdos_lambda.npz"))
    assert np.allclose(btpdos[0], reference["dose"])
    assert np.allclose(btpdos[1], reference["dos"])
    assert np.allclose(btpdos[2], reference["vvdos"])
    assert np.allclose(btpdos[3], reference["cdos"])
    # lambda_to_tau() is also part of the test
    tau = BoltzTraP2.bandlib.lambda_to_tau(vvband, np.ones_like(eband))
    btpdos = BoltzTraP2.bandlib.BTPDOS(eband,
                                       vvband,
                                       cband,
                                       scattering_model=tau)
    assert np.allclose(btpdos[0], reference["dose"])
    assert np.allclose(btpdos[1], reference["dos"])
    assert np.allclose(btpdos[2], reference["vvdos"])
    assert np.allclose(btpdos[3], reference["cdos"])


def test_smoothen_DOS(parabolic_band):
    """smoothen_DOS() should work for the parabolic band model."""
    Ts = (50., 100., 200., 300.)
    # Start with a very noisy estimate of the DOS
    dos = BoltzTraP2.bandlib.DOS(parabolic_band, npts=500)
    de = dos[0][1] - dos[0][0]
    # Analytic reference
    ref = np.sqrt(dos[0])
    ref /= ref.sum() * de
    # Smoothen at several different temperatures and check that normalization
    # is preserved.
    for T in Ts:
        sdos = BoltzTraP2.bandlib.smoothen_DOS(dos[0], dos[1], T)
        assert np.allclose(sdos.sum() * de, 1.0)
    # There can be significant distortions close to the band edges, so
    # only the central part of the smoother dos is tested.
    for T in Ts:
        sdos = BoltzTraP2.bandlib.smoothen_DOS(dos[0], dos[1], T)
        assert np.abs((ref - sdos)[20:-20] * de).max() < 1e-4


def test_calc_N_parabolic(parabolic_band):
    """calc_N() should work for a parabolic band."""
    em = parabolic_band.max()
    dos = BoltzTraP2.bandlib.DOS(parabolic_band)
    # Test for 0 K, with an analytical reference
    mur = np.linspace(0., em, num=100, endpoint=False)[1:]
    r = np.empty_like(mur)
    for imu, mu in enumerate(mur):
        r[imu] = BoltzTraP2.bandlib.calc_N(dos[0],
                                           dos[1],
                                           mu,
                                           0.,
                                           dosweight=1.)
    ref = -(mur / em)**1.5
    assert np.allclose(r, ref, atol=5e-3)
    for imu, mu in enumerate(mur):
        r[imu] = BoltzTraP2.bandlib.calc_N(dos[0], dos[1], mu, 0.)
    ref = -2. * (mur / em)**1.5
    assert np.allclose(r, ref, atol=1e-2)
    # Test that |N| is a monotonically increasing function of T
    mu = .05 * em
    Tr = np.linspace(50., 600., num=100)
    r = np.empty_like(Tr)
    for iT, T in enumerate(Tr):
        r[iT] = BoltzTraP2.bandlib.calc_N(dos[0], dos[1], mu, T)
    assert np.all(np.diff(r) < 0.)


def test_solve_for_mu_parabolic():
    """The intrinsic chemical potential should be at the center of the gap in
    a toy model with two identical parabolic bands, at any temperature.
    """
    Tr = np.arange(0., 700, 50.)
    # Create an analytical DOS for a two-band model with a gap.
    gap = 1e-3
    e = np.linspace(-1e-1, 1e-1, num=1001)
    de = e[1] - e[0]
    dos = np.sqrt(np.maximum(0., np.abs(e) - .5 * gap))
    dos *= 2. / dos.sum() / de
    # Check that mu0 stays at 0 for every temperature
    for T in Tr:
        mu0 = BoltzTraP2.bandlib.solve_for_mu(e,
                                              dos,
                                              1.,
                                              T,
                                              dosweight=1.,
                                              try_center=True)
        assert np.allclose(mu0, 0.)
        # Even if we lift the constraint of mu0 belonging to the histogram
        mu0 = BoltzTraP2.bandlib.solve_for_mu(e,
                                              dos,
                                              1.,
                                              T,
                                              dosweight=1.,
                                              refine=True,
                                              try_center=True)
        assert np.allclose(mu0, 0.)
    # Test that the intrinsic chemical potential moves with temperature as
    # expected when the two bands have different effective masses.
    e = np.linspace(-2e-2, 1e-2, num=1000)
    de = e[1] - e[0]
    dos = np.zeros_like(e)
    dos[e > 0.] += np.sqrt(np.maximum(0., np.abs(e[e > 0.]) - .5 * gap))
    dos[e > 0.] /= dos[e > 0.].sum() * de
    dos[e < 0.] += np.sqrt(np.maximum(0., np.abs(e[e < 0.]) - .5 * gap))
    dos[e < 0.] /= dos[e < 0.].sum() * de
    reference = np.array([
        0.0, -9.00900900901e-05, -0.00018018018018, -0.00027027027027,
        -0.00039039039039, -0.000510510510511, -0.000600600600601,
        -0.000720720720721, -0.000810810810811, -0.000930930930931,
        -0.00102102102102, -0.00111111111111, -0.00123123123123
    ])
    refined_reference = np.array([
        0.0, -8.633429620108184e-05, -0.0001839359725226431,
        -0.0002850740585707863, -0.0003866346009555535, -0.0004957067189869257,
        -0.000604356386707219, -0.0007092502699324298, -0.0008145665937995095,
        -0.0009194604801426412, -0.0010247768008917982, -0.0011259149118833117
    ])
    for r, rr, T in zip(reference, refined_reference, Tr):
        mu0 = BoltzTraP2.bandlib.solve_for_mu(e,
                                              dos,
                                              1.,
                                              T,
                                              dosweight=1.,
                                              try_center=True)
        assert np.allclose(r, mu0, atol=1e-6)
        mu0 = BoltzTraP2.bandlib.solve_for_mu(e,
                                              dos,
                                              1.,
                                              T,
                                              dosweight=1.,
                                              refine=True,
                                              try_center=True)
        assert np.allclose(rr, mu0, atol=1e-6)


def test_calc_cv_parabolic():
    """calc_cv() should work for a single parabolic band."""
    e = np.linspace(0., 1e-1, num=10000)
    de = e[1] - e[0]
    dos = np.sqrt(e)
    dos /= dos.sum() * de
    Tr = np.arange(1., 300., 5.)
    # mu is placed well into the band but far away from its maximum energy,
    # so the result must be linear to a very good approximation
    cv = BoltzTraP2.bandlib.calc_cv(e, dos, [0.01], Tr,
                                    dosweight=1.) / BOLTZMANN_SI
    assert cv.shape == (len(Tr), 1)
    cv = cv[:, 0]
    p = np.polyfit(Tr, cv, 1)
    line = np.polyval(p, Tr)
    assert ((line - cv)**2).sum() / cv.size / cv.var() < 1e-4


def test_calc_cv_Si(si_bands):
    """calc_cv() should reproduce tabulated results for Si."""
    Tr = np.arange(100., 300., 5.)
    mur = np.linspace(0., 0.4, num=100)
    eband, vvband, cband = si_bands
    e, dos = BoltzTraP2.bandlib.DOS(eband.T, npts=2000)
    cv = BoltzTraP2.bandlib.calc_cv(e, dos, mur, Tr,
                                    dosweight=1.) / BOLTZMANN_SI
    assert cv.shape == (len(Tr), len(mur))
    ref = np.load(os.path.join(mydir, "Si_cv.npz"))
    assert np.allclose(cv, ref["cv"])


def test_fermiintegrals_Si(si_bands):
    """fermiintegrals() should reproduce tabulated results for Si."""
    Tr = np.arange(100., 300., 5.)
    mur = np.linspace(0., 0.4, num=100)
    eband, vvband, cband = si_bands
    dose, dos, vvdos, cdos = BoltzTraP2.bandlib.BTPDOS(eband,
                                                       vvband,
                                                       cband=cband)
    N, L0, L1, L2, L11 = BoltzTraP2.bandlib.fermiintegrals(dose,
                                                           dos,
                                                           vvdos,
                                                           mur,
                                                           Tr,
                                                           cdos=cdos)
    ref = np.load(os.path.join(mydir, "Si_fermiintegrals.npz"))
    assert np.allclose(ref["N"], N)
    assert np.allclose(ref["L0"], L0)
    assert np.allclose(ref["L1"], L1)
    assert np.allclose(ref["L2"], L2)
    assert np.allclose(ref["L11"], L11)


def test_calc_Onsager_coefficients(si_bands):
    """calc_Onsager_coefficients() should reproduce known results for Si."""
    Tr = np.arange(100., 300., 5.)
    mur = np.linspace(0., 0.4, num=100)
    data_noder = BoltzTraP2.dft.DFTData(os.path.join(datadir, "Si"))
    fermiintegrals = np.load(os.path.join(mydir, "Si_fermiintegrals.npz"))
    N = fermiintegrals["N"]
    L0 = fermiintegrals["L0"]
    L1 = fermiintegrals["L1"]
    L2 = fermiintegrals["L2"]
    L11 = fermiintegrals["L11"]
    coeffs = BoltzTraP2.bandlib.calc_Onsager_coefficients(
        L0, L1, L2, mur, Tr, data_noder.get_volume(), Lm11=L11)
    ref = np.load(os.path.join(mydir, "Si_Onsager.npz"))
    assert np.allclose(ref["cond"], coeffs[0])
    assert np.allclose(ref["seebeck"], coeffs[1])
    assert np.allclose(ref["kappa"], coeffs[2], atol=50.)
    assert np.allclose(ref["Hall"], coeffs[3])
