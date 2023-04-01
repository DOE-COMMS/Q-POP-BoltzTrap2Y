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

# Test the interpolation and reconstruction code

import os
import os.path
from collections import namedtuple

import ase
import numpy as np
import pytest

import BoltzTraP2
import BoltzTraP2.io
import BoltzTraP2.dft
import BoltzTraP2.fite

mydir = os.path.abspath(os.path.dirname(__file__))
datadir = os.path.join(os.path.dirname(mydir), "data")


@pytest.fixture()
def si_data():
    """Load DFT data for Si and compute the k-point equivalence classes."""
    data_noder = BoltzTraP2.dft.DFTData(os.path.join(datadir, "Si"))
    data_withder = BoltzTraP2.dft.DFTData(
        os.path.join(datadir, "Si"), derivatives=True)
    equivalences = BoltzTraP2.sphere.get_equivalences(
        data_noder.atoms, data_noder.magmom, 2 * len(data_noder.kpoints))
    DataAndEquivalences = namedtuple(
        "DataAndEquivalences", ["data_noder", "data_withder", "equivalences"])
    return DataAndEquivalences(data_noder, data_withder, equivalences)


@pytest.fixture()
def si_coefficients():
    """Load interpolation results for Si."""
    return np.load(os.path.join(mydir, "Si_fitde3D.npz"))


def test_fitde3D_saved_noder(si_data, si_coefficients):
    """fit3D should reproduce known results for Si without derivatives."""
    coeffs_noder = BoltzTraP2.fite.fitde3D(si_data.data_noder,
                                           si_data.equivalences)
    assert np.allclose(si_coefficients["noder"], coeffs_noder)


@pytest.mark.skip(reason="interpolation coefficients can be machine-dependent"
                  " if derivatives are included in the fit")
def test_fitde3D_saved_withder(si_data, si_coefficients):
    """fit3D should reproduce known results for Si with derivatives."""
    coeffs_withder = BoltzTraP2.fite.fitde3D(si_data.data_withder,
                                             si_data.equivalences)
    assert np.allclose(si_coefficients["withder"], coeffs_withder)


def test_getBands_reconstruct(si_data):
    """getBands should be able to reconstruct the original bands."""
    import scipy as sp
    import scipy.linalg as la
    # Low "enhancement factors" can lead to poor-quality interpolation,
    # but a factor of 10 should guarantee that the original bands are
    # reproduced pretty much exactly.
    equivalences = BoltzTraP2.sphere.get_equivalences(
        si_data.data_noder.atoms, si_data.data_noder.magmom,
        10 * len(si_data.data_noder.kpoints))
    coeffs_noder = BoltzTraP2.fite.fitde3D(si_data.data_noder, equivalences)
    ebands = BoltzTraP2.fite.getBands(si_data.data_noder.kpoints, equivalences,
                                      si_data.data_noder.get_lattvec(),
                                      coeffs_noder)[0]
    assert np.allclose(si_data.data_noder.ebands, ebands, atol=1e-8)
    coeffs_withder = BoltzTraP2.fite.fitde3D(si_data.data_withder,
                                             equivalences)
    ebands = BoltzTraP2.fite.getBands(
        si_data.data_withder.kpoints, equivalences,
        si_data.data_withder.get_lattvec(), coeffs_withder)[0]
    assert np.allclose(si_data.data_withder.ebands, ebands, atol=1e-8)


def test_getBTPbands_parallel(si_data, si_coefficients):
    """getBTPBands should work the same on any number of cores."""
    nworkers = (1, 2, 3)
    # Test without derivatives
    results = []
    for i in nworkers:
        results.append(
            BoltzTraP2.fite.getBTPbands(
                si_data.equivalences,
                si_coefficients["withder"],
                si_data.data_noder.get_lattvec(),
                curvature=True,
                nworkers=i))
    for r in results[1:]:
        for e, e0 in zip(r, results[0]):
            assert np.allclose(e, e0)


def test_getBTPbands_None(si_data, si_coefficients):
    """getBTPbands should return None as the curvature by default."""
    r = BoltzTraP2.fite.getBTPbands(
        si_data.equivalences,
        si_coefficients["withder"],
        si_data.data_noder.get_lattvec(),
        curvature=False,
        nworkers=1)
    assert r[2] is None
