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

# Test the public interface of the "sphere" C++ module

import os
import os.path
import xml.etree.ElementTree as etree

import numpy as np
import pytest
import ase

import BoltzTraP2.sphere

mydir = os.path.abspath(os.path.dirname(__file__))
datadir = os.path.join(os.path.dirname(mydir), "data")


def test_compute_bounds_exceptions():
    """compute_bounds should identify invalid arguments."""
    # Lattice vectors with wrong shape
    lattvec = np.empty((4, 3))
    with pytest.raises(ValueError):
        BoltzTraP2.sphere.compute_bounds(lattvec, 1.)
    # Linearly dependent lattice vectors
    lattvec = np.ones((3, 3))
    with pytest.raises(ValueError):
        BoltzTraP2.sphere.compute_bounds(lattvec, 1.)
    # Nonpositive radius
    lattvec = np.eye(3)
    with pytest.raises(ValueError):
        BoltzTraP2.sphere.compute_bounds(lattvec, 0.)


def test_compute_bounds_cube():
    """compute_bounds should identify the cube around a sphere."""
    # Create a simple cubic lattice
    lattvec = np.eye(3)
    for i in range(1, 10):
        assert BoltzTraP2.sphere.compute_bounds(
            lattvec, float(i)).tolist() == [i, i, i]


def test_compute_bounds_Li_and_LiZnSb():
    """compute_bounds should work for realistic lattice vectors."""
    # Test based on Li
    lattvec = np.array([[-3.19204785, 3.19204785,
                         3.19204785], [3.19204785, -3.19204785, 3.19204785],
                        [3.19204785, 3.19204785, -3.19204785]])
    assert BoltzTraP2.sphere.compute_bounds(
        lattvec, 145.47229311).tolist() == [33, 33, 33]
    # Test Based on LZS
    lattvec = np.array([[4.431, -2.2155, 0.], [0., 3.83735835, 0.],
                        [0., 0., 7.157]])
    assert BoltzTraP2.sphere.compute_bounds(
        lattvec, 238.980887299).tolist() == [63, 63, 34]


@pytest.fixture()
def LZS_atoms():
    """Build an ASE atoms object representing LZS."""
    lattvec = np.array([[4.431, -2.2155, 0.], [0., 3.83735835, 0.],
                        [0., 0., 7.157]])
    return ase.Atoms(
        ["Sb", "Sb", "Zn", "Zn", "Li", "Li"],
        cell=lattvec.T,
        scaled_positions=np.array([[0.66666667, 0.33333333, 0.885], [
            0.33333333, 0.66666667, 0.385
        ], [0.66666667, 0.33333333, 0.5], [0.33333333, 0.66666667, 0.],
                                   [0., 0., 0.696], [0., 0., 0.196]]),
        pbc=True)


def test_compute_radius(LZS_atoms):
    """compute_radius should return the right radii."""
    # Test based on Si
    lattvec = 5.45052526 * .5 * (np.ones((3, 3)) - np.eye(3))
    atoms = ase.Atoms("Si", cell=lattvec.T, pbc=True)
    assert np.allclose(
        BoltzTraP2.sphere.compute_radius(atoms, None, 4500), 127.802713321)
    # Test based on LZS
    assert np.allclose(
        BoltzTraP2.sphere.compute_radius(LZS_atoms, None, 19575),
        238.980887299)


def test_get_equivalences(LZS_atoms):
    """"get_equivalences should match know results from BoltzTraP."""
    # Detailed test for Si
    reference = np.load(os.path.join(mydir, "Si_equivalences.npz"))
    reference = [
        frozenset(tuple(i) for i in reference[r].tolist()) for r in reference
    ]
    lattvec = 5.45052526 * .5 * (np.ones((3, 3)) - np.eye(3))
    atoms = ase.Atoms(
        ["Si", "Si"],
        cell=lattvec.T,
        scaled_positions=np.array([[0.125, 0.125, 0.125],
                                   [0.875, 0.875, 0.875]]),
        pbc=True)
    equivalences = [
        frozenset(tuple(i) for i in e.tolist())
        for e in BoltzTraP2.sphere.get_equivalences(atoms, None, 4500)
    ]
    # The ordering of a particular equivalence class is irrelevant, as is
    # the relative order of equivalence classes in the set.
    assert len(reference) == len(equivalences)
    for e in equivalences:
        assert e in reference
    # Test the number of equivalence classes for LZS
    # TODO: Revise whether this is the right magnetic configuration
    equivalences = BoltzTraP2.sphere.get_equivalences(LZS_atoms, None, 19575)
    assert len(equivalences) == 21164


def test_calc_reciprocal_degeneracies():
    """calc_reciprocal_degeneracies should return the same weights as
     determined by VASP reproduce the weights provided by VASP."""
    lattvec = 5.467112115767304 * .5 * (np.ones((3, 3)) - np.eye(3))
    atoms = ase.Atoms(
        ["Si", "Si"],
        cell=lattvec.T,
        scaled_positions=np.array([[0., 0., 0.], [0.25, 0.25, 0.25]]),
        pbc=True)
    kpoints = np.load(os.path.join(mydir, "kpoints.npz"))["kpoints"]
    degeneracies = BoltzTraP2.sphere.calc_reciprocal_degeneracies(
        atoms, None, kpoints)
    degeneracies = degeneracies.astype(np.float64)
    ref = np.array([
        0.00020354, 0.00162833, 0.00162833, 0.00162833, 0.00162833, 0.00162833,
        0.00162833, 0.00162833, 0.00162833, 0.00122125, 0.00488500, 0.00488500,
        0.00488500, 0.00488500, 0.00488500, 0.00488500, 0.00488500, 0.00488500,
        0.00488500, 0.00488500, 0.00488500, 0.00488500, 0.00488500, 0.00488500,
        0.00244250, 0.00122125, 0.00488500, 0.00488500, 0.00488500, 0.00488500,
        0.00488500, 0.00488500, 0.00488500, 0.00488500, 0.00488500, 0.00488500,
        0.00488500, 0.00488500, 0.00244250, 0.00122125, 0.00488500, 0.00488500,
        0.00488500, 0.00488500, 0.00488500, 0.00488500, 0.00488500, 0.00488500,
        0.00488500, 0.00488500, 0.00244250, 0.00122125, 0.00488500, 0.00488500,
        0.00488500, 0.00488500, 0.00488500, 0.00488500, 0.00488500, 0.00488500,
        0.00244250, 0.00122125, 0.00488500, 0.00488500, 0.00488500, 0.00488500,
        0.00488500, 0.00488500, 0.00244250, 0.00122125, 0.00488500, 0.00488500,
        0.00488500, 0.00488500, 0.00244250, 0.00122125, 0.00488500, 0.00488500,
        0.00244250, 0.00122125, 0.00244250, 0.00488500, 0.00977000, 0.00977000,
        0.00977000, 0.00977000, 0.00977000, 0.00977000, 0.00488500, 0.00977000,
        0.00977000, 0.00977000, 0.00977000, 0.00977000, 0.00977000, 0.00977000,
        0.00977000, 0.00977000, 0.00977000, 0.00488500, 0.00488500, 0.00977000,
        0.00977000, 0.00977000, 0.00977000, 0.00977000, 0.00977000, 0.00977000,
        0.00977000, 0.00488500, 0.00488500, 0.00977000, 0.00977000, 0.00977000,
        0.00977000, 0.00977000, 0.00977000, 0.00488500, 0.00488500, 0.00977000,
        0.00977000, 0.00977000, 0.00977000, 0.00488500, 0.00488500, 0.00977000,
        0.00977000, 0.00488500, 0.00488500, 0.00488500, 0.00488500, 0.00977000,
        0.00977000, 0.00977000, 0.00977000, 0.00488500, 0.00977000, 0.00977000,
        0.00977000, 0.00977000, 0.00977000, 0.00977000, 0.00488500, 0.00488500,
        0.00977000, 0.00977000, 0.00977000, 0.00977000, 0.00488500, 0.00488500,
        0.00977000, 0.00977000, 0.00488500, 0.00488500, 0.00488500, 0.00488500,
        0.00977000, 0.00977000, 0.00488500, 0.00977000, 0.00977000, 0.00488500,
        0.00488500, 0.00488500, 0.00488500
    ])
    ref = np.round(ref * degeneracies.sum()).astype(np.int32)
    assert np.array_equal(ref, degeneracies)


def test_calc_reciprocal_stars():
    """calc_reciprocal_stars should expand a set of q points into the correct
    set of stars.
    """
    lattvec = 5.45052526 * .5 * (np.ones((3, 3)) - np.eye(3))
    atoms = ase.Atoms(
        ["Si", "Si"],
        cell=lattvec.T,
        scaled_positions=np.array([[0.125, 0.125, 0.125],
                                   [0.875, 0.875, 0.875]]),
        pbc=True)
    kpoints = np.array([[0., 0., 0.], [.1, .1, .1], [.1, 0., 0.]])
    stars = BoltzTraP2.sphere.calc_reciprocal_stars(atoms, None, kpoints)
    ref = np.zeros((1, 3))
    assert stars[0].shape == ref.shape
    assert np.allclose(stars[0], ref)
    ref = np.array([[-0.1, 0., 0.], [-0.1, -0.1, -0.1], [0., 0., 0.1],
                    [0., 0.1, 0.], [0., -0.1, 0.], [0., 0., -0.1],
                    [0.1, 0.1, 0.1], [0.1, 0., 0.]])
    assert stars[1].shape == ref.shape
    assert np.allclose(stars[1], ref)
    ref = np.array([[-.1, -.1, -.1], [0., 0., .1], [0., .1, 0.], [-.1, 0., 0.],
                    [0., -.1, 0.], [0., 0., -.1], [.1, .1, .1], [.1, 0., 0.]])
    assert stars[2].shape == ref.shape
    assert np.allclose(stars[2], ref)


def test_calc_reciprocal_iksubset():
    """calc_reciprocal_iksubset should be able to reverse the effect of
    calc_reciprocal_stars().
    """
    lattvec = 5.45052526 * .5 * (np.ones((3, 3)) - np.eye(3))
    atoms = ase.Atoms(
        ["Si", "Si"],
        cell=lattvec.T,
        scaled_positions=np.array([[0.125, 0.125, 0.125],
                                   [0.875, 0.875, 0.875]]),
        pbc=True)
    kpoints = np.array([[0., 0., 0.], [.1, .1, .1], [.1, 0., 0.]])
    stars = BoltzTraP2.sphere.calc_reciprocal_stars(atoms, None, kpoints)
    kpoints = []
    for s in stars:
        for k in s:
            kpoints.append(k)
    kpoints = np.array(kpoints)
    indices = BoltzTraP2.sphere.calc_reciprocal_iksubset(atoms, None, kpoints)
    assert len(indices) == 2