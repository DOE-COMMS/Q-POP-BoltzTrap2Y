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

# Test the input and output capabilities of BoltzTraP2

import os
import os.path
import xml.etree.ElementTree as etree

import ase
import ase.io
import ase.io.wien2k
import numpy as np
import pytest
import netCDF4 as nc

import BoltzTraP2
import BoltzTraP2.io
import BoltzTraP2.units

mydir = os.path.abspath(os.path.dirname(__file__))
datadir = os.path.join(os.path.dirname(mydir), "data")


def test_parse_empty_vasp_array():
    """_parse_vasp_array should parse the empty array."""
    assert BoltzTraP2.io._parse_vasp_array("") == []


def test_parse_vasp_array():
    """_parse_vasp_array should parse a non-empty array."""
    a = "1 1.2 1.3e-4 0.9e-5"
    ref = np.array([1., 1.2, 1.3e-4, 0.9e-5])
    assert np.allclose(BoltzTraP2.io._parse_vasp_array(a), ref)


@pytest.fixture()
def si_vasprunxml():
    """Create an xml tree from a vasprun.xml for Si."""
    filename = os.path.join(datadir, "Si.vasp", "vasprun.xml")
    return etree.parse(filename)


@pytest.fixture()
def si_interpolated_vasprunxml():
    """Create an xml tree from a vasprun.xml for Si including interpolation."""
    filename = os.path.join(datadir, "Si.vasp.interp", "vasprun.xml")
    return etree.parse(filename)


def test_parse_vasp_name(si_vasprunxml):
    """_parse_vasp_name should return "unknown system" for Si."""
    assert BoltzTraP2.io._parse_vasp_name(si_vasprunxml) == "unknown system"


def test_parse_vasp_magmom():
    """_parse_vasp_magmom should work for unpolarized, collinear and
    noncollinear spins.
    """
    # Si, unpolarized
    xml = etree.parse(os.path.join(datadir, "Si.vasp", "vasprun.xml"))
    assert BoltzTraP2.io._parse_vasp_magmom(xml) is None
    # PbTe, unpolarized
    xml = etree.parse(
        os.path.join(datadir, "PbTe.vasp.unpolarized", "vasprun.xml"))
    assert BoltzTraP2.io._parse_vasp_magmom(xml) is None
    # PbTe, nocollinear
    xml = etree.parse(
        os.path.join(datadir, "PbTe.vasp.sl", "vasprun.xml"))
    magmom = BoltzTraP2.io._parse_vasp_magmom(xml)
    assert magmom.shape == (2, 3)
    assert np.allclose(magmom, np.ones_like(magmom))
    # LiZnSb, collinear
    xml = etree.parse(os.path.join(datadir, "LiZnSb.vasp", "vasprun.xml"))
    magmom = BoltzTraP2.io._parse_vasp_magmom(xml)
    assert magmom.shape == (6, )
    assert np.allclose(magmom, np.ones_like(magmom))


def test_parse_vasp_fermi(si_vasprunxml):
    """_parse_vasp_fermi should return the Fermi level of Si."""
    assert np.allclose(
        BoltzTraP2.io._parse_vasp_fermi(si_vasprunxml), 5.71756438)


def test_parse_vasp_kinter(si_vasprunxml, si_interpolated_vasprunxml):
    """_parse_vasp_kinter should discriminate between interpolated and
    non-interpolated VASP calculations.
    """
    assert BoltzTraP2.io._parse_vasp_kinter(si_vasprunxml) == 0
    assert BoltzTraP2.io._parse_vasp_kinter(si_interpolated_vasprunxml) == 3


def test_parse_vasp_lvel():
    """_parse_vasp_lvel should discriminate between VASP calculations with
    and without group velocities.
    """
    xml = etree.parse(os.path.join(datadir, "Si.vasp", "vasprun.xml"))
    assert BoltzTraP2.io._parse_vasp_lvel(xml) == True
    xml = etree.parse(os.path.join(datadir, "Si.vasp.noder", "vasprun.xml"))
    assert BoltzTraP2.io._parse_vasp_lvel(xml) == False


def test_detect_vasp_broken_interpolation():
    """_detect_vasp_broken_interpolation should return True only for the 
    Si.vasp.noder.interp.old directory.
    """
    directories = ("Si.vasp", "Si.vasp.noder", "Si.vasp.interp",
                   "Si.vasp.noder.interp", "Si.vasp.interp.old")
    for d in directories:
        xml = etree.parse(os.path.join(datadir, d, "vasprun.xml"))
        assert BoltzTraP2.io._detect_vasp_broken_interpolation(xml) == False
    xml = etree.parse(
        os.path.join(datadir, "Si.vasp.noder.interp.old", "vasprun.xml"))
    assert BoltzTraP2.io._detect_vasp_broken_interpolation(xml) == True


def test_detect_vasp_interpolated_velocities():
    """_detect_vasp_interpolated_velocities should return True only for the
    'interp' data directories.
    """
    directories = ("Si.vasp", "Si.vasp.noder", "Si.vasp.noder.interp",
                   "Si.vasp.noder.interp.old")
    for d in directories:
        xml = etree.parse(os.path.join(datadir, d, "vasprun.xml"))
        assert BoltzTraP2.io._detect_vasp_interpolated_velocities(xml) == False
    directories = ("Si.vasp.interp", "Si.vasp.interp.old")
    for d in directories:
        xml = etree.parse(os.path.join(datadir, d, "vasprun.xml"))
        assert BoltzTraP2.io._detect_vasp_interpolated_velocities(xml) == True


def test_get_vasp_kpoints_path():
    """_get_vasp_kpoints_path should be able to detect all the different
    datasets containing lists of k points.
    """
    reference = {
        "Si.vasp":
        "./kpoints/varray[@name=\"kpointlist\"]",
        "Si.vasp.noder":
        "./kpoints/varray[@name=\"kpointlist\"]",
        "Si.vasp.interp":
        "./calculation/eigenvelocities[@comment=\"interpolated\"]/"
        "kpoints/varray[@name=\"kpointlist\"]",
        "Si.vasp.noder.interp":
        "./calculation/eigenvalues[@comment=\"interpolated\"]/"
        "kpoints/varray[@name=\"kpointlist\"]",
        "Si.vasp.interp.old":
        "./calculation/eigenvalues[@comment=\"interpolated_ibz\"]"
        "/electronvelocities/kpoints/varray[@name=\"kpointlist\"]"
    }
    for d in reference:
        xml = etree.parse(os.path.join(datadir, d, "vasprun.xml"))
        assert BoltzTraP2.io._get_vasp_kpoints_path(xml) == reference[d]


def test_get_vasp_energies_path():
    """_get_vasp_energies_path should be able to detect all the different
    datasets containing lists of energies.
    """
    reference = {
        "Si.vasp":
        "./calculation/eigenvalues/array/set",
        "Si.vasp.noder":
        "./calculation/eigenvalues/array/set",
        "Si.vasp.interp":
        "./calculation/eigenvelocities[@comment=\"interpolated\"]/"
        "eigenvalues/array/set",
        "Si.vasp.noder.interp":
        "./calculation/eigenvalues[@comment=\"interpolated\"]/"
        "eigenvalues/array/set",
        "Si.vasp.interp.old":
        "./calculation/eigenvalues[@comment=\"interpolated_ibz\"]"
        "/electronvelocities/eigenvalues/array/set"
    }
    for d in reference:
        xml = etree.parse(os.path.join(datadir, d, "vasprun.xml"))
        assert BoltzTraP2.io._get_vasp_energies_path(xml) == reference[d]


def test_get_vasp_velocities_path():
    """_get_vasp_energies_path should be able to detect all the different
    datasets containing lists of velocities, or None if the calculation was
    performed with LVEL = F.
    """
    reference = {
        "Si.vasp":
        "./calculation/electronvelocities",
        "Si.vasp.interp":
        "./calculation/eigenvelocities",
        "Si.vasp.interp.old":
        "./calculation/eigenvalues[@comment=\"interpolated_ibz\"]"
        "/electronvelocities"
    }
    for d in reference:
        xml = etree.parse(os.path.join(datadir, d, "vasprun.xml"))
        assert BoltzTraP2.io._get_vasp_velocities_path(xml) == reference[d]
    for d in ("Si.vasp.noder", "Si.vasp.noder.interp"):
        xml = etree.parse(os.path.join(datadir, d, "vasprun.xml"))
        assert BoltzTraP2.io._get_vasp_velocities_path(xml) is None


def test_parse_vasp_structure(si_vasprunxml):
    """_parse_vasp_structure should be able to parse a Si structure."""
    atoms = BoltzTraP2.io._parse_vasp_structure(si_vasprunxml)
    # Test the lattice vectors
    ref = 5.467112115767304 * .5 * (np.ones((3, 3)) - np.eye(3))
    cell = atoms.get_cell()
    assert np.allclose(cell, ref)
    # Test the atomic positions
    ref = np.array([[0., 0., 0.], [.25, .25, .25]])
    positions = atoms.get_scaled_positions()
    assert np.allclose(positions, ref)


def test_parse_vasp_ikpoints(si_vasprunxml):
    """_parse_vasp_ikpoints should return the right k points for Si."""
    kpoints = BoltzTraP2.io._parse_vasp_ikpoints(si_vasprunxml)
    assert kpoints.shape[0] == 165
    ref = np.load(os.path.join(mydir, "kpoints.npz"))["kpoints"]
    assert np.allclose(kpoints, ref)


def test_parse_vasp_eigenvalues(si_vasprunxml):
    """_parse_vasp_eigenvalues should return the right eigenvalues for Si."""
    eigenvalues = BoltzTraP2.io._parse_vasp_eigenvalues(si_vasprunxml)
    assert eigenvalues.shape == (1, 165, 8)
    # Only a few eigenvalues are actually checking
    assert np.allclose(
        eigenvalues[0, 0, :],
        [-6.1962, 5.6258, 5.6258, 5.6258, 8.1852, 8.1852, 8.1852, 8.7682])
    assert np.allclose(
        eigenvalues[0, 20, :],
        [-4.9483, 0.6728, 4.1974, 4.8106, 7.6709, 9.0903, 9.2898, 12.7089])
    assert np.allclose(
        eigenvalues[0, 163, :],
        [-2.1917, -1.7788, 1.6992, 2.1858, 8.7442, 9.6465, 11.3190, 11.4235])


def test_parse_vasp_velocities(si_vasprunxml):
    """_parse_vasp_velocities should return the right velocities for Si."""
    kpoints, velocities = BoltzTraP2.io._parse_vasp_velocities(si_vasprunxml)
    assert kpoints.shape == (4913, 3)
    assert velocities.shape == (1, 4913, 8, 3)
    assert np.allclose(velocities[0, 0, :, :], np.zeros((8, 3)))
    assert np.allclose(velocities[0, 291, :, :],
                       [[-0.4357, 1.3151, 0.4357], [4.3318, -5.7890, -4.3317],
                        [0.0043, -3.6174, -0.0041], [-0.8963, -3.5206, 0.8961],
                        [-1.9016, -2.4715, 1.9016], [2.4795, 4.1444, -2.4795],
                        [-0.0996, 5.1275, 0.0996], [-3.7841, 5.4531, 3.7841]])
    assert np.allclose(
        velocities[0, 1075, :, :],
        [[-1.1381, 3.7815, -0.3560], [2.3511, -2.6098, 2.5102],
         [0.0351, -5.1461, -2.1357], [-0.2462, -2.2663, -2.5013],
         [-5.2723, -2.0726, 2.0552], [8.8141, 5.6665, -0.5984],
         [-9.1441, -2.1866, -8.3536], [0.1183, 6.8992, 8.3927]])


def test_parse_vasprunxml_notfound():
    """parse_vasprunxml should raise a FileNotFoundError if it cannot open the
    file.
    """
    with pytest.raises(FileNotFoundError):
        BoltzTraP2.io.parse_vasprunxml(os.path.join(datadir, "not_there.xml"))


def test_parse_vasprunxml_broken():
    """parse_vasprunxml should raise a ValueError for the old interpolated
    results without derivatives, which are known to be incorrect.
    """
    with pytest.raises(ValueError):
        BoltzTraP2.io.parse_vasprunxml(
            os.path.join(datadir, "Si.vasp.noder.interp.old", "vasprun.xml"))


def test_parse_vasprunxml():
    """parse_vasprunxml should be able to load the vasprun.xml for Si."""
    filename = os.path.join(datadir, "Si.vasp", "vasprun.xml")
    results = BoltzTraP2.io.parse_vasprunxml(filename)
    keys = tuple(sorted(list(results.keys())))
    # Check that all relevant information is in the results and that there are
    # no unknown pieces.
    assert keys == ("E", "atoms", "fermi", "kpoints", "magmom", "name",
                    "nelect", "v")
    assert results["E"].shape == (1, 165, 8)
    assert results["v"].shape == (1, 165, 8, 3)


@pytest.fixture()
def gsrfile():
    """Load the netCDF file for Si."""
    filename = os.path.join(datadir, "Si.abinit", "outsi_DS1_GSR.nc")
    ncf = nc.Dataset(filename, mode="r")
    return ncf


def test_parse_abinit_structure(gsrfile):
    """_parse_abinit_structure should be able to parse a Si structure."""
    atoms = BoltzTraP2.io._parse_abinit_structure(gsrfile)
    # Test the lattice vectors
    ref = 5.16731481286141 * (np.ones(
        (3, 3)) - np.eye(3)) / BoltzTraP2.units.Angstrom
    cell = atoms.get_cell()
    assert np.allclose(cell, ref)
    # Test the atomic positions
    ref = np.array([[0., 0., 0.], [.25, .25, .25]])
    positions = atoms.get_scaled_positions()
    assert np.allclose(positions, ref)


def test_parse_abinitrun():
    """parse_abinitrun should be able to load the GSR.nc file for Si."""
    filename = os.path.join(datadir, "Si.abinit", "outsi_DS1_GSR.nc")
    results = BoltzTraP2.io.parse_abinitrun(filename)
    keys = tuple(sorted(list(results.keys())))
    # Check that all relevant information is in the results and that there are
    # no unknown pieces.
    assert keys == ("E", "atoms", "fermi", "kpoints", "name", "nelect")
    assert results["E"].shape == (1, 408, 10)


def test_W2Kmommat():
    """W2Kmommat should parse .mommat2 files correctly.
    """
    filename = os.path.join(datadir, "Si", "Si.struct")
    dum = ase.io.wien2k.read_struct(filename, ase=False)
    latt = dum[1]
    if latt == "R":
        latt = "P"
    conv = ase.io.wien2k.c2p(latt)
    filename = os.path.join(datadir, "Si", "Si.energy")
    kpoints = BoltzTraP2.io.W2Kene(filename, conv)[0]
    filename = os.path.join(datadir, "Si", "Si.mommat2")
    mommat, nemin, nemax = BoltzTraP2.io.W2Kmommat(filename, kpoints)
    assert nemin == 1
    assert nemax == 6
    ref = np.load(os.path.join(mydir, "Si_mommat_ref.npz"))["mommat"]
    assert np.allclose(mommat, ref)


def test_W2Kfermi():
    """W2Kfermi should return the right value of the Fermi level."""
    filename = os.path.join(datadir, "Si", "Si.scf")
    assert BoltzTraP2.io.W2Kfermi(filename) == .5 * 0.3766817831


def test_read_GENE_struct():
    """read_GENE_struct() should build the right Atoms object."""
    # Li: single-atom system
    filename = os.path.join(datadir, "Li", "Li_BLZTRP.structure")
    atoms = BoltzTraP2.io.read_GENE_struct(filename)
    ref = 0.3192047850E+01 * (np.ones(
        (3, 3)) - 2. * np.eye(3)) / BoltzTraP2.units.Angstrom
    cell = atoms.get_cell()
    assert np.allclose(ref, cell)
    assert atoms.get_chemical_symbols() == ["Li"]
    ref = np.zeros(3)
    assert np.allclose(ref, atoms.get_positions())
    # LZS: multi-atom system
    filename = os.path.join(datadir, "LiZnSb.GENE", "LiZnSb.structure")
    atoms = BoltzTraP2.io.read_GENE_struct(filename)
    ref = np.array([[
        4.4309997546822144, 0.0000000000000000, 0.0000000000000000
    ], [-2.2154998773411063, 3.8373583517174135, 0.0000000000000000], [
        0.0000000000000004, 0.0000000000000008, 7.1570000621175227
    ]]) / BoltzTraP2.units.Angstrom
    cell = atoms.get_cell()
    assert np.allclose(ref, cell)
    assert atoms.get_chemical_symbols() == ["Sb", "Sb", "Zn", "Zn", "Li", "Li"]
    ref = np.array([[
        2.2154998994961068, 1.2791194377812773, 6.3339450549740075
    ], [-0.0000000221549981, 2.5582389139361372, 2.7554450239152461], [
        2.2154998994961064, 1.2791194377812771, 3.5785000310587614
    ], [-0.0000000221549983, 2.5582389139361368, 0.0000000000000000], [
        0.0000000000000003, 0.0000000000000005, 4.9812720432337958
    ], [0.0000000000000001, 0.0000000000000001, 1.4027720121750344]
                    ]) / BoltzTraP2.units.Angstrom
    assert np.allclose(ref, atoms.get_positions())


def test_read_GENE_eneandmat_old():
    """read_GENE_eneandmat() should parse an old-style file correctly."""
    # This old-style file contains only one spin channel, and no momentum
    # matrix elements.
    filename = os.path.join(datadir, "Li", "Li_BLZTRP.energy")
    r = BoltzTraP2.io.read_GENE_eneandmat(filename)
    assert np.allclose(r[0], 0.6940745476E-01 / 2.)
    assert r[1] == 2.0
    assert r[2].shape == (413, 3)
    assert np.allclose(r[2][0], [0., 0., 0.])
    assert np.allclose(r[2][20],
                       [0.3333333333E+00, 0.4166666667E-01, 0.0000000000E+00])
    assert np.allclose(r[2][80],
                       [0.3750000000E+00, 0.8333333333E-01, 0.4166666667E-01])
    assert np.allclose(r[2][391],
                       [-0.4166666667E+00, 0.4166666667E+00, 0.3333333333E+00])
    assert r[3].shape == (3, 413)
    assert np.allclose(
        r[3][:, 20],
        np.array([-0.1963709022E-01, 0.5029658971E+00, 0.8701510047E+00]) / 2.)
    assert r[4] is None


def test_read_GENE_eneandmat_new():
    """read_GENE_eneandmat() should parse a new-style file correctly."""
    # This new-style file contains two spin channels as well as electronic
    # group velocities.
    filename = os.path.join(datadir, "Li.GENE.fromvasp", "Li.energy")
    r = BoltzTraP2.io.read_GENE_eneandmat(filename)
    assert np.allclose(r[0], 0.024951731946)
    assert r[1] == 1.0
    assert r[2].shape == (286, 3)
    assert np.allclose(r[2][0], [0., 0., 0.])
    assert np.allclose(r[2][205], [-0.38095238, 0.47619048, 0.14285714])
    assert r[3].shape == (12, 286)
    assert np.allclose(r[3][:, 205], [
        -1.68363713, 0.04813426, 0.18563185, 0.30062416, 0.47661299,
        0.55504707, -1.68383558, 0.047572, 0.18507326, 0.30008762, 0.47603602,
        0.55446643
    ])
    assert r[4].shape == (286, 12, 3)
    assert np.allclose(
        r[4][205, :, :],
        np.array([[3.32542058e-04, 8.16769966e-05, -1.55575232e-05], [
            1.16308043e-01, -2.65839177e-02, 2.03122912e-02
        ], [-8.43431672e-02, 1.19983508e-01, -1.81672977e-02], [
            -8.65017735e-02, -4.95137622e-02, -8.76724771e-02
        ], [-1.47465873e-02, 1.28310672e-02, 8.57472336e-02], [
            -4.77188129e-02, -1.68221553e-01, 1.73855321e-03
        ], [3.34486748e-04, 8.16769966e-05, -1.55575232e-05], [
            1.16298320e-01, -2.65780836e-02, 2.03103465e-02
        ], [-8.43470565e-02, 1.19973785e-01, -1.81828552e-02], [
            -8.65153863e-02, -4.95254303e-02, -8.76452515e-02
        ], [-1.47485320e-02, 1.28349566e-02, 8.57783487e-02],
                  [-4.76663063e-02, -1.68244890e-01, 1.70160410e-03]]))
