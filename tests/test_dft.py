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

# Test the loading of DFT data and related capabilities

import os
import copy
import tempfile
import os.path
import xml.etree.ElementTree as etree

import numpy as np
import scipy as sp
import scipy.linalg as la
import pytest

import BoltzTraP2
import BoltzTraP2.dft
import BoltzTraP2.units

mydir = os.path.abspath(os.path.dirname(__file__))
datadir = os.path.join(os.path.dirname(mydir), "data")


def test_register_loader():
    """register_loader should append a new loader to the list."""
    n = len(BoltzTraP2.dft.loaders)

    def trivial_function(arg):
        """Do absolutely nothing."""
        pass

    BoltzTraP2.dft.register_loader("TEST", trivial_function)
    assert len(BoltzTraP2.dft.loaders) == n + 1
    assert BoltzTraP2.dft.loaders.pop() == ("TEST", trivial_function)


# Note that the correctness of the individual elements created by the loaders
# is not tested in detail. It is assumed that the io module takes care of that.


def test_VASPLoader_exception():
    """VASPLoader should raise exceptions for non-VASP directories."""
    # Test with a directory that does not exist
    # Note that mktemp is vulnerable to race conditions and should not be used
    # to create actual temporary directories, but is perfect to create a
    # random directory name that does not exist.
    arg = tempfile.mktemp(dir=datadir)
    with pytest.raises(FileNotFoundError):
        BoltzTraP2.dft.VASPLoader(arg)
    # Test with a Wien2k directory
    arg = os.path.join(datadir, "Si")
    with pytest.raises(BoltzTraP2.dft.LoaderError):
        BoltzTraP2.dft.VASPLoader(arg)


def test_VASPLoader():
    """VASPLoader should work correctly for some examples in datadir."""
    # Si, non-spin-polarized, with derivatives
    loader = BoltzTraP2.dft.VASPLoader(os.path.join(datadir, "Si.vasp"))
    lattvec = 5.467112115767304 * 0.5 * (np.ones((3, 3)) - np.eye(3))
    assert np.allclose(lattvec, loader.atoms.get_cell())
    assert loader.atoms.get_chemical_formula() == "Si2"
    assert loader.dosweight == 2.0
    assert loader.kpoints.shape == (165, 3)
    assert np.allclose(loader.fermi, 0.21011661716080218)
    assert loader.ebands.shape == (8, 165)
    assert np.allclose(
        loader.ebands[:, 0],
        np.array([
            -6.1962, 5.6258, 5.6258, 5.6258, 8.1852, 8.1852, 8.1852, 8.7682
        ]) * BoltzTraP2.units.eV,
    )
    assert loader.nelect == 8.0
    assert loader.mommat.shape == (165, 8, 3)
    assert loader.magmom is None
    # Si, non-spin-polarized, without derivatives
    loader = BoltzTraP2.dft.VASPLoader(os.path.join(datadir, "Si.vasp.noder"))
    lattvec = 5.467112115767304 * 0.5 * (np.ones((3, 3)) - np.eye(3))
    assert np.allclose(lattvec, loader.atoms.get_cell())
    assert loader.atoms.get_chemical_formula() == "Si2"
    assert loader.dosweight == 2.0
    assert loader.kpoints.shape == (165, 3)
    assert np.allclose(loader.fermi, 0.21011661716080218)
    assert loader.ebands.shape == (8, 165)
    assert np.allclose(
        loader.ebands[:, 0],
        np.array([
            -6.1962, 5.6258, 5.6258, 5.6258, 8.1852, 8.1852, 8.1852, 8.7682
        ]) * BoltzTraP2.units.eV,
    )
    assert loader.nelect == 8.0
    with pytest.raises(AttributeError):
        loader.mommat
    assert loader.magmom is None
    # LiZnSn, spin-polarized, with derivatives
    loader = BoltzTraP2.dft.VASPLoader(os.path.join(datadir, "LiZnSb.vasp"))
    lattvec = np.array([
        [4.4309997546822144, 0.0000000000000000, 0.0000000000000000],
        [-2.2154998773411063, 3.8373583517174135, 0.0000000000000000],
        [0.0000000000000004, 0.0000000000000008, 7.1570000621175227],
    ])
    assert np.allclose(lattvec, loader.atoms.get_cell())
    assert loader.atoms.get_chemical_formula() == "Li2Sb2Zn2"
    assert loader.dosweight == 1.0
    assert loader.kpoints.shape == (147, 3)
    assert np.allclose(loader.fermi, 5.23774582 * BoltzTraP2.units.eV)
    assert loader.magmom.shape == (6, )
    assert np.allclose(loader.magmom, np.ones_like(loader.magmom))
    # Note that this test involves twice the number of actual bands,
    # combining a spin-up band and a spin-down band.
    assert np.allclose(
        loader.ebands[:, 0],
        np.array([
            -42.4628,
            -42.4618,
            -6.2431,
            -5.1410,
            -2.9092,
            -2.9092,
            -2.8756,
            -2.8756,
            -2.8578,
            -2.7401,
            -2.7401,
            -2.7157,
            -2.7157,
            -2.4807,
            -0.8251,
            3.7239,
            3.7239,
            4.7743,
            4.7743,
            5.0228,
            5.3455,
            5.7265,
            6.8867,
            6.8867,
            7.8775,
            8.9818,
            8.9819,
            9.7743,
            -42.4638,
            -42.4628,
            -6.2423,
            -5.1404,
            -2.9088,
            -2.9088,
            -2.8753,
            -2.8753,
            -2.8574,
            -2.7397,
            -2.7397,
            -2.7153,
            -2.7153,
            -2.4802,
            -0.8239,
            3.7245,
            3.7245,
            4.7749,
            4.7749,
            5.0235,
            5.3472,
            5.7276,
            6.8881,
            6.8881,
            7.8790,
            8.9830,
            8.9830,
            9.7771,
        ]) * BoltzTraP2.units.eV,
    )
    assert loader.nelect == 40.0
    assert loader.mommat.shape == (147, 56, 3)


def test_Wien2kLoader_exception():
    """Wien2kLoader should raise exceptions for non-Wien2k directories."""
    # Test with a directory that does not exist
    arg = tempfile.mktemp(dir=datadir)
    with pytest.raises(FileNotFoundError):
        BoltzTraP2.dft.Wien2kLoader(arg)
    # Test with a VASP directory
    arg = os.path.join(datadir, "Si.vasp")
    with pytest.raises(BoltzTraP2.dft.LoaderError):
        BoltzTraP2.dft.Wien2kLoader(arg)


def test_Wien2kLoader():
    """Wien2kLoader should work correctly for the examples in datadir."""
    # CoSb3, non-spin-polarized, no derivatives
    loader = BoltzTraP2.dft.Wien2kLoader(os.path.join(datadir, "CoSb3"))
    lattvec = 4.5192517 * (np.ones((3, 3)) - 2.0 * np.eye(3))
    assert np.allclose(lattvec, loader.atoms.get_cell())
    assert loader.atoms.get_chemical_formula() == "Co4Sb12"
    assert loader.dosweight == 2.0
    assert loader.kpoints.shape == (1030, 3)
    assert loader.ebands.shape == (174, 1030)
    assert np.allclose(loader.fermi, 0.27735)
    assert np.allclose(
        loader.ebands[:, 0],
        np.array([
            -3.73463841647560,
            -3.73463841642918,
            -3.73463841640897,
            -3.73435523977001,
            -3.73435523974776,
            -3.73435523970237,
            -3.73421620029597,
            -3.73421620026555,
            -3.72480140663210,
            -3.72480140663161,
            -3.72480140663150,
            -3.72458514666301,
            -1.56003537096039,
            -1.55918426079551,
            -1.55918426078980,
            -1.55875742728928,
            -1.55875742728896,
            -1.55875742728873,
            -1.55692219694607,
            -1.55692219693480,
            -1.55528179033583,
            -1.55466315913575,
            -1.55466315913570,
            -1.55466315913564,
            -1.55412373994489,
            -1.55412373992191,
            -1.55412373991452,
            -1.55386947660137,
            -1.55386947660120,
            -1.55386947660096,
            -1.55212412541993,
            -1.55212412541989,
            -1.55212412541960,
            -1.55113168292276,
            -1.55113168290155,
            -1.55113168286976,
            -1.55050063579531,
            -1.55050063579515,
            -1.55050063579500,
            -1.54924987219229,
            -1.54924987217554,
            -1.54924987215810,
            -1.54813339034767,
            -1.54813339031553,
            -1.54727761103978,
            -1.54717993597931,
            -1.54717993597929,
            -1.54717993597913,
            -1.54635477897640,
            -1.54635477897488,
            -1.54635477897309,
            -1.54616747759109,
            -1.54616747759064,
            -1.54599909394886,
            -1.54447752166287,
            -1.54447752165437,
            -1.54447752163775,
            -1.54379884628505,
            -1.54379884623536,
            -1.54379884619455,
            -1.54362744120027,
            -1.54362744120011,
            -1.54362744119979,
            -1.54303137953865,
            -1.54303137953857,
            -1.54303137953813,
            -1.54294855751815,
            -1.54121817996058,
            -1.54121817995473,
            -1.54121817994603,
            -1.54059659340764,
            -1.54059659340643,
            -0.359328350931158,
            -0.182765097587536,
            -0.182765097483322,
            -0.137708845231410,
            -0.137708845228477,
            -0.137708845224276,
            -0.115461600851710,
            -0.115461600848430,
            -0.115461600846927,
            -6.81145418417447e-002,
            -6.81145417813901e-002,
            -6.81145417172872e-002,
            0.200581306115641,
            0.222544948230447,
            0.222544948257286,
            0.293540950236423,
            0.293540950292561,
            0.293540950335412,
            0.306226661180204,
            0.306226661182643,
            0.306226661183737,
            0.312719696603210,
            0.368413337566623,
            0.368413337578068,
            0.368413337582171,
            0.424914501228952,
            0.424914501232549,
            0.425308449073675,
            0.425308449074400,
            0.425308449075034,
            0.439876940454189,
            0.439876940454664,
            0.439876940455189,
            0.448917428333390,
            0.448917428337467,
            0.448917428337561,
            0.466627681747041,
            0.466627681747667,
            0.466627681748022,
            0.468673915777588,
            0.468673915779049,
            0.468673915780460,
            0.472796185743987,
            0.472796185748075,
            0.510010736369430,
            0.510010736372805,
            0.518137080745422,
            0.572289941553569,
            0.593895644389548,
            0.593895644395553,
            0.593895644429100,
            0.624265045316380,
            0.624265045364742,
            0.703737216097031,
            0.703737216097567,
            0.703737216097703,
            0.724230278402419,
            0.724230278459957,
            0.724230278473184,
            0.849200631987119,
            0.849200631987847,
            0.849200631992351,
            0.849496778875068,
            0.849496779061085,
            0.849496779145805,
            0.883663901914265,
            0.883663902135171,
            0.883663902342915,
            0.901040754510086,
            0.982613609948701,
            0.982613609954447,
            0.982613609960709,
            1.00146820979904,
            1.00146821000251,
            1.09428519178085,
            1.09428519212920,
            1.09428519228533,
            1.12259758102634,
            1.12259758103595,
            1.12259758103685,
            1.13729681953289,
            1.22653757458904,
            1.22653757460652,
            1.22653757463700,
            1.24388647847262,
            1.24388647848048,
            1.24388647849231,
            1.25587571662321,
            1.25587571712124,
            1.29723153395052,
            1.29723153432242,
            1.29723153565155,
            1.33023282794837,
            1.33023282795577,
            1.33023282796113,
            1.34373502554745,
            1.34985482959001,
            1.39485949364046,
            1.39485949373952,
            1.41984537743129,
            1.46702146313530,
            1.46702146315164,
        ]) * 0.5,
    )
    with pytest.raises(AttributeError):
        loader.nelect
    with pytest.raises(AttributeError):
        loader.mommat
    # LiZnSb, non-spin-polarized, no derivatives
    loader = BoltzTraP2.dft.Wien2kLoader(os.path.join(datadir, "LiZnSb"))
    lattvec = np.array([
        [4.430999754682214, 0.0, 0.0],
        [-2.2154998773411063, 3.8373583517174135, 0.0],
        [0.0, 0.0, 7.157000062117523],
    ])
    assert np.allclose(lattvec, loader.atoms.get_cell())
    assert loader.atoms.get_chemical_formula() == "Li2Sb2Zn2"
    assert loader.dosweight == 2.0
    assert loader.kpoints.shape == (3915, 3)
    assert loader.ebands.shape == (84, 3915)
    assert np.allclose(loader.fermi, 0.18825)
    assert np.allclose(
        loader.ebands[:, 1933],
        np.array([
            -3.11298118802972,
            -3.11297376497483,
            -1.74629198293860,
            -1.74620289553781,
            -1.74604033012386,
            -1.74600870079737,
            -1.74541740457343,
            -1.74528266884834,
            -1.74455997353295,
            -1.74446776600289,
            -1.74432913797639,
            -1.74427277828611,
            -0.387703638839634,
            -0.380375650672808,
            -0.202799053575902,
            -0.202646299996772,
            -0.201825754216763,
            -0.201516706624145,
            -0.198437416982047,
            -0.197745356769440,
            -0.192632906331954,
            -0.192410467494166,
            -0.189916777885856,
            -0.188218254232355,
            3.028873258369621e-002,
            4.543378546076464e-002,
            0.184186074680079,
            0.199634048735986,
            0.220531078302292,
            0.233551527809091,
            0.476709728671581,
            0.488094554048860,
            0.567559336138779,
            0.596332565590978,
            0.648837043393886,
            0.663113173292134,
            0.706858814167850,
            0.730002835276008,
            0.845855511748313,
            0.857569077614957,
            0.903689670976049,
            0.927260142222797,
            1.01707452238941,
            1.02248838988166,
            1.08010897124484,
            1.11201024315307,
            1.14656793407051,
            1.19164691854229,
            1.24721649114387,
            1.26234246720120,
            1.34714103575420,
            1.36976742361658,
            1.42884709852011,
            1.44934221766529,
            1.50073705647186,
            1.52095442254082,
            1.55075690872352,
            1.56401548433919,
            1.62370924346336,
            1.65244279531598,
            1.72456272592120,
            1.75481203247900,
            1.81304765700583,
            1.81972625856161,
            1.86238761511153,
            1.87840976476296,
            1.93945202759458,
            1.95711229270266,
            2.00901324070075,
            2.02403326427322,
            2.05244635416862,
            2.06211262290772,
            2.10644598691256,
            2.13453196201132,
            2.16825829603857,
            2.18468850958524,
            2.23238665419883,
            2.24133191133065,
            2.27984421018629,
            2.30286126504455,
            2.33993194308868,
            2.35034473080052,
            2.38186216984168,
            2.41888460515334,
        ]) * 0.5,
    )
    with pytest.raises(AttributeError):
        loader.nelect
    with pytest.raises(AttributeError):
        loader.mommat
    # Bi2Te3, spin-polarized, no derivatives
    loader = BoltzTraP2.dft.Wien2kLoader(os.path.join(datadir, "Bi2Te3"))
    lattvec = np.array([
        [2.19300006e00, 1.26612917e00, 1.01656669e01],
        [-2.19300006e00, 1.26612917e00, 1.01656669e01],
        [-4.36662431e-17, -2.53225835e00, 1.01656669e01],
    ])
    assert np.allclose(lattvec, loader.atoms.get_cell())
    assert loader.atoms.get_chemical_formula() == "Bi2Te3"
    assert loader.dosweight == 1.0
    assert loader.kpoints.shape == (4960, 3)
    assert loader.ebands.shape == (174, 4960)
    assert np.allclose(loader.fermi, 0.1726)
    assert np.allclose(
        loader.ebands[:, 0],
        np.array([
            -2.41267929981216,
            -2.41267929981216,
            -2.41251000540029,
            -2.41251000540029,
            -2.41215212391413,
            -2.41215212391413,
            -2.41202548696817,
            -2.41202548696817,
            -2.39757611919898,
            -2.39757611919898,
            -2.39546735788538,
            -2.39546735788538,
            -2.30553269860384,
            -2.30553269860384,
            -2.30544654536280,
            -2.30544654536280,
            -2.30538004584279,
            -2.30538004584279,
            -2.30506284590116,
            -2.30506284590116,
            -2.30474788080467,
            -2.30474788080467,
            -2.30471534186613,
            -2.30471534186613,
            -2.29053115622677,
            -2.29053115622677,
            -2.28971794289034,
            -2.28971794289034,
            -2.28779590431817,
            -2.28779590431817,
            -1.46738917822097,
            -1.46738917822097,
            -1.46717629444052,
            -1.46717629444052,
            -1.46706217371172,
            -1.46706217371172,
            -1.46675464923153,
            -1.46675464923153,
            -1.25115657879234,
            -1.25115657879234,
            -1.25087113975553,
            -1.25087113975553,
            -1.25024501388982,
            -1.25024501388982,
            -1.25017762894675,
            -1.25017762894675,
            -1.25004434989550,
            -1.25004434989550,
            -1.24964229684339,
            -1.24964229684339,
            -0.627717866889758,
            -0.627717866889755,
            -0.573937635684946,
            -0.573937635684944,
            -0.489173180216401,
            -0.489173180216401,
            -0.304110697937364,
            -0.304110697937364,
            -0.268766737230493,
            -0.268766737230493,
            0.150989806000198,
            0.150989806000199,
            0.152868204585908,
            0.152868204585909,
            0.214351744423932,
            0.214351744423934,
            0.230474818393341,
            0.230474818393341,
            0.250282192816133,
            0.250282192816134,
            0.254302393290113,
            0.254302393290113,
            0.260581004300737,
            0.260581004300737,
            0.305981413978614,
            0.305981413978616,
            0.316433350011097,
            0.316433350011098,
            0.358220177532369,
            0.358220177532371,
            0.415344817678491,
            0.415344817678492,
            0.421655967665282,
            0.421655967665282,
            0.439224011789180,
            0.439224011789180,
            0.533154611890972,
            0.533154611890972,
            0.549982656486813,
            0.549982656486814,
            0.708645297161478,
            0.708645297161479,
            0.747590143279449,
            0.747590143279451,
            0.788757005437116,
            0.788757005437116,
            0.793754179408734,
            0.793754179408734,
            0.817667154548906,
            0.817667154548906,
            0.817888305044977,
            0.817888305044980,
            0.827016373423174,
            0.827016373423175,
            0.881014764269847,
            0.881014764269848,
            0.881441548968587,
            0.881441548968587,
            0.945523371224617,
            0.945523371224620,
            0.955161282663239,
            0.955161282663239,
            0.960361198708318,
            0.960361198708319,
            1.04409390922416,
            1.04409390922416,
            1.07326276057588,
            1.07326276057588,
            1.08265928833535,
            1.08265928833535,
            1.11467188239321,
            1.11467188239321,
            1.11569098355775,
            1.11569098355775,
            1.16971166287101,
            1.16971166287102,
            1.18708042299856,
            1.18708042299856,
            1.23283735895467,
            1.23283735895467,
            1.24667127089730,
            1.24667127089730,
            1.25708866439195,
            1.25708866439195,
            1.34949404471641,
            1.34949404471642,
            1.51538887144868,
            1.51538887144868,
            1.51841067520095,
            1.51841067520095,
            1.55318638179925,
            1.55318638179925,
            1.56354593017977,
            1.56354593017977,
            1.60197772182632,
            1.60197772182633,
            1.69377010325688,
            1.69377010325688,
            1.69644555610647,
            1.69644555610647,
            1.73748594308455,
            1.73748594308455,
            1.74792989533757,
            1.74792989533757,
            1.75830703758515,
            1.75830703758515,
            1.76014016811743,
            1.76014016811743,
            1.78725540217984,
            1.78725540217984,
            1.79691530266400,
            1.79691530266400,
            1.80084581086649,
            1.80084581086649,
            1.87951179523627,
            1.87951179523627,
            1.92700196207763,
            1.92700196207763,
            1.93057247357028,
            1.93057247357029,
            1.93842367861019,
            1.93842367861019,
            1.96983622556958,
            1.96983622556958,
        ]) * 0.5,
    )
    with pytest.raises(AttributeError):
        loader.nelect
    with pytest.raises(AttributeError):
        loader.mommat
    # Bi2Te3, non-spin-polarized, with derivatives
    loader = BoltzTraP2.dft.Wien2kLoader(os.path.join(datadir, "Si"))
    lattvec = 2.72526263 * (np.ones((3, 3)) - np.eye(3))
    assert np.allclose(lattvec, loader.atoms.get_cell())
    assert loader.atoms.get_chemical_formula() == "Si2"
    assert loader.dosweight == 2.0
    assert loader.kpoints.shape == (455, 3)
    assert loader.ebands.shape == (6, 455)
    assert np.allclose(loader.fermi, 0.18834089155)
    assert np.allclose(
        loader.ebands[:, 211],
        np.array([
            -0.33037774141413601,
            -9.5387573974942297e-002,
            0.18036284164149730,
            0.27272390344136521,
            0.55698892682197232,
            0.63038485365766084,
        ]) * 0.5,
    )
    with pytest.raises(AttributeError):
        loader.nelect
    assert loader.mommat.shape == (455, 6, 3)


def test_ABINITLoader_exception():
    """ABINITLoader should raise exceptions for non-ABINIT directories."""
    # Test with a directory that does not exist
    # Note that mktemp is vulnerable to race conditions and should not be used
    # to create actual temporary directories, but is perfect to create a
    # random directory name that does not exist.
    arg = tempfile.mktemp(dir=datadir)
    with pytest.raises(FileNotFoundError):
        BoltzTraP2.dft.ABINITLoader(arg)
    # Test with a Wien2k directory
    arg = os.path.join(datadir, "Si")
    with pytest.raises(BoltzTraP2.dft.LoaderError):
        BoltzTraP2.dft.ABINITLoader(arg)


def test_GENELoader_exception():
    """GENELoader should raise exceptions for non-GENE directories."""
    # Test with a directory that does not exist
    arg = tempfile.mktemp(dir=datadir)
    with pytest.raises(FileNotFoundError):
        BoltzTraP2.dft.GENELoader(arg)
    # Test with a Wien2k directory
    arg = os.path.join(datadir, "Si")
    with pytest.raises(BoltzTraP2.dft.LoaderError):
        BoltzTraP2.dft.GENELoader(arg)


def test_GENELoader():
    """GENELoader should work correctly for the two examples in datadir."""
    # Li, non-spin-polarized, without derivatives
    loader = BoltzTraP2.dft.GENELoader(os.path.join(datadir, "Li"))
    lattvec = (0.3192047850e01 * (np.ones(
        (3, 3)) - 2.0 * np.eye(3)) / BoltzTraP2.units.Angstrom)
    assert np.allclose(lattvec, loader.atoms.get_cell())
    assert loader.atoms.get_chemical_formula() == "Li"
    assert loader.dosweight == 2.0
    assert loader.kpoints.shape == (413, 3)
    assert np.allclose(loader.fermi, 0.6940745476e-01 / 2)
    assert loader.ebands.shape == (3, 413)
    assert np.allclose(
        loader.ebands[:, 20],
        np.array([-0.1963709022e-01, 0.5029658971e00, 0.8701510047e00]) / 2.0,
    )
    with pytest.raises(AttributeError):
        loader.nelect
    with pytest.raises(AttributeError):
        loader.mommat
    # Li, spin-polarized, with derivatives (from VASP)
    loader = BoltzTraP2.dft.GENELoader(
        os.path.join(datadir, "Li.GENE.fromvasp"))
    lattvec = (3.231892535921517 * (np.ones(
        (3, 3)) - 2.0 * np.eye(3)) / BoltzTraP2.units.Angstrom)
    assert np.allclose(lattvec, loader.atoms.get_cell())
    assert loader.atoms.get_chemical_formula() == "Li"
    assert loader.dosweight == 1.0
    assert loader.kpoints.shape == (286, 3)
    assert np.allclose(loader.fermi, 0.024951731946)
    assert loader.ebands.shape == (12, 286)
    assert np.allclose(
        loader.ebands[:, 205],
        [
            -1.68363713,
            0.04813426,
            0.18563185,
            0.30062416,
            0.47661299,
            0.55504707,
            -1.68383558,
            0.047572,
            0.18507326,
            0.30008762,
            0.47603602,
            0.55446643,
        ],
    )
    with pytest.raises(AttributeError):
        loader.nelect
    assert loader.mommat.shape == (286, 12, 3)
    assert np.allclose(
        loader.mommat[205, :, :],
        np.array([
            [3.32542058e-04, 8.16769966e-05, -1.55575232e-05],
            [1.16308043e-01, -2.65839177e-02, 2.03122912e-02],
            [-8.43431672e-02, 1.19983508e-01, -1.81672977e-02],
            [-8.65017735e-02, -4.95137622e-02, -8.76724771e-02],
            [-1.47465873e-02, 1.28310672e-02, 8.57472336e-02],
            [-4.77188129e-02, -1.68221553e-01, 1.73855321e-03],
            [3.34486748e-04, 8.16769966e-05, -1.55575232e-05],
            [1.16298320e-01, -2.65780836e-02, 2.03103465e-02],
            [-8.43470565e-02, 1.19973785e-01, -1.81828552e-02],
            [-8.65153863e-02, -4.95254303e-02, -8.76452515e-02],
            [-1.47485320e-02, 1.28349566e-02, 8.57783487e-02],
            [-4.76663063e-02, -1.68244890e-01, 1.70160410e-03],
        ]),
    )


def test_CASTEPLoader_exception():
    """CASTEPLoader should raise exceptions for non-CASTEP directories."""
    # Test with a directory that does not exist
    arg = tempfile.mktemp(dir=datadir)
    with pytest.raises(FileNotFoundError):
        BoltzTraP2.dft.CASTEPLoader(arg)
    # Test with a Wien2k directory
    arg = os.path.join(datadir, "Si")
    with pytest.raises(BoltzTraP2.dft.LoaderError):
        BoltzTraP2.dft.CASTEPLoader(arg)


def test_CASTEPLoader():
    """CASTEPLoader should work correctly for some examples in datadir."""
    # Si, non-spin-polarized
    loader = BoltzTraP2.dft.CASTEPLoader(os.path.join(datadir, "Si.CASTEP"))
    lattvec = 2.733556058 * (np.ones((3, 3)) - np.eye(3))
    assert np.allclose(lattvec, loader.atoms.get_cell())
    assert loader.atoms.get_chemical_formula() == "Si2"
    assert loader.dosweight == 2.0
    assert loader.kpoints.shape == (165, 3)
    assert np.allclose(loader.fermi, 0.166358)
    assert loader.ebands.shape == (23, 165)
    assert np.allclose(
        loader.ebands[:, 0],
        np.array([
            -0.15891922,
            -0.09325013,
            0.00930266,
            0.05187248,
            0.28135290,
            0.29101339,
            0.35553076,
            0.40099755,
            0.44444156,
            0.55167291,
            0.60224349,
            0.65964277,
            0.72946705,
            0.75735457,
            0.77934375,
            0.81626234,
            0.83611211,
            0.88855847,
            0.91090449,
            0.94260295,
            0.98547572,
            1.06291142,
            1.08848286,
        ]),
    )
    assert loader.nelect == 8.0
    with pytest.raises(AttributeError):
        loader.mommat
    assert loader.magmom is None
    # Si, spin-polarized, collinear
    loader = BoltzTraP2.dft.CASTEPLoader(
        os.path.join(datadir, "Si.CASTEP.spinpol"))
    lattvec = 2.733556058 * (np.ones((3, 3)) - np.eye(3))
    assert np.allclose(lattvec, loader.atoms.get_cell())
    assert loader.atoms.get_chemical_formula() == "Si2"
    assert loader.dosweight == 1.0
    assert loader.kpoints.shape == (165, 3)
    assert np.allclose(loader.fermi, 0.166360)
    assert loader.ebands.shape == (52, 165)
    assert np.allclose(
        loader.ebands[:, 0],
        np.array([
            -0.19667437,
            -0.09445701,
            0.11256214,
            0.11256214,
            0.20990416,
            0.27899862,
            0.27899862,
            0.43915862,
            0.52942637,
            0.52942637,
            0.56009783,
            0.58125162,
            0.58125162,
            0.63049844,
            0.83545408,
            0.90556023,
            0.91795071,
            0.91795071,
            0.96945344,
            0.96945345,
            1.07947675,
            1.11640248,
            1.11640248,
            1.14909683,
            1.14909683,
            1.19448599,
            -0.19667350,
            -0.09445607,
            0.11256280,
            0.11256280,
            0.20990466,
            0.27900125,
            0.27900125,
            0.43915659,
            0.52942859,
            0.52942859,
            0.56009440,
            0.58125299,
            0.58125299,
            0.63049790,
            0.83545222,
            0.90556067,
            0.91795138,
            0.91795138,
            0.96945380,
            0.96945380,
            1.07946806,
            1.11640654,
            1.11640654,
            1.14909780,
            1.14909781,
            1.19448921,
        ]),
    )
    assert loader.nelect == 8.0
    with pytest.raises(AttributeError):
        loader.mommat
    assert loader.magmom.shape == (2, )
    assert np.allclose(loader.magmom, np.array([1.0, 2.0]))
    # Si, spin-polarized, non-collinear with spin-orbit coupling
    loader = BoltzTraP2.dft.CASTEPLoader(
        os.path.join(datadir, "Si.CASTEP.noncoll.soc"))
    lattvec = 2.733556058 * (np.ones((3, 3)) - np.eye(3))
    assert np.allclose(lattvec, loader.atoms.get_cell())
    assert loader.atoms.get_chemical_formula() == "Si2"
    assert loader.dosweight == 2.0
    assert loader.kpoints.shape == (2457, 3)
    assert np.allclose(loader.fermi, 0.166359)
    assert loader.ebands.shape == (44, 2457)
    assert np.allclose(
        loader.ebands[:, 0],
        np.array([
            -0.19667408,
            -0.19667348,
            -0.09445698,
            -0.09445639,
            0.11256244,
            0.11256260,
            0.11256289,
            0.11256306,
            0.20990471,
            0.20990557,
            0.27899726,
            0.27899729,
            0.27899841,
            0.27899843,
            0.43915283,
            0.43915571,
            0.52942731,
            0.52942743,
            0.52942790,
            0.52942801,
            0.56009428,
            0.56009608,
            0.58125143,
            0.58125156,
            0.58125241,
            0.58125256,
            0.63049698,
            0.63049835,
            0.83545263,
            0.83545400,
            0.90556129,
            0.90556231,
            0.91795028,
            0.91795050,
            0.91795115,
            0.91795137,
            0.96945298,
            0.96945314,
            0.96945382,
            0.96945396,
            1.07946294,
            1.07946847,
            1.11640619,
            1.11640623,
        ]),
    )
    assert loader.nelect == 8.0
    with pytest.raises(AttributeError):
        loader.mommat
    assert loader.magmom.shape == (2, 3)
    assert np.allclose(loader.magmom,
                       np.array([[1.0, 2.0, 3.0], [2.0, 3.0, 4.0]]))


def test_DFTData_exceptions():
    """DFTData should raise exceptions for invalid directories."""
    # Test with a non-existing directory
    arg = tempfile.mktemp(dir=datadir)
    with pytest.raises(FileNotFoundError):
        BoltzTraP2.dft.DFTData(arg)
    # Test with a directory that does not contain any calculation
    with pytest.raises(ValueError):
        BoltzTraP2.dft.DFTData(datadir)
    # Test with a directory that does not contain derivative information
    with pytest.raises(ValueError):
        BoltzTraP2.dft.DFTData(os.path.join(datadir, "CoSb3"),
                               derivatives=True)


def test_DFTData_VASP():
    """DFTData should process VASP directories correctly."""
    # Tests without derivatives and discarding derivatives
    data = BoltzTraP2.dft.DFTData(os.path.join(datadir, "Si.vasp.noder"))
    assert data.mommat is None
    assert data.source == "VASP"
    data = BoltzTraP2.dft.DFTData(os.path.join(datadir, "Si.vasp"),
                                  derivatives=False)
    assert data.mommat is None
    assert data.source == "VASP"
    # Test with derivatives
    data = BoltzTraP2.dft.DFTData(os.path.join(datadir, "Si.vasp"),
                                  derivatives=True)
    assert data.mommat is not None
    assert data.source == "VASP"
    # Test against data loaded directly
    loaded = BoltzTraP2.dft.VASPLoader(os.path.join(datadir, "Si.vasp"))
    assert data.sysname == loaded.sysname
    assert data.atoms == loaded.atoms
    assert data.dosweight == loaded.dosweight
    assert data.nelect == loaded.nelect
    assert np.array_equal(data.kpoints, loaded.kpoints)
    assert data.fermi == loaded.fermi
    assert np.array_equal(data.ebands, loaded.ebands)
    assert np.array_equal(data.mommat, loaded.mommat)


def test_DFTData_Wien2k():
    """DFTData should process Wien2k directories correctly."""
    # Tests without derivatives and discarding derivatives
    data = BoltzTraP2.dft.DFTData(os.path.join(datadir, "CoSb3"))
    assert data.mommat is None
    assert data.source == "Wien2k"
    data = BoltzTraP2.dft.DFTData(os.path.join(datadir, "Si"),
                                  derivatives=False)
    assert data.mommat is None
    assert data.source == "Wien2k"
    # Test with derivatives
    data = BoltzTraP2.dft.DFTData(os.path.join(datadir, "Si"),
                                  derivatives=True)
    assert data.mommat is not None
    assert data.source == "Wien2k"
    # Test against data loaded directly
    loaded = BoltzTraP2.dft.Wien2kLoader(os.path.join(datadir, "Si"))
    assert data.sysname == loaded.sysname
    assert data.atoms == loaded.atoms
    assert data.dosweight == loaded.dosweight
    assert np.array_equal(data.kpoints, loaded.kpoints)
    assert data.fermi == loaded.fermi
    assert np.array_equal(data.ebands, loaded.ebands)
    assert np.array_equal(data.mommat, loaded.mommat)
    # Check that the number of valence electrons was calculated correctly
    assert data.nelect == 8.0


def test_DFTData_ABINIT():
    """DFTData should process ABINIT directories correctly."""
    # Tests without derivatives and discarding derivatives
    data = BoltzTraP2.dft.DFTData(os.path.join(datadir, "Si.abinit"))
    assert data.mommat is None
    assert data.source == "ABINIT"
    # Test against data loaded directly
    loaded = BoltzTraP2.dft.ABINITLoader(os.path.join(datadir, "Si.abinit"))
    assert data.sysname == loaded.sysname
    assert data.atoms == loaded.atoms
    assert data.dosweight == loaded.dosweight
    assert data.nelect == loaded.nelect
    assert np.array_equal(data.kpoints, loaded.kpoints)
    assert data.fermi == loaded.fermi
    assert np.array_equal(data.ebands, loaded.ebands)
    with pytest.raises(AttributeError):
        assert np.array_equal(data.mommat, loaded.mommat)


def test_DFTData_GENE():
    """DFTData should process GENE directories correctly."""
    # Tests without derivatives and discarding derivatives
    data = BoltzTraP2.dft.DFTData(os.path.join(datadir, "Li"))
    assert data.mommat is None
    assert data.source == "GENE"
    data = BoltzTraP2.dft.DFTData(os.path.join(datadir, "Li.GENE.fromvasp"),
                                  derivatives=False)
    assert data.mommat is None
    assert data.source == "GENE"
    # Test with derivatives
    data = BoltzTraP2.dft.DFTData(os.path.join(datadir, "Li.GENE.fromvasp"),
                                  derivatives=True)
    assert data.mommat is not None
    assert data.source == "GENE"
    # Test against data loaded directly
    loaded = BoltzTraP2.dft.GENELoader(
        os.path.join(datadir, "Li.GENE.fromvasp"))
    assert data.sysname == loaded.sysname
    assert data.atoms == loaded.atoms
    assert data.dosweight == loaded.dosweight
    assert np.array_equal(data.kpoints, loaded.kpoints)
    assert data.fermi == loaded.fermi
    assert np.array_equal(data.ebands, loaded.ebands)
    assert np.array_equal(data.mommat, loaded.mommat)
    # Check that the number of valence electrons was calculated correctly
    assert data.nelect == 3.0


def test_DFTData_CASTEP():
    """DFTData should process CASTEP directories correctly."""
    # MgO test data
    data = BoltzTraP2.dft.DFTData(os.path.join(datadir, "Si.CASTEP"))
    assert data.mommat is None
    assert data.source == "CASTEP"
    # Test against data loaded directly
    loaded = BoltzTraP2.dft.CASTEPLoader(os.path.join(datadir, "Si.CASTEP"))
    assert data.sysname == loaded.sysname
    assert data.atoms == loaded.atoms
    assert data.dosweight == loaded.dosweight
    assert np.array_equal(data.kpoints, loaded.kpoints)
    assert data.fermi == loaded.fermi
    assert np.array_equal(data.ebands, loaded.ebands)
    # Check that the number of valence electrons was calculated correctly
    assert data.nelect == 8.0


def test_DFTData_bandana():
    """DFTData should discard bands correctly on demand."""
    data = BoltzTraP2.dft.DFTData(os.path.join(datadir, "Si"),
                                  derivatives=True)
    assert data.ebands.shape[0] == 6
    assert data.nelect == 8.0
    assert data.mommat.shape[1] == data.ebands.shape[0]
    accepted = data.bandana(-0.02, 0.2)
    assert data.ebands.shape[0] == 3
    assert np.count_nonzero(accepted) == 3
    assert data.nelect == 6.0
    assert data.mommat.shape[1] == data.ebands.shape[0]


def test_DFTData_get_lattvec():
    """DFTData should return the right lattice vectors."""
    data = BoltzTraP2.dft.DFTData(os.path.join(datadir, "CoSb3"))
    assert np.allclose(data.get_lattvec(),
                       data.atoms.get_cell().T * BoltzTraP2.units.Angstrom)


def test_DFTData_get_volume():
    """DFTData should return the right unit cell volume."""
    data = BoltzTraP2.dft.DFTData(os.path.join(datadir, "CoSb3"))
    assert np.allclose(data.get_volume(), np.abs(la.det(data.get_lattvec())))


def test_DFTData_get_formula_count():
    """DFTData should yield the right number of formula units in the
    irreducible cell.
    """
    data = BoltzTraP2.dft.DFTData(os.path.join(datadir, "CoSb3"))
    assert data.get_formula_count() == 4
    data = BoltzTraP2.dft.DFTData(os.path.join(datadir, "Si"))
    assert data.get_formula_count() == 2


### Adding FHI-AIMS tests
def test_AIMSLoader_exception():
    """AIMSLoader should raise exceptions for non-AIMS directories."""
    # Test with a directory that does not exist
    arg = tempfile.mktemp(dir=datadir)
    with pytest.raises(FileNotFoundError):
        BoltzTraP2.dft.AIMSLoader(arg)


def test_AIMSLoader():
    """AIMSLoader should work correctly for the MgO examples."""
    # First test: Unpolarized calculation without SOC for bulk system.
    loader = BoltzTraP2.dft.AIMSLoader(
        os.path.join(datadir, "MgO.AIMS/RHF_ZORA"))
    # Compare the lattice parameter with the value in the geometry.in file.
    # Since "loaders.atoms" is an ASE object, this value is in A.
    lattvec = np.array([
        [3.0097881300000000, 0.0000000000000000, 0.0000000000000000],
        [1.5048940650000004, 2.6065529805888605, 0.0000000000000000],
        [1.5048940650000004, 0.8688509935296203, 2.4574817174618540],
    ])
    assert np.allclose(lattvec, loader.atoms.get_cell())
    # Check that the numbmer of atoms and their positions are correct
    direct = loader.atoms.get_scaled_positions()
    assert direct.shape == (2, 3)
    # Chemical formula, also built from geometry.in
    assert loader.atoms.get_chemical_formula() == "MgO"
    # Non-spin-polarized calculation (from the EIG file)
    assert loader.dosweight == 2.0
    # Fermi Energy (from the .EIG file)
    assert np.isclose(loader.fermi, -0.20414655118011638)
    # Check the shapes of the k-point and energy arrays
    assert loader.kpoints.shape == (1063, 3)
    # Our AIMS loader does not support band derivatives. Check that they have
    # not been loaded.
    with pytest.raises(AttributeError):
        loader.mommat


def test_DFTData_AIMS():
    """DFTData should process AIMS directories correctly."""
    # MgO test data
    data = BoltzTraP2.dft.DFTData(os.path.join(datadir, "MgO.AIMS/RHF_ZORA"))
    assert data.mommat is None
    assert data.source == "AIMS"
    # Test against data loaded directly
    loaded = BoltzTraP2.dft.AIMSLoader(
        os.path.join(datadir, "MgO.AIMS/RHF_ZORA"))
    assert data.sysname == loaded.sysname
    assert data.atoms == loaded.atoms
    assert data.dosweight == loaded.dosweight
    assert np.array_equal(data.kpoints, loaded.kpoints)
    assert data.fermi == loaded.fermi
    assert np.array_equal(data.ebands, loaded.ebands)
    # Check that the number of valence electrons was calculated correctly
    assert data.nelect == 20.0
