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

# Additional tests to try and make sure that the library achieves the right
# reduction of k points as VASP in a few representative cases. These mix the
# sphere and io modules.

import sys
import os
import os.path
import xml.etree.ElementTree as et

import numpy as np
import scipy as sp
import pytest

from BoltzTraP2.io import _parse_vasp_ikpoints
from BoltzTraP2.io import _parse_vasp_velocities
from BoltzTraP2.io import _parse_vasp_array
from BoltzTraP2.io import parse_vasprunxml
from BoltzTraP2.sphere.frontend import calc_reciprocal_degeneracies
from BoltzTraP2.sphere import calc_reciprocal_iksubset

mydir = os.path.abspath(os.path.dirname(__file__))
datadir = os.path.join(os.path.dirname(mydir), "data")


def test_kpoint_reduction_Si():
    """BoltzTraP2 should achieve the same k point reduction as VASP in the case
    of Si with a regular unpolarized calculation and in two additional
    calculations with initial (unphysical) magnetic configurations that lower
    the symmetry, one collinear and other noncollinear.
    """
    for lastpart in ("Si.vasp", "Si.vasp.magmom12", "Si.vasp.noncollinear"):
        dirname = os.path.join(datadir, lastpart)
        filename = os.path.join(dirname, "vasprun.xml")
        xml_tree = et.parse(filename)
        result = _parse_vasp_ikpoints(xml_tree)
        k, v = _parse_vasp_velocities(xml_tree)
        parsed = parse_vasprunxml(filename)
        path = "./kpoints/varray[@name=\"weights\"]"
        xml_weights = xml_tree.find(path)
        weights = []
        for p in xml_weights:
            weights.append(_parse_vasp_array(p.text))
        weights = np.array(weights)
        base = np.round(1. / weights[0])
        weights *= base
        weights = np.round(weights).astype(int).ravel()
        subset = calc_reciprocal_iksubset(parsed["atoms"], parsed["magmom"], k)
        # Assess both the number of k-point equivalence classes and the
        # degeneracy of each k point are assessed.
        assert result.shape[0] == len(subset)
        degeneracies = calc_reciprocal_degeneracies(
            parsed["atoms"], parsed["magmom"], parsed["kpoints"])
        for k, d, w in zip(parsed["kpoints"], degeneracies, weights):
            assert d == w
