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

import os
import os.path
import logging

import numpy as np
import matplotlib.pylab as pl

from environment import data_dir

import BoltzTraP2.dft as BTP
import BoltzTraP2.bandlib as BL
import BoltzTraP2.io as IO
from BoltzTraP2 import sphere
from BoltzTraP2 import fite
from BoltzTraP2 import serialization

logging.basicConfig(
    level=logging.DEBUG, format="{levelname:8s}│ {message:s}", style="{")


def getBands2(kp, equivalences, lattvec, coeffs, atoms, magmom):
    tpii = 2j * np.pi
    phase = np.zeros((len(kp), len(equivalences)), dtype=complex)
    for j, equiv in enumerate(equivalences):
        phase0 = np.exp(tpii * kp @ equiv.T)
        vv = 1j * (lattvec @ equiv.T)
        phase[:, j] = np.sum(phase0, axis=1)
    nstar = np.array([len(equiv) for equiv in equivalences])
    phase /= nstar
    egrid = coeffs @ phase.T
    kstars = sphere.calc_reciprocal_stars(atoms, magmom, kp)
    vvgrid = np.zeros((len(coeffs), 3, len(kp)))
    for ik, kstar in enumerate(kstars):
        phase1 = np.zeros((len(coeffs), 3, len(kstar)))
        for j, equiv in enumerate(equivalences):
            vv = 1j * (lattvec @ equiv.T)
            phase0 = np.exp(tpii * kstar @ equiv.T)
            phase1 += np.outer(coeffs[:, j], np.dot(
                vv, phase0.T)).real.reshape(len(coeffs), 3, -1) / len(equiv)
        vvgrid[:, :, ik] = np.sum((phase1)**2, axis=2) / len(kstar)
    return egrid, vvgrid


def W2Kmat(filename, kpoints):
    """Read the contents of a Wien2k .mat_diag file.

    Args:
        filename: path to the .mat_diag file.
        kpoints: array with the k points of the same system.

    Returns:
        A 3-tuple. The first element is an array with the  matrix elements. The
        second and third are integers delimiting which bands have derivative
        information available.
    """
    nk = len(kpoints)

    fmat = open(filename, 'r')
    matlines = fmat.readlines()
    fmat.close()
    il = 2
    vvr = []
    brk = []
    for ik in range(nk):
        nemin = int(matlines[il].split()[5])
        nemax = int(matlines[il].split()[6])
        brk += [[nemin, nemax]]
        vvk = np.zeros((nemax - nemin + 1, 6))
        for ib in range(nemax - nemin + 1):
            line = matlines[il + 2 + ib]
            for i in range(6):
                vvk[ib, i] = float(line[11 + i * 13:24 + i * 13])
        vvr += [vvk]
        il += (nemax - nemin + 4)
    nemax = np.min(brk, axis=0)[1]
    nemin = np.max(brk, axis=0)[0]
    vvr2 = []
    for ik in range(nk):
        vvr2 += [vvr[ik][nemin - brk[ik][0]:nemax - brk[ik][0] + 1]]
    return np.array(vvr2), nemin, nemax


data = BTP.DFTData(os.path.join(data_dir, "LiZnSb"))
accepted = data.bandana(emin=0, emax=.3)
nemin = accepted.nonzero()[0][0]
nemax = accepted.nonzero()[0][-1] + 1

equivalences = sphere.get_equivalences(data.atoms, data.magmom,
                                       len(data.kpoints) * 5)
print("equivs", len(equivalences), " Should be 21164")

coeffs = fite.fitde3D(data, equivalences)
print("fit done")
metadata = serialization.gen_bt2_metadata(data, data.magmom is None)
serialization.save_calculation("LZS.bt2", data, equivalences, coeffs, metadata)

vvr, nemin2, nemax2 = W2Kmat(
    os.path.join(data_dir, "LiZnSb", "LiZnSb.mat_diag"), data.kpoints)

vvr = vvr[:, nemin - (nemin2 - 1):nemax - (nemin2 - 1)]
lattvec = data.get_lattvec()

eband, vvband, cband = fite.getBTPbands(equivalences, coeffs, lattvec)

nkpt = 51
kp = np.zeros((nkpt, 3))
kp[:, 0] = np.linspace(0, .5, nkpt)
egrid, vvgrid = getBands2(kp, equivalences, lattvec, coeffs, data.atoms,
                          data.magmom)

fig, ax3 = pl.subplots(1, figsize=(6, 4))
fig, ax4 = pl.subplots(1, figsize=(6, 4))
for iband in range(len(data.ebands)):
    ax3.plot(data.kpoints[:400:15, 1], data.ebands[iband, :400:15], "k.")
    ax3.plot(kp[:, 0], egrid[iband])
    ax4.plot(data.kpoints[:400:15, 1], vvr[:400:15, iband, 0], "k.")
    ax4.plot(kp[:, 0], vvgrid[iband, 0])

pl.show()
