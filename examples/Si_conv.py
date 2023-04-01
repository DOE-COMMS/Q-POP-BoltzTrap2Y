# -*- coding: utf-8 -*-
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
import pickle

import numpy as np
import matplotlib.pyplot as pl

from environment import data_dir

from BoltzTraP2 import io as IO
from BoltzTraP2 import dft as BTP
from BoltzTraP2 import sphere
from BoltzTraP2 import fite
from BoltzTraP2 import bandlib as BL
from BoltzTraP2 import serialization

TEMP = np.array([300.])

ids = [
    "kp7", "kp9", "kp11", "kp13", "kp15", "kp17", "kp19", "kp21", "kp23",
    "kp25"
]

Nkp = []
PFmom = []
PFnomom = []
Smom = []
Snomom = []
for icolor, id1 in enumerate(ids):

    data, equivalences, coeffs, metadata = serialization.load_calculation(
        "Si_mom" + id1 + ".bt2")
    data, equivalences, coeffs2, metadata = serialization.load_calculation(
        "Si_nomom" + id1 + ".bt2")

    print(id1, len(data.kpoints))
    Nkp += [len(data.kpoints)]
    lattvec = data.get_lattvec()

    eband, vvband, cband = fite.getBTPbands(equivalences, coeffs, lattvec)
    epsilon, dos, vvdos, cdos = BL.BTPDOS(
        eband, vvband, erange=[-.1, .35], npts=4500)
    mur = epsilon[100:-100]
    N, L0, L1, L2, Lm11 = BL.fermiintegrals(
        epsilon, dos, vvdos, mur=mur, Tr=TEMP)
    N += data.nelect
    UCvol = data.get_volume()
    sigma, seebeck, kappa, Hall = BL.calc_Onsager_coefficients(L0, L1, L2, mur,
                                                               TEMP, UCvol)
    Smom += [seebeck[0, 3067, 0, 0]]
    PFmom += [seebeck[0, 3067, 0, 0]**2 * sigma[0, 3067, 0, 0]]
    eband2, vvband2, cband = fite.getBTPbands(equivalences, coeffs2, lattvec)
    epsilon2, dos2, vvdos2, cdos = BL.BTPDOS(
        eband2, vvband2, erange=[-.1, .35], npts=4500)

    N, L0, L1, L2, Lm11 = BL.fermiintegrals(
        epsilon2, dos2, vvdos2, mur=mur, Tr=TEMP)
    N += data.nelect
    sigma2, seebeck2, kappa, Hall = BL.calc_Onsager_coefficients(
        L0, L1, L2, mur, TEMP, UCvol)

    Snomom += [seebeck2[0, 3067, 0, 0]]
    PFnomom += [seebeck2[0, 3067, 0, 0]**2 * sigma2[0, 3067, 0, 0]]

fig, ax1 = pl.subplots(1, figsize=(8, 4))

ax2 = ax1.twinx()

ax1.plot(np.log(Nkp), np.array(Smom) * 1E6, color="k", label="S")
ax1.plot(np.log(Nkp), np.array(Snomom) * 1E6, color="k", linestyle="--")
ax2.plot(np.log(Nkp), np.array(PFmom) * 1E-10, color="C0", label="S")
ax2.plot(np.log(Nkp), np.array(PFnomom) * 1E-10, color="C0", linestyle="--")
ax1.set_xlabel("#kpts IBZ", fontsize=14)
ax1.set_xticks(np.log(Nkp))
ax1.set_xticklabels(Nkp)
ax1.set_ylabel(r"$S$ [10$^{14}$ $\mu$V/K]", color="k", fontsize=14)
ax1.tick_params(axis='y', colors="k", labelsize=14)
ax1.tick_params(axis='x', colors="k", labelsize=14)
ax2.set_ylabel(
    r"$S^2\sigma/ \tau$ [10$^{14}$ $\mu$W/(cm K$^2$) s]",
    color="C0",
    fontsize=14)
ax2.tick_params(axis='y', colors="C0", labelsize=14)

fig.tight_layout(pad=1.)
fig.savefig("Si_conv.pdf")
pl.show()
