# -*- coding: utf-8 -*
#    BoltzTraP2, a program for interpolating band structures and calculating
#                semi-classical transport coefficients.
#    Copyright (C) 2017-2020 Georg K. H. Madsen <georg.madsen@tuwien.ac.at>
#    Copyright (C) 2017-2020 Jesús Carrete <jesus.carrete.montana@tuwien.ac.at>
#    Copyright (C) 2017-2020 Matthieu J. Verstraete <matthieu.vers1traete@ulg.ac.be>
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

###############################################################################
# Load bt2 file for LiZnSb (from LZS.py) calculate Onsager coeffs and
# Reproduce Fig. 5c from J. Am. Chem. Soc. 128, (2006) 12145
###############################################################################

import os
import os.path
import pickle

import numpy as np
import matplotlib.pylab as pl

from environment import data_dir

from BoltzTraP2 import dft as BTP
from BoltzTraP2 import sphere
from BoltzTraP2 import fite
from BoltzTraP2 import bandlib as BL
from BoltzTraP2 import serialization

data = BTP.DFTData(os.path.join(data_dir, "LiZnSb"))
data.bandana(emin=0, emax=.3)

data, equivalences, coeffs, metadata = serialization.load_calculation(
    "LZS.bt2")

lattvec = data.get_lattvec()
equivalences = sphere.get_equivalences(data.atoms, data.magmom,
                                       len(data.kpoints) * 5)

eband, vvband, cband = fite.getBTPbands(equivalences, coeffs, lattvec)

epsilon, dos, vvdos, cdos = BL.BTPDOS(eband, vvband, npts=6000)

mur = epsilon[100:-100]
TEMP = np.array([500.])
N, L0, L1, L2, Lm11 = BL.fermiintegrals(
    epsilon, dos, vvdos, mur=mur, Tr=[TEMP])
#
UCvol = data.get_volume()
sigma, seebeck, kappa, Hall = BL.calc_Onsager_coefficients(
    L0, L1, L2, mur, TEMP, UCvol)
N0 = N[0] + data.nelect
S = seebeck[0]
sig = sigma[0]
fig, ax2 = pl.subplots(1, figsize=(6, 3))
ax2.set_xlim([-.06, .06])
ax2.set_ylim([0, 120])
ax2.plot(N0, S[:, 0, 0]**2 * sig[:, 0, 0] * 1e-10, "k-", label="xx")
ax2.plot(N0, S[:, 2, 2]**2 * sig[:, 2, 2] * 1e-10, "k--", label="zz")
ax2.set_xlabel("$n$ [e/uc]")
ax2.set_ylabel(r"$S^2\sigma/ \tau$ [10$^{14}$ $\mu$W/(cm K$^2$) s]")
ax2.legend()
fig.tight_layout(pad=1.)

pl.show()
