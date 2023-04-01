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

import os.path

import numpy as np
import matplotlib.pylab as pl

from environment import data_dir

from BoltzTraP2 import io as IO
from BoltzTraP2 import dft as BTP
from BoltzTraP2 import sphere
from BoltzTraP2 import fite
from BoltzTraP2 import bandlib as BL
from BoltzTraP2 import serialization

id1 = "kp9"

data, equivalences, coeffs, metadata = serialization.load_calculation(
    "Si_mom" + id1 + ".bt2")
data, equivalences, coeffs2, metadata = serialization.load_calculation(
    "Si_nomom" + id1 + ".bt2")

id1 = "bands"


class CustomW2kLoader(BTP.GenericWien2kLoader):
    """Custom loader for Wien2k files that will use the right files."""

    def __init__(self, directory):
        BTP.GenericWien2kLoader.__init__(
            self, "Si_" + id1, 2.0, os.path.join(directory, "Si", "Si.scf"),
            os.path.join(directory, "Si", "Si.struct"),
            os.path.join(directory, "Si", "Si_" + id1 + ".energy"),
            os.path.join(directory, "Si", "Si_" + id1 + ".mommat2"))


BTP.register_loader("Wien2k_" + "Si_" + id1, CustomW2kLoader)
dataB = BTP.DFTData(data_dir, derivatives=True)

fig3, ax3 = pl.subplots(1, figsize=(4, 4))

Ef = (np.max(data.ebands[3]) + np.min(data.ebands[4])) / (2 * BL.eV)
nband = len(dataB.ebands)

kp = dataB.kpoints[:24]
lkp = np.linalg.norm(kp, axis=1)
ax3.set_xlim([0, lkp[-1]])
ax3.plot(
    [0, lkp[-1]], [.2167 / BL.eV - Ef, .2167 / BL.eV - Ef], "k:", linewidth=.5)
ax3.set_xticks([0, lkp[-1]])
ax3.tick_params(axis='y', labelsize=16)
ax3.set_xticklabels([r"$\Gamma$", "X"], fontsize=16)
ax3.set_ylim([-1.5, 1.5])
ax3.set_ylabel(r"$\varepsilon$ [eV]", fontsize=16)

for iband in range(nband):
    ax3.plot(lkp, dataB.ebands[iband, :24] / BL.eV - Ef, "k.")

nband = len(data.ebands)
II = np.nonzero((np.abs(data.kpoints[:, 1] - data.kpoints[:, 2]) < .00001) &
                (np.abs(data.kpoints[:, 0]) < .00001))[0]
kp = data.kpoints[II]
lkp = np.linalg.norm(kp, axis=1)
for i in range(nband):
    ax3.plot(lkp, data.ebands[i, II] / BL.eV - Ef, "o", color="C" + str(i))

kp1 = np.outer(np.linspace(0, 1, 101), np.array([0, .5, .5]))
lkp1 = np.linalg.norm(kp1, axis=1)

lattvec = data.get_lattvec()
egrid, vgrid = fite.getBands(kp1, equivalences, lattvec, coeffs)
egrid2, vgrid2 = fite.getBands(kp1, equivalences, lattvec, coeffs2)

for iband in range(nband):
    ax3.plot(lkp1, egrid[iband] / BL.eV - Ef, color="C" + str(iband))
    ax3.plot(
        lkp1,
        egrid2[iband] / BL.eV - Ef,
        color="C" + str(iband),
        linestyle="--")

fig3.tight_layout(pad=1.)
fig3.savefig("Si_pocket.pdf")
pl.show()
