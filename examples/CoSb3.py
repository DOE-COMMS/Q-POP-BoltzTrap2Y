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

###############################################################################
# Load Wien2k results for CoSb3, perform the interpolation and create a band
# diagram and calculate Hall coefficient. Reproduce Fig. 4 from
# Comp. Phys. Comm. 175 (2006) 67–71
###############################################################################

import os
import os.path
import logging

import numpy as np
import matplotlib.pylab as pl
from ase.dft.kpoints import get_special_points
from ase.dft.kpoints import bandpath

from environment import data_dir

from BoltzTraP2 import dft
from BoltzTraP2 import sphere
from BoltzTraP2 import fite
from BoltzTraP2 import bandlib
from BoltzTraP2 import serialization
from BoltzTraP2 import units

logging.basicConfig(
    level=logging.DEBUG, format="{levelname:8s}│ {message:s}", style="{")

dirname = os.path.join(data_dir, "CoSb3")
bt2file = "CoSb3.bt2"

# If a ready-made file with the interpolation results is available, use it
# Otherwise, create the file.
if not os.path.exists(bt2file):
    # Load the input
    data = dft.DFTData(dirname)
    # Select the interesting bands
    data.bandana(emin=data.fermi - .2, emax=data.fermi + .2)
    # Set up a k point grid with roughly five times the density of the input
    equivalences = sphere.get_equivalences(data.atoms, data.magmom,
                                           len(data.kpoints) * 5)
    # Perform the interpolation
    coeffs = fite.fitde3D(data, equivalences)
    # Save the result
    serialization.save_calculation(bt2file, data, equivalences, coeffs,
                                   serialization.gen_bt2_metadata(
                                       data, data.mommat is not None))

# Load the interpolation results from a file
data, equivalences, coeffs, metadata = serialization.load_calculation(bt2file)

# Create a simple path in reciprocal space adapted to the bcc structure
cell = data.get_lattvec()  #-np.eye(3) + .5
points = get_special_points("bcc", cell)
klist = "GHNGP"
PP = [points[k] for k in klist]
kpts, x, X = bandpath(PP, cell, 100)

# Rebuild the bands from the interpolation coefficients
lattvec = data.get_lattvec()
egrid, vgrid = fite.getBands(kpts, equivalences, lattvec, coeffs)

ivbm = int(data.nelect / 2) - 1
fermi = (np.max(egrid[ivbm]) + np.min(egrid[ivbm + 1])) / 2.

# Plot the results

fig1, ax1 = pl.subplots(1, figsize=(6, 4))
ax1.set_ylim([-1, 1])
ax1.set_xlim([x[0], x[-1]])
for iband in range(len(egrid)):
    ax1.plot(x, (egrid[iband] - fermi) / bandlib.eV, "k-")
for l in X:
    ax1.plot([l, l], [-1, 1], "k-")
ax1.set_xticks(X)
ax1.set_xticklabels(klist)
ax1.set_ylabel(r"$\varepsilon - \varepsilon_F$ [eV]")
fig1.tight_layout()

eband, vvband, cband = fite.getBTPbands(
    equivalences, coeffs, lattvec, curvature=True)

epsilon, dos, vvdos, cdos = bandlib.BTPDOS(
    eband, vvband, cband=cband, npts=4000)
mur = epsilon[100:-100]
TEMP = np.array([300.])
N, L0, L1, L2, Lm11 = bandlib.fermiintegrals(
    epsilon, dos, vvdos, mur=mur, Tr=TEMP, cdos=cdos)

# Use the Fermi integrals to obtain the Onsager coefficients
UCvol = data.get_volume()
sigma, seebeck, kappa, Hall = bandlib.calc_Onsager_coefficients(
    L0, L1, L2, mur, TEMP, UCvol, Lm11=Lm11)

# Plot the results
fig2, ax2 = pl.subplots(1, figsize=(4, 4))
nH = 1. / Hall[0, :, 0, 1, 2] * (
    bandlib.Coulomb * data.get_volume() / bandlib.Meter**3)
ax2.plot(nH, (mur - fermi) / bandlib.eV, "k-", label="$1/R_H$")
ax2.plot(N[0] + data.nelect, (mur - fermi) / bandlib.eV, "k--", label="n")
ax2.set_ylim([-1, 1])
ax2.set_xlim([-4, 4])
ax2.set_xlabel("[e/u.c.]")
ax2.set_ylabel(r"$\mu - \varepsilon_F\;\left(\mathrm{eV}\right)$")
ax2.legend()
fig2.tight_layout()

fig3, ax3 = pl.subplots(1, figsize=(4, 4))
ax3.plot(
    nH - (N[0] + data.nelect), (mur - fermi) / bandlib.eV,
    "k-",
    label="$1/R_H-n$")
ax3.set_ylim([-1, 1])
ax3.set_xlim([-1.5, 1.5])
ax3.set_xlabel("[e/u.c.]")
ax3.set_ylabel(r"$\mu - \varepsilon_F\;\left(\mathrm{eV}\right)$")
ax3.legend()
fig3.tight_layout()

pl.show()
