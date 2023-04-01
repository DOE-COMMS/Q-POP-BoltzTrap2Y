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

###############################################################################
# Load Bi2Te3 results from Wien2k, interpolate the bands, compute the Onsager
# transport coefficients and reproduce Fig. 1 from Comp. Phys. Comm. 175
# (2006) 67–71
###############################################################################

import os
import os.path
import pickle
import logging

import numpy as np
import matplotlib.pylab as pl

from environment import data_dir

from BoltzTraP2 import dft as BTP
from BoltzTraP2 import sphere
from BoltzTraP2 import fite
from BoltzTraP2 import bandlib as BL
from BoltzTraP2 import serialization
from BoltzTraP2 import units

logging.basicConfig(
    level=logging.DEBUG, format="{levelname:8s}│ {message:s}", style="{")

dirname = os.path.join(data_dir, "Bi2Te3")
bt2file = "Bi2Te3.bt2"

# If a ready-made file with the interpolation results is available, use it
# Otherwise, create the file.
if not os.path.exists(bt2file):
    # Load the input
    data = BTP.DFTData(os.path.join(data_dir, dirname))
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

# Load the interpolation results
data, equivalences, coeffs, metadata = serialization.load_calculation(bt2file)

# Reconstruct the bands
lattvec = data.get_lattvec()
eband, vvband, cband = fite.getBTPbands(equivalences, coeffs, lattvec)

# Obtain the Fermi integrals for different chemical potentials at
# room temperature.
TEMP = np.array([300.])
epsilon, dos, vvdos, cdos = BL.BTPDOS(eband, vvband, npts=4000)
margin = 9. * units.BOLTZMANN * TEMP.max()
mur_indices = np.logical_and(epsilon > epsilon.min() + margin,
                             epsilon < epsilon.max() - margin)
mur = epsilon[mur_indices]
N, L0, L1, L2, Lm11 = BL.fermiintegrals(
    epsilon, dos, vvdos, mur=mur, Tr=TEMP, dosweight=data.dosweight)

# Compute the Onsager coefficients from those Fermi integrals
UCvol = data.get_volume()
sigma, seebeck, kappa, Hall = BL.calc_Onsager_coefficients(
    L0, L1, L2, mur, TEMP, UCvol)

fermi = BL.solve_for_mu(epsilon, dos, data.nelect, 300, data.dosweight)

# Plot the results
fig1, ax1 = pl.subplots(1, figsize=(6, 3))
ax1.set_xlim([-1, 1])
ax1.set_ylim([-300, 300])
ax1.plot((mur - fermi) / BL.eV, seebeck[0, :, 0, 0] * 1e6, "k-", label="xx")
ax1.plot((mur - fermi) / BL.eV, seebeck[0, :, 2, 2] * 1e6, "k--", label="zz")
ax1.set_xlabel("$\mu$ [eV]")
ax1.set_ylabel("$S$ [$\mu$V/K]")
ax1.legend()
fig1.tight_layout(pad=1.)

fig2, ax2 = pl.subplots(1, figsize=(6, 3))
ax2.set_xlim([-1, 1])
ax2.set_ylim([0, 70])
ax2.plot(
    (mur - fermi) / BL.eV,
    seebeck[0, :, 0, 0]**2 * sigma[0, :, 0, 0] * 1e-10,
    "k-",
    label="xx")
ax2.plot(
    (mur - fermi) / BL.eV,
    seebeck[0, :, 2, 2]**2 * sigma[0, :, 2, 2] * 1e-10,
    "k--",
    label="zz")
ax2.set_xlabel("$\mu$ [eV]")
ax2.set_ylabel(r"$S^2\sigma/ \tau$ [10$^{14}$ $\mu$W/(cm K$^2$) s]")
ax2.legend()
fig2.tight_layout(pad=1.)
pl.show()
