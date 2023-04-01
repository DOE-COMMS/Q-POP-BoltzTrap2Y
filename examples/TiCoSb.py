#!/usr/bin/env python
###############################################################################
# Load Wien2k results for TiCoSb, perform the interpolation, create a band
# diagram and calculate the power factor.
# Reproduce Fig. 16 from J. Chem. Phys. submitted
###############################################################################

import os
import os.path
import logging

import numpy as np
import matplotlib
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
from BoltzTraP2 import fermisurface

#logging.basicConfig(
#    level=logging.DEBUG, format="{levelname:8s}│ {message:s}", style="{")

ninter = 5
dirname = os.path.join(data_dir, "TiCoSb")
bt2file = "TiCoSb_" + str(ninter) + ".bt2"

# If a ready-made file with the interpolation results is available, use it
# Otherwise, create the file.
if not os.path.exists(bt2file):
    data = dft.DFTData(dirname)
    equivalences = sphere.get_equivalences(data.atoms, data.magmom,
                                           len(data.kpoints) * ninter)
    data.bandana(emin=data.fermi - .2, emax=data.fermi + .2)
    coeffs = fite.fitde3D(data, equivalences)
    serialization.save_calculation(
        bt2file, data, equivalences, coeffs,
        serialization.gen_bt2_metadata(data, data.mommat is not None))

data, equivalences, coeffs, metadata = serialization.load_calculation(bt2file)

# Plot band structure
cell = data.get_lattvec()
points = get_special_points(cell, "fcc")
klist = "KLGXWKG"
PP = [points[k] for k in klist]

kpts, x, X = bandpath(PP, cell, 100)
#path = bandpath(PP, cell, 200)
#(x, X, labels) = path.get_linear_kpoint_axis()

# Rebuild the bands from the interpolation coefficients
lattvec = data.get_lattvec()
egrid, vgrid = fite.getBands(kpts, equivalences, lattvec, coeffs)
#egrid, vgrid = fite.getBands(path.kpts, equivalences, lattvec, coeffs)

ivbm = int(data.nelect / 2) - 1
fermi = np.max(egrid[ivbm])

# Plot the results
fig1, (ax1, ax2) = pl.subplots(1,
                               2,
                               figsize=(7, 4),
                               sharey=True,
                               gridspec_kw={'width_ratios': [3, 1]})

try:
    image = pl.imread("FS3.png")
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    axins = inset_axes(ax2, width=1.5, height=1.5, loc=4)
    axins.get_xaxis().set_visible(False)
    axins.get_yaxis().set_visible(False)
    im = axins.imshow(image)
except:
    pass

ymin = -1.
ymax = 0.2
ax1.set_ylim([ymin, ymax])
ax1.set_xlim([x[0], x[-1]])
color_cycle = matplotlib.rcParams["axes.prop_cycle"].by_key()["color"]
color_cycle = list(np.array(color_cycle)[[5, 6, 7, 3, 4, 0, 1, 2, 8, 9]])
for iband in range(ivbm - 1, ivbm - 5, -1):
    ax1.plot(x, (egrid[iband] - fermi) / bandlib.eV,
             color=color_cycle[iband % len(color_cycle)])
iband = ivbm
ax1.plot(x, (egrid[iband] - fermi) / bandlib.eV,
         color=color_cycle[iband % len(color_cycle)])
for l in X:
    ax1.plot([l, l], [ymin, ymax], "k-")
ax1.plot([x[0], x[-1]], [0, 0], "k--")
ax1.set_xticks(X)
klabel = []
for s in klist:
    if (s == "G"):
        klabel += ['$\Gamma$']
    else:
        klabel += ["$" + s + "$"]

ax1.tick_params(axis='x', labelsize=14)
ax1.tick_params(axis='y', labelsize=14)
ax1.set_xticklabels(klabel)
ax1.set_ylabel(r"$\varepsilon - \varepsilon_F$ [eV]", fontsize=14)

eband, vvband, cband = fite.getBTPbands(equivalences,
                                        coeffs,
                                        lattvec,
                                        curvature=False)

epsilon, dos, vvdos, cdos = bandlib.BTPDOS(eband, vvband, npts=4000)

mur_indices = np.logical_and(epsilon > fermi - 2 * bandlib.eV,
                             epsilon < fermi + 1 * bandlib.eV)
mur = epsilon[mur_indices]
TEMP = np.array([300.])
N, L0, L1, L2, Lm11 = bandlib.fermiintegrals(epsilon,
                                             dos,
                                             vvdos,
                                             mur=mur,
                                             Tr=TEMP,
                                             cdos=cdos)

# Use the Fermi integrals to obtain the Onsager coefficients
UCvol = data.get_volume()
sigma, seebeck, kappa, Hall = bandlib.calc_Onsager_coefficients(L0,
                                                                L1,
                                                                L2,
                                                                mur,
                                                                TEMP,
                                                                UCvol,
                                                                Lm11=Lm11)

PF = seebeck[0, :, 0, 0]**2 * sigma[0, :, 0, 0] * 1E-14 * 1E3
ax2.plot(PF, (mur - fermi) / bandlib.eV, "k-")

imax = np.argmax(PF)
l1 = (mur[imax] - fermi) / bandlib.eV
ax2.set_xlim([0, 8])
ax2.plot([0, 8], [0, 0], "k--")
ax1.plot([x[0], x[-1]], [l1, l1], "k:", linewidth=1)
ax2.plot([0, 8], [l1, l1], "k:", linewidth=1)
ax2.set_xticks([0, 4, 8])
ax2.set_ylim([ymin, ymax])
ax2.tick_params(axis='x', labelsize=14)
ax2.set_xlabel("$S^2\sigma$ [mWm$^{-1}$K$^{-2}$]", fontsize=14)
fig1.tight_layout()
#fig1.savefig("TiCoSb_bands.pdf")
pl.show()

fermisurface.plot_fermisurface(data,
                               equivalences,
                               eband,
                               mur[imax] - fermi,
                               color_cycle=color_cycle)

