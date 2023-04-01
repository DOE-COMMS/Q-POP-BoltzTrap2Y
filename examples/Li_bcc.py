import os.path
import copy

import numpy as np
import matplotlib.pylab as plt
import ase
import ase.io

from environment import data_dir

import BoltzTraP2.dft as BTP
import BoltzTraP2.bandlib as BL
import BoltzTraP2.io as IO
from BoltzTraP2 import sphere
from BoltzTraP2 import fite
from BoltzTraP2 import serialization
from BoltzTraP2.misc import ffloat
from BoltzTraP2 import units

# Directory containing the data
data_dir1 = os.path.join(data_dir, "Li")
radix = "Li_BLZTRP"

niter = 60
# If a pregenerated bt2 file with the interpolation exists, read it. Otherwise,
# perform the interpolation and create the file.
bt2filnam = radix + ".bt2"
if os.path.isfile(bt2filnam):
    print("Loading the precalculated results from", bt2filnam)
    data, equivalences, coeffs, metadata = serialization.load_calculation(
        bt2filnam)
    print("done")
else:
    print("No pregenerated bt2 file found; performing a new interpolation")
    data = BTP.DFTData(data_dir1)
    equivalences = sphere.get_equivalences(data.atoms, data.magmom,
                                           niter * len(data.kpoints))
    print("There are", len(equivalences),
          "equivalence classes in the output grid")
    coeffs = fite.fitde3D(data, equivalences)
    serialization.save_calculation(radix + ".bt2", data, equivalences, coeffs,
                                   serialization.gen_bt2_metadata(
                                       data, data.mommat is not None))

lattvec = data.get_lattvec()
eband, vvband, cband = fite.getBTPbands(
    equivalences, coeffs, lattvec, curvature=False)

npts = 4000
Cepsilon, Cdos, Cvvdos, cdos = BL.BTPDOS(eband, vvband, npts=npts)

Tr = np.linspace(200., 600., num=17)
margin = 9. * units.BOLTZMANN * Tr.max()
mur_indices = np.logical_and(Cepsilon > Cepsilon.min() + margin,
                             Cepsilon < Cepsilon.max() - margin)
mur = Cepsilon[mur_indices]

N, L0, L1, L2, Lm11 = BL.fermiintegrals(Cepsilon, Cdos, Cvvdos, mur=mur, Tr=Tr)
Csigma, Cseebeck, kappa, Hall = BL.calc_Onsager_coefficients(
    L0, L1, L2, mur, Tr, data.get_volume())


# Interpolate the relaxation times to the denser grid using the same procedure
# as for the bands themselves.
def read_tauk(filename):
    """Read in data about electron lifetimes on the sparse grids.

    Args:
        filename: path to the file to be read

    Returns:
        An array of scattering rates.
    """
    lines = open(filename, "r", encoding="ascii").readlines()
    # line 1: title string
    # line 2: nk, nspin, Fermi level(Ry)
    linenumber = 1
    tmp = lines[linenumber].split()
    nk, nspin, efermi = int(tmp[0]), int(tmp[1]), float(tmp[2])
    minband = np.infty
    tau = []
    kpoints = []
    for ik in range(nk):
        # k block: line 1 = kx ky kz nband
        linenumber += 1
        tmp = lines[linenumber].split()
        nband = int(tmp[3])
        if nband < minband:
            minband = nband
        kpoints += [np.array(list(map(float, tmp[0:3])))]
        ttau = []
        for ib in range(nband):
            linenumber += 1
            e = ffloat(lines[linenumber].split()[0])
            ttau.append(e)
        tau.append(ttau)
    taus = np.zeros((len(kpoints), minband))
    for i in range(len(kpoints)):
        taus[i] = tau[i][:minband]
    return taus.T


tauDFT = read_tauk(os.path.join(data_dir1, radix + ".tau_k"))
pseudodata = copy.deepcopy(data)
pip = list(range(1, 60)) + list(range(61, 413))
pseudodata.ebands = tauDFT[:, pip]
pseudodata.kpoints = data.kpoints[pip]
pseudocoeffs = fite.fitde3D(pseudodata, equivalences)
tau = fite.getBTPbands(equivalences, pseudocoeffs, lattvec, curvature=False)[0]

epsilon, dos, vvdos, cdos = BL.BTPDOS(
    eband, vvband, npts=npts, scattering_model=tau)
N, L0, L1, L2, Lm11 = BL.fermiintegrals(epsilon, dos, vvdos, mur=mur, Tr=Tr)
sigma, seebeck, kappa, Hall = BL.calc_Onsager_coefficients(
    L0, L1, L2, mur, Tr, data.get_volume())

ii = np.abs(epsilon) < .05
dos = Cdos[ii]
vv = vvdos[0, 0, ii]
vv2 = Cvvdos[0, 0, ii]
coefs = np.polynomial.polynomial.polyfit(vv2 / dos, vv, [1])

Mvvdos = Cvvdos / Cdos * coefs[1]
Mvvdos[np.isnan(Mvvdos)] = 0.
N, L0, L1, L2, Lm11 = BL.fermiintegrals(Cepsilon, Cdos, Mvvdos, mur=mur, Tr=Tr)
Msigma, Mseebeck, kappa, Hall = BL.calc_Onsager_coefficients(
    L0, L1, L2, mur, Tr, data.get_volume())

ctau = np.mean(tau[0])
ii = np.nonzero(N[4] < -1)[0]  # 4 should be 300K
ifermi = ii[0]
Efermi = mur[ifermi]

fig1, ax1 = plt.subplots(1, figsize=(6, 4))
ax1.tick_params(axis='x', labelsize=16)
ax1.tick_params(axis='y', labelsize=16)

ax1.set_xlim([-2.5, 1])
ax1.set_ylim([0, 1250])
ax1.plot([0, 0], [0, 1400], "k:")

dos = BL.smoothen_DOS(epsilon, Cdos, 100.)
ax2 = fig1.add_axes([.58, .25, .3, .3])
ax2.set_xlim([-2.5, 1])
ax2.set_ylim([0, .4])
#ax3.set_yticks([0,.2, .4])
ax2.plot((epsilon - Efermi) / units.eV, dos * units.eV, "k-")
ax2.plot([0, 0], [0, .5], "k:")
ax2.set_ylabel(r"$n$ [eV$^{-1}$]", fontsize="12")

ax1.plot(
    (mur - Efermi) / units.eV, Csigma[4, :, 0, 0] * ctau * 1E-3, label="CRTA")
ax1.plot((mur - Efermi) / units.eV, sigma[4, :, 0, 0] * 1E-3, label="e-ph")
ax1.plot((mur - Efermi) / units.eV, Msigma[4, :, 0, 0] * 1E-3, label="model")
ax1.set_ylabel(r"$\sigma\;\left[\mathrm{kS\,m^{-1}}\right]$", fontsize="16")
ax1.set_xlabel("$\mu$ [eV]", fontsize="16")
ax1.legend(fontsize=16)
fig1.tight_layout(pad=1.)
fig1.savefig("Li_sigma.pdf")
plt.show()

print("CRTA", Cseebeck[4, ifermi, 0, 0] * 1E6)
print("el-ph", seebeck[4, ifermi, 0, 0] * 1E6)
print("Model", Mseebeck[4, ifermi, 0, 0] * 1E6)
