#    plot result in pngfile, coded by Yi Wang 01/29/2020

import os
import logging
import numpy as np
import scipy as sp
import scipy.constants
from BoltzTraP2.units import *
from BoltzTraP2.misc import lexit
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def plotpng(filename, data, metadata="uniform_tau"):
    tt = np.loadtxt(filename, comments="#", usecols=(1,), dtype=np.float)
    tt = sorted(set(tt))
    nT = len(tt)
    mm = np.loadtxt(filename, comments="#", usecols=(0,), dtype=np.float)
    if mm[0] != mm[nT-1]:
        plotfunc (filename, data,  metadata=metadata)
    else:
        plotfunc (filename, data,  abscissa="mu", metadata=metadata, subsample=max(1,nT//2), XX=[-1.2,1.2])
        plotfunc (filename, data,  abscissa="dope", metadata=metadata, subsample=max(1,nT//2), XX=[-1.2e21,1.2e21])


def plotfunc(filename, data, abscissa="T", metadata="uniform_tau", subsample=1, XX=None):
    parse_plot(filename, "ccT", data,  abscissa=abscissa, metadata=metadata, subsample=subsample, XX=XX)
    parse_plot(filename, "cv", data,  abscissa=abscissa, metadata=metadata, subsample=subsample, XX=XX)
    parse_plot(filename, "S", data, abscissa=abscissa, metadata=metadata, subsample=subsample, XX=XX)
    parse_plot(filename, "eT", data, abscissa=abscissa, metadata=metadata, subsample=subsample, XX=XX)
    parse_plot(filename, "ePF", data, abscissa=abscissa, metadata=metadata, subsample=subsample, XX=XX)
    parse_plot(filename, "n", data, abscissa=abscissa, metadata=metadata, subsample=subsample, XX=XX)
    parse_plot(filename, "sigma", data, abscissa=abscissa, metadata=metadata, subsample=subsample, XX=XX)
    parse_plot(filename, "RH", data, abscissa=abscissa, metadata=metadata, subsample=subsample, XX=XX)
    parse_plot(filename, "kappae", data, abscissa=abscissa, metadata=metadata, subsample=subsample, XX=XX)
    parse_plot(filename, "PF", data, abscissa=abscissa, metadata=metadata, subsample=subsample, XX=XX)
    parse_plot(filename, "L", data, abscissa=abscissa, metadata=metadata, subsample=subsample, XX=XX)
    parse_plot(filename, "M", data, abscissa=abscissa, metadata=metadata, subsample=subsample, XX=XX)


def parse_plot(filename, quantity, data, abscissa = "T", metadata="uniform_tau", subsample=1, XX=None):
    """fast plot for scalar quantities"""
    # Try and load the data from the integration step
    # Prepare the abscissas, the third variable and the axis labels
    #colors=('b', 'g', 'r', 'c', 'm', 'y', 'k', 'w')
    colors=('k', 'b', 'r', 'c', 'm', 'g', 'y')
    linestyles = ['-', '--', ':', '-.', '.', ',', 'o', '^', 'v', '<', '>', 's',
              '+', 'x', 'd', '1', '2', '3', '4', 'h', 'p', '|', '_', 'D', 'H']
    n_lines = len(linestyles)
    n_colors = len(colors)

    tt = np.loadtxt(filename, comments="#", usecols=(1,), dtype=np.float)
    nmu = len(tt)
    tt = sorted(set(tt))
    nT = len(tt)
    nmu = nmu//nT

    vuccm3 = data.atoms.get_volume() * 1e-24  # in cm^3

    if abscissa == "T":
        x = tt
        xlabel = r"$T\;\left[\mathrm{K}\right]$"
        dd = np.loadtxt(filename, comments="#", usecols=(2,), dtype=np.float)/vuccm3
        dd = dd[nT//2::nT]
        z = dd[::subsample]
        z = np.where(abs(z) < 1.e15, 0., z)
        zlabel0 = r"$dope = {:.2g}\;\mathrm{{cm^{{-3}}}}$"
    elif abscissa == "dope":
        dd = np.loadtxt(filename, comments="#", usecols=(2,), dtype=np.float)/vuccm3
        dd = np.array(list(map(float,dd))).reshape(nmu,len(dd)//nmu).T
        x = dd[::subsample,:]
        xlabel = "dope"+r"$\left(\mathrm{cm^{-3}}\right)$"
        z = tt[::subsample]
        zlabel0 = r"$T = {:5g}\;\mathrm{{K}}$"
    else:
        # for abscissa == "mu" use
        mm = np.loadtxt(filename, comments="#", usecols=(0,), dtype=np.float)
        mm = np.array(list(map(float,mm))).reshape(nmu,len(mm)//nmu).T
        x = mm[::subsample,:]
        xlabel = r"$\mu - \varepsilon_F\;\left[\mathrm{eV}\right]$"
        z = tt[::subsample]
        zlabel0 = r"$T = {:5g}\;\mathrm{{K}}$"

    ylabel = dict(
        cv=r"$c_v\;\left[\mathrm{J\,mol.atom^{-1}\,K^{-1}}\right]$",
        eT="eT="+r"$\left[(c(\mu)-c(x)\right]/c(\mu)$",
        ePF="ePF="+r"$c(\mu)-c(x)\left[\mathrm{J\,mol.atom^{-1}\,K^{-1}}\right]$",
        ccT=r"$n\;\left[\mathrm{\left|e\right|\,cm^{-3}}\right]$",
        n=r"$n\;\left[\mathrm{\left|e\right|\,cm^{-3}}\right]$",
        DOS=r"$\mathrm{DOS}\;\left[\mathrm{uc^{-1}}\right]$",
        sigma=r"$\sigma^{{\left({}\right)}}/\tau_0\;"
        r"\left[\mathrm{{\Omega^{{-1}}\,m^{{-1}}\,"
        r"s^{{-1}}}}\right]$",
        S=r"$S\;" +
        r"\left[\mu\mathrm{{V\,K^{{-1}}}}\right]$",
        kappae=r"$\kappa_e^{{\left({}\right)}}/\tau_0\;"
        r"\left[\mathrm{{W\,m^{{-1}}\,K^{{-1}}\,s^{{-1}}}}\right]$",
        L=r"$L^{{\left({}\right)}}\;"
        r"\left[\mathrm{{W\,\Omega\,K^{{-2}}}}\right]$",
        PF=r"$\left(S^2\sigma\right)^{{\left({}\right)}}/ \tau_0"
        r"\;\left[\mathrm{{W\,m^{{-1}}\,K^{{-2}}\,s^{{-1}}}}\right]$",
        #M=r"$\sigma/n\,(\tau_0/e)\," +
        M=r"$M\,(\tau_0/e)\," +
        r"[cm^2\,/(V\,s)]$",
        RH=r"$R_H^{{\left({}\right)}}\;" +
        r"\left[\mathrm{{m^3\,C^{{-1}}}}\right]$")
    if metadata == "uniform_lambda":
        ylabel["sigma"] = (r"$\sigma^{{\left({}\right)}}/\lambda_0\;"
                           r"\left[\mathrm{{\Omega^{{-1}}\,m^{{-2}}}}\right]$")
        ylabel["kappae"] = (r"$\kappa_e^{{\left({}\right)}}/\lambda_0\;"
                            r"\left[\mathrm{{W\,m^{{-2}}\,K^{{-1}}}}"
                            r"\right]$")
        ylabel["PF"] = (
            r"$\left(S^2\sigma\right)^{{\left({}\right)}}/ \lambda_0"
            r"\;\left[\mathrm{{W\,m^{{-2}}\,K^{{-2}}}}\right]$")

    ylabel0 = ylabel[quantity]
    # Prepare the ordinate
    y1 = None
    if quantity == "cv":
        y = np.loadtxt(filename, comments="#", usecols=(8,), dtype=np.float)
        y1 = np.loadtxt(filename, comments="#", usecols=(10,), dtype=np.float)
    elif quantity == "eT":
        y = np.loadtxt(filename, comments="#", usecols=(8,), dtype=np.float)
        tmp = np.loadtxt(filename, comments="#", usecols=(10,), dtype=np.float)
        y[y<1.e-16] = 1.e-16
        y = (y-tmp)/y
    elif quantity == "ePF":
        y = np.loadtxt(filename, comments="#", usecols=(8,), dtype=np.float)
        tmp = np.loadtxt(filename, comments="#", usecols=(10,), dtype=np.float)
        y = y-tmp
    elif quantity == "ccT":
        y = np.loadtxt(filename, comments="#", usecols=(12,), dtype=np.float)
    elif quantity == "n":
        #y = np.loadtxt(filename, comments="#", usecols=(2,), dtype=np.float)/vuccm3
        y = np.loadtxt(filename, comments="#", usecols=(12,), dtype=np.float)
    elif quantity == "DOS":
        y = np.loadtxt(filename, comments="#", usecols=(3,), dtype=np.float)
    elif quantity == "S":
        y = np.loadtxt(filename, comments="#", usecols=(4,), dtype=np.float)*1e6
        y1 = np.loadtxt(filename, comments="#", usecols=(11,), dtype=np.float)*1e6
    elif quantity == "sigma":
        y = np.loadtxt(filename, comments="#", usecols=(5,), dtype=np.float)
    elif quantity == "RH":
        y = np.loadtxt(filename, comments="#", usecols=(6,), dtype=np.float)
    elif quantity == "kappae":
        y = np.loadtxt(filename, comments="#", usecols=(7,), dtype=np.float)
    elif quantity == "L":
        y = np.loadtxt(filename, comments="#", usecols=(13,), dtype=np.float)
        t0 = np.loadtxt(filename, comments="#", usecols=(1,), dtype=np.float)
        t1 = np.loadtxt(filename, comments="#", usecols=(5,), dtype=np.float)
        t2 = np.loadtxt(filename, comments="#", usecols=(7,), dtype=np.float)
        y1 = t2/t1/t0
    elif quantity == "M":
        t0 = np.loadtxt(filename, comments="#", usecols=(12,), dtype=np.float)
        t1 = np.loadtxt(filename, comments="#", usecols=(5,), dtype=np.float)
        y = t1/t0*100
    elif quantity == "PF":
        s0 = np.loadtxt(filename, comments="#", usecols=(4,), dtype=np.float)
        s1 = np.loadtxt(filename, comments="#", usecols=(11,), dtype=np.float)
        cond = np.loadtxt(filename, comments="#", usecols=(5,), dtype=np.float)
        y = cond*s0*s0
        y1 = cond*s1*s1

    y = np.array(list(map(float,y))).reshape(nmu,len(y)//nmu).T
    if y1 is not None:
        y1 = np.array(list(map(float,y1))).reshape(nmu,len(y1)//nmu).T

    if abscissa == "T":
        y = y[:, ::subsample]
        if y1 is not None:
            y1 = y1[:, ::subsample]
    else:
        y = y[::subsample, :]
        if y1 is not None:
            y1 = y1[::subsample, :]
        if XX is None:
            XX = [-2.0, 2.0]
        ni = []
        for ix in range(len(x[0,:])):
            if x[:,ix].max() < XX[0] : continue
            if x[:,ix].min() > XX[1] : continue
            ni.append(ix)
        x = x[:,ni]
        y = y[:,ni]
        if y1 is not None:
            y1 = y1[:,ni]

    # Create the plots
    logging.getLogger('matplotlib.font_manager').disabled = True
    pngsizefactor = 2
    plt.rc('font', size=48//pngsizefactor)
    matplotlib.rcParams['lines.linewidth'] = 2
    for c in range(1):
        ylabel = ylabel0
        fig,ax=plt.subplots()
        #fig.set_size_inches(36//pngsizefactor,24//pngsizefactor)
        fig.set_size_inches(24//pngsizefactor,18//pngsizefactor)

        if quantity == "n" or quantity == "ccT":
            #ax.set_yscale('symlog',linthreshy=1e17)
            ax.set_yscale('log')

        if abscissa == "dope":
            ax.set_xscale('symlog',linthreshx=1e17)

        if abscissa == "mu" or abscissa == "dope":
            plt.axvline(x=0., color="#333333", lw=2)

        #plt.figure()
        """
        if quantity == "cv":
            ax.set_ylim(0,2.5*(y-y1).max())
        """

        for iz, zv in enumerate(z):
            zlabel = zlabel0.format(zv)
            if abscissa == "T":
                px = x
                py = y[:,iz]
                if y1 is not None: py1 = y1[:,iz]
            else:
                px = x[iz,:]
                py = y[iz,:]
                if y1 is not None: py1 = y1[iz,:]

            color = iz%n_colors
            if quantity == "cv":
                plt.plot(px, py, linestyle=linestyles[0], color=colors[color], label=zlabel+r",$c(\mu)$")
                plt.plot(px, py1, linestyle=linestyles[1], color=colors[color], label=zlabel+r",$c(x)$")
                #plt.plot(px, py-py1, linestyle=linestyles[2], color=colors[color], label=None)
                #plt.scatter(px[::5], (py-py1)[::5], s=200, marker='o', facecolors='none', color=colors[color], label=zlabel+r",$c(\mu)-c(x)$")
                #plt.scatter(px[::5], (py-py1)[::5], s=200, marker='o', color=colors[color], label=zlabel+r",$c(\mu)-c(x)$")
            elif quantity == "S":
                plt.plot(px, py, linestyle=linestyles[1], color=colors[color], label=zlabel+r",$S_{e}(BTE)$")
                plt.plot(px, py1, linestyle=linestyles[0], color=colors[color], label=zlabel+r",$S_{e}(Th)$")
            elif quantity == "L":
                plt.plot(px, py1, linestyle=linestyles[1], color=colors[color], label=zlabel+r",$S_{e}(BTE)$")
                plt.plot(px, py, linestyle=linestyles[0], color=colors[color], label=zlabel+r",$S_{e}(Th)$")
            elif quantity == "PF":
                plt.plot(px, py, linestyle=linestyles[1], color=colors[color], label=zlabel+r",$PF\left[S_{e}(BTE)\right]$")
                plt.plot(px, py1, linestyle=linestyles[0], color=colors[color], label=zlabel+r",$PF\left[S_{e}(Th)\right]$")
            else:
                plt.plot(px, py, linestyle=linestyles[0], color=colors[color], label=zlabel)

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        legend = plt.legend(loc="best")
        try:
            legend.set_draggable(True)
        except AttributeError:
            pass

        figurefolder = "figures/"
        if not os.path.exists(figurefolder): os.mkdir(figurefolder)
        pfile = figurefolder+abscissa+'-'+quantity + ".png"
        fig.savefig(pfile,bbox_inches='tight')
        plt.close(fig)
