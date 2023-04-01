# -------- energy in eV, temperature in K
# assume every variable starting with a-h  and o-z are real numbers
# common block named comcon
from __future__ import division
import sys
import argparse
import math
import numpy as np
from scipy.constants import physical_constants
from scipy.optimize import brentq, curve_fit
from scipy.integrate import cumtrapz, trapz, simps
from scipy.interpolate import interp1d, splev, splrep, BSpline
from scipy.integrate import quadrature
from scipy.interpolate import UnivariateSpline
from scipy.constants import physical_constants

k_B = physical_constants['Boltzmann constant in eV/K'][0]
Faraday_constant = 96485.3321233100184

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def plot():

    data = np.genfromtxt('thermo.out')
    #zero coupon maturity dates
    x = data[:,0]
    #tenor
    y = data[:,1]
    #rates
    z = 1000.*data[:,2]

    # plotting
    #ax.plot3D(x, y, z, 'green')
    #ax.plot_trisurf(x, y, z, color='white', edgecolors='grey', alpha=0.5)
    #ax.plot_surface(x, y, z, rstride = 2, cstride = 1, cmap = plt.cm.Blues_r)
    xi = np.linspace(min(x), max(x))
    yi = np.linspace(min(y), max(y))
    xi,yi = np.meshgrid(xi,yi)
    zi = griddata((x,y),z,(xi,yi),method='cubic')

    plt.rcParams.update({'font.size': 24})
    plt.rcParams['axes.labelpad'] = 20
    #fig = plt.figure(figsize=(24,24))
    fig = plt.figure()
    fig.set_size_inches(18,18)

    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(xi, yi, zi, rstride=6, cstride=6, cmap=cm.jet,linewidth=3)

    ax.set_title('Chemical potential of electrons')
    ax.set_xlabel('phase fraction of f1')
    ax.set_ylabel('T (K)')
    ax.set_zlabel('$\mu$ (meV)')
    plt.savefig('thermo.png')
    plt.show()

    #S_el, C_el, C_mu, L, Q_el, Q_p, Q_e
    """
    plotcom(1,3,data,com,title="S_el")
    plotcom(1,4,data,com,title="C_el")
    plotcom(1,5,data,com,title="C_mu")
    plotcom(1,6,data,com,title="L")
    plotcom(1,7,data,com,title="Q_el")
    plotcom(1,8,data,com,title="Q_p")
    plotcom(1,9,data,com,title="Q_e")
    """

def plotcom(ix,iy,_data,com,title="xx",ylabel="",filename=""):
    data = _data[_data[:,0]==com]

    x = data[:,ix]
    y = data[:,iy]

    fig,ax = plt.subplots()
    fig.set_size_inches(18,12)
    surf = ax.plot(x, y)

    ax.set_title(title)
    ax.set_xlabel('T (K)')
    if ylabel!="": ax.set_ylabel(ylabel)
    else: ax.set_ylabel(title)
    if filename!="": plt.savefig('{}-{:.3f}.png'.format(filename,com))
    else: plt.savefig('{}-{:.3f}.png'.format(title,com))
    plt.show()

def getdos(xdn, xup, dope, NEDOS, gaussian, edn, eup, vde, ve, tdos): # Line 186
    """

    Parameters
    ----------
    dos : pymatgen.electronic_structure.dos.Dos
        DOS object from pymatgen
    xdn : float
        Minimum energy for integration
    xup : float
        Maximum energy for integration
    dope : float
        Number of electrons to dope (negative means to remove electrons, positive means add electrons)
    dos_grid_size : int
        Number of DOS points have in the energy/density grid
    gaussian_grid_size : int
        Number of Gaussian points to use in the grid mesh around the Fermi energy

    Returns
    -------
    tuple
        Tuple of a (float, float, array, array) of the number of electrons,
        Fermi level shift due to doping, and the arrays of energies and densities on a grid.
    """

    for i,energy in enumerate(ve):
      if energy <-15.0 and tdos[i]==0.0:
        xdn = energy

    n_dos = len(tdos)
    idx = closest(ve,0.0)
    for i in range(idx,n_dos):
        if tdos[i]!=0.0:
            iBoF = i-1
            break

    eBoF = -1.0
    if iBoF>=idx:
      eBoF = ve[iBoF]
      espr = tdos[iBoF+2]-tdos[iBoF+1]
      if espr>0.0:
        espr = tdos[iBoF+1]/espr*vde
        if (espr < vde):
          eBoF = ve[iBoF+1] - espr

    #print("eBoF=", eBoF)
    xdn = max(xdn,edn)
    xup = min(xup,eup)

    e = np.linspace(xdn,xup,NEDOS,dtype=float)
    if gaussian != 0.0:
      e = remesh(xdn, xup, gaussian, 0.0, eBoF, NEDOS)

    dos = refdos(eBoF, 0.0, vde, edn, e, ve, tdos)
    ados = cumtrapz(dos, e, initial=0.0)
    idx = closest(e,0.0)
    for idx1 in range(idx-1, 0, -1):
        if ados[idx1] != ados[idx] : break
    NELECTRONS = ados[idx] - e[idx]/(e[idx1]-e[idx])*(ados[idx1]-ados[idx])

    dF = 0.0
    if dope != 0.0:
        NELECTRONS = NELECTRONS+dope
        idx = closest(ados,NELECTRONS)
        for idx1 in range(idx-1, 0, -1):
            if ados[idx1] != ados[idx] : break
        #if idx == (NEDOS-1) or ados[idx] == ados[idx+1]:
        #print ("NELECTRONS=", NELECTRONS, "idx=", idx, ados[idx], "idx1=", idx1, ados[idx1], "NEDOS=", NEDOS)
        if idx1 <= 0 or idx >= (NEDOS-1) or ados[idx] == ados[idx1]:
            print ("NELECTRONS=", NELECTRONS, "idx=", idx, ados[idx], "idx1=", idx1, ados[idx1], "NEDOS=", NEDOS)
            # we are dopidxng too much
            raise ValueError('Too much doping')
        dF = (NELECTRONS-ados[idx])/(ados[idx1] - ados[idx])*(e[idx1] - e[idx])+e[idx]
                # dF is the shift in the Fermi energy due to doping
        e = e - dF # This is done in a loop (line 289), but I think we can do without

    if gaussian != 0.0 and abs(dope)>0.0001: # why did I do this ***********************
    #if gaussian != 0.0:
      e = remesh(xdn, xup, gaussian, dF, eBoF, NEDOS)

    dos = refdos(eBoF, dF, vde, edn, e, ve, tdos)
    edos = e*dos
    ados = cumtrapz(dos, e, initial=0.0)
    energy = cumtrapz(edos, e, initial=0.0)
    idx = closest(e,0.0)
    NELECTRONS = ados[idx] - e[idx]/(e[idx+1]-e[idx])*(ados[idx+1]-ados[idx])
    E0 = energy[idx] - e[idx]/(e[idx+1]-e[idx])*(energy[idx+1]-energy[idx])

    return NELECTRONS, E0, dF, e, dos, eBoF

def remesh(xdn, xup, gaussian, dF, eBoF, NEDOS):
    """
    refine the dos mesh by using denser mesh around the 0 K Fermi energy in order to decrease the numerical uncertainty
    Parameters
    ----------
    eBoF : Conduction band minimum
    dF : Fermi energy change due to doping
    ve : original e mesh
    NEDOS : original e mesh
    gaussian : parameter used to refine the e mesh near the Fermi energy

    Return
    ------
    e : refined e mesh
    """

    e = np.zeros(NEDOS)
    e[0] = xdn - dF
    xde = 2.0*(xup - xdn)/(NEDOS-1)
    if eBoF>0.0:
        xde = 3.0*(xup - xdn)/(NEDOS-1)
    sigma = -0.5*(gaussian/(xup-xdn))**2
    fac = gaussian/(math.sqrt(2.0*math.pi))
    for i in range(1,NEDOS):
        f1 = 1.0 + fac*math.exp(sigma*(e[i-1])**2)
        if eBoF>0.0:
          if dF < eBoF:
             f1 += fac*math.exp(sigma*((e[i-1]-eBoF+dF))**2)
          else:
             f1 += fac*math.exp(sigma*((e[i-1]+dF))**2)
        e[i] = e[i-1]+xde/f1
    return e

def refdos(eBoF, dF, vde, edn, e, ve, tdos):
    """
    refine the dos mesh by using denser mesh around the 0 K Fermi energy in order to decrease the numerical uncertainty
    Parameter
    ---------
    eBoF : Conduction band minimum
    dF : Fermi energy change due to doping
    e : refined e mesh
    ve : original e mesh
    tdos : original e dos
    Return
    ------
    dos : refined e dos
    """

    dos = np.zeros(len(e))
    n_dos = len(tdos)
    for i in range(0, len(e)):
        tx = e[i] + dF
        kx = int((tx-edn)/vde) # Converted to int, remember the type!
        kx = max([kx,0]) # XXX: is this translated correctly? What is the 1 in fortran?
        kx = min([n_dos-2, kx]) # TODO: the ndos-1 was here before. could be a source of error
        if tdos[kx+1]==0.0 and ve[kx+1]>0.0 and ve[kx+1]<vde:
          # handling near the Top of valence band
          if tx >= 0.0:
            dos[i] = 0.0
          else:
            dos[i] = tdos[kx]*tx/ve[kx]
            #dos[i] = tdos[kx]*(tx/ve[kx])**2
        elif eBoF > 0.0 and tdos[kx]==0.0 and ve[kx+1]-eBoF<vde and ve[kx+1]-eBoF>0.0:
          # handling near the bottom of conduction band
            if tx <= eBoF:
              dos[i] = 0.0
            else:
              dos[i] = tdos[kx+1]*(tx-eBoF)/(ve[kx+1]-eBoF)
        else:
          dos[i] = tdos[kx] + (tdos[kx+1] - tdos[kx])/vde*(tx - ve[kx])
    return dos

def closest(e,val):
    """
    find the index of the band energy which is the close to the energy val
    Parameters
    ----------
    e : float
    array of band energy for the e dos
    val : given value of band energy

    Return
    ------
    index of e that closest to the energy val
    """
    idx = np.abs(e-val).argmin()
    if e[idx] < val:
        idx = idx + 1
    return idx

def gfind(mu_el, pe, pdos, NELECTRONS, Beta, IntegrationFunc=trapz):
    """
    Calculate the number of electron difference from 0K given chemical potential. the purpose is the find the
    chemical potential to make zero of number of electron difference from 0K

    Parameters
    ----------
    mu_el : chemical potential, :math:`\mu`, in the Fermi distribution
    pe : eigenenergies
    pdos : density of states (:math:`n(\varepsilon) \varepsilon`
    NELECTRONS : Total number of electrons in the system at 0K
    Beta : :math:`\frac{1}{T*k_{B}}`

    Returns
    -------
    The number of electron difference from 0K given chemical potential
    """

    tc = Beta*(pe-mu_el)
    tc = tc[np.where(tc<200)]
    k = len(tc)
    fn = pdos[0:k]/(np.exp(tc[0:k])+1.0)
    return IntegrationFunc(fn, pe[0:k])- NELECTRONS


# line 363
def caclf(pe, pdos, NELECTRONS, Beta, mu_ref=0.0, dF=0.0, IntegrationFunc=trapz): #line 363
    """
    Calculate thermal free energy from electronic density of states (e DOS)

    Parameters
    ----------
    pe : band energy array
    pdos : e DOS
    NELECTRONS : total number of electrons
    Beta : 1/(kB*T)

    Returns
    -------
    electron chememical potential, internal energy, entropy, carrier amount, coefficient to cal Seebeck
    """

    #print ("dF=", dF)
    if 1==1:
        mu_el = brentq(gfind, mu_ref-10.0, mu_ref+10.0, args=(pe, pdos, NELECTRONS, Beta, IntegrationFunc), maxiter=10000)
    else:
        t0 = mu_ref
        d0 = gfind(t0, pe, pdos, NELECTRONS, Beta, IntegrationFunc)
        if d0 > 0.0: td = -0.1
        elif d0 <0.0: td = 0.1
        else: return t0
        for i in range(999):
            t1 = t0 + td
            d1 = gfind(t1, pe, pdos, NELECTRONS, Beta, IntegrationFunc)
            if d1*d0 < 0.0: break
            elif d1*d0 == 0.0: break
            t0 = t1
            d0 = d1
            td = td + td
        for i in range(999):
            t2 = (t0 + t1)*0.5
            d2 = gfind(t2, pe, pdos, NELECTRONS, Beta, IntegrationFunc)
            if d2*d0 < 0.0:
                t1 = t2
                d1 = d2
            else:
                t0 = t2
                d0 = d2
            if abs(t1-t0) <1.e-8:
                mu_el = 0.5*(t0+t1)
                break
        #mu_el = brentq(gfind, t0, t1, args=(pe, pdos, NELECTRONS, Beta), maxiter=10000)
        print("xxxxxxx", mu_el,mu_old)
    tc = Beta*(pe-mu_el)
    tc = tc[np.where(tc<200)]
    k1 = len(tc)
    tf = 1.0/(np.exp(tc)+1.0)
    fn = pdos[0:k1]*pe[0:k1]*tf
    u = IntegrationFunc(fn, pe[0:k1])

    k0 = closest(tc,-200)
    tf0 = tf[k0:]
    pdos = pdos[k0:k1]
    pe = pe[k0:k1]
    tf1 = 1.0 - tf0 + 1.e-60 # 1.e-60 is used to avoid log exception
    fn = pdos*(tf0*np.log(tf0)+tf1*np.log(tf1))
    s = IntegrationFunc(fn, pe)

    tf = tf0*(1.0-tf0)
    fn = pdos*tf
    fn2 = pdos*tf*(pe-mu_el)
    Q_el = IntegrationFunc(fn, pe)
    Q_p = IntegrationFunc(fn[pe<=dF], pe[pe<=dF])
    Q_e = IntegrationFunc(fn[pe>dF], pe[pe>dF])
    Y_el = IntegrationFunc(fn2, pe)

    fn = pdos*(pe-mu_el)*tf
    if Q_el!=0.0:
        e_ = IntegrationFunc(fn, pe)/Q_el
        fn = pdos[0:k1]*(pe[0:k1]-mu_el-e_)**2*tf
        cv = IntegrationFunc(fn, pe[0:k1])
    else:
        cv = 0.0
    fn = pdos[0:k1]*(pe[0:k1]-mu_el)**2*tf
    c_mu = IntegrationFunc(fn, pe[0:k1])

#   hole/electron concentration by effective carrier
    tf = tf0*(1.0-tf0)
    fn = pdos*tf
    f2 = interp1d(pe, fn, kind='linear')
    fmu = f2(mu_el)
    x = np.hstack([pe[pe<mu_el],mu_el])
    y = np.hstack([fn[pe<mu_el],fmu])
    W_p = IntegrationFunc(y,x)
    x = np.hstack([mu_el, pe[pe>mu_el]])
    y = np.hstack([fmu, fn[pe>mu_el]])
    W_e = IntegrationFunc(y,x)
    #W_e = IntegrationFunc(fn[pe>mu_el], pe[pe>mu_el])
    #W_e = IntegrationFunc(fn[pe>dF], pe[pe>dF])

#   hole/electron concentration by alternative difination
    fn = pdos*(1.0-tf0)
    f2 = interp1d(pe, fn, kind='linear')
    #fmu = f2(mu_el)
    #x = np.hstack([pe[pe<mu_el],mu_el])
    #y = np.hstack([fn[pe<mu_el],fmu])
    try:
        fmu = f2(dF)
    except:
        fmu = 0.
    x = np.hstack([pe[pe<dF],dF])
    y = np.hstack([fn[pe<dF],fmu])
    Y_p = IntegrationFunc(y,x)

    fn = pdos*tf0
    f2 = interp1d(pe, fn, kind='linear')
    #fmu = f2(mu_el)
    #x = np.hstack([mu_el, pe[pe>mu_el]])
    #y = np.hstack([fmu, fn[pe>mu_el]])
    #print ("mu_el", mu_el, dF)
    try:
        fmu = f2(dF)
    except:
        fmu = 0.
    x = np.hstack([dF, pe[pe>dF]])
    y = np.hstack([fmu, fn[pe>dF]])
    Y_e = IntegrationFunc(y,x)

    return mu_el, u, -s*k_B, cv*k_B*Beta*Beta, Q_el, Y_el, Q_p, Q_e, c_mu*k_B*Beta*Beta, W_p, W_e, Y_p, Y_e



def substr(str1, str2, pos):
  try:
    if str1.index(str2)==pos:
        #print("idx=",str1.index(str2))
        return True
    else:
        return False
  except ValueError:
    return False

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def isint(value):
  try:
    int(value)
    return True
  except ValueError:
    return False


def pregetdos(f): # Line 186
    """

    Parameters
    ----------
    dos_file : the DOSCAR from VASP
    xdn : lower energy to integrate over?
    xup : higher energy to integrate over?
    dope (float): number of electrons to dope (negative means to remove electrons, positive means add electrons)
    eFermi (float): Fermi energy

    Returns
    -------
    DOS (array)
    """
    # read the file
    lines = f.readlines() # read in all lines then determine is it is WIEN2k DOS (in the unit eV) file or VASP DOS file
    # now the first line should be the one with the data, lets remove the one into its own special line
    tmp = lines[0]
    if substr(tmp,"#  BAND", 0):
        tmp = lines[1]
        tmp1 = lines[2]
        if substr(tmp, "#EF=",0) and substr(tmp1, "# ENERGY",0):
            tmp1 = tmp[31:43].replace("NENRG=","")
            if isint(tmp1):
                n_dos = int(tmp1)
                tmp = lines[2]
                lines = lines[3:n_dos+3]
                wienEdos = np.zeros(n_dos)
                ve = np.zeros(n_dos)
                for i, l in enumerate(lines):
                    split_l = l.split(' ')
                    split_l = [k for k in split_l if k != '']
                    ve[i], wienEdos[i] = (float(split_l[0]), float(split_l[1]))
                edn = ve[0]
                eup = ve[n_dos-1]
                ve = np.linspace(edn, eup, n_dos)
                vde = (eup - edn)/(n_dos-1) # This appears to be the change of v per electron, so what is v? Voltage in eV?
                return edn, eup, vde, ve, wienEdos

    tmp = lines[5]
    data_line = tmp[0:32].split(' ') #n_dos >10000, no space left before it in VASP
    data_line.extend(tmp[32:].split(' '))
    # filter out empty spaces
    data_line = [k for k in data_line if k != '']
    #print (data_line)
    eup, edn, n_dos, eFermi = (float(data_line[0]),
                           float(data_line[1]),
                           int(data_line[2]),
                           float(data_line[3])) # we're leaving the last number behind
    lines = lines[6:n_dos+6]

    # vectors
    e = np.linspace(edn, eup, n_dos)
    vaspEdos = np.zeros(n_dos)
    ados = np.zeros(n_dos)

    for i, l in enumerate(lines):
        # why do we need to do this?
        split_l = l.split(' ')
        # filter again
        split_l = [k for k in split_l if k != '']
        if len(split_l)>=5: #spin polarized
            e[i], vaspEdos[i], y, vdos[i], x = (float(split_l[0]), float(split_l[1]), float(split_l[2]), float(split_l[3]), float(split_l[4]))
            vaspEdos[i] += y
            ados[i] += x
        else:
            e[i], vaspEdos[i], ados[i] = (float(split_l[0]), float(split_l[1]), float(split_l[2]))
    return e, vaspEdos, ados, eFermi


def getvol(ff):
    with open (ff, "r") as fp:
        lines = fp.readlines()
        aa = float(lines[1])
        a = [float(x) for x in lines[2].split() if len(x)!=0 ]
        b = [float(x) for x in lines[3].split() if len(x)!=0 ]
        c = [float(x) for x in lines[4].split() if len(x)!=0 ]
        vol = a[0]*b[1]*c[2]+a[1]*b[2]*c[0]+a[2]*b[0]*c[1]-a[2]*b[1]*c[0]-a[0]*b[2]*c[1]-a[1]*b[0]*c[2]
        try:
            nat = [int(x) for x in lines[5].split() if len(x)!=0 ]
        except:
            nat = [int(x) for x in lines[6].split() if len(x)!=0 ]
        return vol, sum(nat)


def dosmix(d0, d1, beta, mixonly=True):
    f0 = d0 + "/DOSCAR"
    f1 = d1 + "/DOSCAR"
    with open (f0, 'r') as fp:
        e0, edos0, ados0, eFermi0 = pregetdos(fp)
    with open (f1, 'r') as fp:
        e1, edos1, ados1, eFermi1 = pregetdos(fp)
    vol0,natom0 = getvol(d0 + "/POSCAR")
    vol1,natom1 = getvol(d1 + "/POSCAR")
    if natom0!=natom1:
        print ("************ERROR in POSCAR file natom0=", natom0, "not equal to natom1=", natom1)
    #print(vol0,vol1)

    e0 -= eFermi0
    e1 -= eFermi1
    ne = len(e0[e0<=0])
    de = -min(e0)/ne
    emin = max(min(e0), min(e1))
    emax = min(max(e0), max(e1))
    e = []
    for i in range(ne+1):
        if i*de <emin: break
        e.append(-i*de)
    e.reverse()
    for i in range(1,len(e0[e0>0])+1,1):
        if i*de >emax: break
        e.append(i*de)

    f2 = interp1d(e0, edos0, kind='linear')
    d0 = f2(e)
    f2 = interp1d(e0, ados0, kind='linear')
    a0 = f2(e)
    f2 = interp1d(e1, edos1, kind='linear')
    d1 = f2(e)
    f2 = interp1d(e1, ados1, kind='linear')
    a1 = f2(e)

    eFermi = eFermi0 + beta*(eFermi1-eFermi0)
    edos = d0 + beta*(d1-d0)
    ados = a0 + beta*(a1-a0)

    # line 197 goes to line 209
    if mixonly: return np.array(e), np.array(edos), np.array(d0), np.array(d1), vol0, vol1, natom0

    eup = e[-1] + eFermi
    edn = e[0] + eFermi
    """
    for j in range(5):
        print('   {}'.format(j))
    print ('{:>16.8}{:>16.8}{:5}{:>16.8}{:>16.8}'.format(eup,edn,len(edos),eFermi,1.0))
    """
    ados = cumtrapz(edos, e, initial=ados[0])
    f2 = interp1d(e, edos, kind='linear')
    NELECTRONS = f2(0.0)
    for j,d in enumerate(e):
        #print('{:>16.8f} {:>16.8e} {:>16.8e}'.format(d+eFermi,edos[j],ados[j]))
        print('{:>16.8f} {:>16.8e} {:>16.8e}'.format(d,edos[j],ados[j]))


def thermo_run(beta, nT, T0, T1, nc, d0, d1, dope=0.0):
    dos_energies, vaspEdos, dosf0, dosf1, vol0, vol1, natom = dosmix(d0, d1, beta)
    toJ = Faraday_constant/natom
    n_dos = len(dos_energies) # number of points in DOS
    edn = dos_energies[0]
    eup= dos_energies[-1]
    vde = (eup - edn)/(n_dos-1) # change in energy per step

    compos = np.linspace(0,1.0, nc)
    T = np.linspace(T0,T1, nT)
    with open ("thermo.out", "w") as fout:
      fout.write ('#com, T, M_el(eV), S_el*toJ, C_el*toJ, C_mu*toJ, L, seebeck_coefficients, Q_el/vol, Y_p/vol, Y_e/vol\n')
      for com in compos:
        vaspEdos = dosf0 + com*(dosf1 - dosf0)
        vol = (vol0 + com*(vol1 - vol0))*1.e-24
        NELECTRONS, E0, dF, e, dos, Eg = getdos(-100, 100, dope, 20000, 2000, edn, eup, vde, dos_energies, vaspEdos)
        #NELECTRONS, E0, dF, e, dos, Eg = getdos(-100, 100, 0, 100000, 10000, edn, eup, vde, dos_energies, vaspEdos)

        for t in T:
            if t==0.0:
                U_el = E0
                F_el = E0
                Q_el = 0
                M_el = 0
                S_el = 0
                C_el = 0
                C_mu = 0
                seebeck_coefficients = 0
                L = 2.443004551768e-08 #1.380649e-23/1.60217662e-19x3.14159265359xx2/3
                Q_el = 0
                Q_p = 0
                Q_e = 0
                Y_p = 0
                Y_e = 0
            else:
                Beta = 1.0e0/(t*k_B)
                M_el,U_el,S_el,C_el,Q_el,Y_el,Q_p,Q_e,C_mu,W_p,W_e,Y_p,Y_e = caclf(e, dos, NELECTRONS, Beta, M_el)
                if Q_el>0.0: seebeck_coefficients = -1.0e6*Y_el/Q_el/t
                F_el = (U_el - t * S_el - E0)
                L = 2.443004551768e-08 #1.380649e-23/1.60217662e-19x3.14159265359xx2/3
                if Q_el > 1.e-16: L = C_el/Q_el*k_B

            fout.write ('{:>8.4f} {:8.1f} {:>12.8f} {} {} {} {} {} {} {} {}\n'.format(com, t, M_el, \
                S_el*toJ, C_el*toJ, C_mu*toJ, L, seebeck_coefficients, Q_el/vol, Y_p/vol, Y_e/vol))
        if com!=compos[-1]:
            fout.write ('\n\n')

    plot()


def thermo_dosmix(args):
    """
    Run thermo_dosmix

    Parameters in args
        beta: initial mixxing coef
        T0: Low temperature limit
        T1: High temperature limit
        nT: Number of temperatures
        nc: Number of compositions
        f0: path of the first DOSCAR file
        f1: path of the second DOSCAR file
    """
    thermo_run(args.beta, args.nT, args.T0, args.T1, args.nc, args.d0, args.d1)

def cmd_parser():
    print("\npydosmix: Code to calculate thermo properties of electrons by dos mixxing")
    print("Copyright \u00a9 to be determined\n")

    """
    thermo_dosmix command
    """
    parser = argparse.ArgumentParser(description='thermo properties of electrons by dos mixxing.')

    subparsers = parser.add_subparsers()

    #SUB-PROCESS: thermo_dosmix
    prun = subparsers.add_parser("thermo_dosmix", help="thermo properties of electrons by dos mixxing.")

    prun.add_argument("-beta", "--beta", dest="beta", nargs="?", type=float, default=0.0,
                      help="initial mixxing coef.\n"
                           "Default: 0.0\n")
    prun.add_argument("-T0", "--T0", dest="T0", nargs="?", type=float, default=0,
                      help="Low temperature limit.\n"
                           "Default: 0\n")
    prun.add_argument("-T1", "--T1", dest="T1", nargs="?", type=float, default=100,
                      help="Hi temperature limit.\n"
                           "Default: 100\n")
    prun.add_argument("-nT", "--nTemperature", dest="nT", nargs="?", type=int, default=101,
                      help="Number of temperatures.\n"
                           "Default: 101\n")
    prun.add_argument("-nC", "--nCom", dest="nc", nargs="?", type=int, default=101,
                      help="Number of compositions.\n"
                           "Default: 101\n")

    prun.add_argument("-d0", "--dir0", dest="d0", type=str,
                      help="path of the first phase.\n")

    prun.add_argument("-d1", "--dir1", dest="d1", type=str,
                      help="path of the second phase.\n")

    prun.set_defaults(func=thermo_dosmix)

    args = parser.parse_args()

    try:
        a = getattr(args, "func")
    except AttributeError:
        parser.print_help()
        sys.exit(0)
    args.func(args)


if __name__ == '__main__':
    cmd_parser()
