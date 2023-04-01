import numpy as np
import math
from scipy.constants import physical_constants
from scipy.optimize import brentq, curve_fit
from scipy.integrate import cumtrapz, trapz, simps
from scipy.interpolate import interp1d, splev, splrep, BSpline
from scipy.integrate import quadrature
from scipy.interpolate import UnivariateSpline
from scipy.constants import physical_constants

k_B = physical_constants['Boltzmann constant in eV/K'][0]


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


