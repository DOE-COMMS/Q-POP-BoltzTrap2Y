import numpy as np
import math
from scipy.interpolate import interp1d, splev, splrep, BSpline
from scipy.integrate import cumtrapz, trapz, simps
from BoltzTraP2.units import Volt


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
    dos : refined dos
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


def pregetdos(f): # Line 186
    """
    to make the code can also handle WIEN2k dos in the unit of eV

    Parameters
    ----------
    f : file descriptor for the DOS file

    Returns
    -------
    xdn : lower energy to integrate over?
    xup : higher energy to integrate over?
    vde : band energy intercal
    e (array): band energy mesh the Fermi energy has been shifted to zero
    DOS (array) : e dos
    """
    # read the file
    lines = f.readlines() # read in all lines then determine is it is WIEN2k DOS (in the unit eV) file or VASP DOS file
    # now the first line should be the one with the data, lets remove the one into its own special line
    tmp = lines[0]
    if tmp.startswith("#  BAND"):
        tmp = lines[1]
        tmp1 = lines[2]
        if tmp.startswith("#EF=") and tmp1.startswith("# ENERGY"):
            tmp1 = tmp[31:43].replace("NENRG=","")
            try:
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
            except:
                pass

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
    # line 197 goes to line 209
    eup = eup - eFermi
    edn = edn - eFermi
    vde = (eup - edn)/(n_dos-1) # This appears to be the change of v per electron, so what is v? Voltage in eV?

    # vectors
    ve = np.linspace(edn, eup, n_dos)
    vaspEdos = np.zeros(n_dos)

    for i, l in enumerate(lines):
        # why do we need to do this?
        split_l = l.split(' ')
        # filter again
        split_l = [k for k in split_l if k != '']
        if len(split_l)>=5: #spin polarized
            t, vaspEdos[i], y, vdos, x = (float(split_l[0]), float(split_l[1]), float(split_l[2]), float(split_l[3]), float(split_l[4]))
            vaspEdos[i] += y
        else:
            t, vaspEdos[i], vdos = (float(split_l[0]), float(split_l[1]), float(split_l[2]))
    return edn, eup, vde, ve, vaspEdos


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


class dosremesh:

    def closest(self, e, val):
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

    def remesh(self, xdn, xup, gaussian, dF, eBoF, NEDOS):
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
        dos : refined dos
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


    def refdos(self, dF, tdos):
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

        dos = np.zeros(len(self.e))
        n_dos = len(self.tdos)
        for i in range(0, len(self.e)):
            tx = self.e[i] + dF
            kx = int((tx-self.xdn)/self.vde) # Converted to int, remember the type!
            kx = max([kx,0]) # XXX: is this translated correctly? What is the 1 in fortran?
            kx = min([n_dos-2, kx]) # TODO: the ndos-1 was here before. could be a source of error
            if self.tdos[kx+1]==0.0 and self.ve[kx+1]>0.0 and self.ve[kx+1]<self.vde:
              # handling near the Top of valence band
              if tx >= 0.0:
                dos[i] = 0.0
              else:
                dos[i] = tdos[kx]*tx/self.ve[kx]
            elif self.eBoF > 0.0 and self.tdos[kx]==0.0 and self.ve[kx+1]-self.eBoF<self.vde and \
                self.ve[kx+1]-self.eBoF>0.0:
                # handling near the bottom of conduction band
                if tx <= self.eBoF:
                  dos[i] = 0.0
                else:
                  dos[i] = tdos[kx+1]*(tx-self.eBoF)/(self.ve[kx+1]-self.eBoF)
            else:
              dos[i] = tdos[kx] + (tdos[kx+1] - tdos[kx])/self.vde*(tx - self.ve[kx])

        try: self.e_orig
        except: self.e_orig = None
        #self.e_orig = None
        if self.e_orig is not None:
            try:
                dos *= self.dos_factor
            except:
                with open("dos_BolzTrap2.dat", "w") as fp:
                    for i,e in enumerate(self.e):
                        fp.write('{} {}\n'.format(e/Volt, dos[i]*Volt))

                f2 = interp1d(self.e_orig, self.dos_orig, kind='linear')
                self.dos_factor = np.ones((len(self.e)), dtype=float)
                e_orig_min = min(self.e_orig)
                e_orig_max = max(self.e_orig)
                for i,e in enumerate(self.e):
                    if dos[i]!=0.0:
                        if e>e_orig_min and e<e_orig_max:
                            self.dos_factor[i] = f2(e)/dos[i]
                dos *= self.dos_factor
                ados = cumtrapz(dos, self.e, initial=0.0)
                idx = self.closest(self.e,0.0)
                for idx1 in range(idx-1, 0, -1):
                    if ados[idx1] != ados[idx] : break
                self.NELECTRONS = ados[idx] - self.e[idx]/(self.e[idx1]-self.e[idx])*(ados[idx1]-ados[idx])
                #print (self.NELECTRONS)
        return dos


    def get_dos_orig(self, dos_orig):
        with open(dos_orig, "r") as fp:
            edn, eup, vde, dos_energies, vaspEdos = pregetdos(fp) # Line 186
        dope = 0.0
        NELECTRONS, E0, dF, self.e_orig, self.dos_orig, Eg = getdos(-100, 100, dope, self.NEDOS, self.gaussian, \
            edn, eup, vde, dos_energies, vaspEdos)
        self.e_orig = self.e_orig*Volt
        #print(self.data.fermi)
        self.dos_orig /= Volt
        with open("dos_orig.dat", "w") as fp:
            for i,e in enumerate(self.e_orig):
                fp.write('{} {}\n'.format(e/Volt, self.dos_orig[i]*Volt))


    def __init__(self, _e, tdos, fermi, dos_orig=None, data=None, NEDOS=2000, gaussian=200):
        ve = _e - fermi

        n_dos = len(tdos)
        idx = self.closest(ve,0.0)
        for i in range(idx,n_dos):
            if tdos[i]!=0.0:
                iBoF = i-1
                break

        self.vde = (ve[-1]-ve[0])/(len(ve)-1)

        eBoF = -1.0
        if iBoF>=idx:
          eBoF = ve[iBoF]
          espr = tdos[iBoF+2]-tdos[iBoF+1]
          if espr>0.0:
            espr = tdos[iBoF+1]/espr*self.vde
            if (espr < self.vde):
              eBoF = ve[iBoF+1] - espr

        self.fermi = fermi
        self.data = data
        self.eBoF = eBoF
        self.ve = ve
        self.xdn = min(ve)
        self.xup = max(ve)

        e = np.linspace(self.xdn,self.xup,NEDOS,dtype=float)
        if gaussian != 0.0:
            e = self.remesh(self.xdn, self.xup, gaussian, 0.0, eBoF, NEDOS)
        self.e = e
        self.NEDOS = NEDOS
        self.gaussian = gaussian
        self.tdos = tdos

        if dos_orig is not None:
            self.get_dos_orig(dos_orig)


    def get_emesh(self):
        return self.e + self.fermi
