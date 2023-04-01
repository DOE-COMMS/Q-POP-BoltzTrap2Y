# -------- energy in eV, temperature in K
# assume every variable starting with a-h  and o-z are real numbers
# common block named comcon
from __future__ import division
import sys
import math
import numpy as np
from scipy.constants import physical_constants
from scipy.optimize import brentq
from scipy.integrate import cumtrapz, trapz, simps
from scipy.interpolate import interp1d
from dfttk.pythelec import getdos, caclf

k_B = physical_constants['Boltzmann constant in eV/K'][0]

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


def dosmix(f0, f1, beta, mixonly=True):
    with open (f0, 'r') as fp:
        e0, edos0, ados0, eFermi0 = pregetdos(fp)
    with open (f1, 'r') as fp:
        e1, edos1, ados1, eFermi1 = pregetdos(fp)

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
    if mixonly: return np.array(e), np.array(edos), np.array(d0), np.array(d1)

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


if __name__ == '__main__':

    # handling the command line option
    # TODO: use proper argparse module for this
    beta = 0.0
    natom = 1
    f0 = ""
    f1 = ""
    count = 1
    nx = 101
    nc = 101
    while (count < len(sys.argv)):
      if (sys.argv[count] == "-f"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        beta = float(sys.argv[count])
      elif (sys.argv[count] == "-nc"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        nc = int(sys.argv[count])
      elif (sys.argv[count] == "-nx"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        nx = int(sys.argv[count])
      elif (sys.argv[count] == "-natom"):
        count = count + 1
        if (count > len(sys.argv)):
          break
        natom = int(sys.argv[count])
      elif (f0 == ""):
        f0 = str(sys.argv[count])
      elif (f1 == ""):
        f1 = str(sys.argv[count])
      count = count + 1

    dos_energies, vaspEdos, dosf0, dosf1 = dosmix(f0, f1, beta)
    n_dos = len(dos_energies) # number of points in DOS
    edn = dos_energies[0]
    eup= dos_energies[-1]
    vde = (eup - edn)/(n_dos-1) # change in energy per step

    compos = np.linspace(0,1.0, nc)
    T = np.linspace(0,100, nx)
    print ('#com, T, M_el(eV)')
    for com in compos:
        vaspEdos = dosf0 + com*(dosf1 - dosf0)
        NELECTRONS, E0, dF, e, dos, Eg = getdos(-100, 100, 0, 20000, 2000, edn, eup, vde, dos_energies, vaspEdos)
        #NELECTRONS, E0, dF, e, dos, Eg = getdos(-100, 100, 0, 100000, 10000, edn, eup, vde, dos_energies, vaspEdos)
    
        for t in T:
            if t==0.0:
                U_el = E0
                F_el = E0
                Q_el = 0
                M_el = 0
            else:
                Beta = 1.0e0/(t*k_B)
                M_el,U_el,S_el,C_el,Q_el,Y_el,Q_p,Q_e,C_mu,W_p,W_e,Y_p,Y_e = caclf(e, dos, NELECTRONS, Beta, M_el)
                if Q_el>0.0: seebeck_coefficients = -1.0e6*Y_el/Q_el/t
                F_el = (U_el - t * S_el - E0)
                if Q_el > 1.e-16: L = C_el/Q_el*k_B
            L = 2.443004551768e-08 #1.380649e-23/1.60217662e-19x3.14159265359xx2/3
            print('{:>8.4f} {:8.1f} {:>12.8f}'.format(com, t, M_el))
        if com!=compos[-1]:
            print ('\n')
    
