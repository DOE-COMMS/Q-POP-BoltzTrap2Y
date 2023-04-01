=========
Changelog
=========

0.0 (2021-01-27)
==================

(Contributor: @YiWang)

- Change List:

 - Changes are made only for constant doping calculations, partially according to Phys. Rev. B, 98 (2018) 224101. The major works include:

   1. revised interface.py 
   2. extended bandlib.py into bandlibEXT.py.
   3. extended io.py into ioEXT.py 
   4. added remesh.py for calculations in low temperature (meshed increased to ~100,000.
   5. improved the computational efficiency on the calculations of the chemical potential of electrons (Fermi levels). See line 332-335 in ``bandlibEXT.py``.
   6. added plotpng.py for figure plot.
 