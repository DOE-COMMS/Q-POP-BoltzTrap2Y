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

# cython: boundscheck=False
# cython: cdivision=True

from BoltzTraP2.sphere.backend cimport Symmetry_analyzer
from BoltzTraP2.sphere.backend cimport Sphere_equivalence_builder
from BoltzTraP2.sphere.backend cimport Degeneracy_counter
from BoltzTraP2.sphere.backend cimport Magmom_type

import numpy as np
import numpy.linalg as la

from libc.stdlib cimport malloc
from libc.stdlib cimport free
from libc.math cimport round
from libcpp.vector cimport vector
cimport numpy as np

np.import_array()


cdef class C_atoms:
    """
    Objects of this class help in translating data from an Atoms object
    into a set of structures compatible with spglib's API.
    """
    # Number of atoms
    cdef int natoms
    # Lattice vectors
    cdef double lattvec[3][3]
    # Atomic positions
    cdef double(*positions)[3]
    # Atomic types
    cdef int[:] types
    # Type of magnetic polarization
    cdef Magmom_type mtype
    # Magnetic moments (if applicable)
    cdef double[:] __magmom
    # Pointer to the beginning of c_magmom or NULL
    cdef double* magmom

    def __cinit__(self, atoms, magmom):
        """Basic constructor taking an ASE Atoms object."""
        self.natoms = len(atoms.numbers)
        lattvec = atoms.get_cell().T
        cdef int i
        cdef int j
        for i in range(3):
            for j in range(3):
                self.lattvec[i][j] = lattvec[i, j]
        positions = atoms.get_scaled_positions()
        self.positions = <double(*)[3]> malloc(
            self.natoms * sizeof(double[3]))
        if self.positions is NULL:
            raise MemoryError
        for i in range(self.natoms):
            for j in range(3):
                self.positions[i][j] = positions[i][j]
        self.types = np.ascontiguousarray(atoms.numbers, dtype=np.intc)
        if magmom is None:
            self.magmom = NULL
        else:
            self.__magmom = np.ascontiguousarray(magmom.ravel())
            self.magmom = &(self.__magmom[0])
        self.mtype = Magmom_type.unpolarized
        if magmom is not None:
            if magmom.ndim == 1:
                self.mtype = Magmom_type.collinear
            elif magmom.ndim == 2:
                self.mtype = Magmom_type.noncollinear

    def __dealloc__(self):
        """Deallocate the dynamically allocated memory."""
        if self.positions != NULL:
            free(self.positions)


def calc_nrotations(atoms, magmom, symprec=1e-4):
    """
    Compute the number of unique rotations in the space group, taking time
    reversal symmetry into account and accounting for the constraints imposed
    by magnetization.
    """
    c_atoms = C_atoms(atoms, magmom)
    analyzer = new Symmetry_analyzer(c_atoms.lattvec,
                                     c_atoms.positions,
                                     &(c_atoms.types[0]),
                                     c_atoms.natoms,
                                     c_atoms.magmom,
                                     c_atoms.mtype,
                                     symprec)
    cdef int nruter = analyzer.get_nrotations()
    return nruter


def calc_tensor_basis(atoms, magmom, symprec=1e-4):
    """
    Return a a set of independent 3x3 tensor that span all possible 3x3
    tensors compatible with the symmetries of the system and its
    magnetization. An order-2 linear response tensor must be a linear
    combination of these.
    """
    c_atoms = C_atoms(atoms, magmom)
    analyzer = new Symmetry_analyzer(c_atoms.lattvec,
                                     c_atoms.positions,
                                     &(c_atoms.types[0]),
                                     c_atoms.natoms,
                                     c_atoms.magmom,
                                     c_atoms.mtype,
                                     symprec)
    cdef int ntensors = analyzer.get_ntensors()
    nruter = np.empty((ntensors, 3, 3))
    cdef double[:, :, :] c_nruter = nruter
    cdef int i
    for i in range(ntensors):
        analyzer.get_tensor(i, &(c_nruter[i, 0, 0]))
    return nruter


def calc_sphere_quotient_set(atoms, magmom, radius, bounds,
                             symprec=1e-4):
    """
    Return a set of equivalence classes for the lattice points in the
    intersection between the supercell based on the lattice of "atoms" with
    the maximum absolute coordinates determined by "bounds", and a sphere of
    radius "radius", with respect to the rotations of the structure in "atoms"
    with the magnetic configuration described by "magmom".
    """
    cdef int i
    c_atoms = C_atoms(atoms, magmom)
    cdef int[:] c_bounds = np.ascontiguousarray(bounds, dtype=np.intc)
    eq_builder = new Sphere_equivalence_builder(c_atoms.lattvec,
                                                c_atoms.positions,
                                                &(c_atoms.types[0]),
                                                c_atoms.natoms,
                                                c_atoms.magmom,
                                                c_atoms.mtype,
                                                radius,
                                                &(c_bounds[0]),
                                                symprec)
    # Get a mapping from the whole point set to a set of irreducible points
    mapping = np.array(eq_builder.get_mapping())
    # And the coordinates of the points in the grid
    grid = np.empty((mapping.size, 3), order="C", dtype=np.intc)
    cdef int[:, :] c_grid = grid
    for i in range(mapping.size):
        eq_builder.get_point(i, &(c_grid[i, 0]))
    # Build the quotient group with the survivors as representatives.
    # Step 1: put equivalent points together.
    indices = mapping.argsort()
    grid = grid[indices, :]
    mapping = mapping[indices]
    # Step 2: find the indices where the mapping changes.
    splitpoints = np.nonzero((mapping[:-1] != mapping[1:]))[0] + 1
    # Step 3: apply the scissors there.
    equivalences = np.split(grid, splitpoints, axis=0)
    return equivalences


def calc_reciprocal_degeneracies(atoms, magmom, kpoints, symprec=1e-4):
    """
    For each reciprocal-space point in the nx3 array kpoints, count the
    number of equivalent reciprocal-space points in a single copy of the
    Brillouin zone. Return an array with those counts.
    """
    cdef int i
    cdef int npoints = kpoints.shape[0]
    c_atoms = C_atoms(atoms, magmom)
    eq_builder = new Degeneracy_counter(c_atoms.lattvec,
                                        c_atoms.positions,
                                        &(c_atoms.types[0]),
                                        c_atoms.natoms,
                                        c_atoms.magmom,
                                        c_atoms.mtype,
                                        symprec)
    # Interrogate the object about the degeneracy of each k point and
    # return the result.
    cdef double[:, :] c_kpoints = np.ascontiguousarray(kpoints)
    nruter = np.empty(npoints, dtype=np.intc)
    for i in range(npoints):
        nruter[i] = eq_builder.count_degeneracy(&(c_kpoints[i, 0]))
    return nruter


def calc_reciprocal_stars(atoms, magmom, kpoints, symprec=1e-4):
    """
    For each reciprocal-space point in the nx3 array kpoints, find all the
    equivalent reciprocal-space points in the first copy of the Brillouin zone.
    Return a list of arrays with those points.
    """
    cdef int i
    cdef int j
    cdef int npoints = kpoints.shape[0]
    c_atoms = C_atoms(atoms, magmom)
    eq_builder = new Degeneracy_counter(c_atoms.lattvec,
                                        c_atoms.positions,
                                        &(c_atoms.types[0]),
                                        c_atoms.natoms,
                                        c_atoms.magmom,
                                        c_atoms.mtype,
                                        symprec)
    # Interrogate the object about the degeneracy of each k point and
    # find a representative of each class in the generated mapping.
    # Each of these classes is a set of points related by translations
    # in reciprocal space by a reciprocal-lattice vector.
    cdef double[:, :] c_kpoints = np.ascontiguousarray(kpoints)
    nruter = []
    cdef double[:, :] c_reprs
    for i in range(npoints):
        eq_builder.count_degeneracy(&(c_kpoints[i, 0]))
        mapping = np.array(eq_builder.get_mapping())
        positions = np.unique(mapping)
        reprs = np.empty((positions.size, 3), order="C")
        c_reprs = reprs
        for j, p in enumerate(positions):
            eq_builder.get_point(p, &(c_reprs[j, 0]))
        nruter.append(reprs)
    return nruter


cpdef calc_reciprocal_iksubset(atoms, magmom, kpoints, symprec=1e-4):
    """
    Return the indices of an irreducible subset of the provided k points, which
    must be passed in as a nx3 array.
    """
    cdef int nstar
    cdef int nstars
    cdef int nsofar
    cdef int i
    cdef int j
    cdef int k
    cdef int l
    cdef double[:, :] stari
    cdef double[:, :] sofar
    stars = calc_reciprocal_stars(atoms, magmom, kpoints, symprec)
    nstars = len(stars)
    sofar = np.empty_like(kpoints, order="C")
    nruter = []
    nsofar = 0
    for i in range(nstars):
        stari = stars[i]
        nstar = len(stars[i])
        found = False
        for j in range(nstar):
            for k in range(nsofar):
                for l in range(3):
                    delta = sofar[k, l] - stari[j, l]
                    delta -= round(delta)
                    if abs(delta) > symprec:
                        break
                else:
                    found = True
                    break
            if found:
                break
        else:
            nruter.append(i)
            sofar[nsofar, :] = stari[0, :]
            nsofar += 1
    return nruter
