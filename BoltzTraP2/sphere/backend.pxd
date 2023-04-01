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

# Wrappers for the C++ code in backend.cpp
from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "backend.hpp":
    cpdef enum Magmom_type:
        unpolarized,
        collinear,
        noncollinear
    cdef cppclass Symmetry_analyzer:
        Symmetry_analyzer(double[3][3],
                          double[][3],
                          int[],
                          int,
                          double *,
                          Magmom_type,
                          double) except +
        int get_nrotations()
        int get_ntensors()
        void get_tensor(int, double[9])
    cdef cppclass Sphere_equivalence_builder:
        Sphere_equivalence_builder(double[3][3],
                                   double[][3],
                                   int[],
                                   int,
                                   double *,
                                   Magmom_type,
                                   double,
                                   int[3],
                                   double) except +
        vector[int] get_mapping()
        void get_point(int, int[3])
    cdef cppclass Degeneracy_counter:
        Degeneracy_counter(double[3][3],
                           double[][3],
                           int[],
                           int,
                           double *,
                           Magmom_type,
                           double) except +
        vector[int] get_mapping()
        int count_degeneracy(double[3])
        void get_point(int, double[3])
