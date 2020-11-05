#cython: language_level=3

'''
Py-PDM - a Python wrapper of PMD
Copyright (C) 2020, Kaiming Cui

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

import numpy as np
cimport numpy as np
np.import_array()
from libc.stdlib cimport free


cdef extern from "PyPDM.c":
    int pdm2(int ne, double datx[], double daty[], double sig[], double f_min, double f_max, double delf, int nbins)
    double* f_array
    double* theta_array
    int nf


def pdm(np.ndarray t, np.ndarray y, np.ndarray s, double f_min, double f_max, double delf, int nbin):
    """
    The core function of the Phase Dispersion Minimization

    Parameters:
        t : np.ndarray
        y : np.ndarray
        s : np.ndarray
            Sigma of each mag, if no sigma, it should be zeros
        f_min : double
            The minima of the frequency range
        f_max : double
            The maxima of the frequency range
        delf : double
            Delta frequency

    Returns:
        frequency_array : np.ndarry
        theta_array : np.ndarray
    """

    if len(t)==len(y)==len(s):
        if t.ndim!=1:
            raise ValueError('Inputs (t, y, s) must be 1-dimensional')
        pass
    else:
        raise ValueError('Dimensions are not same of input parameters')

    if f_min >= f_max:
        raise ValueError('f_min should be less than f_max')
    elif f_max < delf:
        raise ValueError('delf should be less than f_max')

    if nbin > 100:
        raise ValueError('The max allowable bins are 100')
    elif nbin <= 0:
        raise ValueError('The number of bins should be greater than 0')
        
    cdef int data_number = t.size

    t_copy = t.copy()
    y_copy = y.copy()
    sigs_copy = s.copy()

    cdef double [::1] x2 = t_copy
    cdef double [::1] y2 = y_copy
    cdef double [::1] s2 = sigs_copy

    return_code = pdm2(data_number, &x2[0], &y2[0], &s2[0], f_min, f_max, delf, nbin)

    if return_code == 1:
        pass
    elif return_code == -1:
        raise ValueError('pdm: too few points, 100 points at least.')

    freq_np_array = np.asarray(<double[:nf]> f_array)
    theta_np_array = np.asarray(<double[:nf]> theta_array)

    return freq_np_array, theta_np_array

