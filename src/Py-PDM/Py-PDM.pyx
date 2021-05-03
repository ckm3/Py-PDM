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
import warnings
cimport numpy as np
np.import_array()


cdef extern from "PyPDM.c":
    int pdm2(int ne, double datx[], double daty[], double sig[], double f_min, double f_max, double delf, int nbins)
    double* f_array
    double* theta_array
    int nf


def pdm(t, y, s=None, f_min=0, f_max=1, delf=0.1, nbin=10):
    """
    The core function of the Phase Dispersion Minimization

    Parameters
    ----------
    t : array like
        Input time array
    y : array like
        Input value array
    s : array like
        Sigma of each y, if no sigma, it should be zeros
    f_min : double
        The minima of the frequency range
    f_max : double
        The maxima of the frequency range
    delf : double
        Delta frequency
    nbin : int
        The number of bins

    Returns
    -----------
    frequency_array : np.ndarry
    theta_array : np.ndarray
    """

    t = np.ascontiguousarray(t, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)

    if s is None:
        s = np.zeros(len(y))

    s = np.ascontiguousarray(s, dtype=np.float64)

    if len(t)==len(y)==len(s):
        if t.ndim!=1:
            raise ValueError('Inputs (t, y, s) must be 1-dimensional')
        pass
    else:
        raise ValueError('Dimensions are not same of input parameters')

    if f_min >= f_max:
        raise ValueError('f_min should be less than f_max')
    elif f_max <= delf:
        raise ValueError('delf should be less than f_max')

    if nbin > 100:
        raise ValueError('The max allowable bins are 100')
    elif nbin <= 0:
        raise ValueError('The number of bins should be greater than 0')

    stacked_array = np.vstack((t, y, s))

    if ~np.isfinite(stacked_array).all():
        warnings.warn("The input arrays contain some NaNs or Infinites, pdm will ignore those points", RuntimeWarning)
        stacked_array = stacked_array[:, np.isfinite(stacked_array).all(axis=0)]
        
    # Add zero at the begining because the one-based number of C source code
    t_copy = np.insert(stacked_array[0, :], 0, 0).copy()
    y_copy = np.insert(stacked_array[1, :], 0, 0).copy()
    sigs_copy = np.insert(stacked_array[2, :], 0, 0).copy()

    cdef double [::1] x2 = t_copy
    cdef double [::1] y2 = y_copy
    cdef double [::1] s2 = sigs_copy

    cdef int data_number = stacked_array.shape[1]

    return_code = pdm2(data_number, &x2[0], &y2[0], &s2[0], f_min, f_max, delf, nbin)

    if return_code == 1:
        pass
    elif return_code == -1:
        raise ValueError('pdm: too few data points, 100 points at least')
    elif return_code == -2:
        raise ValueError('pdm: too many data points, maximum of 100,000 points allowed')

    freq_np_array = np.asarray(<double[:nf]> f_array)
    theta_np_array = np.asarray(<double[:nf]> theta_array)

    return freq_np_array, theta_np_array


