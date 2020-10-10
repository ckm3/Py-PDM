#cython: language_level=3

import numpy as np
cimport numpy as np
np.import_array()
from libc.stdlib cimport free


cdef extern from "PyPDM.c":
    int pdm2(int ne, double datx[], double daty[], double sig[], double f_min, double f_max, double delf)
    double* f_array
    double* theta_array
    int nf


def pdm(np.ndarray times, np.ndarray mags, np.ndarray sigs, double f_min, double f_max, double delf):
    """
    The core function of the Phase Dispersion Minimization

    Parameters:
        times : np.ndarray
        mags : np.ndarray
        sigs : np.ndarray
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

    if len(times)==len(mags)==len(sigs):
        if times.ndim!=1:
            raise ValueError('Invalid time series dimension')
        pass
    else:
        raise ValueError('Dimensions are not same of input parameters')
        
    cdef int data_number = times.size

    times_copy = times.copy()
    mags_copy = mags.copy()
    sigs_copy = sigs.copy()

    cdef double [::1] x2 = times_copy
    cdef double [::1] y2 = mags_copy
    cdef double [::1] s2 = sigs_copy

    return_code = pdm2(data_number, &x2[0], &y2[0], &s2[0], f_min, f_max, delf)

    if return_code == 1:
        pass
    elif return_code == -1:
        raise ValueError('pdm: too few points, 100 points at least.')

    freq_np_array = np.asarray(<double[:nf]> f_array)
    theta_np_array = np.asarray(<double[:nf]> theta_array)

    return freq_np_array, theta_np_array

