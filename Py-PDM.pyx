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


def pdm(int data_number, np.ndarray x , np.ndarray y, np.ndarray s, double f_min, double f_max, double delf):
    x_copy = x.copy()
    y_copy = y.copy()
    s_copy = s.copy()

    cdef double [:] x2 = x_copy
    cdef double [:] y2 = y_copy
    cdef double [:] s2 = s_copy

    return_code = pdm2(data_number, &x2[0], &y2[0], &s2[0], f_min, f_max, delf)

    if return_code == 1:
        pass
    elif return_code == -1:
        raise ValueError('pdm: too few points, 100 points at least.')

    freq_np_array = np.asarray(<double[:nf]> f_array)
    theta_np_array = np.asarray(<double[:nf]> theta_array)

    return freq_np_array, theta_np_array

