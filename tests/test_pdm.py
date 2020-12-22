from pdmpy import pdm
import numpy as np
import warnings
import pytest

def test_pdm_inputs():
    """
    Test different inputs of the pdm function
    """
    t = np.linspace(0, 20, 199)
    y = np.sin(t)
    s = np.zeros(t.size)

    # - Test points less than 100
    t0 = t[:99]
    y0 = np.sin(t0)
    s0 = np.zeros(t0.size)

    with pytest.raises(ValueError, match="pdm: too few data points, 100 points at least"):
        pdm(t0, y0, s0, 0.01, 1, 0.00001, 10)

    # - Test input ndarrays contain NaNs and Infinites
    t[5] = np.NaN
    y[-5] = np.inf
    s[3] = -np.inf

    with pytest.warns(RuntimeWarning, match="The input arrays contain some NaNs or Infinites, pdm will ignore those points"):
        pdm(t, y, s, 0.01, 1, 0.00001, 10)


    # - Test input data type
    t = 'test_string'
    y = np.zeros(len(t))
    s = y.tolist()

    with pytest.raises(TypeError):
        pdm(t, y, s, 0.01, 1, 0.00001, 10)


    # - Test different length of data


def test_pdm_results():
    # Test the basic function of pdm
    t = np.linspace(0, 20, 199)
    y = np.sin(t)
    s = np.zeros(t.size)
    freq, theta = pdm(t, y, s, 0.01, 1-0.001, 0.00001,10)
    main_freq = freq[np.argmin(theta)]

    assert np.isclose(main_freq, 1/2/np.pi, atol=0.01), "The pdm's main result is wrong"


def test_oversize_nf():
    # Test a bug example from my work
    f_min = 0.0484048455605956
    f_max = 0.2004645119176182
    df = (0.2004645119176182 - 0.0484048455605956)/1e2

    t = np.linspace(0, 20, 199)
    y = np.sin(t)
    s = np.zeros(t.size)
    freq, theta = pdm(t, y, s, f_min, f_max, df, 10)
    
    assert (freq > f_min) & np.isclose(freq[-1], freq[-2], atol=df), "The size of freq is above the range of nf"