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
        pdm(t0, y0, s=s0, f_min=0.01, f_max=1, delf=0.00001, nbin=10)
    
    # - Test points more than 100,000
    t0 = np.linspace(0, 20, 199999)
    y0 = np.sin(t0)
    s0 = np.zeros(t0.size)

    with pytest.raises(ValueError, match="pdm: too many data points, maximum of 100,000 points allowed"):
        pdm(t0, y0, s=s0, f_min=0.01, f_max=1, delf=0.00001, nbin=10)

    # - Test different length of data
    t1 = t[:-2]
    with pytest.raises(ValueError, match='Dimensions are not same of input parameters'):
        pdm(t1, y, s, f_min=0.01, f_max=1, delf=0.00001, nbin=10)

    # - Test different dimension of data
    t2 = t.reshape(199,-1)
    with pytest.raises(ValueError, match=r'.* must be 1-dimensional'):
        pdm(t2, y, s, f_min=0.01, f_max=1, delf=0.00001, nbin=10)
    
    # - Test frequency inputs
    with pytest.raises(ValueError, match='f_min should be less than f_max'):
        pdm(t, y, s, f_min=0.01 + 1, f_max=1, delf=0.00001, nbin=10)

    with pytest.raises(ValueError, match='delf should be less than f_max'):
        pdm(t, y, s, f_min=0.01, f_max=1, delf=1, nbin=10)

    with pytest.raises(ValueError, match='f_min should be less than f_max'):
        pdm(t, y, s, f_min=0.01 + 1, f_max=1, delf=1, nbin=10)
    
    # - Test number of bins
    with pytest.raises(ValueError, match='The max allowable bins are 100'):
        pdm(t, y, s, f_min=0.01, f_max=1, delf=0.1, nbin=101)

    with pytest.raises(ValueError, match='The number of bins should be greater than 0'):
        pdm(t, y, s, f_min=0.01, f_max=1, delf=0.1, nbin=-101)

    # - Test input ndarrays contain NaNs and Infinites
    t[5] = np.NaN
    y[-5] = np.inf
    s[3] = -np.inf

    with pytest.warns(RuntimeWarning, match="The input arrays contain some NaNs or Infinites, pdm will ignore those points"):
        freq, theta = pdm(t, y, s, f_min=0.01, f_max=1, delf=0.00001, nbin=10)
    
    main_freq = freq[np.argmin(theta)]
    assert np.isclose(main_freq, 1/2/np.pi, atol=0.01), "The pdm's main result is wrong"

    # - Test input data type
    t = 'test_string'
    y = np.zeros(len(t))
    s = y.tolist()

    with pytest.raises(ValueError):
        pdm(t, y, s, f_min=0.01, f_max=1, delf=0.00001, nbin=10)


def test_pdm_results():
    # Test the basic function of pdm
    t = np.linspace(0, 20, 199)
    y = np.sin(t)
    s = np.zeros(t.size)
    freq, theta = pdm(t, y, s, f_min=0.01, f_max=1-0.001, delf=0.00001, nbin=10)
    main_freq = freq[np.argmin(theta)]

    assert np.isclose(main_freq, 1/2/np.pi, atol=0.01), "The pdm's main result is wrong"


def test_oversize_nf():
    # Test a bug example from my work
    f_min = 0.0484048455605956
    f_max = 0.2004645119176182
    df = (0.2004645119176182 - 0.0484048455605956)/100

    t = np.linspace(0, 20, 199)
    y = np.sin(t)
    s = np.zeros(t.size)

    freq, theta = pdm(t, y, s, f_min, f_max, df, 10)
    main_freq = freq[np.argmin(theta)]

    assert np.isclose(main_freq, 1/2/np.pi, atol=0.01), "The pdm's main result is wrong"
    assert np.isclose(abs(freq[-1] - freq[-2]), df, rtol=1e-8), "The size of freq is above the range of nf"


def test_large_number_frequencies():
    t = np.linspace(0, 20, 101)
    y = np.sin(t)

    freq, theta = pdm(t, y, None, f_min=0.01, f_max=1-0.001, delf=1e-6, nbin=10)
    main_freq = freq[np.argmin(theta)]

    assert np.isclose(main_freq, 1/2/np.pi, atol=0.01), "The pdm's main result is wrong"

def test_plot():
    import matplotlib.pyplot as plt
    t = np.linspace(0, 20, 199)
    y = np.sin(t)
    freq, theta = pdm(t, y, f_min=0.01, f_max=1, delf=0.001)
    plt.plot(freq, theta)

