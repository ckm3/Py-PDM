from pdmpy import pdm
import numpy as np
# import pytest


def test_pdm_inputs():
    # Test different inputs of the pdm function
    t = np.linspace(0, 20, 99)
    y = np.sin(t)
    s = np.zeros(t.size)
    assert pdm(t, y, s, 0.01, 1-0.001, 0.00001)


def test_pdm_function():
    # Test the basic function of pdm
    t = np.linspace(0, 20, 199)
    y = np.sin(t)
    s = np.zeros(t.size)
    freq, theta = pdm(t, y, s, 0.01, 1-0.001, 0.00001,10)
    main_freq = freq[np.argmin(theta)]
    # primary_freq, secondary_freq = freq[np.argsort(theta)][:2]
    assert np.isclose(main_freq, 1/2/np.pi, atol=0.01), "The pdm's main result is wrong"
    # print(secondary_freq, 1/4/np.pi)
    # assert np.isclose(secondary_freq, 1/4/np.pi, atol=0.01), "The pdm's seoncdary frequency is wrong"