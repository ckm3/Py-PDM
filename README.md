# Py-PDM

[![image](http://img.shields.io/pypi/v/Py-PDM.svg)](https://pypi.python.org/pypi/Py-PDM/) [![Documentation Status](https://readthedocs.org/projects/py-pdm/badge/?version=latest)](https://py-pdm.readthedocs.io/en/latest/?badge=latest) [![CI](https://github.com/ckm3/Py-PDM/actions/workflows/main.yml/badge.svg)](https://github.com/ckm3/Py-PDM/actions/workflows/main.yml) ![GitHub](https://img.shields.io/github/license/ckm3/Py-PDM)
      

A Python wrapper of the Phase Dispersion Minimization (PDM), which is a [C code written by Stellingwerf](https://www.stellingwerf.com/rfs-bin/index.cgi?action=PageView&id=34).

Compared with other Python implementations, with the help of Cython, we can obtain a much faster PDM tool.

![Example result](/docs/source/Py-PDM-example.png)

*The red lines show the true frequency and its n times period.*

# Performance
Compared to the [Pure Python implementation of PDM](https://pyastronomy.readthedocs.io/en/latest/pyTimingDoc/pyPDMDoc/pdm.html) of ``PyAstronomy``:

![Comparison result](/docs/source/Comparison.png)

# Installation
**Before installing, make sure you have already installed the cython and numpy in your Python3 environment.**

To install Py-PDM with pip:

```
pip install py-pdm
```

Alternatively you can install it manually:
```
git clone https://github.com/ckm3/Py-PDM.git
cd Py-PDM
python setup.py install
```

# Usage
```python3
from pdmpy import pdm

freq, theta = pdm(time, y_value, y_sigma, frequency_min, frequency_max, frequency_step, number_of_bins)
```
Please refer to the [documentation](https://py-pdm.readthedocs.io) to see in details.

# Citing
If you find Py-PDM useful in your research, please cite the orginal paper [Stellingwerf, Astrophysical Journal, v224, p953, 1978](https://ui.adsabs.harvard.edu/abs/1978ApJ...224..953S/abstract) and add a footnote of this Github repository.
