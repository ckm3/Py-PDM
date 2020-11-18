# Py-PDM

[![image](http://img.shields.io/pypi/v/Py-PDM.svg)](https://pypi.python.org/pypi/Py-PDM/)

A Python wrapper of the Phase Dispersion Minimization (PDM), which is a [C code written by Stellingwerf](https://www.stellingwerf.com/rfs-bin/index.cgi?action=PageView&id=34).

Compared with other Python implementations, with the help of Cython, we can obtain a much faster PDM tool.

![Example result](/examples/Py-PDM-example.png)

*The red lines show the true frequency and its n times period.*

# Performance
Compared to the [Pure Python implementation of PDM](https://pyastronomy.readthedocs.io/en/latest/pyTimingDoc/pyPDMDoc/pdm.html) of ``PyAstronomy``:

![Comparison result](/examples/Comparison.png)

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
Please refer to the example directory to see in details.
