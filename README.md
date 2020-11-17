# Py-PDM

A Python wrapper of the Phase Dispersion Minimization (PDM), which is a [C code written by Stellingwerf](https://www.stellingwerf.com/rfs-bin/index.cgi?action=PageView&id=34).

Compared with other Python implementations, with the help of Cython, we can obtain a much faster PDM tool.

![Example result](/examples/Py-PDM-example.png)
*The red lines show the true frequency and its harmonics.*

# Installation
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
