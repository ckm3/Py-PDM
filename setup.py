from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np


part_of_readme = '''
                    # Py-PDM

                    A Python wrapper of the Phase Dispersion Minimization (PDM), which is a [C code written by Stellingwerf](https://www.stellingwerf.com/rfs-bin/index.cgi?action=PageView&id=34).

                    Compared with other Python implementations, with the help of Cython, we can obtain a much faster PDM tool.

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
                 '''

extension = Extension(
    name="pdmpy",
    sources=["src/Py-PDM/Py-PDM.pyx"],
    include_dirs=[np.get_include()]
)

setup(
    name="Py-PDM",
    version="0.2",
    author="Kaiming Cui",
    author_email="ckm@nao.cas.cn",
    description="A Python wrapper of the Phase Dispersion Minimization (PDM)",
    long_description=part_of_readme,
    long_description_content_type="text/markdown",
    packages=["Py-PDM"],
    package_dir={"": "src"},
    include_package_data=True,
    url="https://github.com/ckm3/Py-PDM",
    license="GPL-3.0",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics"],
    python_requires='>=3',
    setup_requires=["cython", "numpy"],
    install_requires = ["numpy"],
    ext_modules=cythonize([extension])
)
