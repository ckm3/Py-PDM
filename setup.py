from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np

extension = Extension(
    name="pdmpy",
    sources=["src/Py-PDM/Py-PDM.pyx"],
    include_dirs=[np.get_include()]
)
setup(
    name="Py-PDM",
    version="0.1",
    author="Kaiming Cui",
    author_email="ckm@nao.cas.cn",
    description="A Python wrapper of the Phase Dispersion Minimization (PDM)",
    packages=["Py-PDM"],
    package_dir={"": "src"},
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
    install_requires = ["cython", "numpy"],
    ext_modules=cythonize([extension])
)
