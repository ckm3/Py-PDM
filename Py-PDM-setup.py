from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np

examples_extension = Extension(
    name="pypdm",
    sources=["Py-PDM.pyx"],
    include_dirs=[np.get_include()]
)
setup(
    name="pypdm",
    ext_modules=cythonize([examples_extension])
)