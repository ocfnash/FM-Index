from distutils.core import setup
from Cython.Build import cythonize

setup(name = "FMIndex", ext_modules = cythonize('FMIndex.pyx'))
