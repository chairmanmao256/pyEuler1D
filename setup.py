from setuptools import setup
from Cython.Build import cythonize

setup(
    name='pyEuler1D app',
    ext_modules=cythonize("pyEuler1D.pyx",language_level=3),
    zip_safe=False,
)