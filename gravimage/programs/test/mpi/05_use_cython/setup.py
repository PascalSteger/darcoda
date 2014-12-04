# setup.py - unnecessary if not redistributing the code, see below
import numpy
from distutils.core import setup, Extension
from Cython.Build import cythonize

ext_modules = [
    Extension(
        name="samplemesh",
        sources=["mesh_use.pyx"],
        include_dirs = ["./", numpy.get_include()],
        libraries=["fmesh","fmesh_wrapper"],
        library_dirs=["../build/Linux/bin.release","/usr/local/lib/","/usr/lib"],
        language="c++",)
]

setup(name = 'samplemesh',
      cmdclass = {'build_ext': build_ext},
      ext_modules = ext_modules,
      )
