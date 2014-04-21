from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system

# compile the fortran modules without linking
fortran_mod_comp = 'gfortran rk_adaptive.f90 -c -o rk_adaptive.o -O3 -fPIC'
print fortran_mod_comp
system(fortran_mod_comp)
shared_obj_comp = 'gfortran pyrk_adaptive.f90 -c -o pyrk_adaptive.o -O3 -fPIC'
print shared_obj_comp
system(shared_obj_comp)

ext_modules = [Extension(# module name:
                         'pyrk_adaptive',
                         # source file:
                         ['pyrk_adaptive.pyx'],
                         # other compile args for gcc
                         extra_compile_args=['-fPIC', '-O3'],
                         # other files to link to
                         extra_link_args=['rk_adaptive.o', 'pyrk_adaptive.o'])]

setup(name = 'pyrk_adaptive',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)
