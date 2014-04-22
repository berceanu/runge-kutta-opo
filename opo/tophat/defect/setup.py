from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system

# compile the fortran modules without linking
shared_obj_comp = 'make pyrk_adaptive.o'
print shared_obj_comp
system(shared_obj_comp)

ext_modules = [Extension(# module name:
                         'pyrk_adaptive',
                         # source file:
                         ['pyrk_adaptive.pyx'],
                         # other compile args for gcc
                         #extra_compile_args=[''],
                         # other files to link to
                         extra_link_args=['FFTW3.o', 'nrtype.o', 'ode_path.o', 'global.o', 'rk_adaptive.o', 'pyrk_adaptive.o'])]

setup(name = 'pyrk_adaptive',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)
