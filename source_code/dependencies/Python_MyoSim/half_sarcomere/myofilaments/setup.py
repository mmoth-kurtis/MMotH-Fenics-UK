from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
cymodule = 'kinetics_cython'

setup(
  name='ODE test',
  ext_modules=[Extension(cymodule, [cymodule + '.pyx'],)],
  cmdclass={'build_ext': build_ext},
)
