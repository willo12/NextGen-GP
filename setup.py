#!/usr/bin/env python

# python setup.py build_ext --inplace

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

# there is an option to compile c_nextgen as library: then remove "c_nextgen.c" from sources, and uncomment libraries entry
# in that case, compile library beforehand with make libc_nextgen.so

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("nextgen",
        #                     libraries=["c_nextgen"],   # compile from library libc_nextgen.so. Remove "c_nextgen.c" entry below
                             sources=["nextgen.pyx", "model.c", "c_nextgen.c", "fields.c", "states.c", "basic_ops.c", "my_ops.c", "tree_io.c"],
                             extra_compile_args=["-O3"],
                             include_dirs=[numpy.get_include(),"include"])],
)


#from distutils.core import setup
#from Cython.Build import cythonize

#setup(
#    ext_modules = cythonize("gppar.pyx")
#)

