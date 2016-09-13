"""setup.py for alignment module.

Author: Federico Tomasi
Copyright (c) 2016, Federico Tomasi.
Licensed under the FreeBSD license (see LICENSE.txt).
"""
from distutils.core import setup, Extension

module1 = Extension('string_kernel',
                    sources=['string_kernel.cpp'])

module2 = Extension('sum_string_kernel',
                    sources=['sum_string_kernel.cpp'])

setup(name='string_kernel',
      version='1.0',
      description='This is a demo package',
      ext_modules=[module1, module2])
