#!/usr/bin/python
"""string_kernel setup script.

Author: Federico Tomasi
Copyright (c) 2016, Federico Tomasi.
Licensed under the FreeBSD license (see LICENSE.txt).
"""

from distutils.core import setup, Extension

# Package Version
# sk_module = Extension(
#     'string_kernel.core.src.string_kernel',
#     sources=['string_kernel/core/src/string_kernel.cpp'])
ssk_module = Extension(
    'string_kernel.core.src.sum_string_kernel',
    sources=['string_kernel/core/src/sum_string_kernel.cpp'])
setup(
    name='string_kernel',
    version='0.1a',

    description=(''),
    long_description=open('README.md').read(),
    author='Federico Tomasi',
    author_email='federico.tomasi@dibris.unige.it',
    maintainer='Federico Tomasi',
    maintainer_email='federico.tomasi@dibris.unige.it',
    url='https://github.com/slipguru/icing',
    download_url='https://github.com/slipguru/icing/tarball/0.1a',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Programming Language :: Python',
        'License :: OSI Approved :: BSD License',
        'Topic :: Software Development',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Operating System :: MacOS'
    ],
    license='FreeBSD',
    packages=['string_kernel', 'string_kernel.core', 'string_kernel.core.src',
              ],
    requires=['numpy (>=1.10.1)',
              'scipy (>=0.16.1)',
              'sklearn (>=0.17)',
              'matplotlib (>=1.5.1)',
              'seaborn (>=0.7.0)'],
    # scripts=['scripts/ici_run.py', 'scripts/ici_analysis.py'],
    # ext_modules=[sk_module, ssk_module]
    ext_modules=[ssk_module]
)
