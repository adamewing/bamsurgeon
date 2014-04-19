#!/usr/bin/env python

"""Bamsurgeon Package.
"""

from setuptools import setup, find_packages
from distutils.core import Extension
from distutils.command import build
import os, sys, subprocess

doclines=__doc__.splitlines()

class my_build(build.build): 
    # different order: build_ext *before* build_py, so that 
    # build_py can already use ctypes! 
    sub_commands = [('build_ext', build.build.has_ext_modules), 
        ('build_py', build.build.has_pure_modules), 
        ('build_clib', build.build.has_c_libraries), 
        ('build_scripts', build.build.has_scripts), ] 

setup(name="bs",
    version="0.0.1",
    description=doclines[0],
    long_description="\n".join(doclines[2:]),
    author="Adam Ewing and Li Xia",
    author_email="li.xia@stanford.edu",
    url="https://github.com/adamewing/bamsurgeon",
    license="TBD",
    platforms=["Linux"],
    packages=find_packages(exclude=['ez_setup', 'test', 'doc']),
    include_package_data=True,
    #packages=['lsa'],
    #requires=['numpy(>=1.1)','scipy(>=0.6)','python(>=2.5)','matplotlib(>=0.98)'],
    zip_safe=False,
    install_requires=["python >= 2.7","pysam >= 0.7"],
    provides=['bs'],
    #can insert C++ extension when needed
    #ext_modules = [ Extension('lsa._compcore', 
    #                sources = ['lsa/compcore_wrap.cpp', 'lsa/compcore.cpp'],
    #                                          depends = ['lsa/compcore.hpp'],
    #               )],
    py_modules = ['bs.asmregion','bs.mutableseq','bs.parseamos','bs.replacereads'],
    cmdclass = {'build': my_build},
    data_files = [('',['README','LICENSE.txt'])],
    entry_points = { 
      'console_scripts': [
            'addindel = bs.addindel:run',
            'addsv = bs.addsv:run',
            'addsnv = bs.addsnv:run'
      ]
   },
)
