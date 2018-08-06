"""
banzai_nres - Data reduction pipeline for the LCO NRES instruments using BANZAI
Author
    Curtis McCully (cmccully@lco.global)
    G. Mirek Brandt (gmbrandt@ucsb.edu)
    Timothy D. Brandt (tbrandt@physics.ucsb.edu)

License
    GPL v3.0
July 2018
"""
from setuptools import setup, find_packages


setup(name='banzai_nres',
      author=['Curtis McCully', 'G. Mirek Brandt', 'Timothy D. Brandt'],
      author_email=['cmccully@lco.global', 'gmbrandt@ucsb.edu', 'tbrandt@physics.ucsb.edu'],
      version='0.1.0',
      packages=find_packages(),
      package_dir={'banzai_nres': 'banzai_nres'},
      setup_requires=['pytest-runner'],
      install_requires=['banzai', 'scipy'],
      tests_require=['pytest'],
      entry_points={'console_scripts': ['make_master_bias=banzai_nres.main:make_master_bias_console',
                                        'make_master_dark=banzai_nres.main:make_master_dark_console',
                                        'make_master_trace=banzai_nres.main:make_master_trace_console',
                                        'make_master_trace_blind=banzai_nres.main:make_master_trace_blind_console']})
