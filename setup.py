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


setup(name='lco-banzai-nres',
      author=['Curtis McCully', 'G. Mirek Brandt', 'Timothy D. Brandt'],
      author_email=['cmccully@lco.global', 'gmbrandt@ucsb.edu', 'tbrandt@physics.ucsb.edu'],
      version='0.3.0',
      packages=find_packages(),
      package_dir={'banzai_nres': 'banzai_nres'},
      setup_requires=['pytest-runner'],
      install_requires=['lco-banzai==0.17.0a6', 'numpy>=1.13', 'sphinx', 'coveralls'],
      tests_require=['pytest>=3.5'],
      entry_points = {'console_scripts': ['make_master_bias=banzai_nres.main:make_master_bias',
                                          'make_master_dark=banzai_nres.main:make_master_dark',
                                          'make_master_flat=banzai_nres.main:make_master_flat',
                                          'banzai_nres_reduce_night=banzai_nres.main:reduce_night']})
