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
      version='0.4.1',
      python_requires='>=3.6',
      packages=find_packages(),
      package_dir={'banzai_nres': 'banzai_nres'},
      setup_requires=['pytest-runner'],
      install_requires=['lco-banzai>=0.20.0', 'numpy>=1.14', 'sphinx', 'coveralls'],
      tests_require=['pytest>=3.5'],
      entry_points={'console_scripts': ['banzai_nres_reduce_night=banzai_nres.main:reduce_night']})
