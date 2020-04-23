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
      author=['Curtis McCully', 'G. Mirek Brandt', 'Marshall Johnson', 'Timothy D. Brandt'],
      author_email=['cmccully@lco.global', 'gmbrandt@ucsb.edu', '@lco.global', 'tbrandt@physics.ucsb.edu'],
      version='0.5.0',
      python_requires='>=3.6',
      packages=find_packages(),
      use_scm_version=True,
      package_dir={'banzai_nres': 'banzai_nres'},
      setup_requires=['pytest-runner', 'setuptools_scm'],
      install_requires=['banzai>=0.28.0', 'numpy>=1.13', 'sphinx', 'coveralls', 'sep<=1.0.3',
                        'statsmodels', 'astropy<=3.2.3', 'xwavecal==0.1.9', 'scipy<=1.3.1'],
      tests_require=['pytest>=3.5'],
      entry_points={'console_scripts': ['banzai_nres_reduce_night=banzai_nres.main:reduce_night',
                                        'banzai_nres_run_realtime_pipeline=banzai_nres.main:nres_run_realtime_pipeline',
                                        'banzai_nres_start_stacking_scheduler=banzai_nres.main:nres_start_stacking_scheduler',
                                        'banzai_nres_make_master_calibrations=banzai_nres.main:nres_make_master_calibrations',
                                        'banzai_nres_add_bpm=banzai_nres.main:add_bpm',
                                        'banzai_nres_add_line_list=banzai_nres.main:add_line_list']})
