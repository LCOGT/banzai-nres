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


setup(name='lco-banzai-nres-e2e',
      author=['Curtis McCully', 'G. Mirek Brandt', 'Marshall Johnson', 'Timothy D. Brandt'],
      author_email=['cmccully@lco.global', 'gmbrandt@ucsb.edu', '@lco.global', 'tbrandt@physics.ucsb.edu'],
      version='0.5.1',
      python_requires='>=3.6',
      packages=find_packages(),
      use_scm_version=True,
      package_dir={'banzai_nres': 'banzai_nres'},
      package_data={'banzai_nres': ['data/ThAr_atlas_ESO.txt', 'data/g2v_template.fits']},
      setup_requires=['pytest-runner', 'setuptools_scm'],
      install_requires=['banzai>=1.0.5', 'numpy>=1.13', 'sphinx', 'coveralls', 'sep<=1.10.0',
                        'statsmodels', 'astropy>=4.1rc1', 'xwavecal==0.1.11', 'scipy==1.3.2', 'photutils==0.7.2', 'boto3',
                        'astroquery'],
      tests_require=['pytest>=3.5'],
      entry_points={'console_scripts': ['banzai_nres_reduce_night=banzai_nres.main:reduce_night',
                                        'banzai_nres_run_realtime_pipeline=banzai_nres.main:nres_run_realtime_pipeline',
                                        'banzai_nres_start_stacking_scheduler=banzai_nres.main:nres_start_stacking_scheduler',
                                        'banzai_nres_make_master_calibrations=banzai_nres.main:nres_make_master_calibrations',
                                        'banzai_nres_add_bpm=banzai_nres.main:add_bpm',
                                        'banzai_nres_populate_bpms=banzai_nres.main:add_bpms_from_archive',
                                        'banzai_nres_munge_phoenix=banzai_nres.data.munge_phoenix_models:main',
                                        'banzai_nres_create_db=banzai_nres.main:create_db',
                                        'banzai_nres_populate_phoenix_models=banzai_nres.main:populate_phoenix_models']})
