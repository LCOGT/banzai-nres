import os
import banzai_nres

FRAME_FACTORY = 'banzai_nres.frames.NRESFrameFactory'

CALIBRATION_FRAME_CLASS = 'banzai_nres.frames.NRESCalibrationFrame'

FRAME_SELECTION_CRITERIA = [('type', 'contains', 'NRES')]

ORDERED_STAGES = [
                  'banzai.bpm.BadPixelMaskLoader',
                  'banzai.bias.OverscanSubtractor',
                  'banzai.gain.GainNormalizer',
                  'banzai.trim.Trimmer',
                  'banzai.bias.BiasSubtractor',
                  'banzai.uncertainty.PoissonInitializer',
                  'banzai.dark.DarkSubtractor',
                  'banzai_nres.flats.FlatLoader',
                  # this is turned off because it yields negative fluxes and causes crashing on tracing. See issue #60
                  # 'banzai_nres.background.BackgroundSubtractor',
                  'banzai_nres.wavelength.ArcLoader',
                  'banzai_nres.extract.GetOptimalExtractionWeights',
                  'banzai_nres.extract.WeightedExtract',
                  'banzai_nres.continuum.MaskBlueHookRegion',
                  'banzai_nres.continuum.ContinuumNormalizer',
                  'banzai_nres.continuum.MaskTellurics',
                  'banzai_nres.classify.StellarClassifier',
                  'banzai_nres.rv.RVCalculator'
                  ]

CALIBRATION_MIN_FRAMES = {'BIAS': 5,
                          'DARK': 3,
                          'LAMPFLAT': 5,
                          'DOUBLE': 1,
                          }

CALIBRATION_SET_CRITERIA = {'BIAS': ['binning', 'configuration_mode'],
                            'DARK': ['binning', 'configuration_mode'],
                            'LAMPFLAT': ['binning', 'fiber0_lit', 'fiber1_lit', 'fiber2_lit'],
                            'DOUBLE': ['binning', 'fiber0_lit', 'fiber1_lit', 'fiber2_lit'],
                            'LINELIST': [],
                            }

CALIBRATION_FILENAME_FUNCTIONS = {'BIAS': ('banzai_nres.utils.file_utils.config_to_filename',
                                           'banzai.utils.file_utils.ccdsum_to_filename'),
                                  'DARK': ('banzai_nres.utils.file_utils.config_to_filename',
                                           'banzai.utils.file_utils.ccdsum_to_filename'),
                                  'LAMPFLAT': ('banzai_nres.utils.file_utils.config_to_filename',
                                               'banzai.utils.file_utils.ccdsum_to_filename',
                                               'banzai_nres.fibers.fibers_state_to_filename'),
                                  'DOUBLE': ('banzai_nres.utils.file_utils.config_to_filename',
                                             'banzai.utils.file_utils.ccdsum_to_filename',
                                             'banzai_nres.fibers.fibers_state_to_filename')
                                  }

TELESCOPE_FILENAME_FUNCTION = 'banzai_nres.utils.runtime_utils.get_telescope_filename'

CALIBRATION_IMAGE_TYPES = ['BPM', 'BIAS', 'DARK', 'LAMPFLAT', 'DOUBLE']

LAST_STAGE = {'BIAS': 'banzai.trim.Trimmer',
              'DARK': 'banzai.uncertainty.PoissonInitializer',
              'LAMPFLAT': 'banzai.dark.DarkSubtractor',
              'DOUBLE': 'banzai.dark.DarkSubtractor',
              'TARGET': None,
              }

EXTRA_STAGES = {'BIAS': ['banzai.bias.BiasMasterLevelSubtractor', 'banzai.bias.BiasComparer'],
                'DARK': ['banzai.dark.DarkNormalizer', 'banzai.dark.DarkComparer'],
                'LAMPFLAT': [],
                'DOUBLE': [],
                'TARGET': None,
                }

CALIBRATION_STACKER_STAGES = {'BIAS': ['banzai.bias.BiasMaker'],
                              'DARK': ['banzai.dark.DarkMaker'],
                              'LAMPFLAT': ['banzai_nres.flats.FlatStacker',
                                           'banzai_nres.flats.FlatLoader',
                                           'banzai_nres.traces.TraceInitializer',
                                           'banzai_nres.background.BackgroundSubtractor',
                                           'banzai_nres.traces.TraceRefiner',
                                           'banzai_nres.profile.ProfileFitter'
                                           ],
                              'DOUBLE': ['banzai_nres.wavelength.ArcStacker',  # stack
                                         'banzai_nres.flats.FlatLoader',  # load traces
                                         'banzai_nres.wavelength.ArcLoader',  # load wavelengths, ref_ids, etc...
                                         'banzai_nres.wavelength.LineListLoader',  # load reference lab wavelengths
                                         # 'banzai_nres.background.BackgroundSubtractor',
                                         'banzai_nres.wavelength.IdentifyFeatures',
                                         'banzai_nres.wavelength.WavelengthCalibrate',
                                         'banzai_nres.qc.qc_wavelength.AssessWavelengthSolution'
                                         ]
                              }

# Stack delays are expressed in seconds. 300 is five minutes.
CALIBRATION_STACK_DELAYS = {'BIAS': 300,
                            'DARK': 300,
                            'LAMPFLAT': 300,
                            'DOUBLE': 300,
                            }

SCHEDULE_STACKING_CRON_ENTRIES = {'cpt': {'minute': 0, 'hour': 15},
                                  'tlv': {'minute': 45, 'hour': 8},
                                  'lsc': {'minute': 0, 'hour': 21},
                                  'elp': {'minute': 0, 'hour': 23}}

OBSERVATION_PORTAL_URL = os.getenv('OBSERVATION_PORTAL_URL',
                                   'http://internal-observation-portal.lco.gtn/api/observations/')
CALIBRATE_PROPOSAL_ID = os.getenv('CALIBRATE_PROPOSAL_ID', 'calibrate')

CONFIGDB_URL = os.getenv('CONFIGDB_URL', 'http://configdb.lco.gtn/sites/')

OBSERVATION_REQUEST_TYPES = {'BIAS': 'NRESBIAS', 'DARK': 'NRESDARK', 'DOUBLE': 'ARC'}

# For some extension names, we want to just have corresponding BPM or ERR extensions
EXTENSION_NAMES_TO_CONDENSE = ['SPECTRUM']

# number of days to lookback and stack frames:
CALIBRATION_LOOKBACK = {'BIAS': 2.5, 'DARK': 4.5, 'LAMPFLAT': 0.5, 'DOUBLE': 0.5}

PIPELINE_VERSION = banzai_nres.__version__

# Number of days before proprietary data should become public:
DATA_RELEASE_DELAY = 365

# Proposal ids for data that should be public instantly. Should all be lowercase
PUBLIC_PROPOSALS = ['calibrate', 'standard', '*epo*', 'pointing']

SUPPORTED_FRAME_TYPES = ['BPM', 'BIAS', 'DARK', 'LAMPFLAT', 'TARGET', 'DOUBLE']

# Specify different sources of raw and processed data, if needed.
ARCHIVE_API_ROOT = os.getenv('API_ROOT')
ARCHIVE_AUTH_TOKEN = os.getenv('AUTH_TOKEN')
ARCHIVE_FRAME_URL = f'{ARCHIVE_API_ROOT}frames'
ARCHIVE_AUTH_HEADER = {'Authorization': f'Token {ARCHIVE_AUTH_TOKEN}'}

RAW_DATA_AUTH_TOKEN = os.getenv('RAW_DATA_AUTH_TOKEN', ARCHIVE_AUTH_TOKEN)
RAW_DATA_API_ROOT = os.getenv('RAW_DATA_API_ROOT', ARCHIVE_API_ROOT)
RAW_DATA_FRAME_URL = f'{RAW_DATA_API_ROOT}frames'
RAW_DATA_AUTH_HEADER = {'Authorization': f'Token {RAW_DATA_AUTH_TOKEN}'}

FITS_EXCHANGE = os.getenv('FITS_EXCHANGE', 'archived_fits')

LOSSLESS_EXTENSIONS = ['PROFILE', 'WAVELENGTH']

REDUCED_DATA_EXTENSION_ORDERING = {'BIAS': ['SPECTRUM', 'BPM', 'ERR'],
                                   'DARK': ['SPECTRUM', 'BPM', 'ERR'],
                                   'LAMPFLAT': ['SPECTRUM', 'BPM', 'ERR'],
                                   'DOUBLE': ['SPECTRUM', 'BPM', 'ERR'],
                                   'TARGET': ['SPECTRUM', 'BPM', 'ERR', 'TRACES', 'PROFILE', 'BLAZE', 'WAVELENGTH']}

MASTER_CALIBRATION_EXTENSION_ORDER = {'BIAS': ['SPECTRUM', 'BPM', 'ERR'],
                                      'DARK': ['SPECTRUM', 'BPM', 'ERR'],
                                      'LAMPFLAT': ['SPECTRUM', 'BPM', 'ERR', 'TRACES', 'PROFILE', 'BLAZE'],
                                      'DOUBLE': ['SPECTRUM', 'BPM', 'ERR', 'TRACES', 'PROFILE', 'BLAZE', 'WAVELENGTH',
                                                 'FEATURES']}

REDUCED_DATA_EXTENSION_TYPES = {'ERR': 'float32',
                                'BPM': 'uint8',
                                'SPECTRUM': 'float32',
                                'TRACES': 'int32',  # try uint8
                                'PROFILE': 'float32',
                                'WAVELENGTH': 'float64',
                                }

PHOENIX_MODEL_LOCATION = os.getenv('PHOENIX_FILE_LOCATION', 's3://banzai-nres-phoenix-models-lco-global')

MIN_ORDER_TO_CORRELATE = 77
MAX_ORDER_TO_CORRELATE = 97

GAIA_CLASS = os.getenv('BANZAI_GAIA_CLASS', 'astroquery.gaia.GaiaClass')

SIMBAD_CLASS = os.getenv('BANZAI_SIMBAD', 'astroquery.simbad.Simbad')

# The final trace will be +- this from the center in the y-direction
TRACE_HALF_HEIGHT = 5
