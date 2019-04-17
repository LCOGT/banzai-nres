import operator

from banzai.utils.file_utils import ccdsum_to_filename
from banzai.settings import make_calibration_filename_function
from banzai.utils.instrument_utils import InstrumentCriterion

from banzai import settings

from banzai_nres.fibers import fibers_state_to_filename
from banzai_nres.utils.munge_utils import get_telescope_filename


settings.FRAME_CLASS = 'banzai_nres.images.NRESImage'

settings.FRAME_SELECTION_CRITERIA = [InstrumentCriterion('type', operator.contains, 'NRES')]

settings.ORDERED_STAGES = ['banzai.bpm.BPMUpdater',
                           'banzai.qc.SaturationTest',
                           'banzai.bias.OverscanSubtractor',
                           'banzai.gain.GainNormalizer',
                           'banzai.trim.Trimmer',
                           'banzai.bias.BiasSubtractor',
                           'banzai.dark.DarkSubtractor',
                           'banzai_nres.traces.LoadTrace',
                           'banzai_nres.extract.RectifyTwodSpectrum',
                           'banzai_nres.extract.BoxExtract']

settings.CALIBRATION_MIN_FRAMES = {'BIAS': 5,
                                   'DARK': 3,
                                   'LAMPFLAT': 5,
                                   'TRACE': 1}

# Trace settings
TRACE_FIT_INITIAL_DEGREE_TWO_GUESS = 90  # DO NOT HAPHAZARDLY CHANGE THIS
TRACE_FIT_POLYNOMIAL_ORDER = 4  # DO NOT HAPHAZARDLY CHANGE THIS
TRACE_TABLE_NAME = 'TRACE'
WINDOW_FOR_TRACE_IDENTIFICATION = {'max': 2100, 'min': 2000}  # pixels
MIN_FIBER_TO_FIBER_SPACING = 10  # pixels
MIN_SNR_FOR_TRACE_IDENTIFICATION = 6
# extraction settings
MAX_EXTRACTION_HALF_WINDOW = 10
BOX_EXTRACTION_HALF_WINDOW = 10
#

settings.CALIBRATION_SET_CRITERIA = {'BIAS': ['ccdsum'],
                                     'DARK': ['ccdsum'],
                                     'LAMPFLAT': ['ccdsum', 'fiber0_lit', 'fiber1_lit', 'fiber2_lit'],
                                     'TRACE': ['ccdsum', 'fiber0_lit', 'fiber1_lit', 'fiber2_lit']}

settings.CALIBRATION_FILENAME_FUNCTIONS = {'BIAS': make_calibration_filename_function('BIAS', [ccdsum_to_filename],
                                                                                      get_telescope_filename),
                                           'DARK': make_calibration_filename_function('DARK', [ccdsum_to_filename],
                                                                                      get_telescope_filename),
                                           'LAMPFLAT': make_calibration_filename_function('LAMPFLAT',
                                                                                          [ccdsum_to_filename,
                                                                                           fibers_state_to_filename],
                                                                                          get_telescope_filename),
                                           'TRACE': make_calibration_filename_function('TRACE', [ccdsum_to_filename,
                                                                                       fibers_state_to_filename],
                                                                                       get_telescope_filename)}

settings.CALIBRATION_IMAGE_TYPES = ['BIAS', 'DARK', 'LAMPFLAT', 'TRACE']  # 'WAVELENGTH'

settings.LAST_STAGE = {'BIAS': 'banzai.trim.Trimmer',
                       'DARK': 'banzai.bias.BiasSubtractor',
                       'LAMPFLAT': 'banzai.dark.DarkSubtractor',
                       'TRACE': 'banzai.dark.DarkSubtractor',
                       'DOUBLE': None,
                       'SPECTRUM': None}

settings.EXTRA_STAGES = {'BIAS': ['banzai.bias.BiasMasterLevelSubtractor', 'banzai.bias.BiasComparer'],
                         'DARK': ['banzai.dark.DarkNormalizer', 'banzai.dark.DarkComparer'],
                         'LAMPFLAT': None,
                         'TRACE': None,
                         'DOUBLE': None,
                         'SPECTRUM': None}

settings.CALIBRATION_STACKER_STAGE = {'BIAS': 'banzai.bias.BiasMaker',
                                      'DARK': 'banzai.dark.DarkMaker',
                                      'LAMPFLAT': 'banzai_nres.flats.FlatStacker',
                                      'TRACE': 'banzai_nres.traces.TraceMaker'}  # 'WAVELENGTH':
