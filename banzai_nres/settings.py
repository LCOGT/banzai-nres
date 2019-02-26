import operator

from banzai.settings import Settings
from banzai.utils.file_utils import ccdsum_to_filename
from banzai.calibrations import make_calibration_filename_function
from banzai.context import InstrumentCriterion
from banzai import bias, trim, dark, gain, bpm, qc

from banzai_nres.images import NRESImage
from banzai_nres.fibers import fibers_state_to_filename
from banzai_nres.utils.munge_utils import get_telescope_filename
from banzai_nres import traces

from banzai_nres.flats import FlatStacker


class NRESSettings(Settings):
    FRAME_CLASS = NRESImage

    FRAME_SELECTION_CRITERIA = [InstrumentCriterion('type', operator.contains, 'NRES')]

    ORDERED_STAGES = [bpm.BPMUpdater,
                      qc.SaturationTest,
                      bias.OverscanSubtractor,
                      gain.GainNormalizer,
                      trim.Trimmer,
                      bias.BiasSubtractor,
                      dark.DarkSubtractor,
                      traces.LoadTrace]

    CALIBRATION_MIN_FRAMES = {'BIAS': 5,
                              'DARK': 3,
                              'LAMPFLAT': 5,
                              'TRACE': 1}

    TRACE_FIT_INITIAL_DEGREE_TWO_GUESS = 90  # DO NOT HAPHAZARDLY CHANGE THIS
    TRACE_FIT_POLYNOMIAL_ORDER = 4  # DO NOT HAPHAZARDLY CHANGE THIS
    TRACE_TABLE_NAME = 'TRACE'
    WINDOW_FOR_TRACE_IDENTIFICATION = {'max': 2100, 'min': 2000}  # pixels
    MIN_PEAK_TO_PEAK_SPACING = 10  # pixels
    MIN_SNR_FOR_TRACE_IDENTIFICATION = 6

    CALIBRATION_SET_CRITERIA = {'BIAS': ['ccdsum'],
                                'DARK': ['ccdsum'],
                                'LAMPFLAT': ['ccdsum', 'fiber0_lit', 'fiber1_lit', 'fiber2_lit'],
                                'TRACE': ['ccdsum', 'fiber0_lit', 'fiber1_lit', 'fiber2_lit']}

    CALIBRATION_FILENAME_FUNCTIONS = {'BIAS': make_calibration_filename_function('BIAS', [ccdsum_to_filename],
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

    CALIBRATION_IMAGE_TYPES = ['BIAS', 'DARK', 'LAMPFLAT']

    LAST_STAGE = {'BIAS': trim.Trimmer,
                  'DARK': bias.BiasSubtractor,
                  'LAMPFLAT': dark.DarkSubtractor,
                  'TRACE': dark.DarkSubtractor}

    EXTRA_STAGES = {'BIAS': [bias.BiasMasterLevelSubtractor, bias.BiasComparer],
                    'DARK': [dark.DarkNormalizer, dark.DarkComparer],
                    'LAMPFLAT': [],
                    'TRACE': []}

    CALIBRATION_STACKER_STAGE = {'BIAS': bias.BiasMaker,
                                 'DARK': dark.DarkMaker,
                                 'LAMPFLAT': FlatStacker,
                                 'TRACE': traces.TraceMaker}
