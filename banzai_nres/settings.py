from banzai.settings import Settings
from banzai.utils.file_utils import ccdsum_to_filename
from banzai.calibrations import make_calibration_filename_function
from banzai.context import InstrumentCriterion
from banzai_nres.images import NRESImage
from banzai_nres.fibers import fibers_state_to_filename
from banzai_nres import traces
import operator
from banzai import bias, trim, dark, gain, bpm, qc
from banzai_nres.flats import FlatStacker


def get_telescope_filename(image):
    return image.header.get('TELESCOP', '').replace('nres', 'nrs')


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
                      traces.InitialTraceFit]

    CALIBRATION_MIN_IMAGES = {'BIAS': 5,
                              'DARK': 3,
                              'LAMPFLAT': 5,
                              'TRACE': 1}

    TRACE_FIT_INITIAL_DEGREE_TWO_GUESS = 90  # DO NOT HAPHAZARDLY CHANGE THIS

    CALIBRATION_SET_CRITERIA = {'BIAS': ['ccdsum'],
                                'DARK': ['ccdsum'],
                                'LAMPFLAT': ['ccdsum', 'fiber0_lit', 'fiber1_lit', 'fiber2_lit'],
                                'TRACE': ['ccdsum', 'fiber0_lit', 'fiber1_lit', 'fiber2_lit']}

    CALIBRATION_FILENAME_FUNCTIONS = {'BIAS': make_calibration_filename_function('BIAS', [ccdsum_to_filename], get_telescope_filename),
                                      'DARK': make_calibration_filename_function('DARK', [ccdsum_to_filename], get_telescope_filename),
                                      'LAMPFLAT': make_calibration_filename_function('LAMPFLAT', [ccdsum_to_filename, fibers_state_to_filename], get_telescope_filename),
                                      'TRACE': make_calibration_filename_function('TRACE', [ccdsum_to_filename], get_telescope_filename)}

    CALIBRATION_IMAGE_TYPES = ['BIAS', 'DARK', 'LAMPFLAT']

    LAST_STAGE = {'BIAS': trim.Trimmer,
                  'DARK': bias.BiasSubtractor,
                  'LAMPFLAT': dark.DarkSubtractor,
                  'TRACE': traces.InitialTraceFit}

    EXTRA_STAGES = {'BIAS': [bias.BiasMasterLevelSubtractor, bias.BiasComparer, bias.BiasMaker],
                    'DARK': [dark.DarkNormalizer, dark.DarkComparer, dark.DarkMaker],
                    'LAMPFLAT': [FlatStacker],
                    'TRACE': [traces.TraceMaker]}

