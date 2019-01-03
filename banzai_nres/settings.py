from banzai.settings import Settings
from banzai.utils.file_utils import ccdsum_to_filename
from banzai.context import InstrumentCriterion
from banzai_nres.images import NRESImage
from banzai_nres.fibers import fibers_state_to_filename

import operator
from banzai import bias, trim, dark, gain, bpm, qc
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
                      dark.DarkSubtractor]

    CALIBRATION_MIN_IMAGES = {'BIAS': 5,
                              'DARK': 3,
                              'LAMPFLAT': 5}

    CALIBRATION_SET_CRITERIA = {'BIAS': ['ccdsum'],
                                'DARK': ['ccdsum'],
                                'LAMPFLAT': ['ccdsum', 'fiber1_lit', 'fiber2_lit', 'fiber3_lit']}

    CALIBRATION_FILENAME_FUNCTIONS = {'BIAS': [ccdsum_to_filename],
                                      'DARK': [ccdsum_to_filename],
                                      'LAMPFLAT': [ccdsum_to_filename, fibers_state_to_filename]}

    CALIBRATION_IMAGE_TYPES = ['BIAS', 'DARK', 'LAMPFLAT']

    LAST_STAGE = {'BIAS': trim.Trimmer,
                  'DARK': bias.BiasSubtractor,
                  'LAMPFLAT': dark.DarkSubtractor}

    EXTRA_STAGES = {'BIAS': [bias.BiasMasterLevelSubtractor, bias.BiasComparer, bias.BiasMaker],
                    'DARK': [dark.DarkNormalizer, dark.DarkComparer, dark.DarkMaker],
                    'LAMPFLAT': [FlatStacker]}
