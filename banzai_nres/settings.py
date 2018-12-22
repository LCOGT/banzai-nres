from banzai.settings import Settings
from banzai.context import InstrumentCriterion
from banzai_nres.images import NRESImage

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
                                'LAMPFLAT': ['ccdsum', 'fibers_state']}

    CALIBRATION_IMAGE_TYPES = ['BIAS', 'DARK', 'LAMPFLAT']

    LAST_STAGE = {'BIAS': trim.Trimmer,
                  'DARK': bias.BiasSubtractor,
                  'LAMPFLAT': dark.DarkSubtractor}

    EXTRA_STAGES = {'BIAS': [bias.BiasMasterLevelSubtractor, bias.BiasComparer, bias.BiasMaker],
                    'DARK': [dark.DarkNormalizer, dark.DarkComparer, dark.DarkMaker],
                    'LAMPFLAT': [FlatStacker]}
