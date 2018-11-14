from banzai_nres.images import NRESImage
from banzai import settings
from banzai.context import TelescopeCriterion
import operator
from banzai import bias, trim, dark, gain, bpm, qc
import banzai_nres.bias as nres_bias
import banzai_nres.dark as nres_dark
from banzai_nres.flats import FlatStacker

NRES_CRITERIA = [TelescopeCriterion('camera_type', operator.contains, 'NRES'),
                 TelescopeCriterion('schedulable', operator.eq, True)]

settings.ORDERED_STAGES = [bpm.BPMUpdater,
                           qc.HeaderSanity,
                           qc.ThousandsTest,
                           qc.SaturationTest,
                           bias.OverscanSubtractor,
                           gain.GainNormalizer,
                           trim.Trimmer,
                           bias.BiasSubtractor,
                           dark.DarkSubtractor]

BIAS_IMAGE_TYPES = ['BIAS']
BIAS_LAST_STAGE = trim.Trimmer
BIAS_EXTRA_STAGES = [bias.BiasMasterLevelSubtractor, bias.BiasComparer, nres_bias.BiasMaker]

DARK_IMAGE_TYPES = ['DARK']
DARK_LAST_STAGE = bias.BiasSubtractor
DARK_EXTRA_STAGES = [dark.DarkNormalizer, dark.DarkComparer, nres_dark.DarkMaker]

FLAT_IMAGE_TYPES = ['LAMPFLAT']
FLAT_LAST_STAGE = dark.DarkSubtractor
FLAT_EXTRA_STAGES = [FlatStacker]

IMAGE_CLASS = NRESImage
