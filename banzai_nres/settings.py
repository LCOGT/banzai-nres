from banzai import settings

settings.FRAME_CLASS = 'banzai_nres.images.NRESImage'

settings.FRAME_SELECTION_CRITERIA = [('type', 'contains', 'NRES')]

settings.ORDERED_STAGES = ['banzai.bpm.BPMUpdater',
                           'banzai.qc.SaturationTest',
                           'banzai.bias.OverscanSubtractor',
                           'banzai.gain.GainNormalizer',
                           'banzai.trim.Trimmer',
                           'banzai.bias.BiasSubtractor',
                           'banzai.dark.DarkSubtractor',
                           'banzai_nres.traces.LoadTrace']

settings.CALIBRATION_MIN_FRAMES = {'BIAS': 5,
                                   'DARK': 3,
                                   'LAMPFLAT': 5,
                                   'TRACE': 1}

settings.TRACE_FIT_INITIAL_DEGREE_TWO_GUESS = 90  # DO NOT HAPHAZARDLY CHANGE THIS
settings.TRACE_FIT_POLYNOMIAL_ORDER = 4  # DO NOT HAPHAZARDLY CHANGE THIS
settings.TRACE_TABLE_NAME = 'TRACE'
settings.WINDOW_FOR_TRACE_IDENTIFICATION = {'max': 2100, 'min': 2000}  # pixels
settings.MIN_FIBER_TO_FIBER_SPACING = 10  # pixels
settings.MIN_SNR_FOR_TRACE_IDENTIFICATION = 6

settings.CALIBRATION_SET_CRITERIA = {'BIAS': ['ccdsum'],
                                     'DARK': ['ccdsum'],
                                     'LAMPFLAT': ['ccdsum', 'fiber0_lit', 'fiber1_lit', 'fiber2_lit'],
                                     'TRACE': ['ccdsum', 'fiber0_lit', 'fiber1_lit', 'fiber2_lit']}

settings.CALIBRATION_FILENAME_FUNCTIONS = {'BIAS': ('banzai.utils.file_utils.config_to_filename',
                                                    'banzai.utils.file_utils.ccdsum_to_filename'),
                                           'DARK': ('banzai.utils.file_utils.config_to_filename',
                                                    'banzai.utils.file_utils.ccdsum_to_filename'),
                                           'LAMPFLAT': ('banzai.utils.file_utils.config_to_filename',
                                                        'banzai.utils.file_utils.ccdsum_to_filename',
                                                        'banzai_nres.fibers.fibers_state_to_filename'),
                                           'TRACE': ('banzai.utils.file_utils.config_to_filename',
                                                     'banzai.utils.file_utils.ccdsum_to_filename',
                                                     'banzai_nres.fibers.fibers_state_to_filename')}

settings.TELESCOPE_FILENAME_FUNCTION = 'banzai_nres.utils.runtime_utils.get_telescope_filename'


settings.CALIBRATION_IMAGE_TYPES = ['BIAS', 'DARK', 'LAMPFLAT', 'TRACE']

settings.LAST_STAGE = {'BIAS': 'banzai.trim.Trimmer',
                       'DARK': 'banzai.bias.BiasSubtractor',
                       'LAMPFLAT': 'banzai.dark.DarkSubtractor',
                       'TRACE': 'banzai.dark.DarkSubtractor'}

settings.EXTRA_STAGES = {'BIAS': ['banzai.bias.BiasMasterLevelSubtractor', 'banzai.bias.BiasComparer'],
                         'DARK': ['banzai.dark.DarkNormalizer', 'banzai.dark.DarkComparer'],
                         'LAMPFLAT': [],
                         'TRACE': []}

settings.CALIBRATION_STACKER_STAGE = {'BIAS': 'banzai.bias.BiasMaker',
                                      'DARK': 'banzai.dark.DarkMaker',
                                      'LAMPFLAT': 'banzai_nres.flats.FlatStacker',
                                      'TRACE': 'banzai_nres.traces.TraceMaker'}
