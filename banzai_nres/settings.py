import os

FRAME_FACTORY = 'banzai.images.LCOImageFactory'

FRAME_SELECTION_CRITERIA = [('type', 'contains', 'NRES')]

ORDERED_STAGES = ['banzai.bpm.BadPixelMaskLoader',
                           'banzai.bias.OverscanSubtractor',
                           'banzai.gain.GainNormalizer',
                           'banzai.trim.Trimmer',
                           'banzai.bias.BiasSubtractor',
                           'banzai.dark.DarkSubtractor']

CALIBRATION_MIN_FRAMES = {'BIAS': 5,
                                   'DARK': 3,
                                   'LAMPFLAT': 5}

CALIBRATION_SET_CRITERIA = {'BIAS': ['binning', 'configuration_mode'],
                                     'DARK': ['binning', 'configuration_mode'],
                                     'LAMPFLAT': ['binning', 'fiber0_lit', 'fiber1_lit', 'fiber2_lit']}

CALIBRATION_FILENAME_FUNCTIONS = {'BIAS': ('banzai.utils.file_utils.config_to_filename',
                                                    'banzai.utils.file_utils.ccdsum_to_filename'),
                                           'DARK': ('banzai.utils.file_utils.config_to_filename',
                                                    'banzai.utils.file_utils.ccdsum_to_filename'),
                                           'LAMPFLAT': ('banzai.utils.file_utils.config_to_filename',
                                                        'banzai.utils.file_utils.ccdsum_to_filename',
                                                        'banzai_nres.fibers.fibers_state_to_filename')}

TELESCOPE_FILENAME_FUNCTION = 'banzai_nres.utils.runtime_utils.get_telescope_filename'


CALIBRATION_IMAGE_TYPES = ['BIAS', 'DARK', 'LAMPFLAT']

LAST_STAGE = {'BIAS': 'banzai.trim.Trimmer',
                       'DARK': 'banzai.bias.BiasSubtractor',
                       'LAMPFLAT': 'banzai.dark.DarkSubtractor'}

EXTRA_STAGES = {'BIAS': ['banzai.bias.BiasMasterLevelSubtractor', 'banzai.bias.BiasComparer'],
                         'DARK': ['banzai.dark.DarkNormalizer', 'banzai.dark.DarkComparer'],
                         'LAMPFLAT': []}


CALIBRATION_STACKER_STAGES = {'BIAS': ['banzai.bias.BiasMaker'],
                              'DARK': ['banzai.dark.DarkMaker'],
                               'LAMPFLAT': ['banzai_nres.flats.FlatStacker']}

# Stack delays are expressed in seconds--namely, each is five minutes
CALIBRATION_STACK_DELAYS = {'BIAS': 300,
                            'DARK': 300,
                            'LAMPFLAT': 300}


SCHEDULE_STACKING_CRON_ENTRIES = {'cpt': {'minute': 0, 'hour': 15},
                                  'tlv': {'minute': 45, 'hour': 8},
                                  'lsc': {'minute': 0, 'hour': 21},
                                  'elp': {'minute': 0, 'hour': 23}}

OBSERVATION_PORTAL_URL = os.getenv('OBSERVATION_PORTAL_URL', 'http://internal-observation-portal.lco.gtn/api/observations/')
CALIBRATE_PROPOSAL_ID = os.getenv('CALIBRATE_PROPOSAL_ID', 'calibrate')

CONFIGDB_URL = os.getenv('CONFIGDB_URL', 'http://configdb.lco.gtn/sites/')

OBSERVATION_REQUEST_TYPES = {'BIAS': 'NRESBIAS', 'DARK': 'NRESDARK'}
