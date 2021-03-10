# ----------------------------------------------------------------------------
from ._astropy_init import *   # noqa
# ----------------------------------------------------------------------------

import logging
from banzai.logs import BanzaiLogger

logging.setLoggerClass(BanzaiLogger)
