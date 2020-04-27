from pkg_resources import get_distribution, DistributionNotFound
import logging
from banzai.logs import BanzaiLogger

try:
    __version__ = get_distribution('lco-banzai-nres-e2e').version
except DistributionNotFound:
    # package is not installed
    pass

logging.setLoggerClass(BanzaiLogger)
