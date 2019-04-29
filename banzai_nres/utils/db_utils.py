import os
import logging

from banzai import dbs
from banzai.utils import date_utils, file_utils
from banzai import logs


logger = logging.getLogger(__name__)


def get_raw_path(base_raw_path, runtime_context):
    if base_raw_path is None:
        return None
    timezone = dbs.get_timezone(site=getattr(runtime_context, 'site', None),
                                db_address=getattr(runtime_context, 'db_address', None))
    dayobs = date_utils.get_dayobs(timezone=timezone)
    return os.path.join(base_raw_path, runtime_context.site, runtime_context.instrument_name, dayobs, 'raw')


def post_to_archive(filepath, image=None):
    logger.info('Posting file to the archive', image=image)
    try:
        file_utils.post_to_archive_queue(filepath)
    except Exception:
        logger.error("Could not post to ingester: {error}".format(error=logs.format_exception()), image=image)
