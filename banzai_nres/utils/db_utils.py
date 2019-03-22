import os
import logging

from banzai import dbs
from banzai.utils import date_utils, file_utils
from banzai import logs


logger = logging.getLogger(__name__)


class DataProduct(object):
    def __init__(self, obstype=None, dateobs=None, datecreated=None, is_master=None, is_bad=False,
                 instrument=None, attributes=None, image=None):
        if attributes is None:
            attributes = []
        self.obstype = getattr(image, 'obstype', obstype)
        self.dateobs = getattr(image, 'dateobs', dateobs)
        self.datecreated = getattr(image, 'datecreated', datecreated)
        self.instrument = getattr(image, 'instrument', instrument)
        self.is_master = getattr(image, 'is_master', is_master)
        self.is_bad = getattr(image, 'is_bad', is_bad)
        self.attributes = getattr(image, 'attributes', attributes)
        if hasattr(image, 'attributes'):
            for attribute in image.attributes:
                setattr(self, attribute, getattr(image, attribute))


def get_raw_path(base_raw_path, pipeline_context):
    if base_raw_path is None:
        return None
    timezone = dbs.get_timezone(site=getattr(pipeline_context, 'site', None),
                                db_address=getattr(pipeline_context, 'db_address', None))
    dayobs = date_utils.get_dayobs(timezone=timezone)
    return os.path.join(base_raw_path, pipeline_context.site, pipeline_context.instrument_name, dayobs, 'raw')


def post_to_archive(filepath, image=None):
    logger.info('Posting file to the archive', image=image)
    try:
        file_utils.post_to_archive_queue(filepath)
    except Exception:
        logger.error("Could not post to ingester: {error}".format(error=logs.format_exception()), image=image)
