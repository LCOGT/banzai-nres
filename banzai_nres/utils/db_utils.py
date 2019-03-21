import os

from banzai import dbs
from banzai.utils import date_utils


def get_raw_path(base_raw_path, pipeline_context):
    if base_raw_path is None:
        return None
    timezone = dbs.get_timezone(pipeline_context.site, db_address=pipeline_context.db_address)
    dayobs = date_utils.get_dayobs(timezone=timezone)
    return os.path.join(base_raw_path, pipeline_context.site, pipeline_context.instrument_name, dayobs, 'raw')
