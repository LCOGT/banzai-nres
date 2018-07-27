import os


class TestContext(object):
    """
    Picks out a frame or a set of frames to test. Provides the appropriate
    PipelineContext to pass an NRES test frame to the banzai modules.
    Parameters
    ----------
    filename: None if you want just the path to be included in the context (only frames with OBSTYPE = 'BIAS' are used
    Returns
    -------
    stages_todo: list of banzai.stages.Stage
                 The stages that need to be done
    """
    def __init__(self, filename=None, raw_path='/archive/engineering/lsc/nres01/20180228/raw'):
        self.processed_path = '/tmp'
        self.raw_path = raw_path
        self.filename = filename
        self.post_to_archive = False
        self.db_address = os.environ['DB_URL']
        self.preview_mode = False
        self.rlevel = 0
        self.fpack = True
