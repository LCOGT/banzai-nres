from banzai.dbs import CalibrationImage
from sqlalchemy import Column, SmallInteger


class NRESCalibrationImage(CalibrationImage):
    fibers_state = Column(SmallInteger)
