from banzai.dbs import Base, create_db
from sqlalchemy import Column, String, Integer, Float, Index
import boto3
import banzai.dbs
import os
from glob import glob
import logging


logger = logging.getLogger('banzai')


class PhoenixModel(Base):
    __tablename__ = 'phoenixmodels'
    id = Column(Integer, primary_key=True, autoincrement=True)
    filename = Column(String(100), unique=True)
    location = Column(String(150))
    T_effective = Column(Float)
    log_g = Column(Float)
    metallicity = Column(Float)
    alpha = Column(Float)
    Index('idx_model', 'T_effective', 'log_g', 'metallicity', 'alpha')


def create_nres_db(db_address):
    create_db(db_address)


def populate_phoenix_models(model_location, db_address):
    if 's3' in model_location:
        s3 = boto3.resource('s3')
        model_bucket = s3.Bucket(model_location.replace('s3://', ''))
        model_files = model_bucket.objects.all()
    else:
        # Assume they are on disk
        model_files = glob(os.path.join(model_location, '*.fits'))
    with banzai.dbs.get_session(db_address) as db_session:
        # strip off the s3
        for model_file in model_files:
            if 's3' in model_location:
                filename = model_file.key
                location = model_location
            else:
                filename = os.path.basename(model_file)
                location = os.path.dirname(model_file)

            if 'wave' in filename.lower():
                continue
            equivalence_criteria = {'filename': filename}
            # This naming convention assumes we follow the convention from the munge_phoenix_models.py code.
            _, T_effective, log_g, metallicity, alpha = os.path.splitext(filename)[0].split('_')
            record_attributes = {'filename': filename,
                                 'location': location,
                                 'T_effective': T_effective,
                                 'log_g': log_g,
                                 'metallicity': metallicity,
                                 'alpha': alpha}
            banzai.dbs.add_or_update_record(db_session, PhoenixModel, equivalence_criteria, record_attributes)
            db_session.commit()


def get_phoenix_model_record(db_address, T_effective, log_g, metallicity, alpha):
    with banzai.dbs.get_session(db_address=db_address) as db_session:
        query = PhoenixModel.T_effective == T_effective
        query &= PhoenixModel.log_g == log_g
        query &= PhoenixModel.metallicity == metallicity
        query &= PhoenixModel.alpha == alpha
        model = db_session.query(PhoenixModel).filter(query).first()
        if model is None:
            logger.error('Phoenix model does not exist for these parameters',
                         extra_tags={'T_effective': T_effective, 'log_g': log_g,
                                     'metallicity': metallicity, 'alpha': alpha})
            raise ValueError('Phoenix Model Missing')
    return model
