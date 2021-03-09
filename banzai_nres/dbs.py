from banzai.dbs import Base, create_db, add_or_update_record
from sqlalchemy import Column, String, Integer, Float, Index
import boto3
import banzai.dbs
import os
from glob import glob
import logging
from sqlalchemy import func, desc
from sqlalchemy.ext.hybrid import hybrid_method
import numpy as np

from banzai_nres.utils.phoenix_utils import parse_phoenix_header

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
    radius = Column(Float)
    luminosity = Column(Float)
    mass = Column(Float)
    Index('idx_model', 'T_effective', 'log_g', 'metallicity', 'alpha')
    Index('idx_observer', 'T_effective', 'luminosity')

    @hybrid_method
    def diff_T(self, value):
        return abs(self.T_effective - value)

    @diff_T.expression
    def diff_T(self, value):
        return func.abs(self.T_effective - value)

    @hybrid_method
    def diff_log_g(self, value):
        return abs(self.log_g - value)

    @diff_log_g.expression
    def diff_log_g(self, value):
        return func.abs(self.log_g - value)

    @hybrid_method
    def diff_alpha(self, value):
        return abs(self.alpha - value)

    @diff_alpha.expression
    def diff_alpha(self, value):
        return func.abs(self.alpha - value)

    @hybrid_method
    def diff_metallicity(self, value):
        return abs(self.metallicity - value)

    @diff_metallicity.expression
    def diff_metallicity(self, value):
        return func.abs(self.metallicity - value)

    @hybrid_method
    def diff_luminosity(self, value):
        return abs(self.luminosity - value)

    @diff_luminosity.expression
    def diff_luminosity(self, value):
        return func.abs(self.luminosity - value)


# We define the great circle distance here instead of using astropy because we need it to work inside the db.
def cos_great_circle_distance(sin_ra1, cos_ra1, sin_dec1, cos_dec1, sin_ra2, cos_ra2, sin_dec2, cos_dec2):
    """

    :param sin_ra1: sin(ra1)
    :param cos_ra1: cos(ra1)
    :param sin_dec1: sin(dec1)
    :param cos_dec1: cos(dec1)
    :param sin_ra2: sin(ra2)
    :param cos_ra2: cos(ra2)
    :param sin_dec2: sin(dec2)
    :param cos_dec2: cos(dec2)
    :return: cos(D) where D is the great circle distance

    This is the standard great circle distance from e.g. https://mathworld.wolfram.com/GreatCircle.html
    The only difference is we also use the identity for cos(x1 - x2)
    (e.g. https://mathworld.wolfram.com/TrigonometricAdditionFormulas.html)
    so that we can calculate the sin and cos terms ahead of time.
    """
    cos_distance = sin_dec1 * sin_dec2 + cos_dec1 * cos_dec2 * (cos_ra1 * cos_ra2 + sin_ra1 * sin_ra2)
    return cos_distance


class Classification(Base):
    __tablename__ = 'classifications'
    id = Column(Integer, primary_key=True, autoincrement=True)
    ra = Column(Float)
    dec = Column(Float)
    sin_ra = Column(Float)
    cos_ra = Column(Float)
    sin_dec = Column(Float)
    cos_dec = Column(Float)
    T_effective = Column(Float)
    log_g = Column(Float)
    metallicity = Column(Float)
    alpha = Column(Float)
    Index('idx_radec', "ra", "dec")
    Index('idx_radectrig', "sin_ra", "cos_ra", "sin_dec", "cos_dec")

    @hybrid_method
    def cos_distance(self, sin_ra, cos_ra, sin_dec, cos_dec):
        return cos_great_circle_distance(sin_ra, cos_ra, sin_dec, cos_dec,
                                         self.sin_ra, self.cos_ra, self.sin_dec, self.cos_dec)

    @cos_distance.expression
    def cos_distance(cls, sin_ra, cos_ra, sin_dec, cos_dec):
        return cos_great_circle_distance(sin_ra, cos_ra, sin_dec, cos_dec,
                                         cls.sin_ra, cls.cos_ra, cls.sin_dec, cls.cos_dec)


class ResourceFile(Base):
    __tablename__ = 'resourcefiles'
    id = Column(Integer, primary_key=True, autoincrement=True)
    key = Column(String(100), unique=True, index=True)
    filename = Column(String(100), unique=True)
    location = Column(String(150))


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
        for model_file in model_files:
            # strip off the s3
            if 's3' in model_location:
                filename = model_file.key
                location = model_location
            else:
                filename = os.path.basename(model_file)
                location = os.path.dirname(model_file)

            if 'wave' in filename.lower():
                banzai.dbs.add_or_update_record(db_session, ResourceFile, {'key': 'phoenix_wavelengths'},
                                                {'filename': filename, 'location': location,
                                                 'key': 'phoenix_wavelengths'})
                continue

            # Note that there are 25 header keyword value pairs in a standard phoenix model file all in cgs units
            # FITS headers have 80 character header lines
            if 's3' in model_location:
                header_lines = model_bucket.Object(key=filename).get(Range=f'bytes=0-{80 * 25 - 1}')['Body'].read()
            else:
                with open(os.path.join(location, filename), 'rb') as f:
                    header_lines = f.read(80 * 25)
            T_effective, log_g, metallicity, alpha, radius, luminosity, mass = parse_phoenix_header(header_lines)
            equivalence_criteria = {'filename': filename}
            # This naming convention assumes we follow the convention from the munge_phoenix_models.py code.
            record_attributes = {'filename': filename,
                                 'location': location,
                                 'T_effective': T_effective,
                                 'log_g': log_g,
                                 'metallicity': metallicity,
                                 'alpha': alpha,
                                 'radius': radius,
                                 'luminosity': luminosity,
                                 'mass': mass}
            logger.info('Loading Phoenix Model', extra_tags=record_attributes)
            banzai.dbs.add_or_update_record(db_session, PhoenixModel, equivalence_criteria, record_attributes)
            db_session.commit()


def get_closest_phoenix_models(db_address, T_effective, log_g, metallicity=0.0, alpha=0.0, fixed=None, n=1):
    if fixed is None:
        fixed = []
    with banzai.dbs.get_session(db_address=db_address) as db_session:
        query = []
        for param in fixed:
            query.append(getattr(PhoenixModel, param) == eval(param))
        order = [PhoenixModel.diff_T(T_effective), PhoenixModel.diff_log_g(log_g),
                 PhoenixModel.diff_metallicity(metallicity), PhoenixModel.diff_alpha(alpha)]
        model = db_session.query(PhoenixModel).filter(*query).order_by(*order).limit(n).all()
        if model is None:
            logger.error('Phoenix model does not exist for these parameters',
                         extra_tags={'T_effective': T_effective, 'log_g': log_g,
                                     'metallicity': metallicity, 'alpha': alpha})
            raise ValueError('Phoenix Model Missing')
    return model


def get_closest_HR_phoenix_models(db_address, T_effective, luminosity, metallicity=0.0, alpha=0.0):
    with banzai.dbs.get_session(db_address=db_address) as db_session:
        order = [PhoenixModel.diff_T(T_effective), PhoenixModel.diff_luminosity(luminosity),
                 PhoenixModel.diff_metallicity(metallicity), PhoenixModel.diff_alpha(alpha)]
        model = db_session.query(PhoenixModel).order_by(*order).first()
        if model is None:
            logger.error('Phoenix model does not exist for these parameters',
                         extra_tags={'T_effective': T_effective, 'Luminosity': luminosity,
                                     'metallicity': metallicity, 'alpha': alpha})
            raise ValueError('Phoenix Model Missing')
    return model


def get_resource_file(db_address, key):
    with banzai.dbs.get_session(db_address=db_address) as db_session:
        resource_file = db_session.query(ResourceFile).filter(ResourceFile.key == key).first()
    return resource_file


def get_closest_existing_classification(db_address, ra, dec):
    with banzai.dbs.get_session(db_address=db_address) as db_session:
        # Note the desc here. Because sqlite does not have trig functions, we can't take an arc cos. So we need the
        # value when the cos is maximum (which is theta = minimum)
        order = [desc(Classification.cos_distance(np.sin(np.deg2rad(ra)), np.cos(np.deg2rad(ra)),
                                                  np.sin(np.deg2rad(dec)), np.cos(np.deg2rad(dec))))]
        model = db_session.query(Classification).order_by(*order).first()
    return model


def save_classification(db_address, frame):
    with banzai.dbs.get_session(db_address=db_address) as db_session:
        equivalence_criteria = {'ra': frame.ra, 'dec': frame.dec}
        record_attributes = {'ra': frame.ra, 'dec': frame.dec, 'sin_ra': np.sin(np.deg2rad(frame.ra)),
                             'cos_ra': np.cos(np.deg2rad(frame.ra)), 'sin_dec': np.sin(np.deg2rad(frame.dec)),
                             'cos_dec': np.cos(np.deg2rad(frame.dec)),
                             'T_effective': frame.classification.T_effective,
                             'log_g': frame.classification.log_g, 'metallicity': frame.classification.metallicity,
                             'alpha': frame.classification.alpha}
        add_or_update_record(db_session, Classification, equivalence_criteria, record_attributes)
