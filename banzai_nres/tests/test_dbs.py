import mock
import numpy as np

from banzai_nres.utils import db_utils
from banzai_nres.settings import NRESSettings

from banzai.tests.utils import FakeContext, FakeImage


@mock.patch('banzai_nres.utils.db_utils.dbs.get_timezone', return_value='UTC')
@mock.patch('banzai_nres.utils.db_utils.date_utils.get_dayobs', return_value='42')
def test_get_raw_path(mockday, mockzone):
    fake_context = FakeContext(settings=NRESSettings())
    fake_context.site = 'site'
    fake_context.instrument_name = 'nres'
    new_path = db_utils.get_raw_path(base_raw_path='base/', pipeline_context=fake_context)
    assert new_path == 'base/site/nres/42/raw'
    assert db_utils.get_raw_path(base_raw_path=None, pipeline_context=fake_context) is None


def test_data_product_instantiates_from_image():
    image = FakeImage()
    possible_attributes = ['obstype', 'dateobs', 'datecreated', 'instrument', 'is_master', 'is_bad']
    image.other_attribute = '42'
    image.attributes = ['other_attribute']
    counter = np.arange(len(possible_attributes))
    for attribute, i in zip(possible_attributes, counter):
        setattr(image, attribute, str(i))
    image_db_info = db_utils.DataProduct(image=image)
    possible_attributes += ['attributes']
    for attribute in possible_attributes:
        assert getattr(image, attribute) == getattr(image_db_info, attribute)
    assert image_db_info.other_attribute == image.other_attribute
