from __future__ import absolute_import, division, print_function, unicode_literals
from banzai_nres.bias import BiasMaker
import numpy as np
from astropy.io import fits

from banzai.tests.utils import FakeImage, FakeContext, throws_inhomogeneous_set_exception
from banzai.tests.test_bias_maker import FakeBiasImage

import mock

"""
This tests BiasMaker. This is almost an exact copy of the banzai test.bias.
Functions are only explicitly copied because the @mock.patch needs to point to the correct object type.

The only changes are the error margins, which are now fractional errors,
inside of test_makes_a_sensible_master_bias
"""

def test_min_images():
    bias_maker = BiasMaker(None)
    processed_images = bias_maker.do_stage([])
    assert len(processed_images) == 0


@mock.patch('banzai_nres.bias.Image')
def test_header_master_bias_level_returns_1(mock_image):
    maker = BiasMaker(FakeContext())

    maker.do_stage([FakeBiasImage() for x in range(6)])

    args, kwargs = mock_image.call_args
    header = kwargs['header']
    assert header['BIASLVL'] == 1.0


@mock.patch('banzai_nres.bias.Image')
def test_header_master_bias_level_returns_2(mock_image):
    maker = BiasMaker(FakeContext())

    maker.do_stage([FakeBiasImage(image_multiplier=2.0) for x in range(6)])

    args, kwargs = mock_image.call_args
    header = kwargs['header']
    assert header['BIASLVL'] == 2.0


@mock.patch('banzai_nres.bias.Image')
def test_header_cal_type_bias(mock_image):

    maker = BiasMaker(FakeContext())

    maker.do_stage([FakeBiasImage() for x in range(6)])

    args, kwargs = mock_image.call_args
    header = kwargs['header']
    assert header['OBSTYPE'].upper() == 'BIAS'


@mock.patch('banzai_nres.bias.Image')
def test_raises_an_exception_if_ccdsums_are_different(mock_images):
    throws_inhomogeneous_set_exception(BiasMaker, FakeContext(), 'ccdsum', '1 1')


@mock.patch('banzai_nres.bias.Image')
def test_raises_an_exception_if_epochs_are_different(mock_images):
    throws_inhomogeneous_set_exception(BiasMaker, FakeContext(), 'epoch', '20160102')


@mock.patch('banzai_nres.bias.Image')
def test_raises_an_exception_if_nx_are_different(mock_images):
    throws_inhomogeneous_set_exception(BiasMaker, FakeContext(), 'nx', 105)


@mock.patch('banzai_nres.bias.Image')
def test_raises_an_exception_if_ny_are_different(mock_images):
    throws_inhomogeneous_set_exception(BiasMaker, FakeContext(), 'ny', 107)



@mock.patch('banzai_nres.bias.Image')
def test_makes_a_sensible_master_bias(mock_images):
    """
    Actual test for the bias maker. This assumes a uniform bias across the frame.
    Allowed errors were converted from absolute errors (as in Banzai version)
    to fractional errors in this version.
    """
    nimages = 20
    expected_bias = 1183.0
    expected_readnoise = 15.0

    images = [FakeBiasImage() for x in range(nimages)]
    for image in images:
        image.data = np.random.normal(loc=expected_bias, scale=expected_readnoise,
                                      size=(image.ny, image.nx))

    maker = BiasMaker(FakeContext())
    maker.do_stage(images)

    args, kwargs = mock_images.call_args
    master_bias = kwargs['data']
    assert np.abs(np.mean(master_bias)) < 0.1
    actual_bias = float(kwargs['header']['BIASLVL'])
    assert np.abs((actual_bias - expected_bias)/expected_bias) < 1E-3
    actual_readnoise = np.std(master_bias)
    assert np.abs((actual_readnoise - expected_readnoise / (nimages ** 0.5))/actual_readnoise) < 6E-2
    # expected_readnoise / (nimages ** 0.5) is just the theoretical std after averaging
    # of sqr_root(Var((1/n) \sum_n X_n))
