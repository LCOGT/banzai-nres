from astropy.io import fits
from .trace_utils import get_trace_centroids_from_coefficients

import datetime
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from datetime import datetime

import os
import glob
import numpy as np


class Trace(object):
    """
    Object for storing all the trace related attributes. This gets appended to each Image instance.
    """
    def __init__(self):
        self.coefficients = None
        self.fiber_order = None
        self.high_signal_to_noise_region_bounds = None
        self.has_sufficient_signal_to_noise = None


class Image(object):
    """
    stripped down Banzai image class for testing.
    """
    def __init__(self, pipeline_context, filename=None, data=None, header=None):
        self.filename = filename
        self.header = header
        self.data = data
        if filename is not None:
            self.data = fits.getdata(self.filename)
            self.header = fits.getheader(filename, 1)
        self.pipeline_context = pipeline_context

        """
        new attributes necessary for NRES
        """
        self.trace = Trace()
        self.coordinates = None
        self.ivar = None
        self.fiber_profile = None


class FakeImage(object):
    def __init__(self):
        self.telescope_id = -1
        self.site = 'lsc'
        self.pipeline_context = None
        self.instrument = 'nres01'
        self.filename = 'test.fits'
        self.filter = 'U'
        self.dateobs = datetime(2018, 8, 7)
        self.header = {}
        self.caltype = ''
        self.request_number = '0000331403'
        self.readnoise = 11
        self.block_id = '254478983'
        self.molecule_id = '544562351'
        self.exptime = 30.0
        self.data = np.zeros((500, 500))
        self.trace = Trace()
        self.coordinates = Coordinates()
        self.ny, self.nx = self.data.shape


class CalibrationMaker(object):
    def __init__(self, pipeline_context=None):
        self.pipeline_context = pipeline_context


class Stage(object):
    def __init__(self, pipeline_context=None):
        self.pipeline_context = pipeline_context


class PipelineContext(object):
    def __init__(self, dummy='dummy'):
        self.dummy = dummy


def create_master_calibration_header(images):
    return None


def get_trace_fit_coefficients(image):
    tracecoefficients, fiber_order_tuple = None, None
    if image.header['OBSTYPE'] != 'TRACES':
        observation_date = image.header.get('OBS-DATE')
        fiber_order = image.header.get('FIBRORDR')
        file_path = os.path.join(os.path.dirname(image.filename), 'coefficients/')
        coefficients_file = glob.glob(os.path.join(file_path, 'lscnrs01-fl09-20180726-0025-w00.txt'))[0]
        tracecoefficients = np.genfromtxt(coefficients_file, delimiter="\t", skip_header=2)
        # fiber_order = from the fits file
        fiber_order = '(1,2)'
        fiber_order_tuple = tuple((int(fiber_order[1]), int(fiber_order[3])))

    return tracecoefficients, fiber_order_tuple


def read_images(image_list):
    """
    :param image_list: list of fits image filenames
    :return: list of Image objects which are overscan trimmed
    """
    pipeline_context = PipelineContext()
    images = []
    for filename in image_list:
        image = Image(pipeline_context, filename=filename)
        image.data = overscan_trim(image.data.astype(np.float64))
        images.append(image)
    return images


def overscan_trim(image):
    return image[0:4096, 0:4096] - np.median(image[:, 4096:])


def get_time_stamp(filename):
    """
    :param filename: fits.fz file.
    :return: datetime object, the time stripped from the fits observation date.
    """
    date_string = fits.open(filename)[1].header['DATE-OBS']
    return datetime.datetime.strptime(date_string, "%Y-%m-%dT%H:%M:%S.%f")


def plotter(images):
    for image in images:
        trace_values_versus_xpixel, numtraces, x = get_trace_centroids_from_coefficients(image.trace.coefficients, image)
        plt.figure()
        for i in range(trace_values_versus_xpixel.shape[0]):
            plt.plot(x, trace_values_versus_xpixel[i], 'r')
        minval = 10
        # setting all values less than minval to be equal to minval.
        im_nonneg = image.data * (image.data > minval) + minval * (image.data <= minval)
        plt.imshow(im_nonneg, cmap='gray', norm=LogNorm(vmin=None, vmax=None), interpolation='nearest')
        plt.show()


def load_master_profile():
    return fits.getdata('/home/mbrandt21/Downloads/fiber_profiles/profile_image.fits.fz') # plus header
