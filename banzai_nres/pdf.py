import numpy as np
import matplotlib.pyplot as pl
from astropy import constants
from astropy import units

from banzai.stages import Stage
from banzai_nres import phoenix
import logging
from banzai_nres.qc.qc_science import get_snr

logger = logging.getLogger('banzai')


class MakePDFSummary(Stage):

    def __init__(self, runtime_context):
        super(MakePDFSummary, self).__init__(runtime_context)

    def do_stage(self, image):
        fig1 = pl.figure(figsize=(11, 8.5))

        fiber = image.spectrum.table['fiber']
        wavelength = image.spectrum.table['wavelength'][fiber == image.science_fiber, :]
        flux = image.spectrum.table['normflux'][fiber == image.science_fiber, :]
        order = image.spectrum.table['order'][fiber == image.science_fiber]

        phoenix_loader = phoenix.PhoenixModelLoader(self.runtime_context.db_address)
        template = phoenix_loader.load(image.classification)
        v_over_c_plus_one = 1 + image.meta['RV']/constants.c.to(units.m/units.s).value

        # make the first page showing the spectrum, template, etc.
        pl.subplot(2, 1, 1)
        primary_order = order == 90
        spectrum_line, = pl.plot(np.squeeze(wavelength[primary_order, :]),
                                 np.squeeze(flux[primary_order, :]), color='blue')
        template_line, = pl.plot(template['wavelength']*v_over_c_plus_one, template['flux'], color='red', linewidth=0.5)
        pl.xlim([5140., 5220.])
        pl.ylim([0., np.max(flux[primary_order, :])])
        pl.xlabel('wavelength (Angstroms)')
        pl.ylabel('normalized flux')
        pl.title(image.meta['OBJECT'] + ' on ' + image.meta['DAY-OBS'] + ' from ' + image.meta['TELESCOP'] + '-' +
                 image.meta['SITEID'] + ' for program ' + image.meta['PROPID'])
        pl.legend((spectrum_line, template_line), ('Spectrum', 'Best-Fit Phoenix Template'))

        pl.subplot(2, 3, 4)
        stacked_ccf = np.ones_like(image.ccf['v'][0, :])
        pl.plot([image.meta['RV'] / 1000., image.meta['RV'] / 1000.], [-0.1, 1.1], color='blue')
        for ccf in image.ccf:
            this_ccf = ccf['xcor'] - np.min(ccf['xcor'])
            pl.plot(ccf['v'] + image.meta['BARYCORR'] / 1000., this_ccf / np.max(this_ccf), color='gray', alpha=0.25)
            stacked_ccf *= ccf['xcor']
        stacked_ccf -= np.min(stacked_ccf)
        pl.plot(image.ccf['v'][0, :] + image.meta['BARYCORR'] / 1000., stacked_ccf / np.max(stacked_ccf), color='black')
        pl.xlim([image.meta['RV'] / 1000. - 50, image.meta['RV'] / 1000. + 50])
        pl.ylim([-0.1, 1.1])
        pl.xlabel('barycentric velocity (km/s)')
        pl.ylabel('normalized CCF')
        pl.title('CCF')

        pl.subplot(2, 3, 5)
        snr, snr_wavelength = np.zeros(len(order), dtype=np.float), np.zeros(len(order), dtype=np.float)
        for snr_order in order:
            snr[snr_order - order[0]], snr_wavelength[snr_order - order[0]] = get_snr(image, snr_order)
        pl.plot(snr_wavelength, snr, 'ko')
        rv_orders = (order >= self.runtime_context.MIN_ORDER_TO_CORRELATE) & \
                    (order <= self.runtime_context.MAX_ORDER_TO_CORRELATE)
        pl.plot(snr_wavelength[rv_orders], snr[rv_orders], 'ro')
        pl.xlabel('wavelength (Angstroms)')
        pl.ylabel('peak SNR/pixel in order')
        pl.title('SNR vs wavelength')

        ax = pl.subplot(2, 3, 6)
        ax.set_axis_off()
        pl.title('Summary Information for file')
        line_separation, top_line = 0.065, 0.925
        pl.text(0.1, top_line, image.meta['ORIGNAME'].replace('e00', 'e' + str(self.runtime_context.reduction_level) +
                                                              '-1d') + '.fz')
        pl.text(0.1, top_line - line_separation, 'Teff = {0:1.4g} K'.format(image.meta['TEFF']))
        pl.text(0.1, top_line - line_separation * 2, 'logg = {0:1.2g} (cgs units)'.format(image.meta['LOGG']))
        pl.text(0.1, top_line - line_separation * 3, '[Fe/H] = {0:1.2g}'.format(image.meta['FEH']))
        pl.text(0.1, top_line - line_separation * 4, '[alpha/Fe] = {0:1.2g}'.format(image.meta['ALPHA']))

        pl.text(0.1, top_line - line_separation * 6, 'RV = {0:1.3f} km/s'.format(image.meta['RV'] / 1000.))
        pl.text(0.1, top_line - line_separation * 7, 'RV error = {0:1.3f} km/s'.format(image.meta['RVERR'] / 1000.))
        pl.text(0.1, top_line - line_separation * 8, 'Barycorr = {0:1.3f} km/s'.format(image.meta['BARYCORR'] / 1000.))
        pl.text(0.1, top_line - line_separation * 9, 'BJD_TDB = {0:1.5f}'.format(image.meta['TCORR']))

        pl.text(0.1, top_line - line_separation * 11, 'SNR = {0:1.0f}/pixel @ 5180 Angstroms'.format(image.meta['SNR']))
        pl.text(0.1, top_line - line_separation * 12, 'Exposure time = {0:1.0f} seconds'.format(image.meta['EXPTIME']))

        # Next Page

        # make the plots of individual lines of interest
        line_centers = np.array([3969.63, 4862.764, [5185.14, 5174.22, 5168.84],
                                [5891.68, 5897.65], 6564.73, [6709.73, 6709.88]], dtype=object)
        line_names = np.array(['Ca II H', 'H beta', 'Mg b', 'Na D', 'H alpha', 'Li'])
        line_orders = np.array([117, 96, 90, 79, 71, 70])

        fig2, axes = pl.subplots(nrows=2, ncols=3, figsize=(11, 8.5))

        for ax, line_center, line_name, line_order in zip(axes.flatten(), line_centers, line_names, line_orders):
            make_line_plot(wavelength, flux, order, ax, line_center, line_name,
                           line_order, wavelength_correction=v_over_c_plus_one)

        # Next Page

        # make the plots of the telluric lines
        # NOTE: these are approximate wavelengths, need to look up the actual wavelengths later!
        line_centers = np.array([[6881.81, 6885.73], [7618.23, 7623.09]], dtype=object)
        line_names = np.array(['Telluric B-band', 'Telluric A-band'])
        line_orders = np.array([68, 61])

        fig3, axes = pl.subplots(nrows=2, ncols=1, figsize=(11, 8.5))

        for ax, line_center, line_name, line_order in zip(axes.flatten(), line_centers, line_names, line_orders):
            make_line_plot(wavelength, flux, order, ax, line_center, line_name, line_order)

        # If there is a working exposure meter, a fourth page showing that can go here.
        import pdb; pdb.set_trace()
        image.summary_figures = [fig1, fig2, fig3]
        return image


def make_line_plot(wavelength, flux, order, ax, line_center, line_name, line_order, wavelength_correction=1.0):
    this_order = order == line_order
    non_zero = np.squeeze(wavelength[this_order, :]) != 0.
    ax.plot(np.squeeze(wavelength[this_order, non_zero]), np.squeeze(flux[this_order, non_zero]), color='black')
    if isinstance(line_center, (list, tuple, np.ndarray)):
        for this_center in line_center:
            this_center *= wavelength_correction
            ax.plot([this_center, this_center], [0, 1.1], color='blue')
    else:
        line_center *= wavelength_correction
        ax.plot([line_center, line_center], [0, 1.1], color='blue')
    ax.set_xlabel('wavelength (Angstroms)')
    ax.set_ylabel('normalized flux')
    ax.set_title(line_name)
    ax.set_xlim([np.min(line_center) - 7.5, np.max(line_center) + 7.5])
    ax.set_ylim([0., 1.1])
