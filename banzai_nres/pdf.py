import numpy as np
import matplotlib.pyplot as pl
from matplotlib.backends.backend_pdf import PdfPages

from banzai_nres.qc.qc_science import get_snr

class MakePDFSummary(Stage):
    def do_stage(self, image):
        pl.figure(figsize=(11,8.5)) #inches

        #NOTE: need to change the following to use the e92-summary.pdf format we've agreed on
        pdf_filename = image.meta['RAWIMAGE']+'.pdf'
        pp = PdfPages(pdf_filename)

        wave, flux, order = image.spectrum[image.science_fiber, :]['wave'], image.spectrum[image.science_fiber, :]['normflux'], 
                            image.spectrum[image.science_fiber]['order']

        template = phoenix_loader.load(image.classification)
        ccfs=hdu['CCF'].data

        #make the first page showing the spectrum, template, etc.
        pl.subplot(2, 1, 1)
        primary_order = order == 90
        pl.plot(wave[primary_order,:], flux[primary_order,:], color='blue')
        pl.plot(template['wavelength'], template['flux'], color='red')
        pl.xlim([np.min(wave[primary_order,:]), np.max(wave[primary_order,:])])
        pl.xlabel('wavelength (Angstroms)')
        pl.ylabel('normalized flux')
        pl.title(image.meta['OBJECT']+', '+image.meta['RAWIMAGE'])

        pl.subplot(2, 3, 4)
        stacked_ccf=np.ones_like(ccfs['v'][0,:])
        for ccf in ccfs:
            pl.plot(ccf['v'], ccf['xcor']/np.max(ccf['xcor']), color=gray, alpha=0.25)
            stacked_ccf *= ccf['xcor']
        pl.plot(ccfs['v'][0,:], stacked_ccf/np.max(stacked_ccf), color='black')
        pl.xlim([-50,50])
        pl.xlabel('velocity (km/s)')
        pl.title('CCF')

        pl.subplot(2, 3, 5)
        snr, snr_wave = np.array([len(order)]), np.array([len(order)])
        for snr_order in order:
            snr[snr_order-order[0]], snr_wave[snr_order-order[0]] = get_snr(image, snr_order)
        pl.plot(snr_wave, snr_order)
        pl.xlabel('wavelength (Angstroms)')
        pl.ylabel('peak SNR')
        pl.title('SNR vs wavelength')

        ax=pl.subplot(2, 3, 6)
        ax.set_axis_off()
        pl.title('Estimated Stellar Parameters')
        pl.text(0.1,0.9,'Teff = '+str(image.meta['TEFF'])+' K')
        pl.text(0.1,0.8,'logg = '+str(image.meta['LOGG'])+' (cgs units)')
        pl.text(0.1,0.7,'[Fe/H] = '+str(image.meta['FEH']))
        pl.text(0.1,0.6,'[alpha/Fe] = '+str(image.meta['ALPHA']))

        pl.text(0.1,0.4,'RV = '+str(image.meta['RV']/1000.)+' km/s')
        pl.text(0.1,0.3,'RV error = '+str(image.meta['RVERR']/1000.)+' km/s')
        pl.text(0.1,0.2,'Barycorr = '+str(image.meta['BARYCORR']/1000.)+' km/s')

        pl.text(0.1,0.0,'SNR = '+str(image.meta['SNR']))

        next_page(pp)

        #make the plots of individual lines of interest
        #NOTE: the following wavelengths are currently AIR wavelengths, need to convert to VACUUM
        line_centers = np.array([3968.47, 4861.34, [5183.62, 5172.70, 5167.33], [5889.95, 5895.92], 6562.81, [6707.76, 6707.91]], dtype=object)
        line_names = np.array(['Ca II H', 'H beta', 'Mg b', 'Na D', 'H alpha', 'Li'])
        line_orders = np.array([117, 95, 90, 79, 71, 70])

        fig, axes = pl.subplots(nrows=2,ncols=3)

        for ax, line_center, line_name, line_order in zip(axes.flatten(), line_centers, line_names, line_orders):
            make_line_plot(wave, flux, order, ax, line_center, line_name, line_order)
        next_page(pp)

        #make the plots of the telluric lines
        #NOTE: these are approximate wavelengths, need to look up the actual wavelengths later!
        line_centers = ([[6881.81, 6885.73], [7618.23, 7623.09]], dtype = object)
        line_names = np.array(['Telluric B-band', 'Telluric A-band'])
        line_orders = np.array([68, 61])

        fig, axes = pl.subplots(nrows=2,ncols=1)

        for ax, line_center, line_name, line_order in zip(axes.flatten(), line_centers, line_names, line_orders):
            make_line_plot(wave, flux, order, ax, line_center, line_name, line_order)
        next_page(pp)

        pp.close()



def make_line_plot(wave, flux, order, ax, line_center, line_name, line_order):
    this_order = order == line_order
    ax.plot(wave[this_order,:], flux[this_order,:], color='black')
    for this_center in line_center:
        ax.plot([this_center, this_center], [0, 1], color='blue')
    ax.xlabel('wavelength (Angstroms)')
    ax.ylabel('normalized flux')
    ax.title(line_name)
    ax.xlim([np.min(line_center)-7.5, np.max(line_center)+7.5])

def next_page(pp):
    pl.tight_layout()
    pp.savefig(format='pdf')
    pp.clf()