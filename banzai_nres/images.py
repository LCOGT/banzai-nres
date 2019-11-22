from typing import Union

from banzai_nres.fibers import fiber_states_from_header
from banzai.images import ObservationFrame, LCOObservationFrame, LCOCalibrationFrame, CCDData, Section, HeaderOnly, Data
from banzai.utils import fits_utils
from banzai import dbs
import numpy as np
import logging
import os
from astropy.io import fits

logger = logging.getLogger('banzai')


class NRESCalibrationFrame(LCOCalibrationFrame):
    def __init__(self, hdu_list: list, file_path: str, grouping_criteria: list = None):
        super().__init__(hdu_list, file_path, grouping_criteria)
        self.fiber0_lit, self.fiber1_lit, self.fiber2_lit = fiber_states_from_header(self.meta)

    def num_lit_fibers(self):
        return 1 * self.fiber0_lit + 1 * self.fiber1_lit + 1 * self.fiber2_lit


class DataTable(Data):
    def __init__(self, data, meta, name):
        super().__init__(data, meta, name=name)

    def to_fits(self) -> Union[fits.HDUList, list]:
        return [fits.BinTableHDU(data=self.data, header=self.meta)]


class ArrayData(Data):
    def __init__(self, data, meta, name):
        super().__init__(data, meta, name=name)

    def to_fits(self) -> Union[fits.HDUList, list]:
        return [fits.ImageHDU(data=self.data, header=self.meta)]


class NRESFrameFactory:
    @classmethod
    def open(cls, path, runtime_context) -> ObservationFrame:
        fits_hdu_list = fits_utils.open_fits_file(path)
        hdu_list = []
        for hdu in fits_hdu_list:
            # Move on from the BPM and ERROR arrays
            if 'BPM' in hdu.header.get('EXTNAME', '') or 'ERR' in hdu.header.get('EXTNAME', ''):
                continue
            if hdu.data is None:
                hdu_list.append(HeaderOnly(meta=hdu.header))
            elif isinstance(hdu, fits.BinTableHDU):
                hdu_list.append(DataTable(data=hdu.data, meta=hdu.header, name=hdu.header.get('EXTNAME')))
            elif 'GAIN' in hdu.header:
                condensed_name = hdu.header.get('EXTNAME', '')
                for extension_name_to_condense in runtime_context.EXTENSION_NAMES_TO_CONDENSE:
                    condensed_name = condensed_name.replace(extension_name_to_condense, '')
                if (condensed_name + 'BPM', hdu.header.get('EXTVER')) in fits_hdu_list:
                    bpm_array = fits_hdu_list[condensed_name + 'BPM', hdu.header.get('EXTVER')].data
                else:
                    bpm_array = None

                if (condensed_name + 'ERR', hdu.header.get('EXTVER')) in fits_hdu_list:
                    error_array = fits_hdu_list[condensed_name + 'ERR', hdu.header.get('EXTVER')].data
                else:
                    error_array = None
                hdu_list.append(CCDData(data=hdu.data.astype(np.float32), meta=hdu.header, name=hdu.header.get('EXTNAME'),
                                        mask=bpm_array, uncertainty=error_array))
            else:
                hdu_list.append(ArrayData(data=hdu.data, meta=hdu.header, name=hdu.header.get('EXTNAME')))
        if hdu_list[0].meta.get('OBSTYPE') in runtime_context.CALIBRATION_IMAGE_TYPES:
            grouping = runtime_context.CALIBRATION_SET_CRITERIA.get(hdu_list[0].meta.get('OBSTYPE'), [])
            image = NRESCalibrationFrame(hdu_list, os.path.basename(path), grouping_criteria=grouping)
        else:
            image = LCOObservationFrame(hdu_list, os.path.basename(path))
        image.instrument = cls._get_instrument(image, runtime_context.db_address)

        # TODO: Put all munge code here

        for hdu in image.ccd_hdus:
            if hdu.meta.get('DETSEC', 'UNKNOWN') in ['UNKNOWN', 'N/A']:
                # DETSEC missing?
                binning = hdu.meta.get('CCDSUM', image.primary_hdu.meta.get('CCDSUM', '1 1'))
                data_section = Section.parse_region_keyword(hdu.meta['DATASEC'])
                detector_section = Section(1,
                                           max(data_section.x_start, data_section.x_stop) * int(binning[0]),
                                           1,
                                           max(data_section.y_start, data_section.y_stop) * int(binning[2]))
                hdu.meta['DETSEC'] = detector_section.to_region_keyword()

            # SATURATE Missing?
            def update_saturate(image, hdu, default):
                if hdu.meta.get('SATURATE', 0.0) == 0.0:
                    hdu.meta['SATURATE'] = image.meta.get('SATURATE', 0.0)
                    hdu.meta['MAXLIN'] = image.meta.get('MAXLIN', 0.0)
                if hdu.meta.get('SATURATE', 0.0) == 0.0:
                    hdu.meta['SATURATE'] = (default, '[ADU] Saturation level used')
                    hdu.meta['MAXLIN'] = (default, '[ADU] Non-linearity level')
            if 'sinistro' in image.instrument.type.lower():
                update_saturate(image, hdu, 47500.0)

            elif '1m0' in image.instrument.type:
                # Saturation level from ORAC Pipeline
                update_saturate(image, hdu, 46000.0)
            elif '0m4' in image.instrument.type or '0m8' in image.instrument.type:
                # Measured by Daniel Harbeck
                update_saturate(image, hdu, 64000.0)

            elif 'spectral' in image.instrument.type.lower():
                # These values were given by Joe Tufts on 2016-06-07
                binning = hdu.meta.get('CCDSUM', '1 1')
                n_binned_pixels = int(binning[0]) * int(binning[2])
                update_saturate(image, hdu, 125000.0 * n_binned_pixels / float(hdu.meta['GAIN']))
        return image

    @classmethod
    def _get_instrument(cls, image, db_address):
        site = image.meta.get('SITEID')
        camera = image.meta.get('INSTRUME')
        instrument = dbs.query_for_instrument(db_address, site, camera)
        name = camera
        if instrument is None:
            # if instrument is missing, assume it is an NRES frame and check for the instrument again.
            name = image.meta.get('TELESCOP')
            instrument = dbs.query_for_instrument(db_address, site, camera, name=name)
        if instrument is None:
            msg = 'Instrument is not in the database, Please add it before reducing this data.'
            tags = {'site': site, 'camera': camera, 'telescop': name}
            logger.error(msg, extra_tags=tags)
            raise ValueError('Instrument is missing from the database.')
        return instrument
