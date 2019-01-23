def make_calibration_filename_function(calibration_type, attribute_filename_functions):
    def get_calibration_filename(image):
        name_components = {'site': image.site, 'telescop': image.header.get('TELESCOP', '').replace('nres', 'nrs'),
                           'camera': image.header.get('INSTRUME', ''), 'epoch': image.epoch,
                           'cal_type': calibration_type.lower()}
        cal_file = '{site}{telescop}-{camera}-{epoch}-{cal_type}'.format(**name_components)
        for filename_function in attribute_filename_functions:
            cal_file += '-{}'.format(filename_function(image))
        cal_file += '.fits'
        return cal_file
    return get_calibration_filename
