def config_to_filename(image):
    filename = str(image.configuration_mode)
    filename = filename.replace('nres_full_frame', '')
    filename = filename.replace('full_frame', '')
    filename = filename.replace('default', '')
    filename = filename.replace('central_2k_2x2', 'center')
    return filename
