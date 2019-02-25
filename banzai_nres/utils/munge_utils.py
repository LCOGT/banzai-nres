def get_telescope_filename(image):
    return image.header.get('TELESCOP', '').replace('nres', 'nrs')
