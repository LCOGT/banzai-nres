def get_telescope_filename(image):
    return image.meta.get('TELESCOP', '').replace('nres', 'nrs')
