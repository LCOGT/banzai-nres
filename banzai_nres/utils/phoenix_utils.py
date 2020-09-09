from banzai_nres.utils.fits_utils import fits_bytes_to_header


def plus_minus_to_pm(value, fmt):
    output_str = f'{value:{fmt}}'
    return output_str.replace('+', 'p').replace('-', 'm')


def pm_to_plus_minus(value):
    return value.replace('p', '+').replace('m', '-')


def filename_to_parameters(filename):
    _, T_effective, log_g, metallicity, alpha = filename.replace('.fits', '').split('-')
    T_effective = int(pm_to_plus_minus(T_effective))
    log_g = float(pm_to_plus_minus(log_g))
    metallicity = float(pm_to_plus_minus(metallicity))
    alpha = float(pm_to_plus_minus(alpha))
    return T_effective, log_g, metallicity, alpha


def parameters_to_filename(T_effective, log_g, metallicity, alpha):
    filename = f'phoenix-{plus_minus_to_pm(int(T_effective), "05d")}-{plus_minus_to_pm(log_g, "+0.1f")}'
    filename += f'-{plus_minus_to_pm(metallicity, "+0.1f")}-{plus_minus_to_pm(alpha, "+0.1f")}.fits'
    return filename


def parse_phoenix_header(header_lines):
    # Note all values are in cgs units
    header = fits_bytes_to_header(header_lines)
    T_effective = header['PHXTEFF']
    log_g = header['PHXLOGG']
    metallicity = header['PHXM_H']
    alpha = header['PHXALPHA']
    radius = header['PHXREFF']
    luminosity = header['PHXLUM']
    mass = header['PHXMASS']
    return T_effective, log_g, metallicity, alpha, radius, luminosity, mass