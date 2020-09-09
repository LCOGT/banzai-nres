from astropy.io import fits


def fits_bytes_to_header(header_lines):
    header = fits.Header({})
    for i in range(0, len(header_lines), 80):
        header_line = header_lines[i:i+80]
        keyword = header_line[:8].strip()
        comment_start = header_line.find("\\")
        value = header_line[9:comment_start].strip()
        try:
            value = int(value)
        except:
            try:
                value = float(value)
            except:
                pass
        if value == 'T':
            value = True
        elif value == 'F':
            value = False
        comment = header_line[comment_start + 1:].strip()

        header[keyword] = value, comment
    return header