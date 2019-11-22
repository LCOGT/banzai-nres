def fiber_states_from_header(header):
    parsed_objects = header['OBJECTS'].split('&')
    if len(parsed_objects) < 3:
        return 0, 0, 0
    fiber0_lit = parsed_objects[0].lower() != 'none'
    fiber1_lit = parsed_objects[1].lower() != 'none'
    fiber2_lit = parsed_objects[2].lower() != 'none'
    return int(fiber0_lit), int(fiber1_lit), int(fiber2_lit)


def fibers_state_to_filename(image):
    return str(int(image.fiber0_lit)) + str(int(image.fiber1_lit)) + str(int(image.fiber2_lit))
