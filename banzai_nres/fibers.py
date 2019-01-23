def fiber_states_from_header(header):
    parsed_objects = header['OBJECTS'].split('&')
    fiber0_lit = parsed_objects[0].lower() != 'none'
    fiber1_lit = parsed_objects[1].lower() != 'none'
    fiber2_lit = parsed_objects[2].lower() != 'none'
    return fiber0_lit, fiber1_lit, fiber2_lit


def fibers_state_to_filename(image):
    return str(int(image.fiber0_lit)) + str(int(image.fiber1_lit)) + str(int(image.fiber2_lit))
