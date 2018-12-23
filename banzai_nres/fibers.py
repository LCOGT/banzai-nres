class FibersState(object):
    def __init__(self, fiber0_lit: bool, fiber1_lit: bool, fiber2_lit: bool):
        self.fiber0_lit = fiber0_lit
        self.fiber1_lit = fiber1_lit
        self.fiber2_lit = fiber2_lit

    def fiber_value(self):
        fiber_value = self.fiber0_lit * 1 + self.fiber1_lit * 2 + self.fiber2_lit * 4
        return fiber_value

    def __eq__(self, fiber_state):
        fiber_states_match = fiber_state.fiber0_lit == self.fiber0_lit
        fiber_states_match &= fiber_state.fiber1_lit == self.fiber1_lit
        fiber_states_match &= fiber_state.fiber2_lit == self.fiber2_lit
        return fiber_states_match

    def __gt__(self, fiber_state):
        return self.fiber_value() > fiber_state.fiber_value()

    def __ge__(self, fiber_state):
        return (self.__gt__(fiber_state)) or self.__eq__(fiber_state)

    def __lt__(self, fiber_state):
        return self.fiber_value() < fiber_state.fiber_value()

    def __le__(self, fiber_state):
        return (self.__lt__(fiber_state)) or self.__eq__(fiber_state)

    def __str__(self):
        return str(int(self.fiber0_lit)) + str(int(self.fiber1_lit)) + str(int(self.fiber2_lit))

    @classmethod
    def from_header(cls, header):
        parsed_objects = header['OBJECTS'].split('&')
        fiber0_lit = parsed_objects[0].lower() != 'none'
        fiber1_lit = parsed_objects[1].lower() != 'none'
        fiber2_lit = parsed_objects[2].lower() != 'none'
        return cls(fiber0_lit, fiber1_lit, fiber2_lit)
