def unique_elements_unordered(seq):
    """
    :param seq: list or sequence
    :return: returns unique elements in the list, except the order in which they are found is preserved, e.g.
    inputting [1, 6, 7, 7, 2, 2] will return [1, 6, 7, 2]
    """
    checked = []
    for e in seq:
        if e not in checked:
            checked.append(e)
    return checked
