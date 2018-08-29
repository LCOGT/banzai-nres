def add_nres_attribute(images, attribute_name, class_to_append):
    """
    :param images: banzai image objects
    :param attribute_name: string
    :param class_to_append: A class.
    :return: the same image objects with an instance of the class_to_append appended if it does
    not yet exist or if it is set to None.
    """
    for image in images:
        if not hasattr(image, attribute_name):
            setattr(image, attribute_name, class_to_append())
        else:
            if getattr(image, attribute_name) is None:
                setattr(image, attribute_name, class_to_append())