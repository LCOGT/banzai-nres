from banzai.stages import Stage


class ProfileFitter(Stage):
    def do_stage(self, image):
        # For each order
        # Extract the pixels from the spectrum extension in that order
        # Fit a smooth B-spline to the pixels

        # Evaluate the best fit spline along the center of the trace, normalize, save as BLAZE
        # Evaluate the best fit spline/BLAZE save in an empty array, normalize (PROFILE extension)
        # Divide the original image by best fit spline (Pixel-to-pixel sensitivity)
        return image
