from banzai.stages import Stage


class SpectrumExtractor(Stage):
    def do_stage(self, image):
        # Normalization: Histogram of Mask * Profile**2 / Variance
        # Do a histogram of Mask * Profile * Data / Variance / Normalization
        # Variance = Mask * Profile / Normalization
        return image
