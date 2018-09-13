from banzai_nres import bias


def test_min_images():
    bias_stage = bias.DarkMaker(None)
    assert bias_stage.min_images == 5
