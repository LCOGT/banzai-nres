from banzai_nres import dark


def test_min_images():
    dark_stage = dark.DarkMaker(None)
    assert dark_stage.min_images == 3
