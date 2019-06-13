from banzai_nres.blaze import as_astropy_table
from astropy.table import Table
from banzai.images import DataTable


def test_as_astropy_table():
    output = as_astropy_table(DataTable(data_table={'o': [1]}, name='b'))
    assert isinstance(output, Table)
    assert isinstance(as_astropy_table({'o': [1]}), Table)
