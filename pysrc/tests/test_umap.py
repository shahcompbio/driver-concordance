import numpy as np

from gdan_hcmi_tools.umap import fill_missing

def test_fill_missing():
    data = np.array([
        [1, 1, np.nan],
        [1, 2, 1],
        [1, 2, 2],
    ])
    output = fill_missing(data)
    expected = np.array([
        [1, 1, 1.5],
        [1, 2, 1],
        [1, 2, 2],
    ])
    assert np.all(output == expected), output