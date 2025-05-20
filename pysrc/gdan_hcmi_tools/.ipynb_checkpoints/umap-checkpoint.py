import numpy as np
import umap

from anndata import AnnData
from numpy import ndarray



def fill_missing(data: ndarray) -> ndarray:
    """ Fill missing values with means of non-nan values across rows.

    Deal with missing values by assigning the mean value
    of each bin to missing values of that bin

    Parameters
    ----------
    data : ndarray
        data to fill

    Returns
    -------
    ndarray
        copy of data with missing entries filled
    """

    data = data.copy()

    # Set bins with nan across all cells to 0
    data[:, np.all(np.isnan(data), axis=0)] = 0

    # Mean of each bin, ignoring nan
    bin_means = np.nanmean(data, axis=0)
    bin_means = np.nan_to_num(bin_means, nan=0)
    bin_means = np.tile(bin_means, (data.shape[0], 1))
    data[np.where(np.isnan(data))] = bin_means[np.where(np.isnan(data))]

    return data


def compute_umap(
        adata: AnnData,
        layer_name: str='copy',
        n_components: int=2,
        n_neighbors: int=15,
        min_dist: float=0.1,
        metric: str='euclidean',
    ) -> AnnData:
    """ Cluster cells by copy number.

    Parameters
    ----------
    adata : AnnData
        copy number data
    layer_name : str, optional
        layer with copy number data on which to perform umap dimensionality
        reduction, None for X, by default 'copy'
    n_components : int
        umap n_components param
    n_neighbors : int
        umap n_neighbors param
    min_dist : float
        umap min_dist param

    Returns
    -------
    AnnData
        copy number data with additional `umap_1`, `umap_2` columns
    """

    if layer_name == '_X':
        X = adata.X
    else:
        X = adata.layers[layer_name]

    X = fill_missing(X)

    embedding = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=n_components,
        metric=metric,
        random_state=42,
    ).fit_transform(X)

    adata.obs['UMAP1'] = embedding[:, 0]
    adata.obs['UMAP2'] = embedding[:, 1]

    return adata
