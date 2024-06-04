import numpy as np
from scipy.stats import median_abs_deviation
import anndata


def human_mt_genes_ident(in_adata: anndata.AnnData) -> None:
    """Identify human mitochondrial genes, add the boolean column "mt" to adata.var
    This function also check if the in_adata contains mouse mt genes

    Parameters
    ----------
    in_adata : AnnData
        The raw AnnData object with mitochondrial genes identification

    Returns
    -------
    output : None
    """

    count_true = (in_adata.var_names.str.startswith("mt-")).sum()

    if count_true > 0:
        print("Mouse mt number", count_true)
    else:
        print("No mouse mt genes")
        print("Adding human mt genes to adat.var")
        in_adata.var["mt"] = in_adata.var_names.str.startswith("MT-")


def is_outlier(adata, metric: str, nmads: int) -> pd.Series[bool]:
    """Function that takes a metric, i.e. a column in .obs and the number of MADs (nmad).
    This function is from https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html#filtering-low-quality-reads/,
    which is a part of the best practices for single cell RNA-seq analysis. Please cite the authors if you use this function.

    Parameters
    ----------
    adata : AnnData
        The AnnData object that we want to detect outliers in.

    metric : str
        This is the column that contains the Quality Control (QC) metric we want
        to evaluated the outliers and subsequently filter.

    nmads : int
        The value of MADs (median absolute deviations) to use to define outliers.

    Returns
    -------
    output : pd.Series[bool]
        A boolean series that is True if the value is an outlier and False otherwise.
    """

    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier
