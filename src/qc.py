import numpy as np
from scipy.stats import median_abs_deviation
import anndata

def human_mt_genes_ident(in_adata: anndata.AnnData) -> None:
    """Identify human mitochondrial genes, add the boolean column "mt" to adata.var 
    This function also check if the adata contains mouse mt genes
    NB mitochondrial genes MT- for human data and mt- for mouse data
    
    Parameters
    ----------
    in_adata : AnnData
        The AnnData object to check.
    
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


def is_outlier(adata, metric: str, nmads: int) -> bool:
    """Function that takes a metric, i.e. a column in .obs and the number of MADs (nmad).
    This function is from https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html#filtering-low-quality-reads/,
    which is a part of the best practices for single cell RNA-seq analysis. Please cite the authors if you use this function.
    
    Parameters
    ----------
    adata : AnnData
        The AnnData object to check.
    
    metric : str
        The name of the column that contains the metric to check.
    
    nmads : int
        The number of MADs to use to define outliers.
    
    Returns
    -------
    output : bool
        True if the metric is an outlier, False otherwise.
    """
    
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier