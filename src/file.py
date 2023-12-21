from typing import Union
import fsspec
import anndata
import pandas as pd

def read_file_from_url(url: str) -> Union[anndata.AnnData, pd.DataFrame]:
    """Read a file from a URL and return the appropriate data structure.

    Parameters
    ----------
    url : str
        The URL of the file to read.

    Returns
    -------
    output : Union[anndata.AnnData, pd.DataFrame]
        The data read from the file, which can be either an AnnData object (if ".h5ad" in url)
        or a Pandas DataFrame (if ".csv" in url).
    """

    with fsspec.open(url) as f:
        if ".h5ad" in url:
            output = anndata.read_h5ad(f)
        elif "Summary.csv" in url:
            output = pd.read_csv(f, header=None) #Read the CSV file with the first row as data
        elif ".csv" in url:
            output = pd.read_csv(f)
        else:
            raise ValueError("Unsupported file format. Supported formats are .h5ad and .csv.")
    
    return output