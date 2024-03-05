import anndata
import scanpy as sc
import seaborn as sns


def qc_plots(in_adata: anndata.AnnData, data_source: str, batch: str) -> None:
    """Plot QC plots for the given AnnData object.

    Parameters
    ----------
    in_adata : AnnData
        The AnnData object to check.

    data_source : str
        The data source of the AnnData object.

    batch : str
        The name of the column that contains batch information.
        If batch is None, p5 will not be plotted.

    Returns
    -------
    output : None
    """

    p1 = sns.displot(in_adata.obs["total_counts"], bins=100, kde=False)
    p1.set(xlabel="Distribution of counts per barcodes", ylabel="Counts per barcode")
    p1.fig.savefig(f"displot_{data_source}_total_counts.pdf")

    p2 = sc.pl.violin(
        in_adata,
        keys="total_counts",
        jitter=False,
        save=f"_{data_source}_total_counts.pdf",
    )

    p3 = sc.pl.violin(
        in_adata,
        keys="pct_counts_mt",
        # show=False,
        jitter=False,
        save=f"_{data_source}_pct_mt.pdf",
    )  # jitter=False to remove the dots

    p4 = sc.pl.scatter(
        in_adata,
        x="total_counts",
        y="n_genes_by_counts",
        color="pct_counts_mt",
        # show=False,
        title="#reads x #genes colored by percentage of mitochondiral genes",
        save=f"_{data_source}_counts_n_genes_by_counts.pdf",
    )
    """
    p4.set(
        title="#reads x #genes colored by percentage of mitochondiral genes",
        xlabel="Number of reads per barcode",
        ylabel="Number of expressed genes per barcode",
    )
    """

    ## %mt by batch in violin plots
    if batch is not None:
        violin_params = {
            "adata": in_adata,
            "keys": "pct_counts_mt",
            #'show': False,
            "jitter": False,
            "rotation": 90,
            "title": "Proportion of mitochondrial reads per barcode in function of {batch}",
            "save": f"_{data_source}_pct_mt_by_{batch}.pdf",
        }

        if batch == "emptydrops_IsCell":
            # Create a mapping dictionary
            mapping = {True: "True", False: "False"}
            # Apply the mapping to the emptydrops_IsCell column
            in_adata.obs["emptydrops_IsCell_category"] = in_adata.obs[
                "emptydrops_IsCell"
            ].map(mapping)

            # Update parameters for this case
            violin_params.update({"groupby": "emptydrops_IsCell_category"})
        else:
            # Update parameters for this case
            violin_params.update({"groupby": batch})

        # Plot the violin plot
        p5 = sc.pl.violin(**violin_params)

        """
        p5.set(title="Proportion of mitochondrial reads per barcode per sample in function of emptyDrops cell calling result")
        """

    return
