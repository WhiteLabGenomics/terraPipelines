import anndata
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt


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
        in_adata, keys="pct_counts_mt", jitter=False, save=f"_{data_source}_pct_mt.pdf"
    )  # jitter=False to remove the dots

    p4 = sc.pl.scatter(
        in_adata,
        x="total_counts",
        y="n_genes_by_counts",
        color="pct_counts_mt",
        title="#reads x #genes colored by percentage of mitochondiral genes",
        save=f"_{data_source}_counts_n_genes_by_counts.pdf",
    )

    ## %mt by batch in violin plots
    if batch is not None:
        violin_params = {
            "adata": in_adata,
            "keys": "pct_counts_mt",
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

    return


def pca_variance_ratio_plots(in_adata: anndata.AnnData, nb_pcs: int) -> None:
    """Plot and compute the variance ratio of the PCA.
    Parameters
    ----------
    in_adata : AnnData
        The AnnData object to check.

    nb_pcs : int
        The number of PCs to plot.

    Returns
    -------
    output : None
    """

    # Access the percentage of explained variance for each component
    pca_variance_ratio = in_adata.uns["pca"]["variance_ratio"]
    # Print the percentage of explained variance for each component
    for i, variance in enumerate(pca_variance_ratio):
        print(f"PC{i+1}: {variance * 100}%")

    # Calculate the cumulative variance
    cumulative_variance = pca_variance_ratio.cumsum()
    # Print the cumulative variance
    print("Cumulative Variance:")
    for i, variance in enumerate(cumulative_variance):
        print(f"PC{i+1}: {variance * 100}%")

    # Plot the cumulative variance
    plt.plot(range(1, len(cumulative_variance) + 1), cumulative_variance, marker="o")
    plt.xlabel("Number of PCs")
    plt.ylabel("Cumulative Variance Explained")
    plt.title("Elbow Graph")

    # Add PC number to the graph
    for i, variance in enumerate(cumulative_variance):
        plt.text(i + 1, variance, f"PC{i+1}", ha="center", va="bottom", fontsize=8)

    # plot
    # sc.pl.pca_variance_ratio(adata, log=False, n_pcs=n_pc, save="elbow_plot.png")
    sc.pl.pca_variance_ratio(in_adata, log=False, n_pcs=nb_pcs)

    return
