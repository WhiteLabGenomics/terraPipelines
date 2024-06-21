import pandas as pd


def compare_2_df(df1: pd.DataFrame, df2: pd.DataFrame) -> None:
    """compare the column names in two dataframes and identify any columns with
    identical names, referred to as 'common columns' here. We do not necessarily
    expect to find common columns. If common columns are present, we create
    subsets of the dataframes using these columns and then compare the rows
    within these two subset dataframes.

    The differences raws are printed to the console.

    Parameters
    ----------
    df1 : pd.DataFrame
        The first DataFrame to compare.
    df2 : pd.DataFrame
        The second DataFrame to compare.

    Returns
    -------
    None
    """

    # Identify common columns
    common_columns = df1.columns.intersection(df2.columns)

    if len(common_columns) == 0:
        print("No common features found to compare.")
    else:
        # Select common columns from each DataFrame
        common_df1 = df1[common_columns]
        common_df2 = df2[common_columns]

        # Compare the selected subsets and store the differences
        differences = common_df1.compare(common_df2)

        # Check if there are differences
        if differences.empty:
            print("The common parts of the DataFrames are identical.")
        else:
            print("Differences in the common parts of the DataFrames:")
            print(differences)
