import pandas as pd

def compare_2_df(df1: pd.DataFrame, df2: pd.DataFrame) -> None:
    """Compare two Pandas DataFrames and print the differences.

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

    # Select common columns from each DataFrame
    common_df1 = df1[common_columns]
    common_df2 = df2[common_columns]

    # Compare the selected subsets and store the differences
    differences = common_df1.compare(common_df2)
    print(differences)

    # Check if there are differences
    if differences.empty:
        print("The common parts of the DataFrames are identical.")
    else:
        print("Differences in the common parts of the DataFrames:")
        print(differences)