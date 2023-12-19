from functools import reduce
from pathlib import Path
from typing import List, Set

from adex.type_aliases import Gene
from adex.models import Condition
from polars import DataFrame
import polars as pl


def load_data_per_condition(condition: Condition, path: str) -> List[DataFrame]:
    """
    Loads all the datasets of a certain condition in a list of dataframes
    """

    results = [
        pl.read_parquet(file)
        for file in Path(f"{path}/{condition.name}").glob('*.parquet')
    ]

    if len(results) == 0:
        raise ValueError(f"Possibly wrong path '{path}' provided for files")

    return results


def gene_intersection(dataframes: List[DataFrame]) -> Set[Gene]:
    """
    Returns all the common genes found in a list of dataframes
    """
    common_genes = set()

    for df in dataframes:
        if len(common_genes) == 0:
            common_genes.update(
                df.select("gene").to_series().to_list()
            )
        else:
            common_genes.intersection_update(
                set(df.select("gene").to_series().to_list())
            )

    return common_genes


def common_genes_dataframe(dataframes: List[DataFrame]) -> DataFrame:
    """
    Gives a dataframe with the samples of all the dataframes joined but only for the common genes
    """
    head, *tail = dataframes

    return reduce(
        lambda left, right: left.join(right, on="gene", how="inner"),
        tail,
        head
    )


def high_frequency_genes_dataframe(
    dataframes: List[DataFrame],
    allowed_null_percentage: float = 0.2,
    drop_frequencies_column: bool = True
) -> DataFrame:
    """
    Gives a dataframe with the samples of all the dataframes joined, but only for the genes that appear in a specific
    percent of the samples.

    :param dataframes: dataframes to be joined
    :param allowed_null_percentage: will keep only genes that have a lower than this null percentage across samples
    :param drop_frequencies_column: if frequencies column should be dropped (or kept for exploratory analysis)
    :return:
    """
    head, *tail = dataframes

    outer_joined_df = reduce(
        lambda left, right: left.join(right, on="gene", how="outer"),
        tail,
        head
    )

    filtered_df = (
        outer_joined_df.with_columns(pl.sum_horizontal(pl.all().is_null() / pl.all().count()).alias("Null-Percentage"))
        .filter(pl.col("Null-Percentage") <= allowed_null_percentage)
    )

    if drop_frequencies_column:
        return filtered_df.drop("Null-Percentage")

    return filtered_df


def get_pre_processed_dataset(
    condition: Condition,
    data_path: str,
    metadata_path: str,
    allowed_null_percentage: float = 0.2,
    return_metadata: bool = True
) -> DataFrame:
    """
    :param condition:
    :param data_path: the path where the samples are located
    :param metadata_path: the path where the metadata of the samples is located
    :param allowed_null_percentage: will keep only genes that have a lower than this null percentage across samples
    :param return_metadata: if the metadata columns will be returned as part of the dataframe
    :return: a dataset of a particular condition pre-processed in its final state
    """
    data: List[DataFrame] = load_data_per_condition(condition, data_path)

    # keep only frequent genes between datasets
    data_frequent_genes: DataFrame = high_frequency_genes_dataframe(data, allowed_null_percentage)

    transposed = data_frequent_genes.transpose(include_header=True, header_name='Sample')
    transposed_fixed = (
            transposed
            .rename(transposed.head(1).to_dicts().pop())    # add header
            .slice(1,)                                      # remove first row because it is the header duplicated
            .rename({"gene": "Sample"})                     # fix header
    )

    # join with metadata and keep a sample only if metadata exists for the sample
    transposed_fixed_w_metadata = transposed_fixed.join(
        pl.read_csv(metadata_path, separator="\t").unique(subset=["Sample"]),
        on="Sample",
        how="inner"
    )

    # Make data that comes from different sources use the same value for nulls (e.g. metadata uses 'NA')
    final = transposed_fixed_w_metadata.select(
        pl
            .all()
            .replace(
                mapping={
                    "NA": None
                }
            )
    )

    if return_metadata:
        return final
    else:
        return final.drop(
            columns=[
                "GSE",
                "Experimental Strategy",
                "GPL",
                "Condition",
                "Tissue",
                "Cell Type",
                "Gender",
                "Age",
                "Ethnicity"
            ]
        )