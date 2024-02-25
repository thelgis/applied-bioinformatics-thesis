from dataclasses import dataclass
from functools import reduce
from pathlib import Path
from typing import List, Set, Tuple

from adex.type_aliases import Gene, ConditionName, Color
from adex.models import Condition, METADATA_COLUMNS, DataLoader, ConditionDataLoader, ConditionTissueDataLoader, \
    FileDataLoader
from polars import DataFrame
import polars as pl
import pandas as pd
from pandas.core.series import Series
from matplotlib import pyplot as plt


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
        if len(common_genes) == 0:  # First iteration
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
        lambda left, right: left.join(right, on="gene", how="outer_coalesce"),
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
    data_loader: DataLoader,
    data_path: str,
    metadata_path: str,
    allowed_null_percentage: float = 0.2,
    return_metadata: bool = True
) -> DataFrame:
    """
    :param data_loader: determines the subset of the data that will be loaded
    :param data_path: the path where the samples are located
    :param metadata_path: the path where the metadata of the samples is located
    :param allowed_null_percentage: will keep only genes that have a lower than this null percentage across samples
    :param return_metadata: if the metadata columns will be returned as part of the dataframe
    :return: a dataset of a particular condition pre-processed in its final state
    """

    match data_loader:
        case FileDataLoader(condition, file_name):
            data: List[DataFrame] = [pl.read_parquet(f"{data_path}/{condition.name}/{file_name}")]
        case _:
            data: List[DataFrame] = load_data_per_condition(data_loader.condition, data_path)

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
        pl.read_csv(metadata_path, separator="\t").unique(subset=["Sample"]),  # Filters duplicate rows for a sample in metadata
        on="Sample",
        how="inner"
    )

    # Extra data filtering
    match data_loader:
        case ConditionTissueDataLoader(_, tissue):
            transposed_fixed_w_metadata = transposed_fixed_w_metadata.filter(pl.col("Tissue") == tissue.value)
        case _:
            pass  # nothing to do

    # Make data that comes from different sources use the same value for nulls (e.g. metadata uses 'NA')
    final = transposed_fixed_w_metadata.select(
        pl
        .all()
        .replace(old="NA", new=None)
    )

    if return_metadata:
        return final
    else:
        return final.drop(METADATA_COLUMNS)


@dataclass(frozen=True)
class PlottingColorParameters:
    """
    This class is used to pass the plotting color parameters
    """
    column_that_defines_colors: Series
    target_colors: List[Tuple[ConditionName, Color]]


def plot_condition_2d(
        data_loader: DataLoader,
        method: str,
        x_label: str,
        y_label: str,
        df_to_plot: pd.DataFrame,
        plotting_color_parameters: PlottingColorParameters
) -> None:
    plt.figure()
    plt.figure(figsize=(10, 10))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=14)

    plt.xlabel(x_label, fontsize=20)
    plt.ylabel(y_label, fontsize=20)
    match data_loader:
        case ConditionDataLoader(condition):
            plt.title(f"{method} of {condition.name} Dataset", fontsize=20)
        case ConditionTissueDataLoader(condition, tissue):
            plt.title(f"{method} of {condition.name} Dataset (Tissue: '{tissue.value}')", fontsize=20)
        case FileDataLoader(condition, file_name):
            plt.title(f"{method} of {condition.name} Dataset (File: '{file_name}')", fontsize=20)
        case _:
            raise ValueError(f"DataLoader '{data_loader}' not handled in plotting")

    for target, color in plotting_color_parameters.target_colors:
        indices = plotting_color_parameters.column_that_defines_colors == target
        plt.scatter(
            df_to_plot.loc[indices, x_label],
            df_to_plot.loc[indices, y_label],
            c=color,
            s=50
        )

    targets = [target for target, _ in plotting_color_parameters.target_colors]
    plt.legend(targets, prop={'size': 15})
