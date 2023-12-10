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
