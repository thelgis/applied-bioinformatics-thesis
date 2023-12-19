from adex.helpers import get_pre_processed_dataset
from adex.models import Condition, METADATA_COLUMNS

import logging
import polars as pl
import pandas as pd
import numpy as np

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt


class PcaHelper:

    def __init__(
            self,
            condition: Condition,
            files_path: str,
            metadata_path: str
    ):
        self.condition = condition
        self.dataset: pl.DataFrame = get_pre_processed_dataset(condition, files_path, metadata_path)
        logging.info(f"--- Running PCA for condition '{condition.name}'---")

        samples, genes = self.dataset.shape
        logging.info(f"Loaded dataset for PCA with shape: Samples({samples}), Genes({genes})")

        dataset_only_features: pd.DataFrame = self.dataset.drop("Sample").drop(METADATA_COLUMNS).to_pandas()
        dataset_features_array: np.ndarray = dataset_only_features.values

        dataset_features_normalized = StandardScaler().fit_transform(dataset_features_array)
        logging.info(
            f"Dataset normalised | "
            f"Mean: '{np.nanmean(dataset_features_normalized)}' "
            f"Standard Deviation:'{np.nanstd(dataset_features_normalized)}'"
        )

        self.pca = PCA(n_components=2)
        self.pca_fitted = self.pca.fit_transform(pd.DataFrame(dataset_features_normalized).fillna(0))

        logging.info(f"Explained variation per principal component: {self.pca.explained_variance_ratio_}")
        logging.info(f"----------------------------------------------")

    def as_pandas_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(
            data=self.pca_fitted,
            columns=['PC1', 'PC2']
        )

    def as_polars_dataframe(self) -> pl.DataFrame:
        return pl.DataFrame(
            data=self.pca_fitted,
            schema=['PC1', 'PC2']
        )

    def draw(self) -> None:
        plt.figure()
        plt.figure(figsize=(10, 10))
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=14)

        plt.xlabel('PC1', fontsize=20)
        plt.ylabel('PC2', fontsize=20)
        plt.title(f"Principal Component Analysis of {self.condition.name} Dataset", fontsize=20)

        targets = ['Healthy', self.condition.name]
        colors = ['g', 'r']

        df = self.as_pandas_dataframe()

        for target, color in zip(targets, colors):
            index = self.dataset.to_pandas()['Condition'] == target
            plt.scatter(
                df.loc[index, 'PC1'],
                df.loc[index, 'PC2'],
                c=color,
                s=50
            )

        plt.legend(targets, prop={'size': 15})
