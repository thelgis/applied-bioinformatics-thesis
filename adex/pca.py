from adex.helpers import get_pre_processed_dataset, plot_condition_2d
from adex.models import Condition, METADATA_COLUMNS

import logging
import polars as pl
import pandas as pd
import numpy as np

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


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

    def pca_as_pandas_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(
            data=self.pca_fitted,
            columns=['PC1', 'PC2']
        )

    def pca_as_polars_dataframe(self) -> pl.DataFrame:
        return pl.DataFrame(
            data=self.pca_fitted,
            schema=['PC1', 'PC2']
        )

    def draw(self) -> None:
        plot_condition_2d(
            condition=self.condition,
            method="Principal Component Analysis",
            x_label="PC1",
            y_label="PC2",
            df_to_plot=self.pca_as_pandas_dataframe(),
            condition_column=self.dataset.to_pandas()['Condition'],
        )
