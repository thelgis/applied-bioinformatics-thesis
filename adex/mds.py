from typing import List, Tuple, Optional

from adex.helpers import get_pre_processed_dataset, plot_condition_2d, PlottingColorParameters
from adex.models import METADATA_COLUMNS, ConditionDataLoader, DataLoader, ConditionTissueDataLoader, \
    FileDataLoader, ConditionSequencingTissueDataLoader, DATASET_INFO_COLUMNS
from sklearn.manifold import MDS
from matplotlib import pyplot as plt

import logging
import polars as pl
import pandas as pd

from adex.type_aliases import ConditionName, Color


class MdsHelper:

    def __init__(
            self,
            data_loader: DataLoader,
            files_path: str,
            metadata_path: str,
            datasets_info_path: str,
            allowed_null_percentage: float = 0.2
    ):
        self.random_state = 0
        self.data_loader = data_loader
        self.condition = data_loader.condition
        self.dataset: Optional[pl.DataFrame] = (
            get_pre_processed_dataset(
                data_loader=data_loader,
                data_path=files_path,
                metadata_path=metadata_path,
                datasets_info_path=datasets_info_path,
                allowed_null_percentage=allowed_null_percentage
            )
        )
        self.dataset_only_features: Optional[pd.DataFrame] = None

        if self.dataset is not None and self.dataset.shape[0] > 0:
            match data_loader:
                case ConditionDataLoader(condition):
                    logging.info(f"--- Running MDS for '{condition.name}'---")
                case ConditionTissueDataLoader(condition, tissue):
                    logging.info(f"--- Running MDS for '{condition.name}/{tissue.value}'---")
                case FileDataLoader(condition, file_name):
                    logging.info(f"--- Running MDS for '{condition.name}/{file_name}'---")
                case ConditionSequencingTissueDataLoader(condition, sequencing_technique, tissue):
                    logging.info(f"--- Running MDS for '{sequencing_technique.name}|{condition.name}|{tissue.name}'---")
                case _:
                    raise ValueError(f"DataLoader '{data_loader}' not handled in logging")

            samples, genes = self.dataset.shape
            logging.info(f"Loaded dataset for MDS with shape: Samples({samples}), Genes({genes})")

            self.dataset_only_features = (
                self.dataset
                .drop("Sample")
                .drop(METADATA_COLUMNS)
                .drop(DATASET_INFO_COLUMNS)
                .to_pandas()
                .fillna(0)  # Possibly a problem!!! Although not many are left in this dataset due to pre-processing
            )

    def draw_components_stress_plot(self, max_n: int = 8) -> None:
        def run_mds(n_components: int) -> float:
            mds = MDS(
                n_components=n_components,
                metric=True,
                normalized_stress='auto',
                random_state=self.random_state
            )
            mds.fit_transform(self.dataset_only_features)
            return mds.stress_

        if self.dataset_only_features is not None and self.dataset_only_features.shape[0] > 0:
            stress: List[float] = [
                run_mds(n_components=n)
                for n in range(1, max_n)
            ]

            plt.plot(range(1, max_n), stress)
            plt.xticks(range(1, max_n, 1))
            plt.xlabel('n_components')
            plt.ylabel('stress')
            plt.title = f"n_components/stress for {self.condition.name}"
            plt.show()

    def draw_2d(
        self,
        column_that_defines_colors: str,
        target_colors: List[Tuple[ConditionName, Color]]
    ) -> None:

        if self.dataset_only_features is not None and self.dataset_only_features.shape[0] > 0:
            mds = MDS(
                n_components=2,
                metric=True,
                normalized_stress='auto',
                random_state=self.random_state
            )
            mds_fitted = mds.fit_transform(self.dataset_only_features)
            logging.info(f"MDS Stress: {mds.stress_}")

            plot_condition_2d(
                data_loader=self.data_loader,
                method="MDS",
                x_label="Dim1",
                y_label="Dim2",
                df_to_plot=pd.DataFrame(
                    data=mds_fitted,
                    columns=['Dim1', 'Dim2']
                ),
                plotting_color_parameters=PlottingColorParameters(
                    column_that_defines_colors=self.dataset.to_pandas()[column_that_defines_colors],
                    target_colors=target_colors
                )
            )
