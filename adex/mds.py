from typing import List, Tuple

from adex.helpers import get_pre_processed_dataset, plot_condition_2d, PlottingColorParameters
from adex.models import Condition, METADATA_COLUMNS
from sklearn.manifold import MDS
from matplotlib import pyplot as plt

import logging
import polars as pl
import pandas as pd

from adex.type_aliases import ConditionName, Color


class MdsHelper:

    def __init__(
            self,
            condition: Condition,
            files_path: str,
            metadata_path: str
    ):
        self.condition = condition
        self.dataset: pl.DataFrame = get_pre_processed_dataset(condition, files_path, metadata_path)
        logging.info(f"--- Initializing data to run MDS for condition '{condition.name}'---")

        samples, genes = self.dataset.shape
        logging.info(f"Loaded dataset with shape: Samples({samples}), Genes({genes})")

        self.dataset_only_features = self.dataset.drop("Sample").drop(METADATA_COLUMNS).to_pandas().fillna(0)
        self.random_state = 0

    def draw_components_stress_plot(self, max_n: int = 9) -> None:
        def run_mds(n_components: int) -> float:
            mds = MDS(
                n_components=n_components,
                metric=True,
                normalized_stress='auto',
                random_state=self.random_state
            )
            mds.fit_transform(self.dataset_only_features)
            return mds.stress_

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

        mds = MDS(
            n_components=2,
            metric=True,
            normalized_stress='auto',
            random_state=self.random_state
        )
        mds_fitted = mds.fit_transform(self.dataset_only_features)
        logging.info(f"MDS Stress: {mds.stress_}")

        plot_condition_2d(
            condition=self.condition,
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
