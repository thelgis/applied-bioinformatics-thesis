from dataclasses import dataclass
from enum import Enum
from typing import List, Tuple

from adex.type_aliases import Color, Tissue


class Condition(Enum):
    RA = 1
    T1D = 2
    SSc = 3
    SLE = 4
    SjS = 5


class SequencingTechnique(Enum):
    MICROARRAYS = "Expression profiling by array"
    RNA_SEQ = "Expression profiling by high throughput sequencing"


METADATA_COLUMNS = [
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


class TissueEnum(Enum):
    PERIPHERAL_BLOOD = "Peripheral blood"
    WHOLE_BLOOD = "Whole blood"
    PAROTIC_GLAND = "Parotid gland"
    SALIVARY_GLAND = "Salivary gland"
    SALIVA = "Saliva"
    SKIN = "Skin"
    SYNOVIAL_MEMBRANE = "Synovial membrane"


TISSUE_COLORS: List[Tuple[Tissue, Color]] = [  # TODO build this using the enum with zip
    ("Peripheral blood", "b"),
    ("Whole blood", "g"),
    ("Parotid gland", "r"),
    ("Salivary gland", "c"),
    ("Saliva", "m"),
    ("Skin", "y"),
    ("Synovial membrane", "k")
]


@dataclass(frozen=True)
class DataLoader:
    condition: Condition


@dataclass(frozen=True)
class ConditionDataLoader(DataLoader):
    """
    Loads all data of one condition
    """
    condition: Condition


@dataclass(frozen=True)
class FileDataLoader(DataLoader):
    """
    Loads data from a specific file of a condition
    """
    condition: Condition
    file_name: str


@dataclass(frozen=True)
class ConditionTissueDataLoader(DataLoader):
    """
    Loads all data of a specific condition and tissue
    """
    condition: Condition
    tissue: TissueEnum
