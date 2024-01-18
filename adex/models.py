from enum import Enum
from typing import List, Tuple

from adex.type_aliases import Color, Tissue


class Condition(Enum):
    RA = 1
    T1D = 2
    SSc = 3
    SLE = 4
    SjS = 5


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

TISSUES: List[Tuple[Tissue, Color]] = [
    ("Peripheral blood", "b"),
    ("Whole blood", "g"),
    ("Parotid gland", "r"),
    ("Salivary gland", "c"),
    ("Saliva", "m"),
    ("Skin", "y"),
    ("Synovial membrane", "k")
]

