from enum import Enum


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
