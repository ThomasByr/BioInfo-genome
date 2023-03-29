import os

import string
import random

from ..helper import info

from Bio import Entrez, SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

__all__ = []


def create_data_from_NC(name : str, path: str, NC_list: list[str], region: str) -> None:
  pass
