import os
from dataclasses import dataclass
from datetime import datetime

import string
import pickle
import json

from typing import Any
import requests
import pandas as pd

from ..helper.lib import *

__all__ = ['Tree']
datetime_format = '%Y-%m-%d %H:%M:%S'
non_valid_chars = string.punctuation.join(string.whitespace)
timeouts = (5, 10)


def to_lower(s: str | Any) -> str:
  """
  Convert a string to lowercase.

  ## Parameters
  ```py
  s : str
  ```
  the string to convert
  """
  return str(s).lower()


@dataclass
class Value:
  name: str
  nc: str


class Tree:
  BASE_URL = 'https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/'

  def __init__(self, name: str = None) -> None:
    """
    Build a tree from a given name.\\
    If no name is given, the tree will be built from the overview file.

    ## Parameters
    ```py
    name : str, (optional)
    ```
    the name of the tree to build\\
    defaults to `None`
    """
    self.__name = 'overview' if name is None else name
    self.__data: dict[str, Value] = {}

  def build(self) -> None:
    """
    Build the tree from the given name.\\
    This method should only be called once IF the tree does NOT comes from `union`.
    """
    r = requests.get(f'{self.BASE_URL}{self.__name}.txt', stream=True, timeout=timeouts)
    if r.status_code != 200:
      panic(f'failed to fetch {self.__name}.txt')

    r.raw.decode_content = True

    df: pd.DataFrame
    try:
      df = pd.read_csv(r.raw, sep='\t', low_memory=False)
    except pd.errors.EmptyDataError:
      error(f'failed to parse {self.__name}.txt (file is empty)')
      return

    for _, row in df.iterrows():
      pass

  @staticmethod
  def clean_folders() -> None:
    """
    deletes all data from the `Results` directory\\
    basically leaves the directory empty.
    """
    for root, dirs, files in os.walk('Results', topdown=False):
      for name in files:
        os.remove(os.path.join(root, name))
      for name in dirs:
        os.rmdir(os.path.join(root, name))

  @classmethod
  def update_ids(cls) -> None:
    """
    updates `data/ids` folder with the latest ids\\
    this function will overwrite the existing files.
    """
    url = cls.BASE_URL + 'IDS/'
    # for each file in the ids folder
    for file in os.listdir('data/ids'):
      name = file.split('.')[0]
      # get the file from the server
      r = requests.get(url + file, stream=True)
      # if the file is not found
      if r.status_code == 404:
        error(f'failed to update {name} (file not found)')
        continue
      # if the file is found
      with open(f'data/ids/{file}', 'wb') as f:
        f.write(r.content)
