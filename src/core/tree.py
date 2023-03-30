import os
from dataclasses import dataclass

import string
import pickle
import re
import pandas as pd
import tqdm

import requests
import pandas as pd

from ..helper.lib import info, error, panic

__all__ = ['Tree', 'Value']
datetime_format = '%Y-%m-%d %H:%M:%S'
non_valid_chars = string.punctuation.join(string.whitespace)
timeouts = (5, 10)


@dataclass
class Value:
  name: str
  path: str
  nc: list[str]


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

  def build(self, force_rebuild: bool = False) -> None:
    """
    Build the tree from the given name.\\
    This method should only be called once.

    ## Parameters
    ```py
    >>> force_rebuild : bool, (optional)
    ```
    force the tree to be rebuilt from https request\\
    defaults to `False`
    """
    pickle_path = os.path.join('data', 'tree.pkl')
    if not force_rebuild and os.path.exists(pickle_path):
      info(f'loading tree ({pickle_path}) from pickle')
      with open(pickle_path, 'rb') as f:
        self.__data = pickle.load(f)
      for organism in tqdm.tqdm(self.__data, desc='building tree'):
        path = self.__data[organism].path
        os.makedirs(path, exist_ok=True)
      return

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

    total_rows = len(df.index)
    for _, row in tqdm.tqdm(df.iterrows(), total=total_rows, desc='building tree'):
      # transform non valid chars to underscores
      organism = re.sub(f'[{non_valid_chars}]', '_', row['#Organism/Name'])
      subgroup = re.sub(f'[{non_valid_chars}]', '_', row['SubGroup'])
      group = re.sub(f'[{non_valid_chars}]', '_', row['Group'])
      kingdom = re.sub(f'[{non_valid_chars}]', '_', row['Kingdom'])
      path = os.path.join('Results', kingdom, group, subgroup, organism)

      if organism not in self.__data:
        self.__data[organism] = Value(organism, path, [])

    valid_organisms: set[str] = set()
    self.update_ids()
    ids_files = os.listdir('data/ids')
    for ids in ids_files:
      info(f'updating ids data from ({ids})')
      with open(f'data/ids/{ids}', 'r') as f:
        for line in f.readlines():
          row = line.split('\t')
          if not row[1].startswith('NC'):
            continue
          organism = re.sub(f'[{non_valid_chars}]', '_', row[5])
          if organism in self.__data:
            self.__data[organism].nc.append(row[1])
            valid_organisms.add(organism)

    info(f'found {len(valid_organisms)} valid organisms')
    self.clean_folders()
    for organism in tqdm.tqdm(valid_organisms, desc='creating folders'):
      path = self.__data[organism].path
      os.makedirs(path, exist_ok=True)

    # filter out invalid organisms
    self.__data = {k: v for k, v in self.__data.items() if k in valid_organisms}
    info(f'filtered out {total_rows - len(self.__data)} invalid organisms')

    # save the tree
    with open('data/tree.pkl', 'wb') as f:
      pickle.dump(self.__data, f)
  
  def get_info(self, organism: str) -> Value:
    """
    returns the value of the given organism.

    ## Parameters
    ```py
    >>> organism : str
    ```
    the name of the organism to get the value of
    """
    return self.__data[organism]

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
      info(f'updated {file}')
