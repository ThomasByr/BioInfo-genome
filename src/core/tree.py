import os
from dataclasses import dataclass

import string
import pickle
import re
import pandas as pd
import tqdm

import requests

from dotenv import load_dotenv

from ..helper import info, error, panic, capture

__all__ = ['Tree', 'Value']
datetime_format = '%Y-%m-%d %H:%M:%S'
non_valid_chars = string.punctuation + string.whitespace
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

    load_dotenv()

    __rebuild = os.getenv('REBUILD')
    self.__is_rebuild: bool = False
    if __rebuild is not None:
      if isinstance(__rebuild, str):
        if __rebuild.lower() in {'true', '1'}:
          self.__is_rebuild = True
      elif isinstance(__rebuild, int):
        if __rebuild == 1:
          self.__is_rebuild = True
      elif isinstance(__rebuild, bool):
        self.__is_rebuild = __rebuild

  def build(self, force_rebuild: bool = False, silent: bool = False) -> None:
    """
    Build the tree from the given name.\\
    This method should only be called once.

    ## Parameters
    ```py
    >>> force_rebuild : bool, (optional)
    ```
    force the tree to be rebuilt from https request\\
    defaults to `False`
    ```py
    >>> silent : bool, (optional)
    ```
    suppress all excessively spamming output\\
    defaults to `False`
    """
    pickle_path = os.path.join('data', 'tree.pkl')
    if not (force_rebuild or self.__is_rebuild) and os.path.exists(pickle_path):
      info(f'loading tree ({pickle_path}) from pickle')
      if silent:
        capture.redirect()
      with open(pickle_path, 'rb') as f:
        self.__data.clear()
        self.__data = pickle.load(f)
      for organism in tqdm.tqdm(self.__data, desc='building tree'):
        path = self.__data[organism].path
        if not os.path.exists(path):       # check here to save kernel time
          os.makedirs(path, exist_ok=True) # (useless) extra layer of safety
      if silent:
        capture.stop_redirect()
        info(f'loaded tree ({pickle_path}) from pickle')
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
    self.__data.clear()
    info(f'building tree from online {self.__name}.txt ({total_rows} rows)')
    if silent:
      capture.redirect()
    for _, row in tqdm.tqdm(df.iterrows(), total=total_rows, desc='building tree'):
      # transform non valid chars to underscores
      try:
        organism = re.sub(f'[{non_valid_chars}]', '_', row['#Organism/Name'])
        subgroup = re.sub(f'[{non_valid_chars}]', '_', row['SubGroup'])
        group = re.sub(f'[{non_valid_chars}]', '_', row['Group'])
        kingdom = re.sub(f'[{non_valid_chars}]', '_', row['Kingdom'])
        path = os.path.join('Results', kingdom, group, subgroup, organism)
      except TypeError:
        nan_entries = [
          k for k, v in row.items() if pd.isna(v) and k in {
            '#Organism/Name',
            'SubGroup',
            'Group',
            'Kingdom',
          }
        ]
        error(f'failed to parse row:\n{row}\nnan (needed) entries: {nan_entries}')

      if organism not in self.__data:
        # create new entry
        self.__data[organism] = Value(organism, path, [])

    valid_organisms: set[str] = set()
    if silent:
      capture.stop_redirect()
    self.update_ids()
    ids_files = os.listdir(os.path.join('data', 'ids'))
    for ids in ids_files:
      info(f'updating ids data from {ids}')
      with open(os.path.join('data', 'ids', ids), 'r', encoding='utf-8') as f:
        for line in f.readlines():
          row = line.split('\t')
          if not row[1].startswith('NC'):
            continue
          if (organism := re.sub(f'[{non_valid_chars}]', '_', row[5])) in self.__data:
            self.__data[organism].nc.append(row[1])
            valid_organisms.add(organism)

    info(f'found {len(valid_organisms)} valid organisms')
    self.clean_folders()
    if silent:
      capture.redirect()
    for organism in tqdm.tqdm(valid_organisms, desc='creating folders'):
      path = self.__data[organism].path
      os.makedirs(path, exist_ok=True)

    # filter out invalid organisms
    self.__data = {k: v for k, v in self.__data.items() if k in valid_organisms}
    info(f'filtered out {total_rows - len(self.__data)} invalid organisms')
    if silent:
      capture.stop_redirect()
      info('system file tree reset successfully')

    # save the tree
    with open(pickle_path, 'wb') as f:
      pickle.dump(self.__data, f)

  def get_info(self, organism: str) -> Value:
    """
    returns the value of the given organism.

    ## Parameters
    ```py
    >>> organism : str
    ```
    the name of the organism to get the value of

    ## Returns
    ```py
    Value : Value
    ```
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

  def union(self, o: 'Tree') -> 'Tree':
    """
    returns a new tree that is the union of this tree and the given tree.

    ## Parameters
    ```py
    >>> o : Tree
    ```
    the tree to union with

    ## Returns
    ```py
    Tree : Tree
    ```
    """
    tree = Tree(self.__name)
    tree.__data = {**self.__data, **o.__data}
    return tree
