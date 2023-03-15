import os
from dataclasses import dataclass
from datetime import datetime

import json
from typing import Any
import requests
import pandas as pd

from ..helper.lib import *

__all__ = ['Tree', 'SubGroup', 'Group', 'Kingdom']
datetime_format = '%Y-%m-%d %H:%M:%S'
timeouts = (5, 10)


@dataclass
class SubGroup:
  """
  A subgroup struct that holds :
  - `name` : a string representing the name of the entity
  - `id` : a string representing the genes ids
  """
  name: str
  id: str

  # entries: set[str]
  def __hash__(self) -> int:
    return hash(self.name)

  def to_json(self) -> dict[str, str]:
    return {'name': self.name, 'id': self.id}

  @classmethod
  def from_json(cls, data: dict[str, str]) -> 'SubGroup':
    return cls(data['name'], data['id'])


@dataclass
class Group:
  """
  A group struct that holds :
  - `name` : a string representing the name of the entity
  - `subgroups` : a dictionnary of subgroups names and objects
  """
  name: str
  subgroups: dict[str, SubGroup]

  def __hash__(self) -> int:
    return hash(self.name)

  def to_json(self) -> dict[str, str]:
    return {'name': self.name, 'subgroups': [subgroup.to_json() for subgroup in self.subgroups.values()]}

  @classmethod
  def from_json(cls, data: dict[str, str]) -> 'Group':
    return cls(data['name'], {subgroup['name']: SubGroup(subgroup['name'], subgroup['id']) for subgroup in data['subgroups']})


@dataclass
class Kingdom:
  """
  A group struct that holds :
  - `name` : a string representing the name of the entity
  - `groupe` : a dictionnary of groups names and objects
  """
  name: str
  groups: dict[str, Group]

  def __hash__(self) -> int:
    return hash(self.name)

  def to_json(self) -> dict[str, str]:
    return {'name': self.name, 'groups': [group.to_json() for group in self.groups.values()]}

  @classmethod
  def from_json(cls, data: dict[str, str]) -> 'Kingdom':
    return cls(data['name'], {group['name']: Group(group['name'], {subgroup['name']: SubGroup(subgroup['name'], subgroup['id']) for subgroup in group['subgroups']}) for group in data['groups']})
  
  
def to_lower(s: str | Any) -> str:
  """
  Convert a string to lowercase.

  ## Parameters
  ```py
  s : str
  ```
  the string to convert
  """
  match s:
    case str(s):
      return s.lower()
    case _:
      return str(s).lower()

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
    self.__last_save: str = None  # date
    self.__last_build: str = None  # date
    self.__last_name: str = None  # name
    self.__entries: dict[str, Kingdom] = {}

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
    except pd.errors.ParserError:
      error(f'failed to parse {self.__name}.txt (file is malformed)')
      return

    # Kingdom Group SubGroup
    for _, row in df.iterrows():
      kingdom_name: str
      try:
        kingdom_name = row['Kingdom']
      except KeyError:
        kingdom_name = self.__name
      id_name: str
      try:
        id_name = row['Replicons']
      except KeyError:
        id_name = None
      group_name, subgroup_name = row['Group'], row['SubGroup']
      kingdom_name, group_name, subgroup_name = to_lower(kingdom_name), to_lower(group_name), to_lower(subgroup_name)
      kingdom = self.__entries.get(kingdom_name, Kingdom(kingdom_name, {}))
      group = kingdom.groups.get(group_name, Group(group_name, {}))
      subgroup = group.subgroups.get(subgroup_name, SubGroup(subgroup_name, id_name))

      kingdom.groups[group_name] = group
      group.subgroups[subgroup_name] = subgroup
      self.__entries[kingdom_name] = kingdom

    self.__last_build = datetime.now().strftime(datetime_format)

  @property
  def tree(self) -> dict[str, Kingdom]:
    """
    This getter returns the underlying tree. This property is read-only.\\
    You can browse the tree by using the following code:

    ## Example
    ```py
      from src import Tree
      tree = Tree()
      for kingdom in tree.tree.values():
        print(kingdom.name)
        for group in kingdom.groups.values():
          print(f'  {group.name}')
          for subgroup in group.subgroups.values():
            print(f'    {subgroup.name}')
    ```

    ## Returns
    ```py
      dict[str, Kingdom] : a dictionary of kingdoms
    ```
    The dictionary is indexed by the name of the kingdom.
    """
    return self.__entries

  @property
  def filepath(self) -> str:
    """
    This getter returns the filepath of the tree. This property is read-only.\\
    The filepath is the path to the file that contains the tree.\\
    This attribute is only valid once the tree has been saved.

    ## Returns
    ```py
      str : the filepath
    ```
    """
    return f'data/{self.__name}.json' if self.__last_name is None else self.__last_name

  @property
  def last_save(self) -> str:
    """
    This getter returns the date of the last save. This property is read-only.\\
    The date is formatted as `YYYY-MM-DD HH:MM:SS`.\\
    This attribute only remains valid through the lifetime of the object.

    ## Returns
    ```py
      str : the date
    ```
    """
    return self.__last_save

  @property
  def last_build(self) -> str:
    """
    This getter returns the date of the last build. This property is read-only.\\
    The date is formatted as `YYYY-MM-DD HH:MM:SS`.\\
    This attribute only remains valid through the lifetime of the object.

    ## Returns
    ```py
      str : the date
    ```
    """
    return self.__last_build

  def __repr__(self) -> str:
    return f'<Tree {self.__name}>'

  def __str__(self) -> str:
    return self.__name

  def __len__(self) -> int:
    return len(self.__entries)

  def __iter__(self) -> iter:
    return iter(self.__entries.values())
  
  def __getitem__(self, __k: str |int | slice) -> Kingdom | list[Kingdom]:
    match __k:
      case str(__k):
        return self.__entries[__k]
      case int(__k):
        return list(self.__entries.values())[__k]
      case slice(__k):
        return list(self.__entries.values())[__k]

  def __eq__(self, __o: object) -> bool:
    if not isinstance(__o, Tree):
      return False
    return self.__entries == __o.__entries

  def union(self, other: 'Tree') -> 'Tree':
    """
    This function returns a new tree that contains the union of the two trees.\\
    The new tree will contain all the entries of the two trees.

    ## Parameters
    ```py
      other : Tree
    ```
    The other tree to perform the union with.

    ## Returns
    ```py
      Tree : a new tree
    ```
    """
    new_tree = Tree()
    new_tree.__entries = {**self.__entries, **other.__entries}
    return new_tree

  def save(self, name: str = None) -> None:
    """
    This function saves the tree to a file.\\
    The file will be saved in the `data` directory.

    ## Parameters
    ```py
      name : str, (optional)
    ```
    The name of the file to save.\\
    If no name is given, the name of the tree will be used.\\
    defaults to `None`
    """
    name = self.__name if name is None else name
    try:
      os.mkdir('data')
    except FileExistsError:
      pass
    with open(f'data/{name}.json', 'w', encoding='utf-8') as f:
      f.write(json.dumps(self.__entries, indent=2, default=lambda o: o.to_json()))
    self.__last_save = datetime.now().strftime(datetime_format)
    self.__last_name = f'data/{name}.json'

  @classmethod
  def load(cls, name: str) -> 'Tree':
    """
    This function loads a tree from a file.\\
    The file must be in the `data` directory.

    ## Parameters
    ```py
      name : str
    ```
    The name of the file to load.

    ## Returns
    ```py
      Tree : a new tree
    ```
    """
    try:
      with open(f'data/{name}.json', 'r', encoding='utf-8') as f:
        return cls.from_json(json.load(f), name)
    except FileNotFoundError:
      error(f'failed to load {name}.json (file does not exist)')
      return cls()

  @classmethod
  def from_json(cls, data: dict[str, str], name: str = None) -> 'Tree':
    """
    This function creates a new tree from a JSON object.

    ## Parameters
    ```py
      data : dict[str, str]
    ```
    The JSON object.

    ## Returns
    ```py
      Tree : a new tree
    ```
    """
    tree = cls()
    tree.__name = name
    tree.__entries = {k: Kingdom.from_json(data) for k, data in data.items()}
    return tree
