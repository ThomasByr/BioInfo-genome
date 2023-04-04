import os

import string
import random

from abc import ABCMeta, abstractmethod
from typing import Any, TypeVar

from Bio import Entrez, SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

from ..helper import info, debug


# yapf: disable
class Comparable(metaclass=ABCMeta):
  @abstractmethod
  def __lt__(self, other: Any) -> bool: ...
# yapf: enable

CT = TypeVar('CT', bound=Comparable)

__all__ = ['create_data_from_NC']


def create_data_from_NC(name: str, path: str, NC_list: list[str], region: str) -> int:
  """
  create data from list of NC ids and region\\
  saves data in a file in the given path

  ## Parameters
  ```py
  >>> name : str
  ```
  name of the organism to create data from
  ```py
  >>> path : str
  ```
  path to save the data in (must be a directory, .txt files will be created in it)
  ```py
  >>> NC_list : list[str]
  ```
  list of NC ids `['NC_000001', 'NC_000002', ...]`
  ```py
  >>> region : str
  ```
  region to create data from (e.g. 'CDS', 'tRNA', 'rRNA', 'gene', ...)

  ## Returns
  ```py
  int : no of regions found
  ```
  """
  letters = string.ascii_lowercase
  Entrez.email = ''.join(random.choice(letters) for _ in range(10)) + '@gmail.com'
  NC_i = 1
  no_region_found = 0
  info(f'downloading [{name}]')
  for NC in NC_list:
    info('NC : ' + str(NC_i) + ' / ' + str(len(NC_list)))
    name = name.replace(' ', '_')
    name = name.replace('[', '_')
    name = name.replace(']', '_')
    name = name.replace(':', '_')
    NC_i += 1              # increment here so we do not forget after continue
    debug(f'NC id  = {NC}')
    debug('----------------------------')
    handle_fasta = Entrez.efetch(db='nucleotide', id=NC, rettype='fasta', retmode='text')
    record_fasta = SeqIO.read(handle_fasta, 'fasta')
    debug(record_fasta)
    debug('----------------------------')
    handle_fasta.close()
    handle_text = Entrez.efetch(db='nucleotide', id=NC, retmode='xml')
    record = Entrez.read(handle_text)
    handle_text.close()
    list_file = []
    no_total_features = len(record[0]['GBSeq_feature-table'])
    info(f'\tn° total features : {no_total_features}')
    for i in range(no_total_features):
      feature_location = record[0]['GBSeq_feature-table'][i]['GBFeature_location']
      if (feature_key := record[0]['GBSeq_feature-table'][i]['GBFeature_key']) != region:
        continue
      no_region_found += 1 # we found a region

      # make sure we decrement NC_i because we incremented it at the beginning of the loop
      NC_filename = str(name) + '_' + feature_key + '_NC_' + str(NC_i - 1) + '.txt'

      file_path = os.path.join(path, NC_filename)
      if len(list_file) != 0:
        if NC_filename not in list_file:
          os.remove(file_path)
          list_file.append(NC_filename)
      else:
        try:
          os.remove(file_path)
        except FileNotFoundError:
          pass
        list_file.append(NC_filename)

      with open(file_path, 'a+') as out:
        debug(f'{i + 1}/' + str(len(record[0]['GBSeq_feature-table'])))
        debug(feature_location)
        # todo: tests on regions (part 2.3)
        out.write(feature_key + ' ' + feature_location + '\n')
        if feature_location.find('complement') != -1 and feature_location.find('join') != -1:
          feature_location = feature_location[16:-1]
          x = feature_location.split(',')
          fn = []
          is_valid = True
          for xi in x:
            xi = xi.split('..')
            try:
              (int(xi[0]), int(xi[1]))
            except ValueError:
              is_valid = False
            else:
              if (check_inf_sup(xi[0], xi[1]) == False):
                is_valid = False
              fn.append(FeatureLocation(int(xi[0]), int(xi[1])))
          if not is_valid:
            continue
          f = CompoundLocation(fn)
          debug('COMPLEMENT JOIN')
          debug(f.extract(record_fasta.seq).complement())
          out.write(str(f.extract(record_fasta.seq).complement()))

        elif feature_location.find('complement') != -1:
          feature_location = feature_location[11:-1]
          x = feature_location.split('..')
          try:
            (int(x[0]), int(x[1]))
          except ValueError:
            continue
          else:
            if (check_inf_sup(x[0], x[1]) == False):
              continue
            f = SeqFeature(FeatureLocation(int(x[0]), int(x[1])), type='domain')
            debug('COMPLEMENT')
            debug(f.extract(record_fasta.seq).complement())
            out.write(str(f.extract(record_fasta.seq).complement()))

        elif feature_location.find('join') != -1:
          feature_location = feature_location[5:-1]
          x = feature_location.split(',')
          fn = []
          is_valid = True
          for xi in x:
            xi = xi.split('..')
            try:
              (int(xi[0]), int(xi[1]))
            except ValueError:
              is_valid = False
            else:
              if (check_inf_sup(xi[0], xi[1]) == False):
                is_valid = False
              fn.append(FeatureLocation(int(xi[0]), int(xi[1])))
          if not is_valid:
            continue
          f = CompoundLocation(fn)
          debug('COMPLEMENT JOIN')
          debug(f.extract(record_fasta.seq).complement())
          out.write(str(f.extract(record_fasta.seq).complement()))

        elif feature_location.find('complement') != -1:
          feature_location = feature_location[11:-1]
          x = feature_location.split('..')
          try:
            (int(x[0]), int(x[1]))
          except ValueError:
            continue
          else:
            if (check_inf_sup(x[0], x[1]) == False):
              continue
            f = SeqFeature(FeatureLocation(int(x[0]), int(x[1])), type='domain')
            debug('COMPLEMENT')
            debug(f.extract(record_fasta.seq).complement())
            out.write(str(f.extract(record_fasta.seq).complement()))

        elif feature_location.find('join') != -1:
          feature_location = feature_location[5:-1]
          x = feature_location.split(',')
          fn = []
          is_valid = True
          for xi in x:
            xi = xi.split('..')
            try:
              (int(x[0]), int(x[1]))
            except ValueError:
              is_valid = False
            else:
              if (check_inf_sup(xi[0], xi[1]) == False):
                is_valid = False
              fn.append(FeatureLocation(int(xi[0]), int(xi[1])))
          if not is_valid:
            continue
          f = CompoundLocation(fn)
          debug('JOIN')
          debug(f.extract(record_fasta.seq))
          out.write(str(f.extract(record_fasta.seq)))

        else:
          x = feature_location.split('..')
          try:
            (int(x[0]), int(x[1]))
          except ValueError:
            continue
          else:
            if (check_inf_sup(x[0], x[1]) == False):
              continue
            f = SeqFeature(FeatureLocation(int(x[0]), int(x[1])), type='domain')
            debug('EXTRACT')
            debug(f.extract(record_fasta.seq))
            out.write(str(f.extract(record_fasta.seq)))
        out.write('\n')
  if no_region_found == 0:
    info(f'Selected functional region not found for organism : [{name}]')
    return 0
  info(f'{name} downloaded successfully')
  return no_region_found


# wtf is this atrocious code and why is it here ??
def check_inf_sup(inf: CT, sup: CT) -> bool:
  if (inf <= sup):
    return True
  else:
    return False
