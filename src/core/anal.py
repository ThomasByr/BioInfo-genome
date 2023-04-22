from io import TextIOWrapper
import os

import string
import random

from abc import ABCMeta, abstractmethod
from typing import Any, TypeVar

from Bio import Entrez, SeqIO

from ..helper import info, debug, error

from pathlib import Path
import re


# yapf: disable
class Comparable(metaclass=ABCMeta):
  @abstractmethod
  def __lt__(self, other: Any) -> bool: ...
# yapf: enable

CT = TypeVar('CT', bound=Comparable)

__all__ = ['create_data_from_NC']


def valid_bounds(bounds, seq_length):
  # Vérification nombre pair de bornes
  n_bounds = len(bounds)
  if n_bounds % 2 != 0:
    return 0

  # Vérification que la borne supérieur soit pas plus haute que la longueur de seq
  if int(bounds[n_bounds - 1]) >= seq_length:
    return 0

  # Vérification de l'accroissement des bornes
  bounds_pred = bounds[0]
  for i in range(1, n_bounds):
    if int(bounds[i]) <= int(bounds_pred):
      return 0
    else:
      bounds_pred = bounds[i]
  return 1


def create_data_from_NC(name: str, path: str, NC_list: list[str], region: list[str]) -> int:
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
  >>> region : list[str]
  ```
  list of regions to create data from (e.g. 'CDS', 'tRNA', 'rRNA', 'gene', ...)

  ## Returns
  ```py
  int : no of regions found
  ```
  """
  print(NC_list)
  path = path + "/"
  letters = string.ascii_lowercase
  local_random = random.Random(hash(name))
  Entrez.email = ''.join(local_random.choice(letters) for _ in range(10)) + '@gmail.com'
  debug(f'searching for [{region}] in [{name}] with NC ids {NC_list}')
  NC_i = 1
  info(f'downloading [{name}]')
  for NC in NC_list:
    info(name + ' - NC : ' + str(NC_i) + ' / ' + str(len(NC_list)))

    # increment here so we do not forget after continue
    NC_i += 1
    debug(f'NC id  = {NC}')
    debug('----------------------------')

    try:
      handle_gb = Entrez.efetch(db="nuccore", id=NC, rettype="gbwithparts",
                                retmode="text")                             # ici on récup le fichier .gb
    except Exception as e:
      if '429' in str(e):
        error('429 : too many requests, please wait a few minutes and try again')
        return 0
      error(f'error while fetching NC id {NC} : {e}')
      continue

    try:
      record = SeqIO.read(handle_gb, 'gb')
      debug(record)
      debug('----------------------------')

      last_header = ''
      n_regions = 0

      region_count = []

      bool_introns = False # savoir si les introns sont demandés
      save_intron = False  # savoir si un intron a été écrit dans le fichier introns.txt (effacer fichier si aucun)
      n_introns_total = 0

      # Ajout des CDS si introns
      delete_cds = False
      if "intron" in region and "CDS" not in region:
        region.insert(0, "CDS")
        delete_cds = True
        # vérification au passage si les introns font partie des régions reconnues
      for k in range(0, len(region)):
        if region[k] == "intron":
          bool_introns = True
        region_count.append(0)

      if not os.path.isdir(path):
        os.mkdir(path)
      all_region = [
        "CDS", "centromere", "intron", "mobile_element", "ncRNA", "rRNA", "telomere", "tRNA", "3'UTR", "5'UTR"
      ]

      filename_suffix = '_' + record.annotations.get('organism',
                                                     record.description.split(',')[0]).replace(
                                                       ' ', '_').replace('/', '_')
      if 'mitochondrion' in record.description.lower():
        filename_suffix += '_Mitochondrion'
      if 'chloroplast' in record.description.lower():
        filename_suffix += '_Chloroplast'
      if 'plasmid' in record.description.lower():
        filename_suffix += '_Plasmid'
      filename_suffix += '_' + record.name + '.txt'

      for k in range(0, len(all_region)):
        if os.path.isfile(path + all_region[k] + filename_suffix):
          pass

      filenames = []
      files: list[TextIOWrapper] = []
      for k in range(0, len(region)):
        filenames.append(path + region[k] + filename_suffix)
        files.append(open(filenames[k], 'w'))

      for k in range(0, len(record.features)):
        if record.features[k].type in region:
          n_regions += 1
          region_count[region.index(record.features[k].type)] += 1
          if record.features[k].location.strand == -1:
            prefix = " complement( "
            suffix = ") "
          else:
            prefix = " "
            suffix = ""

          header = record.features[k].type + ' ' + record.annotations.get(
            'organism',
            record.description.split(',')[0]) + ' ' + record.name + ':'
          header2 = prefix + str(record.features[k].location) + suffix
          header2 = header2.replace('{', "( ")
          header2 = header2.replace('}', ') ')
          header2 = header2.replace(':', '..')
          header2 = header2.replace(']', '')
          header2 = header2.replace('(-)', ' ')
          header2 = header2.replace('(+)', ' ')
          header2 = header2.replace('[', '')
          header2 = header2.replace('<', '')
          header2 = header2.replace('>', '')
          if header2 != last_header:

            last_header = header2
            ## Écriture Séquences avec exons et introns
            ## récupération des bornes
            bounds_tmp = re.findall('\d+', header2)
            ## Détection de si complément
            complement = bool(re.search(r'\bcomplement\b', header2))
            join = bool(re.search(r'\bjoin\b', header2)) # pas besoin d'introns si ya pas de join

            bounds = []
            if complement == True:
              for i in range(len(bounds_tmp) - 1, -1, -2):
                bounds.append(bounds_tmp[i - 1])
                bounds.append(bounds_tmp[i])
            else:
              bounds = bounds_tmp

            n_intron = 1
            n_exon = 1

            ## Correction du header pour borne inf
            for i in range(len(bounds) - 1, -1, -2):
              bounds_inf = int(bounds[i - 1])
              bounds_inf += 1
              bounds_sup = int(bounds[i])
              header2 = header2.replace(' ' + bounds_tmp[i - 1] + '.', ' ' + str(bounds_inf) + '.', 1)
              header2 = header2.replace('.' + bounds_tmp[i] + ' ', '.' + str(bounds_sup) + ' ', 1)

            if valid_bounds(bounds, len(record.seq)) == 1:
              n_regions += 1
              region_count[region.index(record.features[k].type)] += 1
              ## CDS termine par TAG TAA TGA

              files[region.index(
                record.features[k].type)].write(header + header2 + '\n' +
                                                str(record.features[k].extract(record.seq)) + '\n')
              # print(header,'\n')
              ## Cherche les Exons et introns
              if (join == True):
                if (bool_introns == True) and (record.features[k].type == "CDS"):
                  save_intron = True
                  # files[region.index("intron")].write(header + header2 + '\n')
                  for i in range(0, len(bounds) - 1):
                    if i % 2 == 0:
                      # exon
                      i_header = header2 + 'Exon ' + str(n_exon)

                      # print(i_header)
                      bounds_inf = int(bounds_tmp[i])
                      bounds_sup = int(bounds_tmp[i + 1])
                      files[region.index(record.features[k].type)].write(header + i_header + '\n')
                      if complement:
                        files[region.index(record.features[k].type)].write(
                          str(record.seq[bounds_inf:bounds_sup].reverse_complement()) + '\n')
                      else:
                        files[region.index(
                          record.features[k].type)].write(str(record.seq[bounds_inf:bounds_sup]) + '\n')

                      n_exon += 1
                    else: # Introns
                      i_header = header2 + 'Intron ' + str(n_intron)

                      if complement:
                        bounds_sup = int(bounds_tmp[i - 1])
                        bounds_inf = int(bounds_tmp[i + 2])
                        # ,':',bounds_inf,'..',bounds_sup)
                        files[region.index("intron")].write(
                          header.replace('CDS', 'intron', 1) + i_header + '\n')
                        files[region.index("intron")].write(
                          str(record.seq[bounds_inf:bounds_sup].reverse_complement()) + '\n')
                      else:
                        bounds_inf = int(bounds_tmp[i])
                        bounds_sup = int(bounds_tmp[i + 1])
                        # ,':',bounds_inf,'..',bounds_sup)
                        files[region.index("intron")].write(
                          header.replace('CDS', 'intron', 1) + i_header + '\n')
                        files[region.index("intron")].write(str(record.seq[bounds_inf:bounds_sup]) + '\n')

                      n_intron += 1
                      n_introns_total += 1
                else:
                  for i in range(0, len(bounds) - 1, 2):
                    i_header = header2 + 'Exon ' + str(n_exon)
                    bounds_inf = int(bounds_tmp[i])
                    bounds_sup = int(bounds_tmp[i + 1])
                    files[region.index(record.features[k].type)].write(header + i_header + '\n')
                    if complement:
                      files[region.index(record.features[k].type)].write(
                        str(record.seq[bounds_inf:bounds_sup].reverse_complement()) + '\n')
                    else:
                      files[region.index(
                        record.features[k].type)].write(str(record.seq[bounds_inf:bounds_sup]) + '\n')

                    n_exon += 1

    except Exception as e:
      error(f'{e}')
      pass

  if n_regions == 0:
    info(f'Selected functional region not found for organism : [{name}]')
    return 0
  info(f'{name} downloaded successfully ({n_regions})')
  print(filenames)
  return n_regions


# wtf is this atrocious code and why is it here ??
def check_inf_sup(inf: CT, sup: CT) -> bool:
  if (int(inf) <= int(sup)):
    return True
  else:
    return False
