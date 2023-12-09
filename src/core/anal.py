import os
import logging

import string
import random
from io import TextIOWrapper

from Bio import Entrez, SeqIO

import re

__all__ = ['create_data_from_stuff']

logger = logging.getLogger()


def valid_bounds(bounds, seq_length) -> bool:
  # check parity
  n_bounds = len(bounds)
  if n_bounds % 2 != 0:
    return False

  # is sup not greater than seq_length ?
  if int(bounds[n_bounds - 1]) >= seq_length:
    return False

  # check if bounds are in order
  bounds_prev = bounds[0]
  for i in range(1, n_bounds):
    if int(bounds[i]) <= int(bounds_prev):
      return False
    else:
      bounds_prev = bounds[i]
  return True


def create_data_from_stuff(name: str, path: str, NC_list: list[str], region: list[str]) -> int:
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
  letters = string.ascii_lowercase
  local_random = random.Random(hash(name))
  Entrez.email = ''.join(local_random.choice(letters) for _ in range(10)) + '@gmail.com'
  logger.debug(f'searching for [{region}] in [{name}] with NC ids {NC_list}')
  NC_i = 1
  logger.info(f'downloading [{name}]')
  for NC in NC_list:
    logger.info(name + ' - NC : ' + str(NC_i) + ' / ' + str(len(NC_list)))

    # increment here so we do not forget after continue
    NC_i += 1
    logger.debug(f'NC id  = {NC}')
    logger.debug('----------------------------')

    try:
      handle_gb = Entrez.efetch(db='nuccore', id=NC, rettype='gbwithparts', retmode='text')
    except Exception as e:
      if '429' in str(e):
        logger.error('429 : too many requests, please wait a few minutes and try again')
        return 0
      logger.error(f'error while fetching NC id {NC} : {e}')
      continue

    try:
      record = SeqIO.read(handle_gb, 'gb')
      logger.debug(record)
      logger.debug('----------------------------')

      last_header = ''
      n_regions = 0

      region_count = []

      bool_introns = False # to know if introns are in region
      save_intron = False  # suppose no intron are saved in file (delete intron file if empty)
      n_introns_total = 0

      # if region is intron, add CDS to region
      delete_cds = False # this is to delete CDS files after if we added it
      if 'intron' in region and 'CDS' not in region:
        region.insert(0, 'CDS')
        delete_cds = True

      # are there introns ?
      for k in range(0, len(region)):
        if region[k] == 'intron':
          bool_introns = True
        region_count.append(0)

      if not os.path.isdir(path):
        os.mkdir(path)

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

      filenames: list[str] = []
      files: list[TextIOWrapper] = []
      for k in range(0, len(region)):
        filenames.append(os.path.join(path, region[k] + filename_suffix))
        files.append(open(filenames[k], 'w'))

      for k in range(0, len(record.features)):
        if record.features[k].type in region:
          n_regions += 1
          region_count[region.index(record.features[k].type)] += 1
          if record.features[k].location.strand == -1:
            prefix = ' complement( '
            suffix = ') '
          else:
            prefix = ' '
            suffix = ''

          header = record.features[k].type + ' ' + record.annotations.get(
            'organism',
            record.description.split(',')[0]) + ' ' + record.name + ':'
          header2 = prefix + str(record.features[k].location) + suffix
          header2 = header2.replace('{', '( ')
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
            # sequences writing with exons and introns
            # bounds recovery
            bounds_tmp: str = re.findall('\d+', header2)
            # complement detection
            complement = bool(re.search(r'\bcomplement\b', header2))
            join = bool(re.search(r'\bjoin\b', header2)) # no need for introns if no join

            bounds: list[str] = []
            if complement:
              for i in range(len(bounds_tmp) - 1, -1, -2):
                bounds.append(bounds_tmp[i - 1])
                bounds.append(bounds_tmp[i])
            else:
              bounds = bounds_tmp

            n_intron = 1
            n_exon = 1

            # header correction for inf bound
            for i in range(len(bounds) - 1, -1, -2):
              bounds_inf = int(bounds[i - 1])
              bounds_inf += 1
              bounds_sup = int(bounds[i])
              header2 = header2.replace(' ' + bounds_tmp[i - 1] + '.', ' ' + str(bounds_inf) + '.', 1)
              header2 = header2.replace('.' + bounds_tmp[i] + ' ', '.' + str(bounds_sup) + ' ', 1)

            if valid_bounds(bounds, len(record.seq)):
              n_regions += 1
              region_count[region.index(record.features[k].type)] += 1
              # CDS ends with w/ TAG TAA TGA

              files[region.index(
                record.features[k].type)].write(header + header2 + '\n' +
                                                str(record.features[k].extract(record.seq)) + '\n')

              # search for exons and introns
              if (join == True):
                if (bool_introns == True) and (record.features[k].type == 'CDS'):
                  save_intron = True
                  # files[region.index('intron')].write(header + header2 + '\n')
                  for i in range(0, len(bounds) - 1):
                    if i % 2 == 0:
                      # exon
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

                    else:
                      # introns
                      i_header = header2 + 'Intron ' + str(n_intron)

                      if complement:
                        bounds_sup = int(bounds_tmp[i - 1])
                        bounds_inf = int(bounds_tmp[i + 2])
                        # ,':',bounds_inf,'..',bounds_sup)
                        files[region.index('intron')].write(
                          header.replace('CDS', 'intron', 1) + i_header + '\n')
                        files[region.index('intron')].write(
                          str(record.seq[bounds_inf:bounds_sup].reverse_complement()) + '\n')
                      else:
                        bounds_inf = int(bounds_tmp[i])
                        bounds_sup = int(bounds_tmp[i + 1])
                        # ,':',bounds_inf,'..',bounds_sup)
                        files[region.index('intron')].write(
                          header.replace('CDS', 'intron', 1) + i_header + '\n')
                        files[region.index('intron')].write(str(record.seq[bounds_inf:bounds_sup]) + '\n')

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
      logger.error(f'{e}')
      pass

  if n_regions == 0:
    logger.info(f'Selected functional region not found for organism : [{name}]')
    return 0
  logger.info(f'{name} downloaded successfully ({n_regions})')
  return n_regions
