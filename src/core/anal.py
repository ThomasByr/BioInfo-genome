import os

import string
import random

from ..helper import info, debug

from Bio import Entrez, SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

__all__ = []


def create_data_from_NC(name: str, path: str, NC_list: list[str], region: str) -> int:
  # todo: add docstring

  letters = string.ascii_lowercase
  Entrez.email = ''.join(random.choice(letters) for _ in range(10)) + '@gmail.com'
  NC_i = 1
  no_region_found = 0
  info()
  info('downloading [' + name + ']')
  for NC in NC_list:
    info('NC : ' + str(NC_i) + ' / ' + str(len(NC_list)))
    name = name.replace(' ', '_')
    name = name.replace('[', '_')
    name = name.replace(']', '_')
    name = name.replace(':', '_')
    NC_i += 1
    debug('NC id  =', NC)
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
    for i in range(len(record[0]['GBSeq_feature-table'])):
      info('\tfeature : ' + str(i + 1) + ' / ' + str(len(record[0]['GBSeq_feature-table'])))
      feature_location = record[0]['GBSeq_feature-table'][i]['GBFeature_location']
      feature_key = record[0]['GBSeq_feature-table'][i]['GBFeature_key']
      if feature_key != region:
        continue
      no_region_found += 1
      NC_filename = str(name) + '_' + feature_key + '_NC_' + str(NC_i) + '.txt'

      file_path = os.path.join(path, name, NC_filename)
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
        debug(i + 1, '/', len(record[0]['GBSeq_feature-table']))
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
  info(f'{name} downloaded')
  return no_region_found


def check_inf_sup(inf, sup) -> bool:
  if (inf <= sup):
    return True
  else:
    return False
