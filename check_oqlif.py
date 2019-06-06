####################################
# Copyright © 2019, University of Texas Southwestern Medical Center. All rights reserved.
# Contributors: Murat Can Cobanoglu
# Department: Lyda Hill Department of Bioinformatics
####################################

import argparse
import os
import numpy as np
import pandas as pd
from scipy.io import mmread

ens2name = dict()
with open('mart_export.txt') as in_f:
  for line in in_f:
    v,k = line.strip().split('\t')
    ens2name[k]=v

def read_expd(folder):
  exp = mmread(folder+"matrix.mtx").todense()
  bc = pd.read_csv(folder+"barcodes.tsv", header=None)
  genes = pd.read_csv(folder+'genes.tsv','\t', header=None)
  genes[1] = [ens2name[k] if k in ens2name else k for _,k in genes[0].iteritems()]
  expd = pd.DataFrame(data=exp, index=genes[1], columns=bc)
  return expd

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Check equality of different oqlif output folders')

  parser.add_argument('folders', nargs='+')
  args = parser.parse_args()
  folders = [f+os.sep if not f.endswith(os.sep) else f for f in args.folders]
  if len(folders) < 2:
    raise ValueError('At least two folders required')

  results = list(map(read_expd, folders))
  for i, f1 in enumerate(folders):
    for j, f2 in enumerate(folders[i+1:]):
      j = (i+1)+j
      print('{2:d}: {0:s}\tvs\t{3:d}: {1:s}'.format(f1, f2, i, j))
      good_shape = (results[i].shape == results[j].shape)
      print('Shape equality: {0:s} | Shapes: {1:s} ; {2:s}'.format(
          str(good_shape), str(results[i].shape), str(results[j].shape)))
      good_index = (np.all(sorted(results[i].index)==sorted(results[j].index)))
      print('Index equality: {0:s}'.format(
          str(good_index) ))
      good_column = (np.all(sorted(results[i].columns)==sorted(results[j].columns)))
      print('Column equality: {0:s}'.format(
          str(good_column) ))
      if good_shape and good_index and good_column:
        e1 = np.floor(results[i].loc[sorted(results[i].index),sorted(results[i].columns)].values).astype(np.uint)
        e2 = np.floor(results[j].loc[sorted(results[j].index),sorted(results[j].columns)].values).astype(np.uint)
        not_eq_idx = np.where(e1!=e2)
        print('Int-cast non-equal entry count: {0:d}'.format(len(not_eq_idx[0])))
        if len(not_eq_idx[0]) > 0:
          good_matrix = False
          outfname = f1.replace(os.sep,'_')+'__vs__'+f2.replace(os.sep,'_')
          with open(outfname,'w') as outf:
            for c1, c2 in zip(*not_eq_idx):
              outf.write('{0:d}, {1:d}, {2:d}, {3:d}\n'.format(c1,c2,e1[c1,c2],e2[c1,c2]))
        else:
          good_matrix = True
      if good_shape and good_index and good_column and good_matrix:
        print('PERFECT MATCH: {0:s}\tvs\t{1:s}\n\n'.format(f1, f2))
      else:
        print('DIFFERENCE DETECTED: {0:s}\tvs\t{1:s}\n\n'.format(f1, f2))
