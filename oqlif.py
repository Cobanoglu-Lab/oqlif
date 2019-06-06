####################################
# Copyright © 2019, University of Texas Southwestern Medical Center. All rights reserved.
# Contributors: Murat Can Cobanoglu
# Department: Lyda Hill Department of Bioinformatics
####################################

import os
import argparse
from collections import defaultdict
import multiprocessing as mp
import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm
from scipy import sparse
from scipy.io import mmwrite


def get_expr_dict(input_tuple):
  """Returns the quantification results as a dict of dicts."""

  fname, contig, start, stop, length, worker_ind = input_tuple
  if worker_ind == None: worker_ind = 0
  if contig is None: 
    desc_str='all_reads'
  else:
    desc_str='chr_'+contig

  region = dict(
    contig=contig,
    start=start,
    stop=stop
  )
  
  expr = {}
  with pysam.AlignmentFile(fname, "rb") as samfile:
    for read in tqdm(samfile.fetch(until_eof=True, multiple_iterators=True, **region), 
                      total=length, position=worker_ind, leave=True, desc=desc_str):
      tags = read.get_tags()
      if tags is not None:
        tags = dict(tags)
        if bcode in tags and label in tags and tags[label] is not None and len(tags[label])>0:
          ens_ids = [tag.split(',')[0] for tag in tags[label].split(';')]
          bc = read.get_tag(bcode)
          if bc not in expr:
            expr[bc] = defaultdict(float)
          for ens in ens_ids:
            expr[bc][ens] += 1./len(ens_ids)
  return expr


def merge_expr_dicts(exprs):
  """Merges a list of expr dicts. Reduce step for parallel exec"""
  expr = {}
  for expr_curr in exprs:
    for bc in expr_curr.keys():
      if bc not in expr:
        expr[bc] = defaultdict(float)
      for g in expr_curr[bc].keys():
        expr[bc][g] += expr_curr[bc][g]
  return expr


def read_ens2name(fname):
  """Returns a dictionary that maps ensembl id to readable names"""
  if args.ens2name != 'n/a':
    with open(args.ens2name) as f:
      ens2name={}
      for line in f:
        try:
          v,k = line.strip().split()
          ens2name[k] = v
        except ValueError:
          print(line)
    return ens2name
  else:
    return None


def get_ind(d, bc, curr): 
  """Faster indexing"""
  try:
    ind = d[bc]
  except KeyError:
    ind = curr
    d[bc] = curr
    curr += 1
  return d, ind, curr


def write(expr, outd, ens2name=None, integer=False):
  """Writes the quantification results to standard format files"""
  
  bc_d = {}
  bc_curr = 1
  
  gene_d = {}
  gene_curr = 1
  all_genes = set()
  tot_entries = 0

  with open(outd+os.sep+'barcodes.tsv','w',newline='') as bc_f:
    for bc in sorted(expr.keys()):
      bc_d, bc_ind, bc_curr = get_ind(bc_d, bc, bc_curr)
      bc_f.write(bc+'\n')
      for g in sorted(expr[bc].keys()):
        all_genes.add(g)
        tot_entries += 1
    all_genes = sorted(list(all_genes))

  with open(outd+os.sep+'genes.tsv','w',newline='') as gene_f:
    if ens2name is not None:
      gene_names = [ens2name[g] if g in ens2name else g for g in all_genes]
      for gene, name in zip(all_genes, gene_names):
        gene_d, gene_ind, gene_curr = get_ind(gene_d, gene, gene_curr)
        gene_f.write(gene+'\t'+name+'\n')
    else:
      for gene in all_genes:
        gene_d, gene_ind, gene_curr = get_ind(gene_d, gene, gene_curr)
        gene_f.write(gene+'\n')

    with open(outd+os.sep+'matrix.mtx','w',newline='') as mtx_f:
      type_str = "integer" if integer else "real" 
      mtx_f.write("%%MatrixMarket matrix coordinate {0:s} general\n%\n".format(type_str))
      mtx_f.write("{0:d} {1:d} {2:d}\n".format(len(all_genes), len(bc_d.keys()), tot_entries))
      for bc in sorted(expr.keys()):
        for g in sorted(expr[bc].keys()):
          if integer:
            count = str(np.floor(expr[bc][g]).astype(np.int))
          else: 
            count = '{0:.3f}'.format(expr[bc][g])
          mtx_f.write(str(gene_d[g])+' '+str(bc_d[bc])+' '+count+'\n')

if __name__ == '__main__':

  parser = argparse.ArgumentParser(description='Quantify reads from BAM files')

  parser.add_argument('bamfile', help='The BAM file that stores the alignment results')

  parser.add_argument('--ens2name', default='n/a', 
    help='A tab-separated text file with only two columns and no header that maps ENS'+\
      ' identifiers to their human readable names. You can download from Ensembl BioMart. '+\
      'Assumes first column is the first column is the human-readable name, second column '+\
      'is the Ensembl identifier. If this argument is missing, only Ensembl identifiers '+\
      'reported in genes.tsv')

  parser.add_argument('--label', default='GX', 
    help='The label for the quantification target. For example, set to TX to quantify transcripts.'+\
      '\n\tDefault: GX')

  parser.add_argument('--bcode', default='CB', 
    help='The label for the tag that holds the barcode information. Default: CB')

  parser.add_argument('--outdir', default='oqlif_gene_bc_matrices', \
    help='The folder where output goes. If folder does not exist, it will be created')

  parser.add_argument('--integer', action='store_true',
    help='If flag is set, counts are cast to integer (by flooring) for discrete counts. '+\
      'If not set, partial counts are output after rounding to 2 decimal places.')

  parser.add_argument('--parallel', action='store_true',
    help='If flag is set, parallelizes execution over system cores. Runs chromosomes in parallel.') 

  args = parser.parse_args()
  label = args.label
  bcode = args.bcode
  outd = args.outdir
  bamfile = args.bamfile
  ens2name = read_ens2name(args.ens2name)

  if not os.path.isdir(outd):
    os.mkdir(outd)

  contigs = []; lengths = []
  with pysam.AlignmentFile(bamfile) as samfile:
    for stat in samfile.get_index_statistics():
      contig = stat.contig
      if contig.isnumeric() or contig=='MT' or contig=='X' or contig=='Y':
        contigs.append(contig)
        lengths.append(stat.total)
  inputs = [(bamfile, contig, None, None, lengths[i], i) for i,contig in enumerate(contigs)]
  if args.parallel:
    with mp.Pool() as pool:
      exprs = pool.map(get_expr_dict, inputs)
  else:
    exprs = list(map(get_expr_dict, inputs))
  expr = merge_expr_dicts(exprs)
  write(expr, outd, ens2name, args.integer)
