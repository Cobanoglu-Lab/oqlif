# oqlif <img src="https://raw.githubusercontent.com/Cobanoglu-Lab/oqlif/master/leaf.png" width="36">
Quantifies scRNA-seq BAM files by taking multi-maps into account.

Takes alignment and single cell barcode annotated BAM files as input, 
and produces MatrixMarket coordinate format barcode-by-gene read tables.
This is a popular output format for scRNA-seq, commonly employed as input
and output for multiple relevant software tools. 

Since `oqlif` quantifies multi-maps, the "read counts" include partial reads.
Hence the default output is real valued, not integer valued. If you prefer
an integer output, pass the `--integer` flag, which floors the output. 

If you have a multi-core system, like most modern machines, the `--parallel`
flag enables parallelized execution. 

    usage: oqlif.py [-h] [--ens2name ENS2NAME] [--label LABEL] [--bcode BCODE]
                    [--outdir OUTDIR] [--integer] [--parallel]
                    bamfile

    Quantify reads from BAM files

    positional arguments:
      bamfile              The BAM file that stores the alignment results

    optional arguments:
      -h, --help           show this help message and exit
      --ens2name ENS2NAME  A tab-separated text file with only two columns and no
                           header that maps ENS identifiers to their human
                           readable names. You can download from Ensembl BioMart.
                           Assumes first column is the first column is the human-
                           readable name, second column is the Ensembl identifier.
                           If this argument is missing, only Ensembl identifiers
                           reported in genes.tsv
      --label LABEL        The label for the quantification target. For example,
                           set to TX to quantify transcripts. Default: GX
      --bcode BCODE        The label for the tag that holds the barcode
                           information. Default: CB
      --outdir OUTDIR      The folder where output goes. If folder does not exist,
                           it will be created
      --integer            If flag is set, counts are cast to integer (by
                           flooring) for discrete counts. If not set, partial
                           counts are output after rounding to 2 decimal places.
      --parallel           If flag is set, parallelizes execution over system
                           cores. Runs chromosomes in parallel.
                           


Contact: Murat Can Cobanoglu, PhD

https://www.utsouthwestern.edu/labs/cobanoglu/
