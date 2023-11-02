#!/usr/bin/env python3

'''MiFish
Command line of MiFish
'''

__author__ = 'Tao Zhu'
__copyright__ = 'Copyright 2023'
__license__ = 'GPL'
__email__ = 'zhutao@edu.k.u-tokyo.ac.jp'
__status__ = 'Production'

import argparse
import subprocess
import sys
from operator import itemgetter
from glob import glob
import os

from Bio.Seq import Seq

from mifish.core import pipeline, version
from packaging.version import Version

parser = argparse.ArgumentParser(description='the command line version of MiFish pipeline. \
    It can also be used with any other eDNA meta-barcoding primers', \
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('seq_dir', help='directory of the amplicon sequencing data file (FASTQ/FASTA)')
parser.add_argument('db', help='amplicon reference database in FASTA format, formatted by BLAST')
parser.add_argument('-d', '--other-data-dir', help='other directory of the amplicon sequencing data file \
    (FASTQ/FASTA). Can specify multiple times. Each directory is considered as a group', action='append')
parser.add_argument('-m', '--min-read-len', help='Minimum read length(bp)', type=int, default=204)
parser.add_argument('-M', '--max-read-len', help='Maximum read length(bp)', type=int, default=254)
parser.add_argument('-f', '--primer-fwd', help='forward sequence of primer (5->3)', default='GTCGGTAAAACTCGTGCCAGC')
parser.add_argument('-r', '--primer-rev', help='reverse sequence of primer (5->3)', default='CATAGTGGGGTATCTAATCCCAGTTTG')
parser.add_argument('-u', '--unoise-min', help='value for the -minsize option in UNOISE3', type=int, default=8)
parser.add_argument('-i', '--blast-min-identity', help='Minimum identity (percentage) for filtering BLASTN results', type=float, default=97.0)
parser.add_argument('-s', '--skip-downstream-analysis', help='Skip abandance statics, phylogenetic and bio-diversity analysis', action='store_true')
parser.add_argument('-k', '--keep-tmp-files', help='Keep temporary files', action='store_true')
parser.add_argument('-o', '--output-dir', help='directory for output', default='.')
parser.add_argument('-t', '--threads', help='number of threads for BLASTN and vsearch', type=int, default=2)
parser.add_argument('-v', '--version', action='version', version='%(prog)s '+version.get())
args = parser.parse_args()

def main():
    # check
    if args.other_data_dir is None:
        data_dir_other_groups = []
    else:
        data_dir_other_groups = args.other_data_dir
    if os.path.isfile(args.db+'.nhr') is False:
        print(f'Error: {args.db} does not seem to be a valid database for NCBI BLAST+', file=sys.stderr)
        exit(1)
    for data_dir in ([args.seq_dir]+data_dir_other_groups):
        if os.path.isdir(data_dir) is False:
            print(f'Error: {data_dir} does not seem to be a valid directory', file=sys.stderr)
            exit(1)
    if os.path.isdir(f'{args.output_dir}/MiFishResult') is True:
        print(f'Warning: Directory {args.output_dir}/MiFishResult has already existed. All the files within it would be deleted', file=sys.stderr)
        os.system(f'rm -rf {args.output_dir}/MiFishResult')
    

    for external_bin in 'Gblocks mafft FastTreeMP blastn fastp seqkit vsearch flash'.split():
        app = os.system(f'which {external_bin} >/dev/null 2>&1')
        if app != 0:
            print(f'Error: no {external_bin} in your system', file=sys.stderr)
            exit(1)
    
    app = subprocess.Popen('vsearch --version', shell=True, stdout=subprocess.PIPE, \
                           stderr=subprocess.STDOUT, encoding='utf-8').communicate()[0].split('\n')[0].split()[1].split('_')[0]
    if Version(app) < Version('2.23.0'):
        print(f'ERROR: The version of vsearch(v{app}) is lower than v2.23.0', file=sys.stderr)
        exit(1)
    if os.path.isdir(args.output_dir) is False:
        os.system(f'mkdir {args.output_dir}')
    if os.path.isdir(f'{args.output_dir}/MiFishResult') is False:
        os.system(f'mkdir {args.output_dir}/MiFishResult')
    os.system(f'mkdir {args.output_dir}/MiFishResult/data')

    rm_p_3 = str(Seq(args.primer_rev).reverse_complement())

    pipeline.runMiFish(data_dir=args.seq_dir, data_dir_other_groups=data_dir_other_groups, \
        min_read_length=args.min_read_len, max_read_length=args.max_read_len, \
    rm_p_5=args.primer_fwd, rm_p_3=rm_p_3, blast_identity=args.blast_min_identity, \
        current_task='', current_message='', debug=True, simple_result=args.skip_downstream_analysis, \
            workdir=args.output_dir, cpu=args.threads, db=args.db, unoise_min=args.unoise_min, keep_tmp=args.keep_tmp_files)
    
    print(f'\n\nPipeline Finished, see {args.output_dir}/MiFishResult', file=sys.stderr)
