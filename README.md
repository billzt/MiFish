# MiFish
This is the command line version of MiFish pipeline. It can also be used with any other eDNA meta-barcoding primers

# External Dependencies
Add these softwares to your system PATH
* [fastp](https://github.com/OpenGene/fastp)
* [FLASH](http://ccb.jhu.edu/software/FLASH/) (v1.2.7)
* [seqkit](https://github.com/shenwei356/seqkit/) (v2.3.0)
* [usearch](https://www.drive5.com/usearch/) (v11.0.667)
* [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) (v2.9.0)
* [MAFFT](https://mafft.cbrc.jp/alignment/software/) (v7.505)
* [Gblocks](https://home.cc.umanitoba.ca/~psgendb/doc/Castresana/Gblocks_documentation.html) (v0.91b)
* [FastTreeMP](http://www.microbesonline.org/fasttree/) (v2.1.11)

# Install
Currently we only support Linux. Please use conda to manage the environment
```
$ conda create -n MiFish python==3.9.13
$ conda activate MiFish
$ pip3 install numpy==1.23.1
$ pip3 install scikit-bio==0.5.6
$ pip3 install PyQt5==5.15.7
$ pip3 install ete3==3.1.2
$ pip3 install duckdb==0.6.1
$ pip3 install XlsxWriter==3.0.3
$ pip3 install cutadapt==4.1
$ pip3 install biopython==1.79
$ git clone https://github.com/billzt/MiFish.git
$ cd MiFish
$ python3 setup.py develop
$ mifish -h
```

# Test
```
$ cd test
$ mifish seq mifishdbv3.83.fa -d seq2
```
There are [six files](#results) in the result directory `MiFishResult`

# Parameters
## Mandatory
```
$ mifish /path/to/your/amplicon/sequencing/directory/ /path/to/your/ref/db/
```
### Directory for amplicon sequencing data (FASTQ/FASTA)
Since MiFish supports multi-sample analysis, amplicon sequencing data in compressed FASTQ/FASTA format should be put in directories. Pass the path of the directory as the first parameter.
Refer to [MiFish's Homepage](http://mitofish.aori.u-tokyo.ac.jp/mifish/) to see the rules of filenames. Here are some examples:
* `MiFish-example-02_S73_L001_R1_001.fastq.gz`
* `MiFish-example-02_S73_L001_R2_001.fastq.gz`
* `DRR126155_1.fastq.bz2`
* `DRR126155_2.fastq.bz2`
* `mydata.1.fq.xz`
* `mydata.2.fq.xz`

### RefDB of your metabarcoding primers
Prepare your RefDB in FASTA format and index it using the `makeblastdb` from [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi). RefDB for an old version of MiFish is in `test/mifishdbv3.83.fa`

The head line of RefDB (FASTA) follows this rule:

`gb|accessionID|species_scientific_name`

Replace blanks with underscores in the species name. Here are examples.

```
>gb|LC021149|Ostorhinchus_angustatus
CACCGCGGTTATACGAGAGGCCCAAGCTGACAATCACCGGCGTAAAGAGTGGTTAATGAC
CCCACAATAATAAAGTCGAACATCTCCAAAGTTGTTGAACACATTCGAAGATATGAAGCT
CTACCACGAAAGTGACTTTACACTCTTTGAACCCACGAAAGCTAGGAAA
>gb|LC579122|Ostorhinchus_angustatus
CACCGCGGTTATACGAGGGGCCCAAGCTGACAATCACCGGCGTAAAGAGTGGTTAATAAC
CCCACAATAATAAAGTCGAACATCTCCAAAGTTGTTGAACACATTCGAAGATATGAAGCT
CTACCACGAAAGTGACTTTACACTCTTTGAACCCACGAAAGCTAGGAAA
>gb|LC717543|Trachidermus_fasciatus
CACCGCGGTTATACGAGAGACTCAAGCTGACAAACACCGGCGTAAAGCGTGGTTAAGCTA
AAAATTTGCTAAAGTCAAACACCTTCAAGACTGTTATACGTACCCGAAGGCAGGAAGCAC
AACCACGAAAGTGACTTTAACTAAGCTGAATCCACGAAAGCTAAGGAA
```

accessionID can be any unique strings. Primers were trimmed off from the sequences.

## Optional (important❗️)
Following optional parameters are designed for MiFish metabarcoding primers. If running with other eDNA primers, change them to satisfy your own primers.

### Length filtering
```
  -m MIN_READ_LEN, --min-read-len MIN_READ_LEN
                        Minimum read length(bp) (default: 204)

  -M MAX_READ_LEN, --max-read-len MAX_READ_LEN
                        Maximum read length(bp) (default: 254)
```
adjust them to satisfy your reads. If uncertain, set `-m 0 -M 99999` to keep all the reads

### Primer sequences
```
  -f PRIMER_FWD, --primer-fwd PRIMER_FWD
                        forward sequence of primer (5->3) (default: GTCGGTAAAACTCGTGCCAGC)
  -r PRIMER_REV, --primer-rev PRIMER_REV
                        reverse sequence of primer (5->3) (default: CATAGTGGGGTATCTAATCCCAGTTTG)
```
change them according to your own primers

## Optional 
Following optional parameters are designed for all metabarcoding primers.

### Group samples
```
  -d OTHER_DATA_DIR, --other-data-dir OTHER_DATA_DIR
                        other directory of the amplicon sequencing data file (FASTQ/FASTA). Can specify multiple times. Each directory is considered as a group (default: None)
```
If your samples are in multiple groups, please arrange them in different directories and use the `-d` parameter for multiple times. e.g. `-d 2nd_group_dir -d 3rd_group_dir`

### Threshold of BLASTN identity
```
  -i BLAST_MIN_IDENTITY, --blast-min-identity BLAST_MIN_IDENTITY
                        Minimum identity (percentage) for filtering BLASTN results (default: 97.0)
```

### Skip downstream analysis
```
  -s, --skip-downstream-analysis
                        Skip abandance statics, phylogenetic and bio-diversity analysis (default: False)
```
Turn on this option if you only want to get taxonomy identification results and do not need other analysis.

### Output directory
```
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        directory for output (default: .)
```
Default is putting `MiFishResult` under your current directory. If you specify another directory `/path/dir/`, it will put results into `/path/dir/MiFishResult`

### Number of threads
```
  -t THREADS, --threads THREADS
                        number of threads for BLASTN and usearch (default: 2)
```
Pass to external programs such as [usearch](https://www.drive5.com/usearch/)

## Results
There are six files in the `MiFishResult` directory.
```
QC.zip
read_stat.xlsx
taxonomy.xlsx
tree.zip (if not using -s)
relative_abandance.json  (if not using -s)
diversity.json  (if not using -s but using -d)
```
The first four files are the same as the web version of [MiFish](http://mitofish.aori.u-tokyo.ac.jp/mifish/).

![QC](https://user-images.githubusercontent.com/1146090/213399391-723d32b7-6aac-4104-a13d-f357ccc20ffd.jpg)
![species](https://user-images.githubusercontent.com/1146090/213399450-384bd293-fb84-46d3-875e-0da5717a2754.jpg)
![tree](https://user-images.githubusercontent.com/1146090/213399478-8265f554-02d7-4882-ba2b-75789ddf732d.jpg)

# An example on using other eDNA primers

See [Riaz](https://github.com/billzt/MiFish/blob/main/Riaz.md)