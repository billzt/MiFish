![favicon](https://user-images.githubusercontent.com/1146090/224640993-186457d3-6820-4d38-97cb-ae8ce4f4c114.png)
# MiFish
This is the command line version of MiFish pipeline. It can also be used with any other eDNA meta-barcoding primers

## Announcement
In the next release (probably in June or July), we will replace [FLASH](http://ccb.jhu.edu/software/FLASH/) and [USEARCH](https://www.drive5.com/usearch/) with [VSEARCH](https://github.com/torognes/vsearch) to suit non-MiFish primers' sequencing data better.

We will also release a docker verison soon.

# References
If you use MiFish Pipeline in your projects, please cite:
* Zhu T, Sato Y, Sado T, Miya M, and Iwasaki W. 2023. MitoFish, MitoAnnotator, and MiFish Pipeline: Updates in ten years.
*Mol Biol Evol*, 40:msad035. https://doi.org/10.1093/molbev/msad035
* Sato Y, Miya M, Fukunaga T, Sado T, Iwasaki W. 2018. MitoFish and MiFish Pipeline: A Mitochondrial Genome Database of Fish with an Analysis Pipeline for Environmental DNA Metabarcoding. *Mol Biol Evol* 35:1553-1555.
* Iwasaki W, Fukunaga T, Isagozawa R, Yamada K, Maeda Y, Satoh TP, Sado T,
Mabuchi K, Takeshima H, Miya M, et al. 2013. MitoFish and MitoAnnotator: a mitochondrial genome database of fish with an accurate and automatic annotation pipeline. *Mol Biol Evol* 30:2531-2540.

If you use MiFish Primers in your projects, please cite:
* Miya M, Sato Y, Fukunaga T, Sado T, Poulsen JY, Sato K, Minamoto T, Yamamoto S, Yamanaka H, Araki H, et al. 2015. MiFish, a set of universal PCR primers for metabarcoding environmental DNA from fishes: detection of more than 230 subtropical marine species. *R Soc Open Sci* 2:150088.

# Install
Currently we only support Linux. Please use conda to manage the environment. If you do not have a Linux OS, or you just want to have a quick look, you can try the [Docker version](https://github.com/billzt/MiFish/blob/main/Docker.md)

## External Dependencies
Add these softwares to your system PATH. You can download all the external executable files [here](http://mitofish.aori.u-tokyo.ac.jp/species/detail/download/?filename=download/external_bin.zip)(except for MAFFT), or compile by yourself.
* [fastp](https://github.com/OpenGene/fastp) (v0.23.2)
* [FLASH](http://ccb.jhu.edu/software/FLASH/) (v1.2.7)
* [seqkit](https://github.com/shenwei356/seqkit/) (v2.3.0)
* [usearch](https://www.drive5.com/usearch/) (v11.0.667)
* [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) (v2.9.0)
* [MAFFT](https://mafft.cbrc.jp/alignment/software/) (v7.505)
* [Gblocks](https://home.cc.umanitoba.ca/~psgendb/doc/Castresana/Gblocks_documentation.html) (v0.91b)
* [FastTreeMP](http://www.microbesonline.org/fasttree/) (v2.1.11)

## Install Steps
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

In Ubuntu, the following library is also needed.
```
sudo apt-get install -y libgl1
```

## Test
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
The range of amplicon lengths (**including primers**). Adjust them to satisfy your own primers. You can estimate the range of  from your reference database file.

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

### Threshold of UNOISE3
```
  -u UNOISE_MIN, --unoise-min UNOISE_MIN
                        value for the -minsize option in UNOISE3 (default: 8)
```
Decrease this value would get higher sensitivity but lower accuracy.

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

### Keep temporary files
```
  -k, --keep-tmp-files  Keep temporary files (default: False)
```
Useful for debug. If you encountered problems, turn it on and share me the `Sample-*` directory in the `MiFishResult` directory.

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
The first four files are the same as the web version of [MiFish](http://mitofish.aori.u-tokyo.ac.jp/mifish/). (Screenshots were from DRR126155 against refDB v3.83)

![QC](https://user-images.githubusercontent.com/1146090/219991073-45f3be06-7610-47fa-9ab3-026711a14114.jpg)
![Species](https://user-images.githubusercontent.com/1146090/219991113-e4dd08ad-ecc2-444a-ad55-d97dd22105de.jpg)
![tree](https://user-images.githubusercontent.com/1146090/213399478-8265f554-02d7-4882-ba2b-75789ddf732d.jpg)

# An example on using other eDNA primers

See [Riaz](https://github.com/billzt/MiFish/blob/main/Riaz.md)

# Tips
1. Please make sure that in a FASTQ/FASTA file, names of reads should start with an identitcal word, such as:
```
@DRR231392.1
@DRR231392.2
@DRR231392.3
```
Otherwise `usearch` cannot work properly.
