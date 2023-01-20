# An example on using other eDNA primers

Here we use the primer [Riaz](https://doi.org/10.1093/nar/gkr732) as an example. 

## Amplicon sequencing data
Amplicon sequencing data were from [Brys et al, 2021](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.15742). We choose pooled samples in timepoint T2 (SRR11479674~SRR11479774, 15 runs) as an example. These 15 accession IDs are:
```
SRR11479674
SRR11479675
SRR11479676
SRR11479677
SRR11479678
SRR11479680
SRR11479681
SRR11479682
SRR11479683
SRR11479769
SRR11479770
SRR11479771
SRR11479772
SRR11479773
SRR11479774
```

Download these data.
```
$ cd test
$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/074/SRR11479674/SRR11479674_1.fastq.gz
$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/074/SRR11479674/SRR11479674_2.fastq.gz
...
...
...
$ mkdir seq_riaz
$ mv SRR*.fastq.gz seq_riaz/
```

## RefDB
In this example, we directly use the local RefDB provided by [Brys et al, 2021](https://zenodo.org/record/3730934). Download it and changed to the MiFish RefDB format. See `test/INBO_riaz.db.fa` as an example.

In other situations, if you do not have a local refDB for your primers, you can follow [CRABS](https://github.com/gjeunen/reference_database_creator) to make a refDB (step 1~6, using MitoFish as original source), then using the awk command to change it to FASTA format.
```
$ awk '{print ">gb|" $1 "|" $9 "\n" $10}' output.tsv >your.db.fa
$ makeblastdb -dbtype nucl -in your.db.fa
```

## Run Pipeline
The primer sequence of Riaz is
```
FWD = 5′-ACTGGGATTAGATACCCC-3′
REV = 5′-TAGAACAGGCTCCTCTAG-3′
```
Try running
```
$ mifish seq_riaz INBO_riaz.db.fa -m 0 -M 99999 -f ACTGGGATTAGATACCCC -r TAGAACAGGCTCCTCTAG -s -t 10
```

It should generate results the same as `test/MiFishResult-example-Riaz`. Species composition are in consistance with the result in Figure 2, [Brys et al, 2021](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.15742)

