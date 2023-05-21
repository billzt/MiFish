# Use the Docker version

## Install Docker

Please refer to [Docker's manual page](https://docs.docker.com/get-docker/).

## Pull the image of MiFish

```bash
docker pull taobioinfo/mifish:latest
```

Old versions (such as v1.0.2) are also [available](https://hub.docker.com/r/taobioinfo/mifish/tags), though not recommended.

## Run

```bash
# Get testing data from GitHub
git clone https://github.com/billzt/MiFish.git
cd MiFish/test

# run
docker run --rm --workdir=/home -v $(pwd):/home docker.io/taobioinfo/mifish seq mifishdbv3.83.fa -d seq2
```

## Limitations

* Running speed is about four times slower than that on native Linux OS.
* Paths of 1) amplicon sequencing data, and 2) RefDB of your metabarcoding primers, must be in **relative path starting from the current working directory**. That is, `./seq` or `./mifishdbv3.83.fa` is OK, but `/home/user/data/seq/` or `/home/user/data/mifishdbv3.83.fa` can not be used.
