# Bootstrap Procedure

These scripts were used to conduct bootstrapping statistics to estimate transposable element repetitiveness in the genome after running the transposable element annotation pipeline. The Python script contains the bootstrapping script, while the shell scripts are example submission scripts.

RepeatMasker must be run on the ANNOTATE.fa databases against original raw sequencing reads (sample of genome) to determine how much of the reads are transposable elements. The amount of bp masked will be used to estimate the % repetitiveness of the read sample. This can be extrapolated to the full genome size. 

## Required
* Python2
* numpy, sys, and argparse Python moduels
* RepeatMasker

## Input 
* RepeatMasker outfile
* FASTA File of sequences

## Parameters
_* means required_

Flag | Default | Description
--------- | -------- | --------


### Run with default parameters to estimate repetitiveness of entire genome 
```
python2.7 bootstrap.py ANNOTATE.fa repeatmasker.out
```

### Run with default parameters to estimate contribution of COPIA to genome
```
python2.7 bootstrap.py <(grep 'COPIA' ANNOTATE.fa | awk '{print $1}') repeatmasker.out
```

### Run with default parameters to estimate contribution of MCL family to genome
```
NUM=1
python2.7 bootstrap.py <(grep 'COPIA' ANNOTATE.fa | awk -F '[>|#]' '$6==$NUM {print $0}') repeatmasker.out
```
