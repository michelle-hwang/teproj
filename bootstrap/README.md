# Bootstrap Procedure
The Python bootsrap.py were used to conduct bootstrapping statistics to estimate transposable element repetitiveness in the genome after running the transposable element annotation pipeline. 

RepeatMasker must be run on the ANNOTATE.fa databases against original raw sequencing reads (sample of genome) to determine how much of the reads are transposable elements. The amount of bp masked will be used to estimate the % repetitiveness of the read sample. This can be extrapolated to the full genome size. 

## Required
* Python2
* numpy, sys, and argparse Python modules
* RepeatMasker

## Input 
* RepeatMasker outfile
* FASTA File of sequences

## Output
Will print to stdout, e.g.:
```
BOOTSRAP RESULTS
Statistic: 0.349015
95 CI: 0.02349, 0.049801
Total BP Masked: 6627238
```

## Parameters

#### Optional Parameters

Flag | Default | Description
--------- | -------- | --------
n | 10000 | Number of iterations
a | 0.05 | Significance level

#### Run with default parameters to estimate repetitiveness of entire genome 
```
python2.7 bootstrap.py ANNOTATE.fa repeatmasker.out
```

#### Run with default parameters to estimate contribution of COPIA to genome
```
python2.7 bootstrap.py <(grep 'COPIA' ANNOTATE.fa | awk '{print $1}') repeatmasker.out
```

#### Run with default parameters to estimate contribution of MCL family to genome
```
NUM=1
python2.7 bootstrap.py <(grep 'COPIA' ANNOTATE.fa | awk -F '[>|#]' '$6==$NUM {print $0}') repeatmasker.out
```
