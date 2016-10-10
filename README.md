# TE Project
### Authors
Michelle Hwang, mchwang@uga.edu
</br> Regina Baucom, rsbaucom@umich.edu
</br>**University of Michigan**
</br>https://sites.lsa.umich.edu/baucom-lab/

### Latest Update
October 2016

## Required Software
* Python 2
* MCL 
* BLAST+

## Included Files
* TE_annotate.sh *(Bash script to run entire pipeline.)*
* annotate-v2.py *(Annotation script. Can be run standalone.)*
* cnv_mcl2na.pl *(Perl script to format MCL output.)*
* sum_bp.py *(Optional script in pipeline to sum base pairs.)*
* get_fasta_seqs.py *(Optional script to grab sequences from a multi-fasta file.)*

## About
This a pipeline used to annotate tranposable elements using a word pattern matching and majority rules method. Accuracy of this method is dependent on the divergence of the species being annotated and the species source of the provided transposable element database. 

The annotation script will sort sequences into transposable element families, determined by the MCL clustering algorithm. Each family will be annotated according to class (I or II), order, and superfamily. 

**** This pipeline assumes host elements are filtered out. ****

#### Annotations currently considered:

Classes | Orders | Superfamilies
------ | ----- | ------ 
1 | LTR, SINE, LINE, RETRO | GYPSY, COPIA, ANGELA, ATHILA, RETROLYC1
2 | TIR, HELITRON | ACDS, ENSPM, CACTA, PIF/HARBRINGER, HAT, MITE, HELITRON

## How to Use
### Run entire pipeline
The bash script file includes a method to run the entire pipeline. You will need to provide the following files:

* File of your sequeneced reads in FASTA format.
* File of your transposable element database in FASTA format.

The pipeline will run through the following steps:

1. BLAST
  * All by All BLAST of your reads.
  * BLAST of reads to NCBI database.
  * BLAST of reads to provided transposable element database.
2. Clustering
  * Clustering of reads via MCL.
  * Formatting of MCL output by "cnv_mcl2na.pl"
3. Annotation of reads by "annoate.py" 
4. Optional output
  * Generates a table of each sequence's annotation in a tab-delimited format
  * "sum_by.py" reports the total base pairs in each family
  
An example command run with default parameters:

```
bash TE_annotate.sh sequenced-reads.fa te-db.fa arabidopsis 1 0 0 0 .2 0.5 5
```

### Run annotate script standalone
The Python annotation script can be run without using the pipeline. One needs, at minimum:

* MCL formatted outfile
* BLAST outfile from custom transposable element database (outfmt=6 + stitle)
* Fasta file of your sequenced reads

The default parameters:

Flag | Default | Description
------------ | ------------ | ------------
-e | 1 | E-value threshold to be considered from BLAST hit
-b | 0 | Bit threshold to be considered from BLAST hit
-p | 0 | Perecent identity threshold to be considered from BLAST hit
-t | 0.2 | Percent of family members annotated to be considered an annotated family
-st | 0.5 | Percent of members in a small family annotated to be considered an annotated family
-ft | 5 | Number of members or below needed to be considered a small family

To run with default parameters:
```
python2.7 annotate.py out.mcl blast.out sequences.fa
```

To run with default parameters and also include BLAST to NCBI results:
```
python2.7 annotate.py out.mcl blast.out sequences.fa -ncbi ncbi-blast.out 
```

For full list of help and parameters:
```
python2.7 annotate.py -h
```

### Output
The annotate.py script will output two files and a report to stdout.
#### FASTA file of all reads and their annotations in the header (ANNOTATE_*.fa)
This file will contain all your sequences with the transposable element annotations in the header. The header will follow the following format separated by "#":

_Sequence ID, Class, Order, Superfamily, MCL family_

#### Text file of family information (FAM_INFO_*.txt)
Summarizes the annotation of each family, one per line, tab-delimited.

#### Stdout report file (RESULTS_*.txt)
This report file will be printed to stdout. To save it, please pipe it into a newfile. It contains a summary of how the script is determining the annotation of each family. 

_Example Case 1: Typical Family_

In this case, we have MCL family 1 with 10 sequences. It passes the annotation threshold, which means that the % of sequences in the family that were annotated passes our designated threshold. Each sequence name is printed out, along with its determined annotation in order of class, order, and superfamily. The last column is the e-value of the BLAST hit that was used to annotate the sequence. If a sequence is determined to be a host element, it is removed. The final line is what the program decides to annotate the family as. 
```
>FAM=1, N=10
Passes annotation threshold.
sequence-1    1    LTR    COPIA    1e-100
sequnece-2    1    LTR    COPIA    0
sequence-3    1    LTR    NA    0.001
sequence-4    1    LTR    COPIA    1e-7
sequence-5    NA    NA    NA    no BLAST hits
sequence-6    HOST    NA    NA    =REMOVED
sequence-7    1    LTR    NA    1e-56
sequence-8    NA    NA    NA    0.053
sequence-9    1    LTR    COPIA    2e-37
sequence-10    1    LTR    COPIA    7e-10
1    LTR    COPIA
```

_Example Case 2: Tie Family_

In this case, we have a tie to annotate this family between two superfamilies, GYPSY and COPIA. Since GYPSY has the higher e-values for the BLAST hit, we choose GYPSY as the annotation.
```
>FAM=22, N=6
Passes annotation threshold.
sequence-1    NA    NA    NA    1e-12
sequence-2    1    LTR    GYPSY    3e-72
sequence-3    1    LTR    GYPSY    8e-68
sequence-4    1    LTR    COPIA    1e-7
sequence-5    1    LTR    COPIA    1e-10
sequence-6    1    LTR    NA    0
1    LTR    GYPSY
```

_Example Case 3: Small Family_

For families that are considered "small" based on our small family threshold parameter, they can follow a different annotation threshold. By default, a small family must have at least 50% of its elements annotated to be considered an annotated family. In this case, only 33% of the sequences are annotated, so this family is considered an NA family. 

```
>FAM=100, N=3
Does not pass annotation threshold.
sequence-1    1    LTR    GYPSY    1e-11
sequence-2    NA    NA    NA    6e-16
sequence-3    NA    NA    NA    5e-10
NA    NA    NA
```

### Helpful Commands

Get original ID of all sequences that are COPIA (or other superfamily): 
```
grep 'COPIA' ANNOTATE.fasta | awk -F '[>|#]' '{print $2}'
```

Get a FASTA file of all sequences that are COPIA (or other superfamily):
```
python2.7 get_fasta_seqs.py ANNOTATE_Dodder.fasta -k <(grep 'COPIA' ANNOTATE_Dodder.fasta | awk -F '[>|#]' '{print $2}')
```

Get a FASTA file of all sequences in a specific MCL family:
```
NUM=6 
python2.7 get_fasta_seqs.py ANNOTATE_Dodder.fasta -k <(awk -F '[>|#]' '$6==$NUM {print $2}' ANNOTATE_Dodder.fasta)
```

Count how many sequences are in a specific MCL familY: 
```
NUM=6
awk -F '[>|#]' '$6==$NUM {print $2}' ANNOTATE_Dodder.fasta | wc -l	
```
