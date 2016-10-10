# TE Project
### Authors
Michelle Hwang, Regina Baucom, **University of Michigan**

### Latest Update
July 2016

## Required Software
* Python 2
* MCL 
* BLAST+

## Included Files
* TE_annotate.sh *(Bash script to run entire pipeline.)*
* annotate-v2.py *(Annotation script. Can be run standalone.)*
* cnv_mcl2na.pl *(Perl script to format MCL output.)*
* sum_bp.py *(Optional script in pipeline to sum base pairs.)*

## About
This a pipeline used to annotate tranposable elements using a word pattern matching and majority rules method. Accuracy of this method is dependent on the divergence of the species being annotated and the species source of the provided transposable element database. 

The annotation script will sort sequences into transposable element families, determined by the MCL clustering algorithm. Each family will be annotated according to class (I or II), order, and superfamily. 

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
  * Generates a table of each family annotation per line
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
The annotate.py script will output 

### Helpful Commands

