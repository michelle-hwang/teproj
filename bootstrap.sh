#!/bin/bash

##### UNTESTED

a='ANNOTATE_Dodder.fasta'
b='ANNOTATE_mglory.fasta'
c='ANNOTATE_Spot.fasta'
d='ANNOTATE_Tobacco.fasta'
e='ANNOTATE_Tomato.fasta'

a2='RM_Dodder.out'
b2='RM_mglory.out'
c2='RM_Spot.out'
d2='RM_Tobacco.out'
e2='RM_Tomato.out'

echo 'GENOME REPETITIVENESS - ALL'
python2.7 bootstrap.py $a $a2
python2.7 bootstrap.py $b $b2
python2.7 bootstrap.py $c $c2
python2.7 bootstrap.py $d $d2
python2.7 bootstrap.py $e $e2


echo `\nGENOME REPETITIVENESS - CLASS I`
python2.7 bootstrap.py <(python2.7 get_fasta_seqs.py $a -k <(awk -F '[>|#]' '$3==1 {print $2}' $a)) $a2
python2.7 bootstrap.py <(python2.7 get_fasta_seqs.py $b -k <(awk -F '[>|#]' '$3==1 {print $2}' $b)) $b2
python2.7 bootstrap.py <(python2.7 get_fasta_seqs.py $c -k <(awk -F '[>|#]' '$3==1 {print $2}' $c)) $c2
python2.7 bootstrap.py <(python2.7 get_fasta_seqs.py $d -k <(awk -F '[>|#]' '$3==1 {print $2}' $d)) $d2
python2.7 bootstrap.py <(python2.7 get_fasta_seqs.py $e -k <(awk -F '[>|#]' '$3==1 {print $2}' $e)) $e2

#python2.7 bootstrap.py <(python2.7 get_fasta_seqs.py ANNOTATE_Dodder.fasta -k <(awk -F '[>|#]' '$3==1 {print $2}' ANNOTATE_Dodder.fasta)) RM_Dodder.out

echo `\nGENOME REPETITIVENESS - CLASS II`
python2.7 bootstrap.py <(python2.7 get_fasta_seqs.py $a -k <(awk -F '[>|#]' '$3==2 {print $2}' $a)) $a2
python2.7 bootstrap.py <(python2.7 get_fasta_seqs.py $b -k <(awk -F '[>|#]' '$3==2 {print $2}' $b)) $b2
python2.7 bootstrap.py <(python2.7 get_fasta_seqs.py $c -k <(awk -F '[>|#]' '$3==2 {print $2}' $c)) $c2
python2.7 bootstrap.py <(python2.7 get_fasta_seqs.py $d -k <(awk -F '[>|#]' '$3==2 {print $2}' $d)) $d2
python2.7 bootstrap.py <(python2.7 get_fasta_seqs.py $e -k <(awk -F '[>|#]' '$3==2 {print $2}' $e)) $e2


echo `\nGENOME REPETITIVENESS - COPIA`
python2.7 bootstrap.py <(grep 'COPIA' $a | awk '{print $1}') $a2 
python2.7 bootstrap.py <(grep 'COPIA' $b | awk '{print $1}') $b2 
python2.7 bootstrap.py <(grep 'COPIA' $c | awk '{print $1}') $c2 
python2.7 bootstrap.py <(grep 'COPIA' $d | awk '{print $1}') $d2 
python2.7 bootstrap.py <(grep 'COPIA' $e | awk '{print $1}') $e2 


echo `\nGENOME REPETITIVENESS - GYPSY`
python2.7 bootstrap.py <(grep 'GYPSY' $a | awk '{print $1}') $a2 
python2.7 bootstrap.py <(grep 'GYPSY' $b | awk '{print $1}') $b2 
python2.7 bootstrap.py <(grep 'GYPSY' $c | awk '{print $1}') $c2 
python2.7 bootstrap.py <(grep 'GYPSY' $d | awk '{print $1}') $d2 
python2.7 bootstrap.py <(grep 'GYPSY' $e | awk '{print $1}') $e2 


echo `\nGENOME REPETITIVENESS - ENSPM`
python2.7 bootstrap.py <(grep 'ENSPM' $a | awk '{print $1}') $a2 
python2.7 bootstrap.py <(grep 'ENSPM' $b | awk '{print $1}') $b2 
python2.7 bootstrap.py <(grep 'ENSPM' $c | awk '{print $1}') $c2 
python2.7 bootstrap.py <(grep 'ENSPM' $d | awk '{print $1}') $d2 
python2.7 bootstrap.py <(grep 'ENSPM' $e | awk '{print $1}') $e2 


echo `\nGENOME REPETITIVENESS - ACDS`
python2.7 bootstrap.py <(grep 'ACDS' $a | awk '{print $1}') $a2 
python2.7 bootstrap.py <(grep 'ACDS' $b | awk '{print $1}') $b2 
python2.7 bootstrap.py <(grep 'ACDS' $c | awk '{print $1}') $c2 
python2.7 bootstrap.py <(grep 'ACDS' $d | awk '{print $1}') $d2 
python2.7 bootstrap.py <(grep 'ACDS' $e | awk '{print $1}') $e2 


echo `\nGENOME REPETITIVENESS - NA FAMILIES`
python2.7 bootstrap.py <(python2.7 get_fasta_seqs.py $a -k <(awk -F '[>|#]' '$3=="NA" {print $2}' $a)) $a2
python2.7 bootstrap.py <(python2.7 get_fasta_seqs.py $b -k <(awk -F '[>|#]' '$3=="NA" {print $2}' $b)) $b2
python2.7 bootstrap.py <(python2.7 get_fasta_seqs.py $c -k <(awk -F '[>|#]' '$3=="NA" {print $2}' $c)) $c2
python2.7 bootstrap.py <(python2.7 get_fasta_seqs.py $d -k <(awk -F '[>|#]' '$3=="NA" {print $2}' $d)) $d2
python2.7 bootstrap.py <(python2.7 get_fasta_seqs.py $e -k <(awk -F '[>|#]' '$3=="NA" {print $2}' $e)) $e2

