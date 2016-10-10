#!/bin/bash

echo 'Dodder'
python2.7 bootstrap.py <(grep 'COPIA' ANNOTATE_Dodder.fasta | awk -F '[>|#]' '$6==1 {print $0}') RM_Dodder.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Dodder.fasta | awk -F '[>|#]' '$6==2 {print $0}') RM_Dodder.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Dodder.fasta | awk -F '[>|#]' '$6==3 {print $0}') RM_Dodder.out
python2.7 bootstrap.py <(grep 'COPIA' ANNOTATE_Dodder.fasta | awk -F '[>|#]' '$6==4 {print $0}') RM_Dodder.out
python2.7 bootstrap.py <(grep 'COPIA' ANNOTATE_Dodder.fasta | awk -F '[>|#]' '$6==5 {print $0}') RM_Dodder.out
python2.7 bootstrap.py <(grep 'COPIA' ANNOTATE_Dodder.fasta | awk -F '[>|#]' '$6==6 {print $0}') RM_Dodder.out
python2.7 bootstrap.py <(grep 'COPIA' ANNOTATE_Dodder.fasta | awk -F '[>|#]' '$6==7 {print $0}') RM_Dodder.out
python2.7 bootstrap.py <(grep 'COPIA' ANNOTATE_Dodder.fasta | awk -F '[>|#]' '$6==8 {print $0}') RM_Dodder.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Dodder.fasta | awk -F '[>|#]' '$6==9 {print $0}') RM_Dodder.out
python2.7 bootstrap.py <(grep 'COPIA' ANNOTATE_Dodder.fasta | awk -F '[>|#]' '$6==10 {print $0}') RM_Dodder.out

echo 'Mglory'
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_mglory.fasta | awk -F '[>|#]' '$6==1 {print $0}') RM_mglory.out
python2.7 bootstrap.py <(grep 'COPIA' ANNOTATE_mglory.fasta | awk -F '[>|#]' '$6==2 {print $0}') RM_mglory.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_mglory.fasta | awk -F '[>|#]' '$6==3 {print $0}') RM_mglory.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_mglory.fasta | awk -F '[>|#]' '$6==4 {print $0}') RM_mglory.out
python2.7 bootstrap.py <(grep 'COPIA' ANNOTATE_mglory.fasta | awk -F '[>|#]' '$6==5 {print $0}') RM_mglory.out
python2.7 bootstrap.py <(grep 'ACDS' ANNOTATE_mglory.fasta | awk -F '[>|#]' '$6==6 {print $0}') RM_mglory.out
python2.7 bootstrap.py <(grep 'ACDS' ANNOTATE_mglory.fasta | awk -F '[>|#]' '$6==7 {print $0}') RM_mglory.out
python2.7 bootstrap.py <(grep 'ENSPM' ANNOTATE_mglory.fasta | awk -F '[>|#]' '$6==8 {print $0}') RM_mglory.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_mglory.fasta | awk -F '[>|#]' '$6==9 {print $0}') RM_mglory.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_mglory.fasta | awk -F '[>|#]' '$6==10 {print $0}') RM_mglory.out


echo 'Spot'
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Spot.fasta | awk -F '[>|#]' '$6==1 {print $0}') RM_Spot.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Spot.fasta | awk -F '[>|#]' '$6==2 {print $0}') RM_Spot.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Spot.fasta | awk -F '[>|#]' '$6==3 {print $0}') RM_Spot.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Spot.fasta | awk -F '[>|#]' '$6==4 {print $0}') RM_Spot.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Spot.fasta | awk -F '[>|#]' '$6==5 {print $0}') RM_Spot.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Spot.fasta | awk -F '[>|#]' '$6==6 {print $0}') RM_Spot.out
python2.7 bootstrap.py <(grep 'ACDS' ANNOTATE_Spot.fasta | awk -F '[>|#]' '$6==7 {print $0}') RM_Spot.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Spot.fasta | awk -F '[>|#]' '$6==8 {print $0}') RM_Spot.out
python2.7 bootstrap.py <(grep 'COPIA' ANNOTATE_Spot.fasta | awk -F '[>|#]' '$6==9 {print $0}') RM_Spot.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Spot.fasta | awk -F '[>|#]' '$6==10 {print $0}') RM_Spot.out


echo 'Tobacco' 
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Tobacco.fasta | awk -F '[>|#]' '$6==1 {print $0}') RM_Tobacco.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Tobacco.fasta | awk -F '[>|#]' '$6==2 {print $0}') RM_Tobacco.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Tobacco.fasta | awk -F '[>|#]' '$6==3 {print $0}') RM_Tobacco.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Tobacco.fasta | awk -F '[>|#]' '$6==4 {print $0}') RM_Tobacco.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Tobacco.fasta | awk -F '[>|#]' '$6==5 {print $0}') RM_Tobacco.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Tobacco.fasta | awk -F '[>|#]' '$6==6 {print $0}') RM_Tobacco.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Tobacco.fasta | awk -F '[>|#]' '$6==7 {print $0}') RM_Tobacco.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Tobacco.fasta | awk -F '[>|#]' '$6==8 {print $0}') RM_Tobacco.out
python2.7 bootstrap.py <(grep 'COPIA' ANNOTATE_Tobacco.fasta | awk -F '[>|#]' '$6==9 {print $0}') RM_Tobacco.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Tobacco.fasta | awk -F '[>|#]' '$6==10 {print $0}') RM_Tobacco.out


echo 'Tomato'
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Tomato.fasta | awk -F '[>|#]' '$6==1 {print $0}') RM_Tomato.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Tomato.fasta | awk -F '[>|#]' '$6==2 {print $0}') RM_Tomato.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Tomato.fasta | awk -F '[>|#]' '$6==3 {print $0}') RM_Tomato.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Tomato.fasta | awk -F '[>|#]' '$6==4 {print $0}') RM_Tomato.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Tomato.fasta | awk -F '[>|#]' '$6==5 {print $0}') RM_Tomato.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Tomato.fasta | awk -F '[>|#]' '$6==6 {print $0}') RM_Tomato.out
python2.7 bootstrap.py <(grep 'NA' ANNOTATE_Tomato.fasta | awk -F '[>|#]' '$6==7 {print $0}') RM_Tomato.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Tomato.fasta | awk -F '[>|#]' '$6==8 {print $0}') RM_Tomato.out
python2.7 bootstrap.py <(grep 'GYPSY' ANNOTATE_Tomato.fasta | awk -F '[>|#]' '$6==9 {print $0}') RM_Tomato.out
python2.7 bootstrap.py <(grep 'NA' ANNOTATE_Tomato.fasta | awk -F '[>|#]' '$6==10 {print $0}') RM_Tomato.out






