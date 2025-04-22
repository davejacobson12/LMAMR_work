#!/bin/sh
source activate qiime1

#analyze ecological resilience metrics using a suite of rscripts

cat /Users/dave/Desktop/reference_based_mapping/ddig.txt | while read pathway name; do /Users/dave/Desktop/reference_based_mapping/run_pathRich.sh $pathway $name; done
cat /Users/dave/Desktop/reference_based_mapping/ddig.txt | while read pathway name; do /Users/dave/Desktop/reference_based_mapping/faithPD.sh $pathway $name; done
cat /Users/dave/Desktop/reference_based_mapping/ddig.txt | while read pathway name; do /Users/dave/Desktop/reference_based_mapping/unclass.sh $pathway $name; done


Rscript /Users/dave/Desktop/reference_based_mapping/ecol_plot.r 

rm /Users/dave/Desktop/reference_based_mapping/delete.txt
rm /Users/dave/Desktop/reference_based_mapping/tempRich.txt
rm /Users/dave/Desktop/reference_based_mapping/*R.tre
rm /Users/dave/Desktop/reference_based_mapping/*align.fasta
rm /Users/dave/Desktop/reference_based_mapping/*taxa.fasta

source deactivate qiime1

