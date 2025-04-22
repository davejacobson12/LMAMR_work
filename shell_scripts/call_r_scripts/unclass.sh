#!/bin/sh

awk -F "\t" 'NR == 1 || $1 ~ /g__/ || $1 ~ /unclass/' $1\_humann2.txt > $1\_processed.txt

sed 's/# Pathway/Pathway|Species/' $1\_processed.txt | tr "|" "\t" > $1\_forR.txt

Rscript unclassAbund.R $1\_forR.txt $1 $2 
