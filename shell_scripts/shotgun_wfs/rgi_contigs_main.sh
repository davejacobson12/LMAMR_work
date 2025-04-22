#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
# out and err files are in your home folder. Organize according to your preferences. 
#$ -e /home/dave/rgi-5.1.0/rgi_err
#$ -o /home/dave/rgi-5.1.0/rgi_out
#$ -q Single.q
unset SGE_ROOT
#echo $JOB_ID
#NSLOTS=4

home/dave/rgi-5.1.0/rgi main -i /home/dave/rgi-5.1.0/$1\_contigs.fasta.gz -o /state/partition1/$1\_rgi -n 8 --include_loose --local --clean --low_quality 

scp -C -l 500000 /state/partition1/$1\_rgi/* 10.1.1.1:/home/dave/rgi-5.1.0/$1\_rgi/

rm -r /state/partition1/$1\_rgi/
