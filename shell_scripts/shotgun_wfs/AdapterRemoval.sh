#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
# out and err files are in your home folder. Organize according to your preferences.
#$ -e /home/dave/lmamr47_AR.err
#$ -o /home/dave/lmamr47_AR.out
#$ -q Single.q
unset SGE_ROOT
#echo $JOB_ID
#NSLOTS=2


/share/apps/bin/AdapterRemoval --file1 /home/dave/lmamr47_fastqs/$1\_R1_001.fastq.gz  --file2 /home/dave/lmamr47_fastqs/$1\_R2_001.fastq.gz  --collapse --threads 2 --trimns --maxns 0 --minquality 30 --minlength 30 --trimqualities --gzip --basename /state/partition1/$1

scp -C -l 500000 /state/partition1/$1* 10.1.1.1:/home/dave/lmamr47_analysis_ready

rm /state/partition1/$1*

