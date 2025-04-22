#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
# out and err files are in your home folder. Organize according to your preferences. 
#$ -e /home/dave/MBIO5623/usearch/usearch_out/usearch.err
#$ -o /home/dave/MBIO5623/usearch/usearch_out/usearch.out
#$ -q Single.q
unset SGE_ROOT
#echo $JOB_ID
#NSLOTS=1
### provide these arguments again
# other arguments for AdapterRemoval should be used as we discussed in class. Q score cutoff, length trimming, etc.

#Create OTU tables from shotgun reads?


## SCRIPT BEGINS ##
## Give the path to the folder containing your raw reads files
path=/home/dave/MBIO5623/community_compositions/10million_ez_mapping/

/home/dave/MBIO5623/usearch/usearch10.0.240_i86linux32 -closed_ref $path/$1\_renamed.fasta -db /home/dave/MBIO5623/usearch/eztaxon_qiime_full.fasta --stepwords 20 -wordlength 12 -maxaccepts 0 -maxrejects 0 -strand both -uc /state/partition1/$1.uc

scp -C -l 500000 /state/partition1/$1.* 10.1.1.1:/home/dave/MBIO5623/usearch/usearch_out/
rm /state/partition1/$1.*
