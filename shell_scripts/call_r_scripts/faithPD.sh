 #!/bin/sh



sed 's/g__.*s__//g' $1com_matrix.txt | awk -F "\t" '{print$1}' | while read line; do grep "$line" -m 1 eztaxon_id_taxonomy.txt | awk -v a="$line" '{print$1"\t"a}' ; done | while read number taxa ; do grep "$number" -A 1 eztaxon_qiime_full.fasta | sed "s/${number}/${taxa}/g" ;done  > $1\_taxa.fasta

sed 's/g__.*s__//' $1com_matrix.txt > $1edited_comm.txt

mafft $1\_taxa.fasta > $1\_align.fasta

make_phylogeny.py -i $1\_align.fasta -o $1\_R.tre

Rscript faithPD.R $1\_R.tre $1 $1edited_comm.txt $2


