#This lays out the code to run the analysis - it would need some adaptation to run on HPC


### I found the required SRA names for different datasets. For Matses I only downloaded the samples that were reported as having treponema, from the baboon datasets I downloaded the four datasets with the highest reported amount of Treponema succinofaciens, and for Hadza I downloaded metagenome reads from 4 women as women were reported as having higher treponema abundance so I picked four women
# Used this general format to download the respective SRA datasets
../programs/sratoolkit.2.8.2-1-mac64/bin/fastq-dump --split-3 sra_file

### Downloaded treponema genomes from NCBI and concatenated into one file (treponema.fasta).Then I realized that I downloaded non-complete genomes, so I had to figure out which genomes were complete circular and concatenate those into one file.

grep "complete\ genome" treponema.fasta  | awk '{print$3}' | while read line; do cat t_$line.fna >> trep_complete.fasta; done

##build a bowtie database with the trep_complete, so now the treponema_ref base name will be used for the reference mapping

../../shotgun_analysis/new/ReferenceBasedMapping/bowtie_dir/bowtie2-build trep_complete.fasta treponema_ref

### for each of the fastq files retrieved from the datasets, I did a quality check

ls *.fastq | awk -F”_” ‘{print$1}’| sort | uniq | while read line; do fastqc $line\_1.fastq; done

ls *.fastq | awk -F”_” ‘{print$1}’| sort | uniq | while read line; do fastqc $line\_2.fastq; done

### the quality looks good for each of my files, all are q>30

### Next is quality filtering and read merging because each of my file samples are from paired end sequencing - I’ll use adapter removal

ls *.fastq | awk -F"_" '{print$1}'| sort | uniq | while read line; do AdapterRemoval --file1 $line\_1.fastq --file2 $line\_2.fastq --maxns 0 --trimns --minquality 30 --trimqualities --collapse --basename $line.AR --minlength 30 1>$line.out 2>$line.err ; done


# map to treponema reference 

ls *.err| sed 's/.err//' | while read line; do ../../shotgun_analysis/new/ReferenceBasedMapping/bowtie_dir/bowtie2 -x ../treponema_genomes/treponema_ref -1 $line.AR.pair1.truncated -2 $line.AR.pair2.truncated -U $line.AR.collapsed -S $line.sam --local --no-unal -p 2 1>$line.sam.out 2>$line.sam.err; done

### make bam files

ls *.sam | awk -F"." '{print$1}' | while read line; do ../../shotgun_analysis/new/ReferenceBasedMapping/samtools_dir/samtools view -S -b -o $line.bam $line.sam; done

### sort bam files

ls *.sam | awk -F"." '{print$1}' | while read line; do ../../shotgun_analysis/new/ReferenceBasedMapping/samtools_dir/samtools sort $line.bam -o sorted_$line.bam; done

### To get the proportion of reads that are for each reference from each SRR sample

ls *.sam | sed 's/.sam//' | while read line; do grep -v "^@" $line.sam | awk '{print$1"\t"$3}' | sort | uniq | awk '{print$2}' | sort | uniq -c > $line\_reads.txt ; done

### to get the data into a tabular format that I could easily edit in excel

 ls *.sam | sed 's/.sam//' | while read line; do echo $line | tr "\n" "\t"; awk '{print$1}' $line\_reads.txt | tr "\n" "\t" ; echo; done > hadza_proportion.txt

#### Because I have mostly treponema succinofaciens reads, I will build VCF files with succinofacines as the reference 

### first make the sussinifaciens ref

../../shotgun_analysis/new/ReferenceBasedMapping/bowtie_dir/bowtie2-build t_succinifaciens.fna succinifaciens_ref

### create sam files by mapping to succinifaciens 

ls *.AR.collapsed | sed 's/.AR.collapsed//' |while read line; do ../../shotgun_analysis/new/ReferenceBasedMapping/bowtie_dir/bowtie2 -x succinifaciens_ref -1 $line.AR.pair1.truncated -2 $line.AR.pair2.truncated -U $line.AR.collapsed -S $line.sam --local --no-unal -p 2 1>$line.sam.out 2>$line.sam.err; done

### create bam files 

ls *.sam | awk -F"." '{print$1}' | while read line; do ../../shotgun_analysis/new/ReferenceBasedMapping/samtools_dir/samtools view -S -b -o $line.bam $line.sam; done

## sorted bam files

ls *.sam | awk -F"." '{print$1}' | while read line; do ../../shotgun_analysis/new/ReferenceBasedMapping/samtools_dir/samtools sort $line.bam -o sorted_$line.bam; done

### created t succinifaciens reference for vcf file building

../../shotgun_analysis/new/ReferenceBasedMapping/samtools_dir/samtools faidx t_succinifaciens.fna

### created vcf files

ls *.sam | sed 's/.sam//' | while read line; do ../../shotgun_analysis/new/ReferenceBasedMapping/samtools_dir/samtools mpileup -uv -A -d 10000 --skip-indels -B -q 20 -Q 30 -f t_succinifaciens.fna sorted_$line.bam > $line.vcf ; done

### used the KEGG database to determine the start and end sequence of the desired gene (gyrB etc)

### use the multiFasta.py from class to convert the succinifaciens genome file to a single line fasta.

multiFasta.py ../final_project/succinofaciens_mapping/t_succinifaciens.fna  > succinifaciens_single.fna 

### Because this file has both the genome and the plasmid, i used grep to pull out the fasta line that corresponds to the genome

grep "NC_015385.1" -A 1  succinifaciens_single.fna > succinifaciens_genome.fna

### then ran the vcf to fasta script

### this is the script (vcf2fasta.py) that takes the vcf output and then converts it to a fasta, using the reference genome and depth of 5

ls *.vcf | sed 's/.vcf//' | while read line; do python3 ../vcf_test.py $line.vcf succinifaciens_genome.fna $line\_test.fna; done

### wrote a python script to pull out gene sequence from each mapped fasta file for the given marker genes

### for rpoB
ls *.vcf | sed 's/.vcf//' | while read line; do python3 ../read_reference.py $line\_mapped.fna 2617208 2621692 $line\_rpoB.fna ; done

### for gyrB
ls *.vcf | sed 's/.vcf//' | while read line; do python3 ../read_reference.py $line\_mapped.fna 2729846 2731771 $line\_gyrB.fna ; done

### concat into all one file

cat *_gyrB.fna >> gyrB_all.fna

### use muscle to align

muscle -in gyrB_all.fna -physout gyrB.align.phy

### use phyML to build tree

### open gyrB.align.phy_phyml_tree.txt in mega display newick tree

### used borrelia recurrent as an out group, so I downloaded the filed and gunzipped it

### Then pulled out just the genome

/Users/dave/Mbio5810/python_work/multiFasta.py b_recurrentis_comlete.fna  > b_recur_single.fna

grep "^>NC_011244.1" -A 1 b_recur_single.fna > b_recur_full_genome.fna

## got the sequences for gyrB and rpoB for the trep succinifaciens reference genome and borrelia

python3 ../read_reference.py succinifaciens_genome.fna 2617208 2621692 succ_rpoB.fna 

python3 ../read_reference.py succinifaciens_genome.fna 2729846 2731771 succ_gyrB.fna 

python3 ../read_reference.py b_recur_full_genome.fna 472832 474751  borrelia_gyrB.fna

python3 ../read_reference.py b_recur_full_genome.fna 422058 425525  borrelia_rpoB.fna

#### Instead of doing muscle and alignment with physout, do with just -out. This will give multiline Fasta -convert to single line Fasta.

muscle -in gyrB_all.fna -out gyrB_align.fna

muscle -in rpoB_all.fna -out rpoB_align.fna

python3 /Users/dave/Mbio5810/python_work/multiFasta.py rpoB_align.fna > rpoB_align_single.fasta

grep "^>" gyrB_align_single.fasta | while read line; do echo "$line"; grep -h -A 1 -w "$line" *single.fasta | grep -v "^>" | grep -v "^--$"| tr -d "\n" ; echo ; done > markers_concat.fna

### added four more markers from that should refine the phylogenetic relationships according to the literature. Found in both trep succinifacines and borrelia

# ribosomal protein S2
# ribosomal protein L1
# ribosomal protein L11
# ribosomal protein L5

ls *.vcf | sed 's/.vcf//' | while read line; do python3 ../read_reference.py $line\_mapped.fna 1294644 1295510 $line\_rpsB.fna ; done
ls *.vcf | sed 's/.vcf//' | while read line; do python3 ../read_reference.py $line\_mapped.fna 2626623 2627300 $line\_rplA.fna ; done
ls *.vcf | sed 's/.vcf//' | while read line; do python3 ../read_reference.py $line\_mapped.fna 2627300 2627734 $line\_rplK.fna ; done
ls *.vcf | sed 's/.vcf//' | while read line; do python3 ../read_reference.py $line\_mapped.fna 2607798 2608349 $line\_rplE.fna ; done


python3 ../read_reference.py succinifaciens_genome.fna 1294644 1295510  succ_rpsB.fna
python3 ../read_reference.py succinifaciens_genome.fna 2626623 2627300  succ_rplA.fna
python3 ../read_reference.py succinifaciens_genome.fna 2627300 2627734  succ_rplK.fna
python3 ../read_reference.py succinifaciens_genome.fna 2607798 2608349  succ_rplE.fna

grep "^>" rplA_align.single.fasta | while read line; do echo "$line"; grep -h -A 1 -w "$line" *_align.single.fasta | grep -v "^>" | grep -v "^--$" | tr -d "\n" ; echo ; done > final_markers_concat.fna


muscle -in final_markers_concat.fna -physout final_markers.aligned.phy




















