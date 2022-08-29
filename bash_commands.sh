## 1.download RNA-seq raw data from NCBI database using SRA-toolkit
module load sra-toolkit/2.8.1
cd /rds/general/user/ys2121/home/Bumblebee/sra/fastq
# batch download according to a list
prefetch --option-file SRP172774_run_acc_list.txt
# unzip the raw .sra files into .fastq files
fastq-dump *.sra --outdir /rds/general/user/ys2121/home/Bumblebee/sra/fastq


## 2.quality assessment of the .fastq files using fastqc
module load anaconda3/personal
module load fastqc
for file in /rds/general/user/ys2121/home/bumblebee/fastq/*.fastq
do
fastqc -f fastq ${file} -t 20 -o /rds/general/user/ys2121/home/bumblebee/multiqc -d /rds/general/user/ys2121/home/bumblebee/multiqc
done


## 3.align the fastq files to the reference genome iyBomTerr1.2
module load gcc/11.2.0
module load hisat2/2.0.4
cd /rds/general/user/ys2121/home/bumblebee/hisat2/iyBomTerr1.2

# extract exon and splice sites from the reference genome using HISAT2 provided python scripts
module load anaconda3/personal
cd /rds/general/user/ys2121/home/bumblebee/reference/iyBomTerr1.2
python hisat2_extract_splice_sites.py iyBomTerr1.2_genomic.gtf > iyBomTerr1.2.splice_sites.txt
python hisat2_extract_exons.py iyBomTerr1.2_genomic.gtf > iyBomTerr1.2.exons.txt
#build a exon and splice site reference 
hisat2-build -p 48 --ss iyBomTerr1.2.splice_sites.txt --exon iyBomTerr1.2.exons.txt\
 iyBomTerr1.2_genomic.fna iyBomTerr1.2 > iyBomTerr1.2.log 2> iyBomTerr1.2.err

# perform splice-aware alignment to all the .fastq files
treatments="CON CLO IMI"
castes="Q W"
# need to manually change the file names to "accension_colony_treatment_4_caste.fastq"
# example: SRR8288022_C02_CLO_4_Q.fastq
for treatment in $treatments; do
  for caste in $castes; do  
    echo "running with ${treatment}"
    echo "running with ${caste}"
    fastqs=`ls /rds/general/user/ys2121/home/bumblebee/fastq/*_${treatment}_4_${caste}.fastq` #link the two sequencing files for one sample
    colonies=$(ls $fastqs | cut -d '_' -f 2 | sort -u) #get colony information from file names
  
    for colony in $colonies; do
      echo "running with ${colony}"
      fastqs_colony=`ls /rds/general/user/ys2121/home/bumblebee/fastq/*_${treatment}_4_${caste}.fastq | grep "${colony}" | tr '\n' ','`
      echo "$fastqs_colony"
      output=/rds/general/user/ys2121/home/bumblebee/hisat2/iyBomTerr1.2/${treatment}_${colony}_${caste}
      echo "hisat2 -x/rds/general/user/ys2121/home/bumblebee/reference/iyBomTerr1.2/iyBomTerr1.2 \
       -U ${fastqs_colony} --known-splicesite-infile /rds/general/user/ys2121/home/bumblebee/reference/iyBomTerr1.2.splice_sites.txt\
       --threads=48 -S ${output}.sam" >> hisat2_iyBomTerr1.2.sh
    done 
  done
done
sh hisat2_iyBomTerr1.2.sh > hisat2.iyBomTerr1.2.sh.log 2> hisat2.iyBomTerr1.2.sh.err

# rename aligned files to .sam 
cd /rds/general/user/ys2121/home/bumblebee/sam/iyBomTerr1.2
find . -depth -name "*.fastq.sam" -exec sh -c 'f="{}"; mv -- "$f" "${f%.fastq.sam}.sam"' \;


## 4.sort .sam files and turn them into .bam format using samtools
module load samtools/1.3.1
cd /rds/general/user/ys2121/home/bumblebee/bam/iyBomTerr1.2
for file in /rds/general/user/ys2121/home/bumblebee/sam/iyBomTerr1.2/*.sam; do
  name=$(echo $(basename $file .sam))
  samtools view -bS "$file" | samtools sort --threads 48 -o $name.sorted.bam
done


## 5.check the quality of alignments using qualimap
## output = qualimapReport.html
module load anaconda3/personal
module load qualimap/2.2.1
cd /rds/general/user/ys2121/home/bumblebee/multiqc
#batch quality check
qualimap multi-bamqc -r -d quali_path.txt -outdir /rds/general/user/ys2121/home/bumblebee/bam/iyBomTerr1.2/qualimap
multiqc /rds/general/user/ys2121/home/bumblebee/bam/iyBomTerr1.2/qualimap --filename iyBomTerr1.2_bam_qc_report.html


## 6.count transcript reads for each genes using HTseq 
module load anaconda3/personal
source activate HTseq
cd /rds/general/user/ys2121/home/bumblebee/htseq/iyBomTerr1.2
htseq-count -f bam -r pos -s reverse -i gene_id \
 /rds/general/user/ys2121/home/bumblebee/bam/iyBomTerr1.2/*.bam \
 /rds/general/user/ys2121/home/bumblebee/reference/iyBomTerr1.2/iyBomTerr1.2_genomic.gtf \
 1>/rds/general/user/ys2121/home/bumblebee/htseq/iyBomTerr1.2/htseq_count_iyBomTerr1.2.txt
 2>/rds/general/user/ys2121/home/bumblebee/htseq/iyBomTerr1.2/htseq_summary_iyBomTerr1.2.summary

