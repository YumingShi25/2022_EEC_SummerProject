# 2022_EEC_SummerProject
Scripts used in the Imperial College MRes Ecology, Evolution and Conservation Programme summer project,
Identifying hub genes by weighted gene co-expression network analysis to better understand molecular functional responses to insecticide exposure in bumblebee (Bombus terrestris)

a summary of bash scripts submitted to HPC for bioinformatic analysis can be found in bash_commands.sh
section includes:
1. download RNA-seq raw data from NCBI database using SRA-toolkit
2. quality assessment of the .fastq files using fastqc
3. align the fastq files to the reference genome iyBomTerr1.2
4. sort .sam files and turn them into .bam format using samtools
5. check the quality of alignments using qualimap
6. count transcript reads for each genes using HTseq

R codes can be found in R_codes.R, including the following sections:
1. normalize htseq count matrix using DESeq2
2. perform WGCNA using WGCNA
3. calculate module-wise relations in WGCNA
4. modules-treatments correlation in WGCNA
5. export the network to cytoscape
6. identify hub genes using dgha
7. gene ontology enrichment analysis using TopGO, including
        GO enrichment analysis for modules
        GO enrichment analysis for hub genes
