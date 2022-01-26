		# STEP 1 - QUALITY CONTROL AND TRIMMING
# Tools: 
#	- FastQC
#	- Trimmomatic
# Input files:
#	- ./Fastq/*.fastq - raw data files
#	- ./adapters.fa - Tru_Seq adapter sequences
# Output files:
#	- ./Fastq_trimmed/*.trim.fastq.gz - trimmed fastq files
#	- ./Fastq_trimmed/*.trimG.fastq.gz - trimmed 5'G-clipped fastq files
#	- ./Fastq/*.fastqc.html, ./Fastq_trimmed/*.fastqc.html - quality control files

mkdir ./Fastq_trimmed
path_in='./Fastq'
path_out='./Fastq_trimmed'
for i in $path_in/*.fastq.gz; do {
res_file=${i%.fastq*}.trim.fastq.gz
res_file=${res_file/$path_in/}
res_file=$path_out$res_file
# quality control of raw data
fastqc $i -o ./Fastq;
# adapter trimming and deletion of low quality reads
java -jar /usr/share/java/trimmomatic.jar SE -threads 8 $i $res_file ILLUMINACLIP:./adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35;
# single 5'-Guanine clipping
zcat $res_file | awk '{if(NR%4==2) { if(substr($1, 1, 1) == "G") { print substr($1,2); getline; print $1; getline; print substr($1,2)} else {print $1; getline; print $1; getline; print $1}} else { print $1}}' | gzip > ${res_file%.trim.fastq.gz*}.trimG.fastq.gz;
# quality control of trimmed data
fastqc ${res_file%.trim.fastq.gz*}.trimG.fastq.gz -o ./Fastq_trimmed/;
} ; done


		# STEP 2 - ALIGNMENT
# Tools:
#	- STAR
#	- Samtools
# Input files:
# 	- ./gencode.v33.primary_assembly.genome.fa - genome fasta file, source - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh38.primary_assembly.genome.fa.gz
# 	- ./gencode.v33.primary_assembly.annotation.gtf - annotation file, source -  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.primary_assembly.annotation.gtf.gz
#	- ./Fastq_trimmed/*.trimG.fastq.gz - trimmed 5'G-clipped fastq files
# Output files:
#	- ./bam/*.bam - indexed sorted bam files

# build genome index
STAR --runMode genomeGenerate \
	--runThreadN 8 \
	--genomeDir ./GRCh38.primary_assembly.genome.gencode.STARindex/  \
	--genomeFastaFiles ./GRCh38.primary_assembly.genome.fa \
	--sjdbGTFfile ./gencode.v33.primary_assembly.annotation.gtf \
	--sjdbOverhang 75
# align trimmed fastq
mkdir ./bam
for i in ./Fastq_trimmed/*.trimG.fastq.gz; do {
STAR --runMode alignReads \
--runThreadN 8 \
--genomeDir ./GRCh38.primary_assembly.genome.gencode.STARindex/ \
--readFilesIn $i \
--readFilesCommand zcat \
--outFileNamePrefix ${i%.trimG.fastq*} \
--outSAMtype BAM SortedByCoordinate \
--outFilterMismatchNoverReadLmax 0.06 \
--outSAMstrandField intronMotif \
--outFilterMultimapNmax 1000000;
} ; done
# index bam
ls ./bam/*.bam | xargs -ISMPL bash -c "samtools index -@ 8 SMPL"

		# STEP 3 - DELETION OF READS FROM rRNA LOCUSES AND mtDNA
# Tools:
#	- R, R packages (rtracklayer)
#	- bedtools
#	- RSeQC (split_bam.py)
# Input files:
# 	- ./GRCh38_rRNA.bed - RSEQC rRNA bed file, source - http://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/GRCh38_rRNA.bed.gz/download
# 	- ./GCF_000001405.39_GRCh38.p13_genomic.gtf - Refseq GTF file, source - https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109.20200815/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz
# 	- ./GCF_000001405.39_GRCh38.p13_assembly_report.txt - Refseq genome assembly report file, source - https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109.20200815/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt
#	- ./bam/*.bam - indexed sorted bam files
# Output files:
#	- ./RSEQC_refseq_rRNA_mt.sorted.bed - sorted bed12 containing rRNA locuses and MT genome
#	- ./bam/*.ex.bam - filtered bam, reads of interest
#	- ./bam/*.in.bam - rRNA and MT genome reads
#	- ./bam/*.junk.bam - unmapped reads

# generating rRNA/chrM bed12 file from Refseq GTF and rRNA bed12 provided RSEQC
R
library(rtracklayer)
gtf=readGFF('GCF_000001405.39_GRCh38.p13_genomic.gtf')
genome=read.table('GCF_000001405.39_GRCh38.p13_assembly_report.txt',sep='\t',header=F)
rRNA_bed=gtf[gtf$type=='gene'&gtf$gene_biotype=='rRNA',c(1,4,5,9,6,7)]
rseqc_bed=read.table('GRCh38_rRNA.bed',sep='\t',header=F)
# convert chromosome and scaffold names from Refseq to Gencode style
rRNA_bed_chr=merge(rRNA_bed,genome[genome$V2=='assembled-molecule',c(7,10)],by.x='seqid',by.y='V7')
rRNA_bed_chr=rRNA_bed_chr[,c(7,2:6)]
rRNA_bed_scaf=merge(rRNA_bed,genome[genome$V2!='assembled-molecule',c(5,7)],by.x='seqid',by.y='V7')
rRNA_bed_scaf=rRNA_bed_scaf[,c(7,2:6)]
colnames(rRNA_bed_chr)=colnames(rRNA_bed_scaf)
rRNA_bed_new=rbind(rRNA_bed_chr,rRNA_bed_scaf)
# convert chromosome and scaffold names to Gencode style for GRCh38_rRNA.bed
rseqc_bed_chr=merge(rseqc_bed,genome[genome$V2=='assembled-molecule',c(7,10)],by.x='V1',by.y='V10')
rseqc_bed_chr=rseqc_bed_chr[,1:6]
rseqc_bed_scaf=merge(rseqc_bed,genome[genome$V2!='assembled-molecule',c(5,10)],by.x='V1',by.y='V10')
rseqc_bed_scaf=rseqc_bed_scaf[,c(13,2:6)]
colnames(rseqc_bed_chr)=colnames(rseqc_bed_scaf)
rseqc_bed_new=rbind(rseqc_bed_chr,rseqc_bed_scaf)
# generating chrM raws
mt1=cbind(as.character(genome[genome$V10=='chrM',10]),1,genome[genome$V10=='chrM',9],'MT',NA,'+')
mt2=cbind(as.character(genome[genome$V10=='chrM',10]),1,genome[genome$V10=='chrM',9],'MT',NA,'-')
mt=as.data.frame(rbind(mt1,mt2))
# generating common bed for rRNA locuces and mtDNA
colnames(rRNA_bed_new)=colnames(rseqc_bed_new)=colnames(mt)
rseqc_refseq_rrna_mt_bed=rbind(rRNA_bed_new,rseqc_bed_new,mt)
# convert bed to bed12 format
rseqc_refseq_rrna_mt_bed12=cbind(rseqc_refseq_rrna_mt_bed,rseqc_refseq_rrna_mt_bed[,2:3],0,1,(as.numeric(rseqc_refseq_rrna_mt_bed$V3)-as.numeric(rseqc_refseq_rrna_mt_bed$V2)),0)
write.table(rseqc_refseq_rrna_mt_bed12,'RSEQC_refseq_rRNA_mt.bed',quote=F,sep='\t',row.names=F,col.names=F)
q()
n
# sort bed
bedtools sort -i ./RSEQC_refseq_rRNA_mt.bed > ./RSEQC_refseq_rRNA_mt.sorted.bed
rm ./RSEQC_refseq_rRNA_mt.bed 
# split bam into .ex.bam (reads of interest), .in.bam (rRNA and MT reads) and .junk.bam (unmapped reads)
ls ./bam/*.bam | xargs -ISMPL bash -c "split_bam.py -i SMPL -r ./RSEQC_refseq_rRNA_mt.sorted.bed -o SMPL"



		# STEP 4 - GENERATING CAGE TSS CLUSTERS
# Tools: 
#	- samtools
#	- DPI (identify_tss_peaks.sh), source - https://github.com/hkawaji/dpi1
#	- PromoterPipeline_20150516 (level1.py, depends on pysam), source - http://genome.gsc.riken.jp/plessy-20150516/PromoterPipeline_20150516.tar.gz
# Input files:
#	- ./bam/*.ex.bam - filtered bam, reads of interest
# Output files:
#	- ./osc/*.osc - coverage files
#	- ./osc/*.CTSS.gz - coverage files dor DPI input
#	- ./DPI/ - DPI output, including robust CAGE TSS clusters (./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.bed.gz)

mkdir ./osc
mkdir ./DPI
# index *.ex.bam
ls ./bam/*.ex.bam | xargs -ISMPL bash -c "samtools index -@ 8 SMPL"
# generate *.CTSS.txt for DPI input
path_in='./bam'
path_out='./osc'
for i in $path_in/*.ex.bam; do {
res_file=${i%.ex.bam*}.osc
res_file=${res_file/$path_in/}
res_file=$path_out$res_file
python2 ./PromoterPipeline_20150516/level1.py -q 30 -o $res_file $i;
sed '/^#/ d' $res_file | tail -n +2 | awk 'BEGIN {OFS="\t"};{print $2, $3,$4, $1, $6, $5}'  > ${res_file%.osc*}.CTSS.txt;
LC_ALL=C sort -k 1,1 -k2,2n ${res_file%.osc*}.CTSS.txt >  ${res_file%.osc*}.sorted.CTSS.txt;
gzip ${res_file%.osc*}.sorted.CTSS.txt
} ; done
# generate CAGE TSS clusters
./dpi1-master/identify_tss_peaks.sh -g ./GRCh38.primary_assembly.genome.gencode.STARindex/chrNameLength.txt -i './osc/*.CTSS.txt.gz' -o ./DPI



		# STEP 5 - CALCULATE READ COUNTS AND TPM VALUES FOR EACH ROBUST CAGE TSS CLUSTER
# Tools: 
#	- bedtools
#	- R, R packages (rlist,stringr)
#	- PromoterPipeline_20150516 (level1.py, depends on pysam), source - http://genome.gsc.riken.jp/plessy-20150516/PromoterPipeline_20150516.tar.gz
# Input files:
#	- ./osc/*.CTSS.gz - coverage files
#	- ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.bed.gz - robust CAGE TSS clusters
#	- ./bam/*.ex.bam - filtered bam, reads of interest
# Output files:
#	- ./osc/*.counts - tables of read counts by cluster for each sample
#	- ./Expression/total_tags.txt - library size by sample table
#	- ./Expression/counts.txt - summary table of read counts by sample for each cluster
#	- ./Expression/tpm.txt - summary table of normalized read counts (TPM - tags per million) by sample for each cluster


mkdir ./Expression
gzip -d ./osc/*.gz
gzip -d ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.bed.gz
# generate files containing read counts for CAGE TSS clusters 
path_in='./osc'
for i in $path_in/*.txt; do {
res_file=${i%.Aligned.sortedByCoord.out.bam.CTSS.txt *}.bed
res_file=${res_file/$path_in/}
res_file=$path_in$res_file
awk -v OFS='\t' '{$3 = $2; print}' $i > $res_file;
awk -F'\t' '$2 != "0"' $res_file > ${res_file%.bed*}.nonzero.bed
intersectBed -s -wo -a ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.bed -b ${res_file%.bed*}.nonzero.bed > ${res_file%.bed*}.int.bed;
awk -F '\t' '{a[$4] += $17} END{for (i in a) print i, a[i]}' ${res_file%.bed*}.int.bed > ${res_file%.bed*}.counts;
rm ${res_file%.bed*}.int.bed;
rm ${res_file%.bed*}.nonzero.bed;
rm $res_file;
} ; done

# calculate library size (number of reads with MAPQ>=30) for each sample
ls ./bam/*.ex.bam|awk -F"/" '{print $NF}'> ./Expression/bam_names.txt
ls *.ex.bam | xargs -ISMPL bash -c "samtools view -@ 8 -q 30 -c SMPL" > ./Expression/total_tags_temp.txt
paste -d'\t' ./Expression/bam_names.txt ./Expression/total_tags_temp.txt > ./Expression/total_tags.txt
rm ./Expression/total_tags_temp.txt
rm ./Expression/bam_names.txt

# combine counts by samples for each CAGE TSS cluster to single table
R
library(rlist)
library(stringr)
dpi=read.table('./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.bed')
dpi=as.data.frame(dpi$V4)
colnames(dpi)='V1'
files=list.files(path = './osc/', pattern = '.counts',full.names=T)
f=function(x) {merge(dpi,read.table(x,sep=' '),by='V1',all.x=T)}
m=apply(as.data.frame(files),1,f)
l=list.cbind(m)
rownames(l)=l$V1
leven=l[ ,!c(TRUE,FALSE)]
leven[is.na(leven)] <- 0
names=list.files(path = './osc/', pattern = '.counts',full.names=F)
colnames(leven)=apply(as.data.frame(names),1,function(x) {substr(x, 1, 4)})
# calculate TPM by samples for each CAGE TSS in single table
lib_size=read.table('./Expression/total_tags.txt',sep='\t',row.names=1)
rownames(lib_size)=substr(rownames(lib_size), 1, 4)
lib_size=lib_size[colnames(leven),]
tpm=t(apply(leven,1,function(x){x*1000000/lib_size}))
write.table(leven,'./Expression/counts.txt',quote=F,sep='\t',col.names=T,row.names=T)
write.table(tpm,'./Expression/tpm.txt',quote=F,sep='\t',col.names=T,row.names=T)
q()
n



		# STEP 6 - ESTIMATE TSS FOR MUSCLE TISSUE USING RNA-SEQ DATA
# Tools: 
#	- SEASTAR, source - https://github.com/Xinglab/SEASTAR
# Input files:
#	- ./GSE120862/*.merged.bam - 12 merged by condition bam files from RNA-Seq experiment (only strand-specific data) - /https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120862
#	- ./Homo_sapiens.GRCh38.101.gtf - Ensembl annotation file, source - http://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
#	- ./index/ - bowtie index files, source - https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
#	- ./chr_size.txt - genome size file
# Output files:
#	- /seastar/tmp/tsgtf/nrtss.annotation - TSS annotation file, see more about other output in https://github.com/Xinglab/SEASTAR
mkdir ./seastar
SEASTAR.sh \
-A ./GSE120862/D0.bam,./GSE120862/D2.bam,./GSE120862/D3.bam,./GSE120862/DR0.bam,./GSE120862/DR2.bam,./GSE120862/DR3.bam \
-B ./GSE120862/P0.bam,./GSE120862/P2.bam,./GSE120862/P3.bam,./GSE120862/PR0.bam,./GSE120862/PR2.bam,./GSE120862/PR3.bam \
-o ./seastar/ -g ./Homo_sapiens.GRCh38.101.gtf -i ./chr_size.txt -s ./index/ -t U -S s -c 0.1 -p 8 -b U


		# STEP 7 - ANNOTATE CAGE TSS CLUSTERS
# Tools: 
#	- R, R packages (rtracklayer)
#	- bedtools
#	- ./add_utrs_to_gff.py, source - https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/add_utrs_to_gff/add_utrs_to_gff.py
# Input files:
#	- ./Homo_sapiens.GRCh38.101.gtf - Ensembl annotation file, source - http://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
#	- ./GCF_000001405.39_GRCh38.p13_genomic.gff - Refseq annotation GFF, source -  - https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109.20200815/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz
#	- ./GCF_000001405.39_GRCh38.p13_genomic.gtf - Refseq annotation GTF, source - https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109.20200815/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz
#	- ./GCF_000001405.39_GRCh38.p13_assembly_report.txt - Refseq assembly report file, source - https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109.20200815/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt
#	- ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.bed - robust CAGE TSS clusters
#	- ./Expression/tpm.txt - summary table of normalized read counts (TPM - tags per million) by sample for each cluster
#	- ./seastar/tmp/tsgtf/nrtss.annotation - SEASTAR output TSS annotation file
# Output files:
#	- ./Annotation/annotated_CAGE_TSS_clusters.txt - annotated by genes and genomic features CAGE TSS clusters
mkdir ./Annotation
mkdir ./Annotation/temp
Rscript ./Annotate_peaks.R


		# STEP 8 - DIFFERENTIAL GENE EXPRESSION AND DIFFERENTIAL CAGE TSS USAGE
# Tools: 
#	- R, R packages (DESeq,DESeq2,DRIMSeq)
# Input files:
#	- ./Expression/counts.txt - summary table of read counts by sample for each cluster
#	- ./Annotation/annotated_CAGE_TSS_clusters.txt - annotated by genes and genomic features CAGE TSS clusters
#	- ./Expression/design.txt - design table of experiment (samples in rows, conditions in columns)
# Output files:
#	- ./Dif_gene_expression_and_tss_usage/P1-D_genes.txt, - ./Dif_gene_expression_and_tss_usage/P3-D_genes.txt, - ./Dif_gene_expression_and_tss_usage/P6-D_genes.txt - Differential gene expression for each condition
#	- ./Dif_gene_expression_and_tss_usage/Dif_gene_expression.txt - summary of differential gene expression
#	- ./Expression/DESeq2_normalized_counts.txt - DESeq2 normalized counts for CAGE TSS clusters
#	- ./Dif_gene_expression_and_tss_usage/P1-D_dif_tss_usage.txt, ./Dif_gene_expression_and_tss_usage/P1-D_dif_tss_usage.txt, ./Dif_gene_expression_and_tss_usage/P1-D_dif_tss_usage.txt - Differential TSS usage by gene for each condition
#	- ./Dif_gene_expression_and_tss_usage/Dif_tss_usage.txt - summary of differential TSS usage
mkdir ./Dif_gene_expression_and_tss_usage
Rscript ./Dif_gene_expression_and_tss_usage.R


		# STEP 9 - PROMOTER SELECTION, DIFFERENTIAL PROMOTER EXPRESSION AND DIFFERENTIAL USAGE
# Tools: 
#	- R, R packages (DESeq,DESeq2,DRIMSeq)
#	- bedtools
#	- bigBedToBed (kentUtils) - source https://github.com/ENCODE-DCC/kentUtils/blob/master/src/utils/bigBedToBed/
# Input files:
#	- ./Expression/counts.txt - summary table of read counts by sample for each cluster
#	- ./Annotation/annotated_CAGE_TSS_clusters.txt - annotated by genes and genomic features CAGE TSS clusters
#	- ./Expression/design.txt - design table of experiment (samples in rows, conditions in columns)
#	- ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.bed - robust CAGE TSS clusters
#	- ./Expression/DESeq2_normalized_counts.txt - DESeq2 normalized counts for CAGE TSS clusters
#	- ./Promoters/Chip_seq_data/*.bb - bigBed files from GTRD database, Chip-Seq meta-clusters dor each transcription factor, source - http://gtrd20-06.biouml.org/downloads/current/gtrdHub/hg38/bigBed/ChIP-seq/Meta-clusters_by_TF/
#	- ./Promoters/ATAC_seq_data/*.bed.gz - bed files from ENCODE database, ATAC-Seq data for human gastrocnemius medialis tissue (4 persons)
#	- ./Promoters/DNase_seq_data/*.bed.gz - bed files from ENCODE database, DNase-Seq data for human gastrocnemius medialis tissue (4 persons)
#	- ./chrNameLength.txt - genome size file
# Output files:
#	- ./Promoters/OCR.bed - bed with open chromatin regions based on ATAC-Seq and DNase-Seq data
#	- ./Promoters/Promoters.txt - summarized data for all promoters determined in experiment
#	- ./Promoters/Dif_expression/P1-D.txt, - ./Promoters/Dif_expression/P3-D.txt, - ./Promoters/Dif_expression/P6-D.txt - Differential promoter expression for each condition
#	- ./Promoters/Dif_expression/Dif_promoter_usage.txt - summary of differential promoter expression
#	- ./Promoters/Dif_expression/P1-DEP.txt, - ./Promoters/Dif_expression/P3-D_promoter_usage.txt, - ./Promoters/Dif_expression/P6-D_promoter_usage.txt - Differential promoter usage for each condition
#	- ./Promoters/Dif_expression/Dif_promoter_usage.txt - summary of differential promoter usage

mkdir ./Promoters
mkdir ./Promoters/temp
mkdir ./Promoters/Chip_seq_data
mkdir ./Promoters/ATAC_seq_data
mkdir ./Promoters/DNase_seq_data
mkdir ./Promoters/Dif_expression

	# Download Chip-seq Meta-clusters for each TF
wget -r -l1 -t1 -nd -N -np -A.bb -erobots=off --directory-prefix ./Promoters/Chip_seq_data/ http://gtrd20-06.biouml.org/downloads/current/gtrdHub/hg38/bigBed/ChIP-seq/Meta-clusters_by_TF/ 
	# Convert bigBed to bed
ls ./Promoters/Chip_seq_data/*.bb | xargs -ISMPL bash -c "bigBedToBed SMPL SMPL.bed"
	# Combine all meta-clusters to one sorted bed
cat ./Promoters/Chip_seq_data/*.bed > ./Promoters/Chip_seq.bed
sort -k1,1 -k2,2 < ./Promoters/Chip_seq_data/Chip_seq.bed > ./Promoters/Chip_seq_data/Chip_seq.sorted.bed
	# Download ATAC-seq IDR thresholded peaks for each human sample (gastrocnemius medialis tissue) from ENCODE
#ENCSR689SDA
wget --directory-prefix ./Promoters/ATAC_seq_data/ https://www.encodeproject.org/files/ENCFF666UUK/@@download/ENCFF666UUK.bed.gz
#ENCSR308HPZ
wget --directory-prefix ./Promoters/ATAC_seq_data/ https://www.encodeproject.org/files/ENCFF461RJQ/@@download/ENCFF461RJQ.bed.gz
#ENCSR258JCL
wget --directory-prefix ./Promoters/ATAC_seq_data/ https://www.encodeproject.org/files/ENCFF034WCY/@@download/ENCFF034WCY.bed.gz
#ENCSR823ZCR
wget --directory-prefix ./Promoters/ATAC_seq_data/ https://www.encodeproject.org/files/ENCFF953QCT/@@download/ENCFF953QCT.bed.gz
gunzip ./Promoters/ATAC_seq_data/*.gz 
	# Download DNase-seq bed narrow peaks for each human sample (gastrocnemius medialis tissue) from ENCODE
#ENCSR686WJL
wget --directory-prefix ./Promoters/DNase_seq_data/ https://www.encodeproject.org/files/ENCFF780TDG/@@download/ENCFF780TDG.bed.gz
#ENCSR520BAD
wget --directory-prefix ./Promoters/DNase_seq_data/ https://www.encodeproject.org/files/ENCFF033EQB/@@download/ENCFF033EQB.bed.gz
#ENCSR791BHE
wget --directory-prefix ./Promoters/DNase_seq_data/ https://www.encodeproject.org/files/ENCFF884XQO/@@download/ENCFF884XQO.bed.gz
#ENCSR856XLJ
wget --directory-prefix ./Promoters/DNase_seq_data/ https://www.encodeproject.org/files/ENCFF909ISK/@@download/ENCFF909ISK.bed.gz
gunzip ./Promoters/DNase_seq_data/*.gz 
	# Combine ATAC-seq data to one sorted bed and merge open chromatin intervals
cat ./Promoters/ATAC_seq_data/*.bed >> ./Promoters/ATAC_seq_data/ATAC-seq.bed
sort -k 1,1 -k2,2n ./Promoters/ATAC_seq_data/ATAC-seq.bed > ./Promoters/ATAC_seq_data/ATAC-seq.sorted.bed
cat ./Promoters/ATAC_seq_data/ATAC-seq.sorted.bed | awk '{print $1 "\t" $2 "\t" $3}' > ./Promoters/ATAC_seq_data/ATAC_seq.sorted2.bed
bedtools merge -i ./Promoters/ATAC_seq_data/ATAC_seq.sorted2.bed > ./Promoters/ATAC_seq_data/ATAC_seq.sorted2.merged.bed
	# Combine DNase-seq data to one sorted bed and merge open chromatin intervals
cat ./Promoters/DNase_seq_data/*.bed >> ./Promoters/DNase_seq_data/DNase-seq/DNase-seq.bed
sort -k 1,1 -k2,2n ./Promoters/DNase_seq_data/DNase-seq/DNase-seq.bed > ./Promoters/DNase_seq_data/DNase-seq/DNase-seq.sorted.bed
cat ./Promoters/DNase_seq_data/DNase-seq.sorted.bed | awk '{print $1 "\t" $2 "\t" $3}' > ./Promoters/DNase_seq_data/DNase_seq.sorted2.bed
bedtools merge -i ./Promoters/DNase_seq_data/DNase_seq.sorted2.bed > ./Promoters/DNase_seq_data/DNase_seq.sorted2.merged.bed
	# Get intersected regions for ATAC-seq and DNase-seq data
bedtools intersect -a ./Promoters/ATAC_seq_data/ATAC_seq.sorted2.merged.bed -b ./Promoters/DNase_seq_data/DNase_seq.sorted2.merged.bed >./Promoters/OCR.bed

Rscript ./Promoters.R
