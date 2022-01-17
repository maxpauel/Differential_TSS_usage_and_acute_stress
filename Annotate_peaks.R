		# Add UTRs to Refseq GFF annotation and generate new GFF
system('python3 ./add_utrs_to_gff.py ./GCF_000001405.39_GRCh38.p13_genomic.gff > ./GCF_000001405.39_GRCh38.p13_genomic_with_utr.gff')

		# Edit Refseq GTF annotation to change chromosome names to Ensembl nomenclature
refseqGTF=read.table('./GCF_000001405.39_GRCh38.p13_genomic.gtf',sep='\t')
refseqGenome=read.table('./GCF_000001405.39_GRCh38.p13_assembly_report.txt',sep='\t')
refseqGenome=refseqGenome[1:24,]
refseqGTF_changed=merge(refseqGTF,refseqGenome[,c(3,7)],by.x='V1',by.y='V7')
refseqGTF_changed=refseqGTF_changed[,c(10,2:9)]
write.table(refseqGTF_changed,'./GCF_000001405.39_GRCh38.p13_genomic_red.gtf',quote=F,col.names=F,row.names=F,sep='\t')
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects
gc()

		# Prepare Refseq UTRs table for annotation
library('rtracklayer')
ref=readGFF('./GCF_000001405.39_GRCh38.p13_genomic_with_utr.gff')
ref2=readGFF('./GCF_000001405.39_GRCh38.p13_genomic_red.gtf')
refseqGenome=read.table('./GCF_000001405.39_GRCh38.p13_assembly_report.txt',sep='\t')
refseqGenome=refseqGenome[1:24,]
ref3=unique(ref[ref$type=='five_prime_UTR'|ref$type=='three_prime_UTR',c(1:8,22)])
ref4=merge(ref3,refseqGenome[,c(3,7)],by.x='seqid',by.y='V7')
ref4=ref4[,c(10,2:9)]
ref5=unique(ref2[,c(9:10)])
ref6=merge(ref4,ref5,by='transcript_id',all.x=T)
write.table(ref6,'./Annotation/temp/refseq_utrs.txt',quote=F,sep='\t',row.names=F)
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects
gc()

		# Load annotation files and CAGE TSS clusters .bed
library('rtracklayer')
TC=read.table('./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.chr.bed')
en=readGFF('./Homo_sapiens.GRCh38.101.gtf')
ref=readGFF('./GCF_000001405.39_GRCh38.p13_genomic_red.gtf')
utr=read.table('./Annotation/temp/refseq_utrs.txt',header=T,sep='\t')

		# Prepare .bed with annotated by Ensembl TSS
ent=unique(en[(en$type=='exon'& en$exon_number=='1')|en$type=='gene',])
entp=unique(ent[ent$strand=='+',c(1,4,4,9,6,7)])
entm=unique(ent[ent$strand=='-',c(1,5,5,9,6,7)])
colnames(entp)=colnames(entm)
ensTSS=unique(rbind(entp,entm))
write.table(ensTSS,'./Annotation/temp/ensTSS.bed',sep='\t',row.names=F,col.names=F,quote=F)

		# Prepare .bed with annotated by Refseq TSS
reft=unique(ref[(ref$type=='exon'& ref$exon_number=='1')|ref$type=='gene',])
reftp=unique(reft[reft$strand=='+',c(1,4,4,9,6,7)])
reftm=unique(reft[reft$strand=='-',c(1,5,5,9,6,7)])
colnames(reftp)=colnames(reftm)
refsTSS=unique(rbind(reftp,reftm))
write.table(refsTSS,'./Annotation/temp/refsTSS.bed',sep='\t',row.names=F,col.names=F,quote=F)

		# Find closest upstream and downstream distances between CAGE TSS and annotated TSS
system('bedtools sort -i ./Annotation/temp/ensTSS.bed >./Annotation/temp/ensTSS.sorted.bed;
bedtools sort -i ./Annotation/temp/refsTSS.bed >./Annotation/temp/refsTSS.sorted.bed;
bedtools closest -s -iu -D a -a ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.chr.bed -b ./Annotation/temp/ensTSS.sorted.bed > ./Annotation/temp/ensTSS_up.bed;
bedtools closest -s -id -D a -a ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.chr.bed -b ./Annotation/temp/ensTSS.sorted.bed > ./Annotation/temp/ensTSS_down.bed;
bedtools closest -s -iu -D a -a ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.chr.bed -b ./Annotation/temp/refsTSS.sorted.bed > ./Annotation/temp/refsTSS_up.bed;
bedtools closest -s -id -D a -a ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.chr.bed -b ./Annotation/temp/refsTSS.sorted.bed > ./Annotation/temp/refsTSS_down.bed')

		# Create temporary work tables with distances
tssenu=read.table('./Annotation/temp/ensTSS_up.bed')
tssenu=tssenu[tssenu$V16!='.',]
tssend=read.table('./Annotation/temp/ensTSS_down.bed')
tssend=tssend[tssend$V16!='.',]
tssreu=read.table('./Annotation/temp/refsTSS_up.bed')
tssreu=tssreu[tssreu$V16!='.',]
tssred=read.table('./Annotation/temp/refsTSS_down.bed')
tssred=tssred[tssred$V16!='.',]
tssen=merge(tssenu,tssend,by='V4',all=T)
tssre=merge(tssreu,tssred,by='V4',all=T)
tssen=unique(tssen[,c(1,16,19,34,37)])
tssre=unique(tssre[,c(1,16,19,34,37)])
tss_dist=merge(tssen, tssre, by='V4')
colnames(tss_dist)=c('name','up_ens_g','up_ens_d','down_ens_g','down_ens_d','up_ref_g','up_ref_d','down_ref_g','down_ref_d')
write.table(tss_dist,'./Annotation/temp/tss_distance.txt',sep='\t',row.names=F,col.names=F,quote=F)

		# Determine closest annotated Ensembl TSS between upstream and downstream 
tssen$V19.y=tssen$V19.y*(-1)
tssen$eq=tssen$V19.x<tssen$V19.y
tssen_up=unique(tssen[tssen$eq=='TRUE',c(1,3)])
tssen_up$V19.x=tssen_up$V19.x*(-1)
tssen_d=unique(tssen[tssen$eq=='FALSE',c(1,5)])
colnames(tssen_d)=colnames(tssen_up)

		# Determine closest annotated Refseq TSS between upstream and downstream 
tssre$V19.y=tssre$V19.y*(-1)
tssre$eq=tssre$V19.x<tssre$V19.y
tssre_up=unique(tssre[tssre$eq=='TRUE',c(1,3)])
tssre_up$V19.x=tssre_up$V19.x*(-1)
tssre_d=unique(tssre[tssre$eq=='FALSE',c(1,5)])
colnames(tssre_d)=colnames(tssre_up)

		# Annotate CAGE TSS clusters by Ensembl TSS (distance < 51 bp)
ens_tss=unique(tssen[tssen$V19.x < 51|tssen$V19.y < 51,])
ens_tss1=unique(ens_tss[ens_tss$eq=='FALSE',c(1,4)])
ens_tss2=unique(ens_tss[ens_tss$eq=='TRUE',c(1,2)])
colnames(ens_tss1)=colnames(ens_tss2)=c('name','gene')
TSSens=unique(rbind(ens_tss1,ens_tss2))
TSSens=TSSens[complete.cases(TSSens$name),]

		# Annotate remaining CAGE TSS clusters by Refseq TSS (distance < 51 bp)
ost1=as.data.frame(setdiff(TC$V4,TSSens$name))
colnames(ost1)='V4'
ref_tss=merge(tssre,ost1,by='V4')
ref_tss=unique(ref_tss[ref_tss$V19.x < 51|ref_tss$V19.y < 51,])
ref_tss1=unique(ref_tss[ref_tss$eq=='FALSE',c(1,4)])
ref_tss2=unique(ref_tss[ref_tss$eq=='TRUE',c(1,2)])
colnames(ref_tss1)=colnames(ref_tss2)=c('name','gene')
TSSref=unique(rbind(ref_tss1,ref_tss2))
TSSref=TSSref[complete.cases(TSSref$name),]

		#  Prepare .bed with annotated UTRs
utr=read.table('./Annotation/temp/refseq_utrs.txt',header=T,sep='\t')
u5e=unique(en[en$type=='five_prime_utr',c(1,4,5,9,6,7)])
u3e=unique(en[en$type=='three_prime_utr',c(1,4,5,9,6,7)])
u5r=unique(utr[utr$type=='five_prime_UTR',c(2,5,6,10,7,8)])
u3r=unique(utr[utr$type=='three_prime_UTR',c(2,5,6,10,7,8)])
write.table(u5e,'./Annotation/temp/ensU5.bed',sep='\t',row.names=F,col.names=F,quote=F)
write.table(u3e,'./Annotation/temp/ensU3.bed',sep='\t',row.names=F,col.names=F,quote=F)
write.table(u5r,'./Annotation/temp/refU5.bed',sep='\t',row.names=F,col.names=F,quote=F)
write.table(u3r,'./Annotation/temp/refU3.bed',sep='\t',row.names=F,col.names=F,quote=F)

		# Intersect CAGE TSS with annotated UTRs
system('bedtools sort -i ./Annotation/temp/ensU5.bed >./Annotation/temp/ensU5.sorted.bed;
bedtools sort -i ./Annotation/temp/ensU3.bed >./Annotation/temp/ensU3.sorted.bed;
bedtools sort -i ./Annotation/temp/refU5.bed >./Annotation/temp/refU5.sorted.bed;
bedtools sort -i ./Annotation/temp/refU3.bed >./Annotation/temp/refU3.sorted.bed;
bedtools intersect -s -wo -a ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.chr.bed -b ./Annotation/temp/ensU5.sorted.bed > ./Annotation/temp/ensU5.txt;
bedtools intersect -s -wo -a ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.chr.bed -b ./Annotation/temp/ensU3.sorted.bed > ./Annotation/temp/ensU3.txt;
bedtools intersect -s -wo -a ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.chr.bed -b ./Annotation/temp/refU5.sorted.bed > ./Annotation/temp/refU5.txt;
bedtools intersect -s -wo -a ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.chr.bed -b ./Annotation/temp/refU3.sorted.bed > ./Annotation/temp/refU3.txt')

		# Annotate remaining CAGE TSS clusters by Ensembl 5'UTR
ost2=as.data.frame(setdiff(TC$V4,union(TSSens$name,TSSref$name)))
colnames(ost2)='V4'
u5e=read.table('./Annotation/temp/ensU5.txt')
u5e=unique(u5e[,c(4,16)])
u5e=merge(u5e,ost2,by='V4')

		# Annotate remaining CAGE TSS clusters by Refseq 5'UTR
ost3=as.data.frame(setdiff(ost2$V4,u5e$V4))
colnames(ost3)='V4'
u5r=read.table('./Annotation/temp/refU5.txt')
u5r=unique(u5r[,c(4,16)])
u5r=merge(u5r,ost3,by='V4')

		# Annotate remaining CAGE TSS clusters by Ensembl 3'UTR
ost4=as.data.frame(setdiff(ost3$V4,u5r$V4))
colnames(ost4)='V4'
u3e=read.table('./Annotation/temp/ensU3.txt')
u3e=unique(u3e[,c(4,16)])
u3e=merge(u3e,ost4,by='V4')

		# Annotate remaining CAGE TSS clusters by Refseq 3'UTR
ost5=as.data.frame(setdiff(ost4$V4,u3e$V4))
colnames(ost5)='V4'
u3r=read.table('./Annotation/temp/refU3.txt')
u3r=unique(u3r[,c(4,16)])
u3r=merge(u3r,ost5,by='V4')

		#  Prepare .bed with annotated CDSs
CDSe=unique(en[en$type=='CDS',c(1,4,5,9,6,7)])
CDSr=unique(ref[ref$type=='CDS',c(1,4,5,9,6,7)])
write.table(CDSe,'./Annotation/temp/ensCDS.bed',sep='\t',row.names=F,col.names=F,quote=F)
write.table(CDSr,'./Annotation/temp/refCDS.bed',sep='\t',row.names=F,col.names=F,quote=F)

		# Intersect CAGE TSS with annotated CDSs
system('bedtools sort -i ./Annotation/temp/ensCDS.bed >./Annotation/temp/ensCDS.sorted.bed;
bedtools sort -i ./Annotation/temp/refCDS.bed >./Annotation/temp/refCDS.sorted.bed;
bedtools intersect -s -wo -a ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.chr.bed -b ./Annotation/temp/ensCDS.sorted.bed > ./Annotation/temp/ensCDS.txt;
bedtools intersect -s -wo -a ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.chr.bed -b ./Annotation/temp/refCDS.sorted.bed > ./Annotation/temp/refCDS.txt')

		# Annotate remaining CAGE TSS clusters by Ensembl CDS
ost6=as.data.frame(setdiff(ost5$V4,u3r$V4))
colnames(ost6)='V4'
CDSe=read.table('./Annotation/temp/ensCDS.txt')
CDSe=unique(CDSe[,c(4,16)])
CDSe=merge(CDSe,ost6,by='V4')

		# Annotate remaining CAGE TSS clusters by Refseq CDS
ost7=as.data.frame(setdiff(ost6$V4,CDSe$V4))
colnames(ost7)='V4'
CDSr=read.table('./Annotation/temp/refCDS.txt')
CDSr=unique(CDSr[,c(4,16)])
CDSr=merge(CDSr,ost7,by='V4')

		#  Prepare .bed with annotated exons of non protein-coding transcripts
exe=unique(en[en$type=='exon',c(1,4,5,9,6,7)])
exr=unique(ref[ref$type=='exon',c(1,4,5,9,6,7)])
write.table(exe,'./Annotation/temp/exe.bed',sep='\t',row.names=F,col.names=F,quote=F)
write.table(exr,'./Annotation/temp/exr.bed',sep='\t',row.names=F,col.names=F,quote=F)

		# Intersect CAGE TSS with annotated non protein-coding transcripts
system('bedtools sort -i ./Annotation/temp/exe.bed >./Annotation/temp/exe.sorted.bed;
bedtools sort -i ./Annotation/temp/exr.bed >./Annotation/temp/exr.sorted.bed;
bedtools intersect -s -wo -a ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.chr.bed -b ./Annotation/temp/exe.sorted.bed > ./Annotation/temp/exe.txt;
bedtools intersect -s -wo -a ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.chr.bed -b ./Annotation/temp/exr.sorted.bed > ./Annotation/temp/exr.txt')

		# Annotate remaining CAGE TSS clusters by Ensembl non protein-coding transcripts
ost8=as.data.frame(setdiff(ost7$V4,CDSr$V4))
colnames(ost8)='V4'
exe=read.table('./Annotation/temp/exe.txt')
exe=unique(exe[,c(4,16)])
exe=merge(exe,ost8,by='V4')

		# Annotate remaining CAGE TSS clusters by Refseq non protein-coding transcripts
ost9=as.data.frame(setdiff(ost8$V4,exe$V4))
colnames(ost9)='V4'
exr=read.table('./Annotation/temp/exr.txt')
exr=unique(exr[,c(4,16)])
exr=merge(exr,ost9,by='V4')

		#  Prepare .bed with annotated genes
ge=unique(en[en$type=='gene',c(1,4,5,9,6,7)])
gr=unique(ref[ref$type=='gene',c(1,4,5,9,6,7)])
write.table(ge,'./Annotation/temp/ge.bed',sep='\t',row.names=F,col.names=F,quote=F)
write.table(gr,'./Annotation/temp/gr.bed',sep='\t',row.names=F,col.names=F,quote=F)

		# Intersect CAGE TSS with annotated genes
system('bedtools sort -i ./Annotation/temp/ge.bed >./Annotation/temp/ge.sorted.bed;
bedtools sort -i ./Annotation/temp/gr.bed >./Annotation/temp/gr.sorted.bed;
bedtools intersect -s -wo -a ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.chr.bed -b ./Annotation/temp/ge.sorted.bed > ./Annotation/temp/ge.txt;
bedtools intersect -s -wo -a ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.chr.bed -b ./Annotation/temp/gr.sorted.bed > ./Annotation/temp/gr.txt')

		# Annotate remaining CAGE TSS clusters by Ensembl genes (corresponding to introns)
ost10=as.data.frame(setdiff(ost9$V4,exr$V4))
colnames(ost10)='V4'
ge=read.table('./Annotation/temp/ge.txt')
ge=unique(ge[,c(4,16)])
ge=merge(ge,ost10,by='V4')

		# Annotate remaining CAGE TSS clusters by Refseq genes (corresponding to introns)
ost11=as.data.frame(setdiff(ost10$V4,ge$V4))
colnames(ost11)='V4'
gr=read.table('./Annotation/temp/gr.txt')
gr=unique(gr[,c(4,16)])
gr=merge(gr,ost11,by='V4')

		# Get remaining unannotated CAGE TSS clusters (corresponding to intergenic regions)
ost12=as.data.frame(setdiff(ost11$V4,gr$V4))
colnames(ost12)='V4'

		# Write full table with annotated CAGE TSS clusters with Ensembl IDs and gene symbols
TSSref=cbind(TSSref,'TSS_Refseq')
TSSens=cbind(TSSens,'TSS_Ensembl')
u5e=cbind(u5e,'5UTR_Ensembl')
u5r=cbind(u5r,'5UTR_Refseq')
u3e=cbind(u3e,'3UTR_Ensembl')
u3r=cbind(u3r,'3UTR_Refseq')
CDSe=cbind(CDSe,'CDS_Ensembl')
CDSr=cbind(CDSr,'CDS_Refseq')
exe=cbind(exe,'exon_non_cod_Ensembl')
exr=cbind(exr,'exon_non_cod_Refseq')
ge=cbind(ge,'intron_Ensembl')
gr=cbind(gr,'intron_Refseq')
ost12=cbind(ost12,'NA','intergenic')
colnames(TSSref)=
colnames(u5e)=
colnames(u5r)=
colnames(u3e)=
colnames(u3r)=
colnames(CDSe)=
colnames(CDSr)=
colnames(exe)=
colnames(exr)=
colnames(ge)=
colnames(gr)=
colnames(TSSens)=colnames(ost12)=c('name','gene','annotation')
ann=rbind(TSSens,u5e,u3e,CDSe,exe,ge,
TSSref,u5r,u3r,CDSr,exr,gr,ost12)
dist=tss_dist[,c(1,3,5,7,9)]
ann=merge(ann,dist,by='name',all.x=T)
ann=unique(ann)

ref_ann=ann[ann$annotation=='TSS_Refseq'|
ann$annotation=='5UTR_Refseq'|
ann$annotation=='3UTR_Refseq'|
ann$annotation=='CDS_Refseq'|
ann$annotation=='exon_non_cod_Refseq'|
ann$annotation=='intron_Refseq',]

ens_ann=ann[ann$annotation=='TSS_Ensembl'|
ann$annotation=='5UTR_Ensembl'|
ann$annotation=='3UTR_Ensembl'|
ann$annotation=='CDS_Ensembl'|
ann$annotation=='exon_non_cod_Ensembl'|
ann$annotation=='intron_Ensembl',]

int=ann[ann$annotation=='intergenic',]

ens_names=en[en$type=='gene',c(9,11)]
ens_ann2=merge(ens_ann,ens_names,by.x='gene',by.y='gene_id',all.x=T)
ref_ann2=merge(ref_ann,ens_names,by.x='gene',by.y='gene_name',all.x=T)
ens_ann2=ens_ann2[,c(8,1:7)]
ref_ann2=ref_ann2[,c(1,8,2:7)]
colnames(ref_ann2)=colnames(ens_ann2)
ann2=rbind(ens_ann2,ref_ann2)
intron=ann2[ann2$annotation=='intron_Refseq'|ann2$annotation=='intron_Ensembl',]
intron=cbind(NA,NA,intron)
int=cbind(NA,NA,NA,NA,int[,-2])
colnames(intron)=colnames(int)
nonann=rbind(intron,int)
colnames(nonann)[1:4]=c('gene_name','gene_id','intron_gene_name','intron_gene_id')
ann3=cbind(ann2[,1:2],NA,NA,ann2[,3:ncol(ann2)])
colnames(ann3)=colnames(nonann)

		# Annotate "intergenic" and "intron" CAGE TSS clusters using TSS annotation from RNA-Seq data
	# read TSS annotation from SEASTAR output
cuff=read.table('./seastar/tmp/tsgtf/nrtss.annotation',header=T,sep='\t')

	# filter TSS annotation by class_code (see more http://cole-trapnell-lab.github.io/cufflinks/file_formats/)
cuff_filt=cuff[cuff$class_code=='='|cuff$class_code=='o'|cuff$class_code=='j',]

	# create bed file with extended borders (+50 bp for upstream and downstream) for TSS annotation
cuff_bed=cuff_filt[,c(7:9,1,14,10)]
cuff_bed$tss_ss=cuff_bed$tss_ss-50
cuff_bed$tss_es=cuff_bed$tss_es+50
write.table(cuff_bed,'./Annotation/temp/nrtss.bed',sep='\t',row.names=F,col.names=F,quote=F)
	
	# intersect CAGE TSS clusters with extended TSS regions from SEASTAR output
system('bedtools sort -i ./Annotation/temp/nrtss.bed >./Annotation/temp/nrtss.sorted.bed;
bedtools intersect -s -wo -a ./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.chr.bed -b ./Annotation/temp/nrtss.sorted.bed > ./Annotation/temp/novelTSS.txt')

	# Add data with new annotated CAGE TSS clusters to annotation table
cuff_int=read.table('./Annotation/temp/novelTSS.txt',sep='\t')
cuff_int_red=cuff_int[,c(4,16)]
cuff_red=cuff[,c(1,4)]
cuff_int_ann=merge(cuff_red,cuff_int_red,by.x='tss_id',by.y='V16')
ntss=merge(nonann,cuff_int_ann,by.x='name',by.y='V4')
ntsse=merge(ens_names,ntss,by.x='gene_name',by.y='gene_short_name',all.y=T)
ntsse=unique(ntsse[,c(1,2,6,7,3,8:12)])
colnames(ntsse)=colnames(ann3)
ntsse$intron_gene_name=ntsse$intron_gene_id=NA
ntsse$annotation='Verified_by_RNA-seq'
nonann_red=nonann[nonann$name %in% setdiff(nonann$name,ntsse$name),]
ann4=ann3[ann3$annotation!='intron_Refseq'&ann3$annotation!='intron_Ensembl',]
full_ann1=unique(rbind(nonann_red,ann4,ntsse))

	# Add mean (for each condition) TPM values to annotation table
TPMs=read.table('./Expression/tpm.txt',row.names=1)
design=read.table('./Expression/design.txt',header=T,sep='\t')
D=design[design$condition=='D',1]
P0=design[design$condition=='P0',1]
P1=design[design$condition=='P1',1]
P3=design[design$condition=='P3',1]
P6=design[design$condition=='P6',1]
D=apply(TPMs[,D],1,mean)
P0=apply(TPMs[,P0],1,mean)
P1=apply(TPMs[,P1],1,mean)
P3=apply(TPMs[,P3],1,mean)
P6=apply(TPMs[,P6],1,mean)
mean=cbind(D,P0,P1,P3,P6)
full_ann2=merge(full_ann1,mean,by.x='name',by.y='row.names')
TCann=full_ann2[complete.cases(full_ann2$gene_name),c(1,2,11:15)]

	# Determine low expressed CAGE TSS clusters
max1=tapply(TCann$D,TCann$gene_name,max)
max2=tapply(TCann$P0,TCann$gene_name,max)
max3=tapply(TCann$P1,TCann$gene_name,max)
max4=tapply(TCann$P3,TCann$gene_name,max)
max5=tapply(TCann$P6,TCann$gene_name,max)
max=cbind(max1,max2,max3,max4,max5)
TCann_max=merge(TCann,max,by.x='gene_name',by.y='row.names')

for (i in 1:nrow(TCann_max)) {
	maj=ifelse(TCann_max[i,3]>0.1*TCann_max[i,8]|
	TCann_max[i,4]>0.1*TCann_max[i,9]|
	TCann_max[i,5]>0.1*TCann_max[i,10]|
	TCann_max[i,6]>0.1*TCann_max[i,11]|
	TCann_max[i,7]>0.1*TCann_max[i,12], 'major','minor')
	TCann_max[i,13]=maj
}

full_ann=unique(merge(full_ann2,TCann_max[,c(1,2,13)],all.x=T,by=c('gene_name','name')))
names(full_ann)[ncol(full_ann)]='major'

write.table(full_ann,'./Annotation/annotated_CAGE_TSS_clusters.txt',quote=F,sep='\t',row.names=F)













