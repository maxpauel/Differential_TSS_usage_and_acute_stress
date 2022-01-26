	# Generate .bed with major TSS clusters
tc_ann=read.table('./Annotation/annotated_CAGE_TSS_clusters.txt',header=T)
tc_ann_maj=tc_ann[tc_ann$major=='major',c(1,2)]
tc=read.table('./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.bed')
tc$V7=tc$V7-2000
tc$V8=tc$V8+2000
pm=merge(tc,tc_ann_maj,by.x='V4',by.y='name')
write.table(unique(pm[,c(2,3,4,1,5,6)]),'./Promoters/temp/TSS_maj.bed',quote=F,col.names=F,sep='\t',row.names=F)
	# Generate .bed with extended promoter regions around "major" TSS clusters (-2000..+2000bp around "major" CAGE TSS cluster)
write.table(unique(pm[,c(2,7,8,1,5,6)]),'./Promoters/temp/promoter_maj_extended.bed',quote=F,col.names=F,sep='\t',row.names=F)
system('sort -k 1,1 -k2,2n ./Promoters/temp/TSS_maj.bed > ./Promoters/temp/TSS_maj.sorted.bed')
system('bedtools sort -i ./Promoters/temp/promoter_maj_extended.bed > ./Promoters/temp/promoter_maj_extended.sorted.bed')
system('bedtools merge -s -c 6 -o distinct  -i ./Promoters/temp/promoter_maj_extended.sorted.bed > ./Promoters/temp/promoter_maj_extended.sorted.merged.bed')
	# Get intersection between extended promoter regions and open chromatin region (OCR)
system('bedtools intersect -a ./Promoters/temp/promoter_maj_extended.sorted.merged.bed -b ./Promoters/OCR.bed > ./Promoters/temp/promoter_chromatin_intersect.bed')
system('sort -k 1,1 -k2,2n ./Promoters/temp/promoter_chromatin_intersect.bed >  ./Promoters/temp/promoter_chromatin_intersect.sorted.bed')
	# Merge OCR intervals by distance 70bp
system('bedtools merge -d 70 -i ./Promoters/temp/promoter_chromatin_intersect.sorted.bed > ./Promoters/temp/promoter_chromatin_intersect.sorted.merge70.bed')
	# Get distance from CAGE TSS clusters to OCR intervals
system('bedtools closest -d -a ./Promoters/temp/TSS_maj.sorted.bed -b ./Promoters/temp/promoter_chromatin_intersect.sorted.merge70.bed > ./Promoters/temp/TSS_to_openr_chromatin_distance.bed')
	# Get mean expression for each CAGE TSS cluster
mean=apply(tc_ann[,11:15],1,mean)
am=cbind(tc_ann,mean)
n=unique(am[,c(2,17)])
tp=read.table('./Promoters/temp/promoter_maj_extended.sorted.bed')
tpn=merge(tp,n,by.x='V4',by.y='name')
tpn$V5=tpn$mean
tpn=tpn[,c(2,3,4,1,5,6)]
write.table(tpn,'./Promoters/temp/prom_maj_TSS_TPM.bed',quote=F,col.names=F,sep='\t',row.names=F)
	# Get open chromatin regions asocciated with major TSS clusters (distance < 201bp)
a=read.table('./Promoters/temp/TSS_to_openr_chromatin_distance.bed')
a200=unique(a[a$V10<201,])
tpn=read.table('./Promoters/temp/prom_maj_TSS_TPM.bed')
tpnap=merge(tpn,a200,by='V4')
tpnap=tpnap[,c(12,13,14,1,5,6)]
write.table(tpnap,'./Promoters/temp/prom_open.bed',quote=F,col.names=F,sep='\t',row.names=F)
	# Merge open chromatin regions asocciated with major TSS clusters and get maximal TFBS coverage on OCR interval
system('sort -k 1,1 -k2,2n ./Promoters/temp/prom_open.bed > ./Promoters/temp/prom_open.sorted.bed')
system('bedtools merge -s -c 6,5 -o distinct,sum -i ./Promoters/temp/prom_open.sorted.bed > ./Promoters/temp/prom_open.sorted.merged.bed')
temp1=read.table('./Promoters/temp/prom_open.sorted.merged.bed',sep='\t')
temp2=cbind(temp1[,1:3],1:nrow(temp1),temp1$V5,temp1$V4)
write.table(temp2,'./Promoters/temp/prom_open.sorted.mergedN.bed',quote=F,col.names=F,sep='\t',row.names=F)
system('bedtools intersect -u -a ./Promoters/Chip_seq_data/Chip_seq.sorted.bed -b ./Promoters/temp/prom_open.sorted.mergedN.bed > ./Promoters/temp/prom_open.sorted.merged.chip.int.bed')
system('bedtools genomecov -bga -i ./Promoters/temp/prom_open.sorted.merged.chip.int.bed -g ./chrNameLength.txt > ./Promoters/temp/prom_open.sorted.merged.chip.int.BedGraph')
temp1=read.table('./Promoters/temp/prom_open.sorted.merged.chip.int.BedGraph',sep='\t')
temp2=cbind(temp1[,1:3],1:nrow(temp1),temp1$V4)
write.table(temp2,'./Promoters/temp/prom_open.sorted.merged.chip.int.BedGraph.bed',quote=F,col.names=F,sep='\t',row.names=F)
system('bedtools intersect -wao -a ./Promoters/temp/prom_open.sorted.merged.chip.int.BedGraph.bed -b ./Promoters/temp/prom_open.sorted.mergedN.bed > ./Promoters/temp/prom_open.sorted.intersect2.txt')
chi=read.table('./Promoters/temp/prom_open.sorted.intersect2.txt',header=F,sep='\t')
max=tapply(chi$V5,chi$V9,max)
max=as.data.frame(max)
pr=read.table('./Promoters/temp/prom_open.sorted.mergedN.bed',sep='\t')
prm=merge(max,pr,by.x='row.names',by.y='V4')
write.table(prm,'./Promoters/temp/prom_open.max_per_prom.txt',quote=F,row.names=F,sep='\t')
	# Generate final table with individual promoters associated with open chromatin
ap=merge(a200,prm,by.x=c('V8','V9','V6'),by.y=c('V2','V3','V6'))
ap=ap[,c(9,1,2,7,10,11,12,14)]
apt=merge(ap,tpn,by='V4')	#!
n=tapply(apt$V5,apt$Row.names,max)
n=as.data.frame(n)
apn=merge(apt,n,by.x=c('Row.names','V5'),by.y=c('row.names','n'))
tc=read.table('./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.bed')
tc=tc[,c(4,7,8)]
apntc=merge(apn,tc,by='V4')
write.table(apntc,'./Promoters/temp/open_chr.length.txt',quote=F,row.names=F,sep='\t')
	# Get all CAGE TSS clusters associated with OCR
system('bedtools intersect -s -wao -a ./Promoters/temp/prom_open.sorted.bed -b ./Promoters/temp/prom_open.sorted.mergedN.bed > ./Promoters/temp/all_TC_open.txt')
	# Get all CAGE TSS without OCR association (>200bp distance to OCR)
an200=a[a$V10>200,c(1:4,10,6)]
write.table(an200,'./Promoters/temp/nonchr_TC.bed',quote=F,row.names=F,sep='\t',col.names=F)
system('bedtools merge -s -d 200 -c 5,6 -o min,distinct -i ./Promoters/temp/nonchr_TC.bed > ./Promoters/temp/nonchr_TC.merged.bed')
temp1=read.table('./Promoters/temp/nonchr_TC.merged.bed',sep='\t')
temp2=read.table('./Promoters/temp/prom_open.sorted.mergedN.bed',sep='\t')
temp3=cbind(temp1[,1:3],(nrow(temp2)+1):(nrow(temp2)+nrow(temp1)),temp1[,4:5])
write.table(temp3,'./Promoters/temp/nonchr_TC.mergedN.bed',quote=F,row.names=F,sep='\t',col.names=F)
system('bedtools intersect -s -wao -a ./Promoters/temp/nonchr_TC.bed -b ./Promoters/temp/nonchr_TC.mergedN.bed >./Promoters/temp/all_TC_closed.txt')
	# Combine all CAGE TSS clusters annotated by promoters_id
cTC=read.table('./Promoters/temp/all_TC_closed.txt')
oTC=read.table('./Promoters/temp/all_TC_open.txt')
cTCr=cTC[,1:10]
oTCr=oTC[,1:10]
TCr=rbind(oTCr,cTCr)
write.table(TCr,'./Promoters/temp/all_TC.bed',quote=F,row.names=F,col.names=F,sep='\t')
	# Get table of promoters without OCR
n=unique(am[,c(2,17)])
bn=merge(cTC,n,by.x='V4',by.y='name')
bnmax=tapply(bn$mean,bn$V10,max)
bnmax=as.data.frame(bnmax)
bm1=merge(bn,bnmax,by.x=c('V10','mean'),by.y=c('row.names','bnmax'))
bnm=tapply(bn$mean,bn$V10,sum)
bnm=as.data.frame(bnm)
bm2=merge(bm1,bnm,by.x=c('V10'),by.y=c('row.names'))
tc=read.table('./DPI/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.bed',sep='\t')
tc=tc[,c(4,7,8)]
m=merge(bm2,tc,by='V4')
m=m[,c(1:6,12,8:ncol(m))]
write.table(m,'./Promoters/temp/promoters_closed.txt',quote=F,row.names=F,sep='\t')
	# Get chip-seq coverage for promoters without OCR
cl_bord=cbind(as.data.frame(m$V1),m$V7.y-2000,m$V8.y+2000,m$V10,m$bnm,as.data.frame(m$V6))
write.table(cl_bord,'./Promoters/temp/closed_border.bed',quote=F,row.names=F,sep='\t',col.names=F)
system('sort -k 1,1 -k2,2n ./Promoters/temp/closed_border.bed > ./Promoters/temp/closed_border.sorted.bed')
system('bedtools intersect -u -a ./Promoters/Chip_seq_data/Chip_seq.sorted.bed -b ./Promoters/temp/closed_border.sorted.bed > ./Promoters/temp/closed_border.sorted.int.chip.bed')
system('bedtools genomecov -bga -i ./Promoters/temp/closed_border.sorted.int.chip.bed -g ./chrNameLength.txt > ./Promoters/temp/closed_border.sorted.int.chip.BedGraph')
temp1=read.table('./Promoters/temp/closed_border.sorted.int.chip.BedGraph',sep='\t')
temp2=cbind(temp1[,1:3],1:nrow(temp1),temp1$V4)
write.table(temp2,'./Promoters/temp/closed_border.sorted.int.chip.BedGraph.bed',quote=F,col.names=F,sep='\t',row.names=F)
system('bedtools intersect -wao -a ./Promoters/temp/closed_border.sorted.int.chip.BedGraph.bed -b ./Promoters/temp/closed_border.sorted.bed > ./Promoters/temp/closed_border.sorted.chip.inter.txt')
cl_chip=read.table('./Promoters/temp/closed_border.sorted.chip.inter.txt',header=F,sep='\t')
cl=read.table('./Promoters/temp/promoters_closed.txt',header=T,sep='\t')
max=tapply(cl_chip$V5,cl_chip$V9,max)
max=as.data.frame(max)
closed=merge(cl,max,by.x='V10',by.y='row.names')
write.table(closed,'./Promoters/temp/promoters_closed_final.txt',quote=F,row.names=F,sep='\t')
	# Get full table of promoters
open=read.table('./Promoters/temp/open_chr.length.txt',header=T,sep='\t')
closed=read.table('./Promoters/temp/promoters_closed_final.txt',header=T,sep='\t')
closed=closed[,-12]
closed=closed[,c(2,1,3:7,17,14,9,15,16,8,15,16)]
colnames(closed)=colnames(open)
all=rbind(open,closed)
	# Add information about bidirect promoters to full table of promoters
system('bedtools intersect -S -a ./Promoters/temp/prom_open.sorted.mergedN.bed -b ./Promoters/temp/prom_open.sorted.mergedN.bed > ./Promoters/temp/bidirect.bed')
###add bidirect
bid=read.table('./Promoters/temp/bidirect.bed',header=F)
bid=bid$V4
bid=as.data.frame(bid)
bid=cbind(bid,'bidirect')
all=merge(all,bid,by.x='Row.names',by.y='bid',all=T)
	# Estimate differential expression of promoters
counts=read.table('./Expression/counts.txt',header=T,sep='\t')
TC=read.table('./Promoters/temp/all_TC.bed')
TCr=TC[,c(4,10)]
TC_counts=merge(TCr,counts,by.x='V4',by.y='row.names')
prom_exp=aggregate(TC_counts[,3:52],by=list(TC_counts$V10), sum)
data=prom_exp[,2:51]
library(DESeq)
library(DESeq2)
design=read.table('./Expression/design.txt',row.names=1,header=T)
design$person=as.factor(design$person)
design1=design[design$condition=='D'|design$condition=='P1',]
design2=design[design$condition=='D'|design$condition=='P3',]
design3=design[design$condition=='D'|design$condition=='P6',]
data1=data[,as.numeric(rownames(design1))]
data2=data[,as.numeric(rownames(design2))]
data3=data[,as.numeric(rownames(design3))]

data1 = newCountDataSet( data1, design1 )
dds=DESeqDataSetFromMatrix(countData=counts(data1),
colData=design1,
design=~person+condition)
dds <- DESeq(dds)
res1=results(dds,contrast=c('condition','P1','D'))
write.table(res1,'./Promoters/Dif_expression/P1-D.txt',quote=F,row.names=T,sep='\t')

data2 = newCountDataSet( data2, design2 )
dds=DESeqDataSetFromMatrix(countData=counts(data2),
colData=design2,
design=~person+condition)
dds <- DESeq(dds)
res2=results(dds,contrast=c('condition','P3','D'))
write.table(res2,'./Promoters/Dif_expression/P3-D.txt',quote=F,row.names=T,sep='\t')

data3 = newCountDataSet( data3, design3 )
dds=DESeqDataSetFromMatrix(countData=counts(data3),
colData=design3,
design=~person+condition)
dds <- DESeq(dds)
res3=results(dds,contrast=c('condition','P6','D'))
write.table(res3,'./Promoters/Dif_expression/P6-D.txt',quote=F,row.names=T,sep='\t')

res1=res1[,c(2,6)]
res2=res2[,c(2,6)]
res3=res3[,c(2,6)]
colnames(res1)=c('log2FoldChange_1h','Padj_1h')
colnames(res2)=c('log2FoldChange_3h','Padj_3h')
colnames(res3)=c('log2FoldChange_6h','Padj_6h')
res=cbind(res1,res2,res3)
write.table(res,'./Promoters/Dif_expression/DEP.txt',quote=F,row.names=T,sep='\t')
	# Add information about DE promoters to full table of promoters
promoters_1=merge(res,all,by.x='row.names',by.y='Row.names')
write.table(promoters_1,'./Promoters/temp/promoters_1.txt',quote=F,row.names=F,sep='\t')
	# Annotate full table of promoters
tc_ann_red=unique(tc_ann[,c(1,2,3)])
promoters_2=merge(promoters_1,tc_ann_red,by.y='name',by.x='V4',all.x=T)
write.table(promoters_2,'./Promoters/temp/promoters_2.txt',quote=F,row.names=F,sep='\t')

tss=read.table('./Promoters/temp/all_TC.bed',header=F)
tss=tss[,c(4,10)]
tsa=merge(tss,tc_ann_red,by.y='name',by.x='V4',all.x=T)
tsa=unique(tsa[,2:4])
promoters_3=merge(promoters_1,tsa,by.y='V10',by.x='Row.names',all.x=T)
write.table(promoters_3,'./Promoters/temp/promoters_3.txt',quote=F,row.names=F,sep='\t')
	# Estimate differential usage of promoters
promoters_3r=promoters_3[,c(1,23)]
colnames(promoters_3r)=c('feature_id','gene_id')
counts=read.table('./Expression/DESeq2_normalized_counts.txt',header=T)
tss_counts=merge(tss,counts,by.x='V4',by.y='row.names')
ex=aggregate(tss_counts[,3:52],by=list(tss_counts$V10), sum)
exa=merge(promoters_3r,ex,by.x='feature_id',by.y='Group.1')
exa=unique(exa)
exa$feature_id=1:nrow(exa)

library(DRIMSeq)
data=exa
design=read.table('./Expression/design.txt',row.names=1,header=T)
design1=design[design$condition=='D'|design$condition=='P1',]
design2=design[design$condition=='D'|design$condition=='P3',]
design3=design[design$condition=='D'|design$condition=='P6',]
data1=data[,c(1,2,(as.numeric(rownames(design1))+2))]
data2=data[,c(1,2,(as.numeric(rownames(design2))+2))]
data3=data[,c(1,2,(as.numeric(rownames(design3))+2))]

metadata=as.data.frame(cbind(colnames(data1)[3:ncol(data1)],design1$condition))
colnames(metadata)=c('SampleName', 'condition')
samples <- data.frame(sample_id = metadata$SampleName,group = metadata$condition)
data1$feature_id=as.character(data1$feature_id)
d <- dmDSdata(counts =data1, samples =samples)
design_full <- model.matrix(~ group, data = samples)
set.seed(123)
d <- dmPrecision(d, design = design_full)
common_precision(d)
d <- dmFit(d, design = design_full, verbose = 1)
d <- dmTest(d, coef = "group2", verbose = 1)
design(d)
res <- results(d)
res1 <- res[order(res$pvalue, decreasing = FALSE), ]

metadata=as.data.frame(cbind(colnames(data2)[3:ncol(data2)],design2$condition))
colnames(metadata)=c('SampleName', 'condition')
samples <- data.frame(sample_id = metadata$SampleName,group = metadata$condition)
data2$feature_id=as.character(data2$feature_id)
d <- dmDSdata(counts =data2, samples =samples)
design_full <- model.matrix(~ group, data = samples)
set.seed(123)
d <- dmPrecision(d, design = design_full)
common_precision(d)
d <- dmFit(d, design = design_full, verbose = 1)
d <- dmTest(d, coef = "group2", verbose = 1)
design(d)
res <- results(d)
res2 <- res[order(res$pvalue, decreasing = FALSE), ]

metadata=as.data.frame(cbind(colnames(data3)[3:ncol(data3)],design3$condition))
colnames(metadata)=c('SampleName', 'condition')
samples <- data.frame(sample_id = metadata$SampleName,group = metadata$condition)
data3$feature_id=as.character(data3$feature_id)
d <- dmDSdata(counts =data3, samples =samples)
design_full <- model.matrix(~ group, data = samples)
set.seed(123)
d <- dmPrecision(d, design = design_full)
common_precision(d)
d <- dmFit(d, design = design_full, verbose = 1)
d <- dmTest(d, coef = "group2", verbose = 1)
design(d)
res <- results(d)
res3 <- res[order(res$pvalue, decreasing = FALSE), ]

write.table(res1,'./Promoters/Dif_expression/P1-D_dif_promoter_usage.txt',quote=F,sep='\t',row.names=F)
write.table(res2,'./Promoters/Dif_expression/P3-D_dif_promoter_usage.txt',quote=F,sep='\t',row.names=F)
write.table(res3,'./Promoters/Dif_expression/P6-D_dif_promoter_usage.txt',quote=F,sep='\t',row.names=F)

res1=res1[,c(1,5)]
res2=res2[,c(1,5)]
res3=res3[,c(1,5)]
colnames(res1)[2]='Dif_usage_adj_pvalue_1h'
colnames(res2)[2]='Dif_usage_adj_pvalue_3h'
colnames(res3)[2]='Dif_usage_adj_pvalue_6h'

l=list(res1,res2,res3)
res=Reduce(merge, lapply(l, function(x) data.frame(x, rn = x$gene_id)))
write.table(res[,-2],'./Promoters/Dif_expression/Dif_promoter_usage.txt',quote=F,sep='\t',row.names=F)
	# Add information about differential usage of promoters to full table of promoters
promoters_4=merge(res[,-2],promoters_3,by.x='gene_id',by.y='gene_name',all.y=T)
colnames(promoters_4)[c(5,1,27,14,15,16,23,12,17,19,18,26)]=c('promoter_id',
'gene_name',
'gene_id',
'chr',
'OCR_start',
'OCR_end',
'strand',
'main_tss_cluster',
'distance_tss_cluster_to_open_chromatin',
'mean_promoter_expression',
'max_tfbs_density',
'bidirect')
promoters_4=promoters_4[order(as.numeric(promoters_4$promoter_id)),c(5,1,27,14,15,16,23,12,17,19,18,26,6:11,2:4)]
write.table(promoters_4,'./Promoters/Promoters.txt',quote=F,sep='\t',row.names=F)
