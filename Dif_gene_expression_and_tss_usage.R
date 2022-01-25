		# Prepare count table for annotated CAGE TSS clusters
counts=read.table('./Expression/counts.txt',sep='\t')
tc_ann=read.table('./Annotation/annotated_CAGE_TSS_clusters.txt',sep='\t',header=T)
tc_ann_red=unique(tc_ann[tc_ann$annotation!='intron_Ensembl'&tc_ann$annotation!='intron_Refseq'&tc_ann$annotation!='intergenic',1:2])
counts_ann=merge(tc_ann_red,counts,by.x='name',by.y='row.names')

		# Prepare count table by gene
data_genes=na.omit(aggregate(counts_ann[,3:52],by=list(counts_ann$gene_name), sum))
counts_genes=data_genes[,2:51]
rownames(counts_genes)=data_genes$Group.1

		# Estimate differential gene expression by conditions
library(DESeq)
library(DESeq2)
design=read.table('./Expression/design.txt',row.names=1,header=T)
design$person=as.factor(design$person)
design1=design[design$condition=='D'|design$condition=='P1',]
design2=design[design$condition=='D'|design$condition=='P3',]
design3=design[design$condition=='D'|design$condition=='P6',]
data1=counts_genes[,as.numeric(rownames(design1))]
data2=counts_genes[,as.numeric(rownames(design2))]
data3=counts_genes[,as.numeric(rownames(design3))]

data1 = newCountDataSet( data1, design1 )
dds=DESeqDataSetFromMatrix(countData=counts(data1),
colData=design1,
design=~person+condition)
dds <- DESeq(dds)
res1=results(dds,contrast=c('condition','P1','D'))
write.table(res1,'./Dif_gene_expression_and_tss_usage/P1-D_genes.txt',quote=F,row.names=T,sep='\t')

data2 = newCountDataSet( data2, design2 )
dds=DESeqDataSetFromMatrix(countData=counts(data2),
colData=design2,
design=~person+condition)
dds <- DESeq(dds)
res2=results(dds,contrast=c('condition','P3','D'))
write.table(res2,'./Dif_gene_expression_and_tss_usage/P3-D_genes.txt',quote=F,row.names=T,sep='\t')

data3 = newCountDataSet( data3, design3 )
dds=DESeqDataSetFromMatrix(countData=counts(data3),
colData=design3,
design=~person+condition)
dds <- DESeq(dds)
res3=results(dds,contrast=c('condition','P6','D'))
write.table(res3,'./Dif_gene_expression_and_tss_usage/P6-D_genes.txt',quote=F,row.names=T,sep='\t')

res1=res1[,c(2,6)]
res2=res2[,c(2,6)]
res3=res3[,c(2,6)]
colnames(res1)=c('log2FoldChange_1h','Padj_1h')
colnames(res2)=c('log2FoldChange_3h','Padj_3h')
colnames(res3)=c('log2FoldChange_6h','Padj_6h')
res1$DEG_1h=(res1$Padj_1h<0.01)&(res1$log2FoldChange_1h< -0.3219281)|(res1$log2FoldChange_1h> 0.3219281)
res2$DEG_3h=(res2$Padj_3h<0.01)&(res2$log2FoldChange_3h< -0.3219281)|(res2$log2FoldChange_3h> 0.3219281)
res3$DEG_6h=(res3$Padj_6h<0.01)&(res3$log2FoldChange_6h< -0.3219281)|(res3$log2FoldChange_6h> 0.3219281)
l=list(res1,res2,res3)
res=Reduce(merge, lapply(l, function(x) data.frame(x, rn = row.names(x))))
colnames(res)[1]='gene_name'
write.table(res,'./Dif_gene_expression_and_tss_usage/Dif_gene_expression.txt',quote=F,sep='\t',row.names=F)

		# Get DESeq2 normalized counts for CAGE TSS clusters
data = newCountDataSet( counts, design )
dds=DESeqDataSetFromMatrix(countData=counts(data),
colData=design,
design=~person+condition)
dds <- DESeq(dds)
counts_norm=counts(dds, normalized=TRUE)
write.table(counts_norm,'./Expression/DESeq2_normalized_counts.txt',quote=F,row.names=T,sep='\t')

		# Estimate differential CAGE TSS usage by conditions for each gene
counts_norm_ann=merge(tc_ann_red,counts_norm,by.x='name',by.y='row.names')
counts_norm_ann$name=rownames(counts_norm_ann)
colnames(counts_norm_ann)[1:2] = c('feature_id','gene_id')
library(DRIMSeq)
data=na.omit(counts_norm_ann)
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

write.table(res1,'./Dif_gene_expression_and_tss_usage/P1-D_dif_tss_usage.txt',quote=F,sep='\t',row.names=F)
write.table(res2,'./Dif_gene_expression_and_tss_usage/P3-D_dif_tss_usage.txt',quote=F,sep='\t',row.names=F)
write.table(res3,'./Dif_gene_expression_and_tss_usage/P6-D_dif_tss_usage.txt',quote=F,sep='\t',row.names=F)

res1=res1[,c(1,5)]
res2=res2[,c(1,5)]
res3=res3[,c(1,5)]
colnames(res1)[2]='adj_pvalue_1h'
colnames(res2)[2]='adj_pvalue_3h'
colnames(res3)[2]='adj_pvalue_6h'
res1$Dif_TSS_Usage_1h=res1$adj_pvalue_1h<0.05
res2$Dif_TSS_Usage_3h=res2$adj_pvalue_3h<0.05
res3$Dif_TSS_Usage_6h=res3$adj_pvalue_6h<0.05
l=list(res1,res2,res3)
res=Reduce(merge, lapply(l, function(x) data.frame(x, rn = row.names(x))))
write.table(res[,-2],'./Dif_gene_expression_and_tss_usage/Dif_tss_usage.txt',quote=F,sep='\t',row.names=F)

