	# Generate table with DESeq2_normalized counts for promoters
counts=read.table('./Expression/DESeq2_normalized_counts.txt',header=T,sep='\t')
tc=read.table('./Promoters/temp/all_TC.bed',sep='\t')
tc=tc[,c(4,10)]
tc_counts=merge(tc,counts,by.x='V4',by.y='row.names')
ex=aggregate(abe[,3:52],by=list(abe$V10), sum)
data=ex[,2:51]
rownames(data)=ex$Group.1
	# Get lists of differentially expressed promoters
a=read.table('./Promoters/Dif_expression/DEP.txt',header=T)
d1=unique(a[a$log2FoldChange_1h< -0.3219281&a$Padj_1h<0.01,])
write.table(rownames(d1),'./CRC_clustering/temp/d1.txt',row.names=F,quote=F,col.names=F)
d3=unique(a[a$log2FoldChange_3h< -0.3219281&a$Padj_3h<0.01,])
write.table(rownames(d3),'./CRC_clustering/temp/d3.txt',row.names=F,quote=F,col.names=F)
d6=unique(a[a$log2FoldChange_6h< -0.3219281&a$Padj_6h<0.01,])
write.table(rownames(d6),'./CRC_clustering/temp/d6.txt',row.names=F,quote=F,col.names=F)
u1=unique(a[a$log2FoldChange_1h>0.3219281&a$Padj_1h<0.01,])
write.table(rownames(u1),'./CRC_clustering/temp/u1.txt',row.names=F,quote=F,col.names=F)
u3=unique(a[a$log2FoldChange_3h>0.3219281&a$Padj_3h<0.01,])
write.table(rownames(u3),'./CRC_clustering/temp/u3.txt',row.names=F,quote=F,col.names=F)
u6=unique(a[a$log2FoldChange_6h>0.3219281&a$Padj_6h<0.01,])
write.table(rownames(u6),'./CRC_clustering/temp/u6.txt',row.names=F,quote=F,col.names=F)
all=unique(c(rownames(d1),rownames(d3),rownames(d6),rownames(u1),rownames(u3),rownames(u6)))
write.table(all,'./CRC_clustering/temp/all.txt',row.names=F,quote=F,col.names=F)
all=read.table('./CRC_clustering/temp/all.txt',header=F)
	# Rank DESeq2_normalized expression for promoters (from 1 to 4 for each person)
tab_ranked=t(rbind(apply(data[,c(1,21,31,41)],1,rank),
apply(data[,c(2,22,32,42)],1,rank),
apply(data[,c(3,23,33,43)],1,rank),
apply(data[,c(4,24,34,44)],1,rank),
apply(data[,c(5,25,35,45)],1,rank),
apply(data[,c(6,26,36,46)],1,rank),
apply(data[,c(7,27,37,47)],1,rank),
apply(data[,c(8,28,38,48)],1,rank),
apply(data[,c(9,29,39,49)],1,rank),
apply(data[,c(10,30,40,50)],1,rank)))
	# Extract subset of differentially expressed promoters from table with ranked expression
dep=tab_ranked[all$V1,]
write.table(tab_ranked,'./CRC_clustering/tab_ranked.txt',quote=F,sep='\t')
write.table(dep,'./CRC_clustering/tab_ranked_dep.txt',quote=F,sep='\t')

