library(Seurat)
source("heatmap.R")
source("niceDraw.R")
library(dplyr)
library(tidyr)
library(ggplot2)
library(NMF)
library(reshape)
library(cowplot)

genFigures<-function(seur,outdir="outdir")
{
backup=seur

system(paste("mkdir",outdir))

seur<-SetAllIdent(seur,"CellType")

print("Generate TSNEs")
pdf(paste(outdir,"TSNE_CellType.pdf",sep="/"))
TSNEPlot(seur,)
dev.off();

pdf(paste(outdir,"TSNE_Lane.pdf",sep="/"))
TSNEPlot(SetAllIdent(seur,"lane"),)
dev.off();


print("Generate Cell type bar plot")
tmp=cast(seur@data.info,CellType~lane,value="nGene")
tmp=data.frame(tmp)
colnames(tmp)=c("CellType","lane 3","lane 4")
mlt<-melt(tmp,id="CellType")
tmp=mlt
colnames(tmp)=c("CellType","lane","Number")
p=ggplot(tmp,aes(x=CellType,y=Number,fill=lane))+geom_bar(stat="identity",position="dodge")+xlab("")+ylab("Number of Cells")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_y_continuous(expand = c(0,0))
pdf(paste(outdir,"BarPlot_celltype.pdf",sep="/"))
plot(p)
dev.off()

print("Generate Cell type bar plot")
tmp=cast(seur@data.info,CellType~lane,value="nGene")
tmp=data.frame(tmp)
colnames(tmp)=c("CellType","lane 3","lane 4")
mlt<-melt(tmp,id="CellType")
tmp=mlt
colnames(tmp)=c("CellType","lane","Number")
p=ggplot(tmp,aes(x=CellType,y=Number,fill=lane))+geom_bar(stat="identity",position="dodge")+xlab("")+ylab("Number of Cells")+coord_flip()+scale_y_continuous(expand = c(0,0))
pdf(paste(outdir,"BarPlot_celltype_flip.pdf",sep="/"))
plot(p)
dev.off()








load("markers.Robj")

print("Generate Feature Plots")
system(paste("mkdir",paste(outdir,"featurePlots",sep="/")))
tab=tab[tab[,3]%in%rownames(seur@data),]
lst=tab[,2]

if(FALSE)#for(l in lst)
{
p=niceFeaturePlot(seur,sub(" ","",paste("hg19_",l,sep="")),low="blue")+ggtitle(l)
pdf(paste(outdir,"/featurePlots/",l,".pdf",sep=""))
plot(p)
dev.off()

}

p=niceFeaturePlot(seur,"G1S",datainfo=T,low="blue")+ggtitle("G1S score")
pdf(paste(outdir,"/G1S.pdf",sep=""))
plot(p)
dev.off()



p=niceFeaturePlot(seur,"G2M",datainfo=T,low="blue")+ggtitle("G2M score")
pdf(paste(outdir,"/G2M.pdf",sep=""))
plot(p)
dev.off()


seur<-SetAllIdent(seur,"CellType_ext2")
uniq=unique(as.character(seur@ident))
print(uniq)
uniq=uniq[uniq!="Radial Glia_3"]
print(uniq)

seur=SubsetData(seur,WhichCells(seur,uniq))

print("Generate Gene by Cluster Heatmap")

tab=tab[!duplicated(tab[,2]),]
inter=tab[,3]
ave=AverageExpression(seur,genes.use=inter)
rownames(ave)=sub("hg19_","",rownames(ave))
tab=data.frame(tab[,1])


mn=apply(ave,1,max)
print(mn)
ave=ave[mn>3,]
tab=tab[mn>3,]
tab=data.frame(tab)

print(dim(tab))
print(dim(ave))

ave=data.frame(ave)[order(tab[,1],decreasing=T),]
tab=tab[order(tab[,1],decreasing=T),]

tab=data.frame(tab)

colnames(tab)="Cell Type"
tab["Cell Type"]=as.character(tab[,1])

filename=paste(outdir,"genes_celltype.pdf",sep="/")
rowNames=c("Neuroepithelial","Radial Glia","Astrocytes","Inhibitory Neurons","Neurons: lane 4","Neurons: lane 3")
ave=ave[,c(3,6,1,2,5,4)]
print(colnames(ave))

colnames(ave)=rowNames
aheatmap(ave,Rowv=seq(nrow(ave),1),Colv=seq(ncol(ave),1),annRow=tab,filename=filename)


rowNames=c("Neuroepithelial","Radial Glia","Astrocytes","Inhibitory Neurons","Neurons: lane 4","Neurons: lane 3")


print("Compare to Pollen!")
load("Pollen.Robj")
COR=getCOR(seur,Pollen)
COR=data.frame(COR)[,rev(c(4,5,2,1,6,3))]
COR=data.frame(COR)[c(4,2,1,3),]
ann=data.frame(factor(rev(c(3,4,4,4,4,4))))
colnames(ann)="Lane"
aheatmap(t(COR),Rowv=seq(ncol(COR),1),Colv=seq(nrow(COR),1),filename=paste(outdir,"Pollen.pdf",sep="/"),annRow=ann , labRow=rowNames)


print("Compare to new")
load("devo_keep.Robj")
COR=getCOR(seur,devo_keep)
COR=data.frame(COR)[,rev(c(4,5,2,1,6,3))]
ann=data.frame(factor(rev(c(3,4,4,4,4,4))))
colnames(ann)="Lane"
load("list_new.Robj")
tab=c()
type=c("Excitatory","Newborn Excitatory","Inhibitory","Newborn Inhibitory","IPC","Astrocytes","Radial Glia")
for(i in 1:length(lst)){tab=c(tab,rep(type[i],length(lst[[i]])))}
tab=rev(tab)
tab=data.frame(tab)
colnames(tab)="Cell Type"
COR=COR[rev(unlist(lst)),]
aheatmap(t(COR),Rowv=seq(ncol(COR),1),Colv=seq(nrow(COR),1),filename=paste(outdir,"new_Pollen.pdf",sep="/"),annRow=ann , labRow=rowNames,annCol=tab,width=14)
aheatmap(t(COR),Rowv=seq(ncol(COR),1),filename=paste(outdir,"new_Pollen_clustered.pdf",sep="/"),annRow=ann ,labRow=rowNames)

print("Compare to cell type new")
devo_keep<-SetAllIdent(devo_keep,"CellType")
COR=getCOR(seur,devo_keep)
COR=data.frame(COR)[,rev(c(4,5,2,1,6,3))]
COR=COR[c(7,1,4,6,3,5,2),]
aheatmap(t(COR),Rowv=seq(ncol(COR),1),Colv=seq(nrow(COR),1),filename=paste(outdir,"new_Pollen_cellType.pdf",sep="/"),annRow=ann , labRow=rowNames)


print("Compare to Dronc")
load("dronc.Robj")

dronc=SubsetData(dronc,WhichCells(dronc,c("Astrocytes","Excitatory","Inhibitory")))
COR=getCOR(seur,dronc)

COR=data.frame(COR)[,rev(c(4,5,2,1,6,3))]
aheatmap(t(COR),Rowv=seq(ncol(COR),1),filename=paste(outdir,"Dronc.pdf",sep="/"),annRow=ann , labRow=rowNames)

print("Compare to Organoids!")
load("org11a.Robj")
org11a<-SetAllIdent(org11a,"tree.ident")
COR=getCOR(seur,org11a)

COR=data.frame(COR)[,rev(c(4,5,2,1,6,3))]
aheatmap(t(COR),Rowv=seq(ncol(COR),1),filename=paste(outdir,"org11a.pdf",sep="/"),annRow=ann ,labRow=rowNames)


print("Compare to Organoids, Cell Type Names!")
load("sub_org.Robj")
COR=getCOR(seur,sub_org)

COR=data.frame(COR)[,rev(c(4,5,2,1,6,3))]
COR=data.frame(COR)[c(4,5,1,2,3),]
print(max(COR))
print(min(COR))
aheatmap(t(COR),Rowv=seq(ncol(COR),1),Colv=seq(nrow(COR),1),filename=paste(outdir,"suborg11a.pdf",sep="/"),annRow=ann ,labRow=rowNames,breaks=seq(.2,.8,length=100))


print("Compare to Organoids, Forebrain!")
load("fore.Robj")
COR=getCOR(seur,fore)
print(max(COR))
print(min(COR))

COR=data.frame(COR)[,rev(c(4,5,2,1,6,3))]
COR=data.frame(COR)[rev(rownames(COR)),]
aheatmap(t(COR),Rowv=seq(ncol(COR),1),Colv=seq(nrow(COR),1),filename=paste(outdir,"forebrain.pdf",sep="/"),annRow=ann ,labRow=rowNames,breaks=seq(.2,.8,length=100))

seur=backup

seur<-SetAllIdent(seur,"Clusts_nn_louvain_100_PC_13")


print("Compare to Pollen!")
COR=getCOR(seur,Pollen)

aheatmap(t(COR),Rowv=seq(ncol(COR),1),filename=paste(outdir,"Pollen_clust.pdf",sep="/"))

print("Compare to new")
load("devo_keep.Robj")
COR=getCOR(seur,devo_keep)
load("list_new.Robj")
tab=c()
type=c("Excitatory","Newborn Excitatory","Inhibitory","Newborn Inhibitory","IPC","Astrocytes","Radial Glia")
for(i in 1:length(lst)){tab=c(tab,rep(type[i],length(lst[[i]])))}
tab=rev(tab)
tab=data.frame(tab)
colnames(tab)="Cell Type"
COR=COR[rev(unlist(lst)),]
save(COR,file="temp1.Robj")
save(tab,file="temp2.Robj")
aheatmap(t(COR),Rowv=seq(ncol(COR),1),Colv=seq(nrow(COR),1),filename=paste(outdir,"new_Pollen_clust.pdf",sep="/") ,annCol=tab,width=14)











print("Compare to Dronc")

dronc=SubsetData(dronc,WhichCells(dronc,c("Astrocytes","Excitatory","Inhibitory")))
COR=getCOR(seur,dronc)

aheatmap(t(COR),Rowv=seq(ncol(COR),1),filename=paste(outdir,"Dronc_clust.pdf",sep="/"))

print("Compare to Organoids!")
org11a<-SetAllIdent(org11a,"tree.ident")
COR=getCOR(seur,org11a)

aheatmap(t(COR),Rowv=seq(ncol(COR),1),filename=paste(outdir,"org11a_clust.pdf",sep="/"))


print("Compare to Organoids, Cell Type Names!")
COR=getCOR(seur,sub_org)
print(max(COR))
print(min(COR))
aheatmap(t(COR),Rowv=seq(ncol(COR),1),filename=paste(outdir,"suborg11a_clust.pdf",sep="/"),breaks=seq(.2,.8,length=100))


print("Compare to Organoids, Forebrain!")
COR=getCOR(seur,fore)
print(max(COR))
print(min(COR))
aheatmap(t(COR),Rowv=seq(ncol(COR),1),filename=paste(outdir,"forebrain_clust.pdf",sep="/"),breaks=seq(.2,.8,length=100))


print("Generate TSNEs")
pdf(paste(outdir,"TSNE_CellType_clust.pdf",sep="/"))
TSNEPlot(seur,T)
dev.off();


print("Generate Gene by Cluster Heatmap")

load("markers.Robj")

print("Generate Feature Plots")
tab=tab[tab[,3]%in%rownames(seur@data),]


tab=tab[!duplicated(tab[,2]),]
inter=tab[,3]



ave=AverageExpression(seur,genes.use=inter)
rownames(ave)=sub("hg19_","",rownames(ave))
tab=data.frame(tab[,1])

mn=apply(ave,1,max)

ave=ave[mn>4,]
tab=tab[mn>4,]
tab=data.frame(tab)


ave=data.frame(ave)[order(tab[,1],decreasing=T),]
tab=tab[order(tab[,1],decreasing=T),]

tab=data.frame(tab)


colnames(tab)="Cell Type"
tab["Cell Type"]=as.character(tab[,1])

colnames(ave)=sub("X","",colnames(ave))

filename=paste(outdir,"genes_celltype_clust.pdf",sep="/")
aheatmap(ave,Rowv=seq(nrow(ave),1),annRow=tab,filename=filename)






}

if(!interactive())
{
load("lane34.Robj")
genFigures(lane34,outdir="outdir_newest2")
}



