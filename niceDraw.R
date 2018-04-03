library(ggplot2)
library(cowplot)
library(Seurat)

niceFeaturePlot<-function(seur,gene_name,datainfo=F,low="gray",high="red",size=.5,usePCA=F,useScale=F)
{

dat=seur@tsne.rot
if(useScale)
{
seur<-RegressOut(seur,"nGene",genes.regress=c(gene_name,seur@var.genes[1:10]))
}

if(usePCA)
{
dat=seur@pca.rot
}
print(head(dat))
if(!datainfo)
{dat["Gene"]=as.numeric(seur@data[gene_name,])}
else{dat["Gene"]=as.numeric(seur@data.info[,gene_name])}
if(useScale)
{
print("hi")
print(corner(seur@scale.data))
dat["Gene"]=as.numeric(seur@scale.data[gene_name,])
}
print(head(dat))
p=c()
dat=dat[order(dat$Gene),]
if(usePCA)
{

p=ggplot(dat,aes(x=PC1,y=PC2,color=Gene))
}
else{
p=ggplot(dat,aes(x=tSNE_1,y=tSNE_2,color=Gene))
}
p=p+geom_point(size=size)+scale_colour_gradient(low=low,high=high)+ggtitle(gene_name)
return(p)

}





