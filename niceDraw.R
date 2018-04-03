library(ggplot2)
library(cowplot)
library(Seurat)

###
##This function reads in a Seurat object and a gene name and creates a Featureplot of the gene.
##It reorders the cells so that highly expressed cells are drawn on top of lowly expressed cells (helpful when there is high drop out and lots of cells on top of one another)
##Other parameters:
##low: Specifies the color used for lowly expressed cells
##high: Specifies the color used for lowly expressed cells
##size: Specifies the size of the points in the feature plot
##datainfo: Specifies if the feature being plotted is stored in the meta data (in the data.info dataframe)
##usePCA: Specified if one wants to plot using pca instead of TSNE
##useScale: If true, uses data in seur@scale.data, else uses data in seur@data
##
##Returns:
##p: A ggplot2 plot of the feature plot of interest
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





