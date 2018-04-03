library(Seurat)

##Given two seurat objects, takes the average of the clusters in both 
##and looks at the correlation using the genes specified by the genes argument.

getCOR<-function(seur1,seur2,useAve=T,genes=c())
{
lst=c()
for(i in seur2@var.genes){lst<-c(lst,paste("hg19_",i,sep=""))}
if(length(genes)>0)
{
lst=c()
for(i in genes){lst<-c(lst,paste("hg19_",i,sep=""))}
}
inter=intersect(lst,rownames(seur1@data))
print(length(inter))



ave=AverageExpression(seur2,genes.use=sub("hg19_","",inter))
ave2=AverageExpression(seur1,genes.use=inter)
if(!useAve)
{
ave2=as.matrix(seur1@data[inter,])
}

rownames(ave)=sub("hg19_","",inter)

COR=cor(ave,ave2)

return(COR)

}
