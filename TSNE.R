library(useful)
library(reshape)
library(Seurat)
#source("~/SingleCell/niceDraw.R")
print("Loaded!")
library(Rtsne)


getTSNE<-function(seur,numPC,perp=30,lstPC=c())
{
X=seur@pca.rot[1:numPC]
d=seur@pca.obj[[1]]$d[1:numPC]
#d=d/max(d)

for(i in 1:numPC){X[i]=d[i]*X[i]}
if(length(lstPC)>0){X=X[lstPC]}


tsn=Rtsne(X,verbose=T,pca=F,perplexity=perp)

tsn=data.frame(tsn$Y)
rownames(tsn)=rownames(X)
colnames(tsn)=c("tSNE_1","tSNE_2")

seur@tsne.rot=tsn

return(seur)

}


Read10X_nice<-function(lst)
{

dat=Read10X(lst)

cols=colnames(dat)
new_cols=c()
for(cur in cols){
c_new=cur
if(substr(cur,1,1)%in%c("A","T","G","C")){c_new=paste("1_",cur,sep="")}
new_cols<-c(new_cols,c_new)
}

colnames(dat)=new_cols

return(dat)
}

