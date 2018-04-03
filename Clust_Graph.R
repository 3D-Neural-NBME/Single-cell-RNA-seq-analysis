##
##Modified version of the code from:
##https://github.com/broadinstitute/BipolarCell2016
##
##All credit to the Authors of that package (K Shekhar)
##Also included is the license they published under
##

library(Matrix)
library(reshape)
library(igraph)
library(RANN)


##Feed in a matrix, X, with rows corresponding to cells and columns to PCs
##the number of nearest neighbors to use (nn), and specify the type (louvain or infomap)
##Returns the associated clustering.
clust_Graph<-function(X,nn,getGraph=FALSE,type="louvain")
{

nearest=nn2(X,X,k=nn+1, treetype = "bd", searchtype="priority")
print("Found nearest neighbors")
nearest$nn.idx = nearest$nn.idx[,-1]
nearest$nn.dists = nearest$nn.dists[,-1] #Convert to a similarity score
nearest$nn.sim = 1*(nearest$nn.dists >= 0 )
    
    
edges = melt(t(nearest$nn.idx)); colnames(edges) = c("B", "A", "C"); edges = edges[,c("A","B","C")]
edges$B = edges$C; edges$C=1
    

edges = unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))

        
NN = nearest$nn.idx
jaccard_dist = apply(edges, 1, function(x) length(intersect(NN[x[1], ],NN[x[2], ]))/length(union(NN[x[1], ], NN[x[2], ])) )
        
edges$C = jaccard_dist
edges = subset(edges, C != 0)
edges$C = edges$C/max(edges$C)


colnames(edges) = c("First_Node", "Second_Node", "Weight")

print("Make Sparse Adj!")

adj_sparse<-sparseMatrix(i=edges$First_Node,j=edges$Second_Node,x=edges$Weight,symmetric=TRUE)

print("Make graph!")

g=graph_from_adjacency_matrix(adj_sparse,mode="undirected",weighted=TRUE)

if(getGraph){return(g)}

print("cluster!")
if(type=="louvain")
{
clust = cluster_louvain(g)
}
else
{

clust = cluster_infomap(g)
}
return(clust)

}


##Same as above, except takes in a seurat object, and number of PCs (numPC) to use, and returns a seurat object with the clustering added to the data.info dataframe and in the ident variable.
clustGraph_Seurat<-function(seur,nn=50,numPC=7,type="louvain",lstPC=c())
{

    X=seur@pca.rot[1:numPC]
    d=seur@pca.obj[[1]]$d[1:numPC]
    for(i in 1:numPC){X[i]=d[i]*X[i]}
    if(length(lstPC)>0){X=X[,lstPC]}
    clust<-clust_Graph(X,nn,type=type)
    seur@data.info[paste("Clusts_nn",type,toString(nn),"PC",toString(numPC),sep="_")]=clust$membership;
    seur<-SetAllIdent(seur,paste("Clusts_nn",type,toString(nn),"PC",toString(numPC),sep="_"))
    return(seur);
    
}
