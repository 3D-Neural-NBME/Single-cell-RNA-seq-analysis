# Single-cell-RNA-seq-analysis

This directory contains code used for the analysis of the single cell data in our manuscript. For a run through of the analysis (from expression matrix to clustering to plotting) see Analysis.Rmd and Analysis.html. Note that our analysis requires an earlier versions of Seurat (prior to version 2, we used 1.3 to test). All analysis prior to the code provided was performed with CellRanger, as specified in the manuscript.

The code makeFigs.R was the code used to generate the figures in the paper (in addition to other figures). Note that it relies on data not included in this github repositories (due to space restrictions), but can be downloaded from the respective publications. We do include the data from our publication (as a sparse matrix of log counts per million) and the data from Pollen et al 2015 as a Seurat object.

In addition to the analysis data, we included some scripts of general use for Single Cell analysis:
1) TSNE.R: Performs TSNE (corrects for some issues with how older versions of Seurat normalize the pca data)
2) Clust_Graph.R: Contains the code for clustering (louvain or infomap). Note this is modified version of the code used in Karthik et al 2016.
3) niceDraw.R: Functions to draw feature plots, etc, in a (hopefully) pleasing way.
