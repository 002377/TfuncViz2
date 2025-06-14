% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TfuncViz.R
\name{clean_genes}
\alias{clean_genes}
\title{Clean gene list by removing genes not in Seurat object}
\usage{
clean_genes(gene_list, seurat_obj = NULL)
}
\arguments{
\item{gene_list}{A vector of gene names}

\item{seurat_obj}{Seurat object}
}
\value{
A vector of genes present in the Seurat object
}
\description{
Clean gene list by removing genes not in Seurat object
}
