% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FunGenPipe.R
\name{annTFbyFA}
\alias{annTFbyFA}
\title{Depricated: Get annotations for PGSEA from TRANSFAC using factor names (factor.FA)}
\usage{
annTFbyFA(fac2drDF, eset = NULL, patternExcludeFactorFA = NULL)
}
\arguments{
\item{fac2drDF}{Data frame transcription factor}

\item{eset}{Expression set with anotations ENTREZID within fData(); these are used to limit TFs to ENTREZIDs; default NULL}

\item{patternExcludeFactorFA}{Regular expression pattern for exclusion of TFs based on their names (i.e., factor.FA); default NULL;}
}
\value{
A data frame with annotations and collapsed gene Symbols and ENTREZIDs
}
\description{
Import TRANSFAC data.frame fac2drDF and filter using ENTREZIDs (may be collapsed by ";") from an ExpressionSet.
Data frame includes columns: factor.AC site.AC factor.BSq gene.AC gene.DRac2 gene.DRdbName2 factorURL factor.FA gene.SD gene.DE
TFs may be excluded based on a pattern matching their names, for example:
\code{"[^(hsa)]-miR"} will exclude non-human microRNAs
\code{"(-miR)|(isoform)"} will exclude all microRNAs and isoforms
Old TODO: consider creating new factors for fac2drDF_entrez_sym or even more general, e.g. using fac2drDF or  or fac2drDF_osXX
}
