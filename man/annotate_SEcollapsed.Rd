% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FunGenPipe.R
\name{annotate_SEcollapsed}
\alias{annotate_SEcollapsed}
\title{Add collapsed annotations SYMBOLSs, ENTREZIDs (;-separated) to eset using an AnnotationDbi package}
\usage{
annotate_SEcollapsed(eset, annPackage)
}
\arguments{
\item{eset}{ExpressionSet}

\item{annPackage}{Character AnnotationDbi package name, e.g. "clariomsmousetranscriptcluster.db"}
}
\value{
ExpressionSet eset with collapsed annotations in slot fData(eset) and annotations(eset) set to DB name (annPackage)
}
\description{
Sets fData(eset) and annotation(eset); keeps probes w/o annotations.
Warning: SYMBOLSs and ENTREZIDs are unique and individually sorted; thus, they are not matches by their position within the collapsed strings.
}
