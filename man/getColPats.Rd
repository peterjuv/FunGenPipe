% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FunGenPipe.R
\name{getColPats}
\alias{getColPats}
\title{Depricated: Get color patterns with names from another column from ExpressionSet slot phenoData.}
\usage{
getColPats(eset)
}
\arguments{
\item{eset}{ExpressionSet}
}
\value{
List of colors, each of length equal to the number of rows in phenoData
}
\description{
Returns columns from phenoData that start with "color_" and have their suffix in common with another column;
Names of colors are set from the associated column, e.g., \code{list(color_Sex = stats::setNames(color_Sex, Sex), ...)}
}
\section{Old implementation}{

\code{
getColPats <- function(eset) {
    assertthat::has_attr(eset, "phenoData")
    colPats <- list()
    aNames <- colnames(Biobase::pData(eset))
    cNames <- aNames[grep("^color_", aNames, perl=TRUE)]
    fNames <- sub("^color_", "", cNames)
    for (i in seq_len(length(fNames))) {
        colPats[[fNames[i]]] <- stats::setNames(Biobase::pData(eset)[[cNames[i]]], Biobase::pData(eset)[[fNames[i]]])
    }
    colPats
}
}
}

\examples{
\dontrun{
getColPats(Biobase::pData(dataRaw/eset))
getColPats(Biobase::pData(eset))
}
}
