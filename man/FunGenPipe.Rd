% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FunGenPipe.R
\docType{package}
\name{FunGenPipe}
\alias{FunGenPipe}
\title{Functional genomics functions and pipes}
\description{
Tools for building function genomics pipelines
}
\section{Genomics}{

Functions that operate with ranges and sequences
dechiperMeltDNAtime
}

\section{Tagets meta-data}{

Functions for manipulating colors in targets.
getColPats - depricated
getColPats2
}

\section{Plot distributions, heatmap, PCA & MDS using targets meta-data}{

Note: the order od data columns should (always) be kept in the order of targets.
plotDistrTargets2
pheatmapTargets2
plotPCAtargets2
plotPCAtargets - depricated
plotMDStargets2
plotMDStargets2_MASS - depricated
}

\section{Plot distributions using limma-style data}{

plotDistr_RGList
boxplotRLE
}

\section{From KBlag_doxy_Cla.R For EKocar_fibIT.R}{

writeExprsGEO
}

\section{From Osijek_HFHSD For KBlag_doxy_Cla.R}{

esetAnnSEs
getWriteHeatmap_probesTTcontrasts
getWriteHeatmap_probesTTcomparisons
getWriteHeatmap_PgseaTTcontrasts
getWriteHeatmap_PgseaTTcomparisons
annKeggEntrez
gscKeggEntrez
}

\author{
Peter Juvan, \email{peter.juvan@gmail.com}
}