% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MethyLiution.R
\name{removeSnpProbesAndImprintingGenes}
\alias{removeSnpProbesAndImprintingGenes}
\title{Step 2: Filter SNP probes and imprinting genes}
\usage{
removeSnpProbesAndImprintingGenes(begin, mset_file, qc_probe_plot_file,
  logfile, array, col.name = NULL, mset.lumi = NULL)
}
\arguments{
\item{begin}{Step the pipeline begins at}

\item{mset_file}{Path to file in which dataset object is stored between steps (auto-filled by runPipeline() function)}

\item{qc_probe_plot_file}{Path to file where probe QC information will be written (auto-filled by runPipeline() function)}

\item{logfile}{Path to log file (auto-filled by runPipeline() function)}

\item{array}{Microarray probeset used: valid options are "450K" and "EPIC" (default "450K")}

\item{col.name}{List of column names, generated by setColumnNames() (auto-filled by runPipeline() function)}

\item{mset.lumi}{Microarray dataset object (auto-filled by runPipeline() function)}
}
\value{
mset.lumi dataset, with SNP probes and imprinting genes filtered out
}
\description{
This probably isn't something you want to run directly. Just use the runPipeline() function.
}
\keyword{internal}
