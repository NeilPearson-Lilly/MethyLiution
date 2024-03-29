% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MethyLiution.R
\name{checkReportedSex}
\alias{checkReportedSex}
\title{Step 3: Check that metadata for samples matches inferences}
\usage{
checkReportedSex(begin, mset_file, logfile, gender_stringent,
  sex_label_cluster_file, outdir, col.name = NULL, mset.lumi = NULL)
}
\arguments{
\item{begin}{Step the pipeline begins at}

\item{mset_file}{Path to file in which dataset object is stored between steps (auto-filled by runPipeline() function)}

\item{logfile}{Path to log file (auto-filled by runPipeline() function)}

\item{gender_stringent}{Should the pipeline discard samples where inferred gender is inconsistent with metadata?}

\item{sex_label_cluster_file}{Path to file where probe QC information will be written (auto-filled by runPipeline() function)}

\item{outdir}{Path to output directory}

\item{col.name}{List of column names, generated by setColumnNames() (auto-filled by runPipeline() function)}

\item{mset.lumi}{Microarray dataset object (auto-filled by runPipeline() function)}
}
\value{
mset.lumi dataset
}
\description{
This probably isn't something you want to run directly. Just use the runPipeline() function.
}
\keyword{internal}
