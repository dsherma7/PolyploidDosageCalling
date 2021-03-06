% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LinearModel.R
\name{GetDosage}
\alias{GetDosage}
\title{Linear Regression on PCR Variables}
\usage{
GetDosage(AllData, index = 1:length(AllData), Controls = NULL,
  Target = "Assay1")
}
\arguments{
\item{AllData}{A serial data frame of all the relevant values including Job, Plate, RQ, CT, delta CT, and Sample Name.}

\item{Controls}{data.frame containing the Control values for this experiment. Must have column names of Actual.Name, Name, Dosage; where Actual.Name is what is found exactly in the data sets, Name is the desired output name (often the same as Actual.Name), and Dosage is the Dosage for each control.}

\item{Results.Only}{FALSE means that the AllData is actually a three-dimensonal list with the third dimension being each sheet in the data sets. If AllData only contains the "Results' section then set this to TRUE.}

\item{Index}{The range of the rows to be used in the fit of the AllData data frame.}
}
\description{
Uses a regression technique to model log(Dosage) as a function of delta CT values.
}
