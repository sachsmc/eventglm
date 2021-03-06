% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{colon}
\alias{colon}
\title{Chemotherapy for Stage B/C colon cancer}
\format{
A data frame with 929 rows and 17 variables:
\describe{

\item{id}{id}
\item{study}{1 for all patients}
\item{rx}{Treatment - Obs(ervation), Lev(amisole), Lev(amisole)+5-FU}
\item{sex}{1=male}
\item{age}{in years}
\item{obstruct}{obstruction of colon by tumour}
\item{perfor}{perforation of colon}
\item{adhere}{adherence to nearby organs}
\item{nodes}{number of lymph nodes with detectable cancer}
\item{time}{days until death or censoring}
\item{status}{censoring status}
\item{differ}{differentiation of tumour (1=well, 2=moderate, 3=poor)}
\item{extent}{Extent of local spread (1=submucosa, 2=muscle, 3=serosa, 4=contiguous structures)}
\item{surg}{time from surgery to registration (0=short, 1=long)}
\item{node4}{more than 4 positive lymph nodes}
\item{etype}{event type: 1=recurrence,2=death}
\item{event}{event indicator: censored, death}
}
}
\usage{
colon
}
\description{
These are data from one of the first successful trials of adjuvant
chemotherapy for colon cancer. Levamisole is a low-toxicity compound
previously used to treat worm infestations in animals; 5-FU is a moderately
toxic (as these things go) chemotherapy agent. There are only one record per patient for the death outcome (or censoring). This is redistributed from
the survival package, with a small modification to include only the death outcome.
}
\seealso{
\link[survival]{colon}
}
\keyword{datasets}
