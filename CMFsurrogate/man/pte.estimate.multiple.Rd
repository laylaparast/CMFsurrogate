\name{pte.estimate.multiple}
\alias{pte.estimate.multiple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates the proportion of treatment effect explained by the optimal combination of multiple surrogate markers using  a calibrated model fusion approach
}
\description{
Estimates the proportion of treatment effect explained by the optimal combination of multiple surrogate markers using  a calibrated model fusion approach}
\usage{
pte.estimate.multiple(sob, yob, aob, var = TRUE, rep = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sob}{surrogates
}
  \item{yob}{
primary outcome, y
}
  \item{aob}{
treatment indicator
}
  \item{var}{
TRUE or FALSE, if variance/SE of PTE is being requested
}
  \item{rep}{
if var is TRUE, number of resampled draws to use for bootstrap
}
}
\value{
  \item{pte.es }{Estimate of the proportion of treatment effect explained (PTE)}
  \item{pte.se }{if var = TRUE, estimate of the standard error of the PTE}

}
\references{
Wang X, Parast L, Tian L, Han L, Cai T. Robust Approach to Combining Multiple Markers to Improve Surrogacy. (Under Review)
}
\examples{
data(example.data)
out=pte.estimate.multiple(sob=example.data$sob, yob=example.data$yob, aob=example.data$aob, var = FALSE)
out
}