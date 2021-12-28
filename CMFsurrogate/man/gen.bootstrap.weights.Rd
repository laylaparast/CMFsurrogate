\name{gen.bootstrap.weights}
\alias{gen.bootstrap.weights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate bootstrap sample
}
\description{
Generate bootstrap sample
}
\usage{
gen.bootstrap.weights(n, num.perturb = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
sample size
}
  \item{num.perturb}{
number of replicates/resamples
}
}
\value{
matrix with n rows and num.perturb columns of indeces
}
