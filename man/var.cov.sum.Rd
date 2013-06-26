\name{var.cov.sum}

\title{Estimated sum of autocovariances from time series
}
\description{For each time series in the columns of the data matrix, \code{X}, selects an optimal ARMA model (according to an information criteria); then, fits such model and analyses the corresponding residuals. If all the ARMA models are suitable, returns a vector containing the corresponding sums the autocovariances. If some ARMA model is not suitable, returns the word "error".
}
\usage{
var.cov.sum(X = 1:100, lag.max = 50, p.max = 3, q.max = 3, ic = "BIC", 
alpha = 0.05, num.lb = 10)
}

\keyword{internal}

