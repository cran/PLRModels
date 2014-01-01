\name{var.cov.matrix}

\title{Estimated variance-covariance matrix from time series
}
\description{
Selects an optimal ARMA model (according to an information criteria) for the time series in the data vector, \code{x}; then, fits such model and analyses the corresponding residuals. If the ARMA model is suitable, returns the \code{n x n} variance-covariance matrix corresponding to \code{n} consecutive variables in the ARMA process. If the ARMA model is not suitable, it informs the user with a message.
}
\usage{
var.cov.matrix(x = 1:100, n = 4, p.max = 3, q.max = 3, ic = "BIC",  p.arima=NULL, 
q.arima=NULL, alpha = 0.05, num.lb = 10)
}

\keyword{internal}

