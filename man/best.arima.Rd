\name{best.arima}
\alias{best.arima}
\title{Best Arima model according some information criterion
}
\description{Obtains the orders p, q, P and Q of the best ARIMA (p, d, q)x(P, D, Q)_s according one of the following information criteria: AIC, AICC, BIC.
It allows ARIMAs with both constant term and differentiation (d+D!=0), besides the ordinary cases.
}
\usage{
best.arima(x = x, order.max = c(0, 0, 0), seasonal = list(order.max = c(0, 0, 0), 
period = 1), include.mean = NULL, criterio = NULL, dist.max.crit = NULL, 
method = NULL)
}
\keyword{internal}

