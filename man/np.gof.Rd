\name{np.gof}
\alias{np.gof}

\title{
Goodness-of-Fit tests in nonparametric regression models
}
\description{
This routine tests the equality of a nonparametric regression curve, \eqn{m}, and a given function, \eqn{m_0}, from a sample \eqn{{(Y_i, t_i): i=1,...,n}}, where:
\deqn{Y_i= m(t_i) + \epsilon_i.}
The unknown function \eqn{m} is smooth, fixed equally spaced design is considered, and the random errors, \eqn{{\epsilon_i}}, are allowed to be time series. The test statistic used for testing the null hypothesis, \eqn{H0: m = m_0}, derives from a Cramer-von-Mises-type functional distance between a nonparametric estimator of \eqn{m} and \eqn{m_0}.
}
\usage{
np.gof(data = data, m0 = NULL, h.seq = NULL, w = NULL, 
estimator = "NW", kernel = "quadratic", time.series = FALSE, 
Tau.eps = NULL, h0 = NULL, lag.max = 50, p.max = 3, 
q.max = 3, ic = "BIC", num.lb = 10, alpha = 0.05)
}
\arguments{
 \item{data}{
\code{data[, 1]} contains the values of the response variable, \eqn{Y};

\code{data[, 2]} contains the values of the explanatory variable, \eqn{t}.
}
  \item{m0}{the considered function in the null hypothesis. If \code{NULL} (the default), the zero function is considered.}
 \item{h.seq}{the statistic test is performed using each bandwidth in the vector \code{h.seq}. If \code{NULL} (the default), 10 equidistant values between zero and a quarter of the range of \eqn{{t_i}} are considered.}
  \item{w}{support interval of the weigth function in the test statistic. If \code{NULL} (the default), \eqn{(q_{0.1}, q_{0.9})} is considered, where \eqn{q_p} denotes the quantile of order \eqn{p} of \eqn{{t_i}}.}
  \item{estimator}{allows us the choice between \dQuote{NW} (Nadaraya-Watson) or \dQuote{LLP} (Local Linear Polynomial). The default is \dQuote{NW}.}
  \item{kernel}{allows us the choice between \dQuote{gaussian}, \dQuote{quadratic} (Epanechnikov kernel), \dQuote{triweight} or \dQuote{uniform} kernel. The default is \dQuote{quadratic}.}
  \item{time.series}{it denotes whether the data are independent (FALSE) or if data is a time series (TRUE). The default is FALSE.}
  \item{Tau.eps}{it contains the sum of autocovariances associated to the random errors of the regression model. If NULL (the default), the function tries to estimate it: it fits an ARMA model (selected according to an information criterium) to the residuals from the fitted nonparametric regression model and, then, it obtains the sum of the autocovariances of such ARMA model.}
  \item{h0}{if \code{Tau.eps=NULL}, \code{h0} contains the pilot bandwidth used for obtaining the residuals to construct the default for \code{Tau.eps}. If \code{NULL} (the default), a quarter of the range of \eqn{{t_i}} is considered.}
  \item{lag.max}{if \code{Tau.eps=NULL}, \code{lag.max} contains the maximum delay used to construct the default for \code{Tau.eps}. The default is 50.}
  \item{p.max}{if \code{Tau.eps=NULL}, the ARMA model is selected between the models ARMA(p,q) with 0<=p<=\code{p.max} and 0<=q<=\code{q.max}. The default is 3.}
  \item{q.max}{if \code{Tau.eps=NULL}, the ARMA model is selected between the models ARMA(p,q) with 0<=p<=\code{p.max} and 0<=q<=\code{q.max}. The default is 3.}
  \item{ic}{if \code{Tau.eps=NULL}, \code{ic} contains the information criterion used to suggest the ARMA model. It allows us to choose between: "AIC", "AICC" or "BIC" (the default).}
  \item{num.lb}{if \code{Tau.eps=NULL}, it checks the suitability of the selected ARMA model according to the Ljung-Box test and the t-test. It uses up to \code{num.lb} delays in the Ljung-Box test. The default is 10.}
  \item{alpha}{if \code{Tau.eps=NULL}, \code{alpha} contains the significance level which the ARMA model is checked. The default is 0.05.}
}
\details{
A weight function (specifically, the indicator function \bold{1}\eqn{_{[w[1] , w[2]]}}) is introduced in the test statistic to allow elimination (or at least significant reduction) of boundary effects from the estimate of \eqn{m(t_i)}.

If \code{Tau.eps=NULL} and the routine is not able to suggest an approximation for \code{Tau.eps}, it warns the user with a message saying that the model could be not appropriate and then it shows the results. In order to construct \code{Tau.eps}, the procedures suggested in Muller and Stadmuller (1988) and Herrmann \emph{et al.} (1992) can be followed.

The implemented statistic test particularizes that one in Gonzalez Manteiga and Vilar Fernandez (1995) to the case where the considered class in the null hypothesis has only one element.
}
\value{A list with a dataframe containing:
  \item{h.seq}{sequence of bandwidths used in the test statistic.}
  \item{Q.m}{values of the test statistic (one for each bandwidth in \code{h.seq}).}
  \item{Q.m.normalised}{normalised value of Q.m.}
  \item{p.value}{p-values of the corresponding statistic tests (one for each bandwidth in \code{h.seq}).}
  Moreover, if \code{data} is a time series and \code{Tau.eps} is not especified:
  \item{pv.Box.test}{p-values of the Ljung-Box test for the model fitted to the residuals.}
  \item{pv.t.test}{p-values of the t.test for the model fitted to the residuals.}
  \item{ar.ma}{ARMA orders for the model fitted to the residuals.}
}
\references{
Biedermann, S. and Dette, H. (2000) Testing linearity of regression models with dependent errors by kernel based methods. \emph{Test} \bold{9}, 417-438.

Gonzalez-Manteiga, W. and Aneiros-Perez, G. (2003) Testing in partial linear regression models with dependent errors. \emph{J. Nonparametr. Statist.} \bold{15}, 93-111.

Gonzalez-Manteiga, W. and Cao, R. (1993) Testing the hypothesis of a general linear model using nonparametric regression estimation. \emph{Test} \bold{2}, 161-188.

Gonzalez Manteiga, W. and Vilar Fernandez, J. M. (1995) Testing linear regression models using non-parametric regression estimators when errors are non-independent. \emph{Comput. Statist. Data Anal.} \bold{20}, 521-541.

Herrmann, E., Gasser, T. and Kneip, A. (1992) Choice of bandwidth for kernel regression when residuals are correlated. \emph{Biometrika} \bold{79}, 783-795

Muller, H.G. and Stadmuller, U. (1988) Detecting dependencies in smooth regression models. \emph{Biometrika} \bold{75}, 639-650
}
\author{German Aneiros Perez \email{ganeiros@udc.es}

Ana Lopez Cheda \email{ana.lopez.cheda@udc.es}}
\seealso{
Other related functions are \code{\link{np.est}}, \code{\link{par.gof}} and \code{\link{plrm.gof}}.
}
\examples{
# EXAMPLE 1: REAL DATA
data <- matrix(10,120,2)
data(barnacles1)
barnacles1 <- as.matrix(barnacles1)
data[,1] <- barnacles1[,1]
data <- diff(data, 12)
data[,2] <- 1:nrow(data)

np.gof(data)



# EXAMPLE 2: SIMULATED DATA
## Example 2a: dependent data

set.seed(1234)
# We generate the data
n <- 100
t <- ((1:n)-0.5)/n
m <- function(t) {0.25*t*(1-t)}
f <- m(t)
f.function <- function(u) {0.25*u*(1-u)}

epsilon <- arima.sim(list(order = c(1,0,0), ar=0.7), sd = 0.01, n = n)
y <-  f + epsilon
data <- cbind(y,t)

## Example 2a.1: true null hypothesis
np.gof(data, m0=f.function, time.series=TRUE)

## Example 2a.2: false null hypothesis
np.gof(data, time.series=TRUE) 

}

\keyword{Statistical Inference}
\keyword{Regression}
\keyword{Time Series}
\keyword{Nonparametric Statistics}

