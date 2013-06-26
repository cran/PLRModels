
###########################################################################
# It fits an optimal ARMA (according to an information criteria) to the data (an only series, that is, a vector).
# It analyses the residual and it returns the variances-covariances matrix of such ARMA model.
# If the ARMA model is not suitable, it returns "error".
###########################################################################

var.cov.matrix <- function(x=1:100, n=4, p.max=3, q.max=3, ic="BIC", alpha=0.05, num.lb=10) {
# n x n is the dimension of the output matrix


p.q <- best.arima(x=x, order.max=c(p.max,0,q.max), include.mean=FALSE, criterio=ic)[1,]

fitted.model <- arima(x=x, order=c(p.q[1,1],0,p.q[1,2]), include.mean=FALSE)

fitdf <- sum(fitted.model$arma[1:2])


pv.lb.t <- c(rep(0,num.lb+1))


for (i in 1:num.lb)
pv.lb.t[i] <- Box.test(x=residuals(fitted.model), lag = fitdf+i, type = "Ljung-Box", fitdf = fitdf)$p.value


pv.lb.t[num.lb + 1] <- t.test(residuals(fitted.model), mu=0)$p.value




if (min(pv.lb.t)<alpha) {
	cat("The ", ic, " model does not overcome the residual analysis", "\n")
	return("error")
		}

else {
#if ((length(fitted.model$model$phi)+length(fitted.model$model$theta))==0) Var.Cov.x <- diag(1,n)
if (sum(fitted.model$arma[1:2])==0) Var.Cov.x <- diag(1,n)

	else
	Var.Cov.x <- toeplitz(ARMAacf(ar=fitted.model$model$phi, ma=fitted.model$model$theta, lag.max=n-1))
  
v.x <- var(x)
v.x <- as.numeric(v.x)

Var.Cov.x <- v.x * Var.Cov.x                 
return (Var.Cov.x)                                                        
}
}