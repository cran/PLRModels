
plrm.gcv <- function(data=data, b.equal.h=TRUE, num.b=NULL, num.h=NULL, 
                     estimator="NW", kernel="quadratic")
{

if (!is.matrix(data))  stop("data must be a matrix")
if (ncol(data)<3)  stop("data must have at least 3 columns: y, x, t")

if (!is.logical(b.equal.h)) stop("b.equal.h must be logical")

if ( (!is.null(num.b)) && (length(num.b) !=1) ) stop ("num.b must be an only value") 
if ( (!is.null(num.b)) && (!is.numeric(num.b)) )  stop ("num.b must be numeric") 
if ( (!is.null(num.b)) && (num.b<=0) ) stop ("num.b must be a positive value") 

if ( (!is.null(num.h)) && (length(num.h) !=1) ) stop ("num.h must be an only value") 
if ( (!is.null(num.h)) && (!is.numeric(num.h)) )  stop ("num.h must be numeric") 
if ( (!is.null(num.h)) && (num.h<=0) ) stop ("num.h must be a positive value") 

if ((estimator != "NW") & (estimator != "LLP"))  stop("estimator=NW or estimator=LLP is required")

if ((kernel != "quadratic") & (kernel != "Epanechnikov") & (kernel != "triweight") & (kernel != "gaussian") & (kernel != "uniform"))  stop("kernel must be one of the following: quadratic, Epanechnikov, triweight, gaussian or uniform")



kernel.function <- get(kernel)       
n <- nrow(data)
p <- ncol(data)-2
x <- data[, 2:(1+p)]
y <- data[, 1]
t <- data[, p+2]

if (!is.matrix(x))  x <- as.matrix(x)


if ((b.equal.h==TRUE) & ((!is.null(num.b)) | (!is.null(num.h))) ) {num.b <- max(num.b, num.h); num.h <- num.b}
else if ((b.equal.h==TRUE) & (is.null(num.b)) & (is.null(num.h)) ) {num.b <- 50; num.h <- num.b}
else if ((b.equal.h==FALSE) & ((!is.null(num.b)) & (is.null(num.h))) ) num.h <- num.b
else if ((b.equal.h==FALSE) & ((is.null(num.b)) & (!is.null(num.h))) ) num.b <- num.h
else if ((b.equal.h==FALSE) & (is.null(num.b)) & (is.null(num.h)) ) {num.b <- 50; num.h <- num.b}


if (b.equal.h) GCV <- matrix(0,num.b, 1)
else GCV <- matrix(0,num.b, num.h)

b.h.GCV <- data.frame(cbind(b=0,h=0))

  
  
a <- as.matrix(abs(outer(t, t,"-")))
for (i in 1:n) {a[i,i] <- -1000}
a <- as.vector(a[a!=-1000])
  
b.min <- quantile(a,0.05)
b.max <- (max(t)-min(t))*0.25
  
b.seq <- seq(b.min,b.max,length.out=num.b)
h.seq <- seq(b.min,b.max,length.out=num.h)
    


W.g <- function(t=t, g=NULL, estimator=estimator, kernel.function=kernel.function)
{

if (estimator=="NW") {
  
  Zmat <- outer(t, t, "-")
  Umat <- Zmat/g

  Kmat <- kernel.function(Umat)
  Kmat[(Kmat<0)] <- 0 

  S0 <- apply(Kmat, 1, sum)
  Kmat <- Kmat/S0

}      
  
else if (estimator=="LLP")
        
  Zmat <- outer(t, t, "-")
  Umat <- Zmat/g

  Kmat <- kernel.function(Umat)
  Kmat[(Kmat<0)] <- 0

  S0 <- apply(Kmat, 1, sum)
  S1 <- apply(Kmat*Zmat, 1, sum)
  S2 <- apply(Kmat*Zmat^2, 1, sum)

  Kmat <- Kmat * (S2 - Zmat*S1)/(S0*S2 - S1^2)

}



if (b.equal.h) {
  
  for (i in 1:num.b) {
    
    W.b <- W.g(t=t, g=b.seq[i], estimator=estimator, kernel.function=kernel.function)
    
    XX.b <- (diag(n) - W.b) %*% x
    
    XI <- solve(t(XX.b)%*%XX.b) %*% t(XX.b) %*% (diag(n) - W.b)
    
 
    A.bh <- W.b + XX.b %*% XI
      
    RSS.bh <- sum(((diag(n) - A.bh) %*% y)^2)/n
      
    GCV[i,] <- RSS.bh / (1 - sum(diag(A.bh))/n )^2
      
  } # for i
} # if

else {
  
  for (i in 1:num.b) {
  
    W.b <- W.g(t=t, g=b.seq[i], estimator=estimator, kernel.function=kernel.function)

    XX.b <- (diag(n) - W.b) %*% x
  
    XI <- solve(t(XX.b)%*%XX.b) %*% t(XX.b) %*% (diag(n) - W.b)
  
  
    for (j in 1:num.h) {

      W.h <- W.g(t=t, g=h.seq[j], estimator=estimator, kernel.function=kernel.function)
    
      XX.h <- (diag(n) - W.h) %*% x
    

      A.bh <- W.h + XX.h %*% XI
  
      RSS.bh <- sum(((diag(n) - A.bh) %*% y)^2)/n
    
      GCV[i,j] <- RSS.bh / (1 - sum(diag(A.bh))/n )^2

    } # for i
  } # for j
} # else



index.GCV <- order(GCV)[1]
  
  if (b.equal.h==FALSE) {
      index.h <- 1+trunc((index.GCV-1)/num.b)
      index.b <- index.GCV - (index.h-1)*num.b
      b.h.GCV <- c(b.seq[index.b], h.seq[index.h])
      GCV.opt<- GCV[index.b, index.h]
  }
    
  
  else if (b.equal.h==TRUE) {
      b.h.GCV <- c(b.seq[index.GCV], h.seq[index.GCV])
      GCV.opt <- GCV[index.GCV]
  }


list(bh.opt=b.h.GCV, GCV.opt=GCV.opt, GCV=GCV, b.seq=b.seq, h.seq=h.seq)

}

