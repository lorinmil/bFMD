#' Analyze Balanced Modules
#'
#' Find subset of variables composing a balanced fuctional module within high dimensional data
#'
#'
#'
#' @param x n x p data matrix of variables containing a balanced module as a subset
#' @param y n X 1 vector of outcome measurements
#' @param maxSize maximum module size to be considered (to limit computation time)
#' @param conv the algorithm converges if the set of k variables selected remains the same for conv iterations
#' @param decline if the balance metric declines for decline straight choices of k > kmax, kmax variables are chosen
#'
#'
#' @return list of module variables, the module size, balance metric for the module, x, y
#'
#' @references Sigg, C. D. and Buhmann, J. M. (2008) Expectation-Maximization
#' for Sparse and Non-Negative PCA. In \emph{Proceedings of the 25th
#' International Conference on Machine Learning} (pp. 960-967).
#'
#' @export
#'
#' @examples
#' dat <- simDatHub()    # generate small test dataset
#' FMD(x=dat$x, y=dat$y)

library(irr)

FMD <- function(x, y, maxSize=ncol(x), maxiter=200, conv=5, decline=20) {
  x <- scale(x)
  y <- scale(y)
  maxBal <- -1000
  maxK <- 0
  counter <- 0
  covs <- as.vector(cor(x, y))
  w <- x %*% diag(covs)
  # idx <- 0
  # w <- apply(x, 2, function(x_col){
  #   idx <- idx + 1
  #   x_col*covs[idx]
  # })
  
  
  # w <- diag(as.vector(cor(x,y))) %*% cor(x) %*% diag(as.vector(cor(x,y)))
  
  
  for(j in 2:maxSize) {
    if(counter > decline) break
    counter <- counter + 1
    prout <- nnspc(Xp=w, k=j, conv=conv, maxiter=maxiter)$sel
    fmdout <- prout > 0
    bal <- balance(y, x, fmdout)
    
    if(bal >= maxBal) {
      counter <- 0
      maxBal <- bal
      maxK <- j
    }
    maxK <- max(maxK, 2)
  }
  prout_nnsp <- nnspc(Xp=w, k=maxK, conv=conv, maxiter=maxiter)
  prout <- prout_nnsp$sel
  fmdpos <- prout > 0
  fmdout <- list(select=fmdpos,
                 w=w,
                 loadings=prout_nnsp$w,
                 size=maxK, 
                 balance=balanceKappa(fmdpos, x=x, y=y)$kappa$value, 
                 x=x, 
                 y=y)
  #fmdout <- list(select=fmdpos,size=maxK, x=x, y=y)
  class(fmdout) <- "fmd"
  fmdout
}

balance <- function(y, xmat, chosen) {
  if(sum(chosen) < 2) return(0)
  cormat <- cor(xmat)[chosen, chosen]
  # cormat <- cor(xmat[,chosen])
  #  print(cormat)
  xycormat <-cor(xmat, y)[chosen]
  w <- diag(xycormat) %*% cormat %*% diag(xycormat)
  balout <- (sum(w) - sum(diag(w)))/sum(abs(w))
  balout
}

balancePlot <- function(x, y, select)  {
  selx <- x[, select]
  covxy <- cov(selx, y)
  ordxy <- order(covxy)
  ordmatrix <- selx[,ordxy]
  heatmap.2(cov(ordmatrix),Rowv=NA, Colv=NA, dendrogram="none",density.info = "none")
}

plot.fmd <- function(x) balancePlot(x$x, x$y, x$select)

balanceKappa <- function(x, y, select){
  xsel <- x[, select]
  ysel <- y
  covxy <- cov(xsel, y)
  covx <- cov(xsel)
  covsAgree <- NULL
  varsAgree <- NULL
  for (i in 1:(ncol(xsel)-1)) {
    for (j in (i+1):ncol(xsel)) {
      covsAgree <- c(covsAgree, as.numeric(covxy[i] * covxy[j] > 0))
      varsAgree <- c(varsAgree, as.numeric(covx[i,j] > 0))
    }
  }
  kappaout <- list(kappa=kappa2(data.frame(covsAgree, varsAgree)), covsAgree=sum(covsAgree)/length(covsAgree), varsAgree=sum(varsAgree)/length(varsAgree))
  kappaout
}

kappa.fmd <- function(x){
  balanceKappa(x=x$x, y=x$y, select=x$select)
}



nnspc <- function(Xp, k, conv, maxiter, verbose=TRUE) {
  d <- ncol(Xp); n <- nrow(Xp)
  #Randomly generate an initial loading vector
  w <- abs(rnorm(d))
  #w <- w/normv(d)
  #Scale the initialized loading vector
  w <- w/sqrt(sum(w^2))
  #Initialize everything as being selected
  prevsel <- rep(1,d)
  #Initialize the iteration counter
  ii <- 0
  #Initialize the convergence measure
  counter <- 0
  #Loop through until convergence (or maximum number of iterations)
  while ((ii < maxiter) & (counter < conv)){  
    #Multiply the data by the loading vector
    z <- Xp %*% w
    #So like t(Xp) * Xp * w ?
    w_star <- as.vector(t(Xp) %*% z)
    w_star[w_star < 0] <- 0
    if(k < d) {
      s <- sort(w_star, decreasing=TRUE)
      w <- pmax(w_star - s[k + 1], 0)* sign(w_star)
    } else {
      w <- w_star
    }
    w <- w/sqrt(sum(w^2))
    ii <- ii + 1
    sel <- as.integer(w > 0)
    if (sum(abs(sel - prevsel)) == 0) {
      counter= counter + 1
    } else{
      counter <- 0
      prevsel <- sel
    }
  }
  if(verbose > 0 & ii == maxiter ) print("convergence failed")
  #print(paste("k= ", k, "iterations= ", ii, "conv= ", counter))
  return(list(sel=sel,
              w=w))
}






