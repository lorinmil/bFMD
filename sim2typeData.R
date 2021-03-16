#' Generate data set for a sequence of 2 functional modules
#'
#'
#' \code{simPathData} generates a data set representing a latent pathway through 2 blocks data
#' y <- x <- g where X and are latent factors for different types of data. Other structure is also generated for realism.
#' The two blocks are X and G, and x and g are corresponding latent variables of that type which compose the pathway. In addition other
#' latent variables are generated: s is a subset of X which is independent of G and functional for y; h is a subset of G
#' whish is independent of X and functional for y; g is functional for y through the pathway and also directly; x' and g'
#' are nonfunctional and correlated with each other. The variance of g and h is set to 1.0.
#' The variance of s is set to equal the variance of x. The equation error for y = b*g is 1.0. All other parameters for the latent model are fixed
#' by specifying alpha and Rsq. The observed indicators of the latent factors have correlation rind with the latent factor.
#'
#' @param alpha Relative path importance, 0 <= alpha <= 1. It is the ratio of the effect of g on y through the pathway to the
#' total effect of g
#' @param Rsq R^2 for G determining X and for X and G type variables determining Y. This describes the ammount of variance and covariance determined by the latent factors.
#' @param nsample Sample size
#' @param rind Correlation of X indicators with X latent factor, likewise for G indicators. It also describes the correlations amongst the indicators and hence the clustering of varaiables.
#' @param nx Number of X indicators in pathway
#' @param nxPrime number of nonfunctional X indicators correlated with G
#' @param ns Number of functional X variables not in pathway
#' @param ntranscriptNoise Number of random noise variables of type X
#' @param ng Number of G variables in pathway
#' @param ngPrime Number of nonfunctional G variables correlated with X.
#' @param nh Number of functional G variables not on pathway
#' @param ngeneNoise Number of random noise variables of type G
#'
#' @return{list  ordered as y, transcripts=c(x, x', s, independent noise),
#' genes=c(g, g', h, independent noise), true transcripts identifying pathway transcripts,
#' true genes identifying pathway genes, data generation parameters}
#'
#' @references Zhang F, Miecznikowski JC, Tritchler DL.
#' Identification of supervised and sparse functional genomic pathways.
#' \emph{Stat Appl Genet Mol Biol.} 2020 issue-1.
#'
#' @export
simPathData <- function(alpha = 0.8,
                           Rsq = 0.7,
                           nsample = 1000,
                           rind = 0.8,
                           nx = 60,
                           nxPrime = 50,
                           ns = 50,
                           ntranscriptNoise = 50,
                           ng = 50,
                           ngPrime = 50,
                           nh = 50,
                           ngeneNoise = 50)
{
  cat("generating obs...")

  noiseSD <- 1
  input <- list(alpha = alpha,
                Rsq = Rsq,
                nsample = nsample,
                rind = rind,
                nx = nx,
                nxPrime = nxPrime,
                ns = ns,
                ntranscriptNoise = ntranscriptNoise,
                ng = ng,
                ngPrime = ngPrime,
                nh = nh,
                ngeneNoise = ngeneNoise,
                noiseSD = noiseSD)


  input <- c(input, Gcolmn = input$ng + input$ngPrime + input$nh + input$ngeneNoise,
             actvGcolmn = input$ng,
             Xcolmn = input$nx + input$nxPrime + input$ns + input$ntranscriptNoise,
             actvXcolmn = input$nx)


  obs <- genIndicators(alpha = input$alpha,
                       Rsq = input$Rsq,
                       nsample = input$nsample,
                       rind = input$rind,
                       nx = input$nx,
                       nxPrime = input$nxPrime,
                       ns = input$ns,
                       ntranscriptNoise = input$ntranscriptNoise,
                       ng = input$ng,
                       ngPrime = input$ngPrime,
                       nh = input$nh,
                       ngeneNoise = input$ngeneNoise,
                       noiseSD = input$noiseSD)

  cat("Done\n")
  return(list(Obs = obs, Input = input))
}


genIndicators <- function(alpha, Rsq, nsample,
                          rind, nx, nxPrime, ns, ntranscriptNoise,
                          ng, ngPrime, nh, ngeneNoise, noiseSD){

  set.seed(NULL)

  set.seed(sample(x = seq(1,20000),size = 1))
  tauInd <- (1 - rind^2)/rind^2
  dfg <- NULL
  dfy <- NULL
  dfx <- NULL
  for (i in 1:nsample){

    sigma <- calSigma(setPar(alpha, Rsq))
    latent <- genLatent(sigma)

    y <- rnorm(1, latent["Yprime"], sqrt(tauInd * sigma["Yprime", "Yprime"]))


    x <- c(rnorm(nx, latent["X"], sqrt(tauInd * sigma["X", "X"])),
           rnorm(nxPrime, latent["Xprime"], sqrt(tauInd * sigma["Xprime", "Xprime"])),
           rnorm(ns, latent["S"], sqrt(tauInd * sigma["S", "S"])),
           rnorm(ntranscriptNoise, 0, noiseSD * sqrt(tauInd * sigma["Yprime", "Yprime"])))


    g <- c(rnorm(ng, latent["G"], sqrt(tauInd * sigma["G", "G"])),
           rnorm(ngPrime, latent["Gprime"], sqrt(tauInd * sigma["Gprime", "Gprime"])),
           rnorm(nh, latent["H"], sqrt(tauInd * sigma["H", "H"] )),
           rnorm(ngeneNoise, 0, noiseSD * sqrt(tauInd * sigma["Yprime", "Yprime"])))


    dfy <- rbind(dfy, y)
    dfx <- rbind(dfx, x)
    dfg <- rbind(dfg, g)
  }
  truegenes <- c(rep(1, ng), rep(0, ngPrime + nh + ngeneNoise) )
  truetranscripts <- c(rep(1, nx), rep(0, nxPrime + ns + ntranscriptNoise) )
  df <- list(y = dfy - mean(y),
             transcripts = scale(dfx, scale=FALSE), # only centered
             genes = scale(dfg, scale=FALSE), 			# only centered
             truetranscripts = truetranscripts,
             truegenes = truegenes)

  dimnames(df$transcripts) <- NULL
  dimnames(df$genes) <- NULL
  dimnames(df$y) <- NULL
  set.seed(NULL)
  return(df)
}

library(MASS)

genLatent <- function(sigma2prime){

  latent <- mvrnorm(1, rep(0, 7), sigma2prime)
  return(latent)
}


calSigma <- function(sysPar){

  cc <- unlist(sysPar)
  for(i in 1:length(cc)){assign(names(cc)[i],cc[i])}

  Sigmaprime <- matrix(c(syxgsh + sxg * alpha + sg + alpha^2 * ss + alpha^2 * sh, # 1,1
                         sxg * sqrt(alpha) + sg * sqrt(alpha), # 1,2
                         sg, # 1,3
                         alpha * ss, #1,4
                         alpha * sh, #1,5
                         sxg * sqrt(alpha) + sg * sqrt(alpha), #2,1
                         sxg + sg * alpha, #2,2
                         sg * sqrt(alpha), #2,3
                         0, #2,4
                         0, #2,5
                         sg, #3,1
                         sg * sqrt(alpha), #3,2
                         sg, #3,3
                         0, #3,4
                         0, #3,5
                         alpha * ss, #4,1
                         0, #4,2
                         0, #4,3
                         ss, #4,4
                         0, #4,5
                         alpha * sh, #5,1
                         0, #5,2
                         0, #5,3
                         0, #5,4
                         sh ), 5,5)
  SigInd <- Sigmaprime[2:3, 2:3]
  Sig2Prime <- matrix(0, nrow=7, ncol=7)
  Sig2Prime[1:5, 1:5] <- Sigmaprime
  Sig2Prime[6:7, 6:7] <- SigInd
  rownames(Sig2Prime) <- c("Yprime", "X", "G", "S", "H", "Xprime", "Gprime")
  colnames(Sig2Prime) <- c("Yprime", "X", "G", "S", "H", "Xprime", "Gprime")
  return(Sig2Prime)
}


setPar <- function(alpha, Rsq)
{
  sg <- 1

  sh <- sg
  bxg <- sqrt(alpha)
  byx.g <- sqrt(alpha)
  byg.x <- 1 - alpha
  bys <- alpha
  byh <- alpha
  tau <- (1 - Rsq)/Rsq
  sxg <- alpha * tau
  ss <- sxg + sg * alpha

  syxgsh <- tau + (alpha ^ 2) * tau * (1 + tau) * (1 + alpha)
  list(alpha=alpha, Rsq=Rsq, sg=sg, ss=ss, sh=sh,
       bxg=bxg, 	byx.g=byx.g, byg.x=byg.x, bys=bys, byh=byh,
       tau=tau, sxg=sxg, syxgsh=syxgsh)
}
