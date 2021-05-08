###############################################################
### Script to replicate the results presented in the manuscript
### "An R Package for the Analysis of Geostatistical Count Data
###  Using Gaussian Copulas"
### Authors: Zifei Han (hanzifei1@gmail.com)
###          Victor De Oliveira (victor.deoliveira@utsa.edu)
###############################################################

##################################################
### Section 3.1. Data simulation and visualization
##################################################

library("gcKrig")

# Simulate ten random fields with negative binomial marginal and exponential correlation function
grid <- seq(0.05, 0.95, by = 0.1)
xloc <- expand.grid(x = grid, y = grid)[, 1]
yloc <- expand.grid(x = grid, y = grid)[, 2]
set.seed(12345)
sim1 <- simgc(locs = cbind(xloc,yloc), sim.n = 10,
              marginal = negbin.gc(mu = 5, od = 1),
              corr = matern.gc(range = 0.3, kappa = 0.5, nugget = 0.1))

# Visualize one of the simulated data
plot(sim1, index = 1, plottype = "Text", col = 2)
plot(sim1, index = 1, plottype = "Dot", col = 2)
plot(sim1, index = 1, plottype = "3D", col = 2)

###############################################################
### Section 3.6. Computation of Frechet–Hoeffding upper bounds
###############################################################

# The Frechet Hoeffding upper bounds for the correlations between negative binomial variables
NBmu <- seq(0.01, 15, by = 0.02)
fhub1 <- vector("numeric", length(NBmu))
for (i in 1:length(NBmu)) {
     fhub1[i] <- FHUBdiscrete(marg1 = "nb", marg2 = "nb", mu1 = 10, od1 = 0.2,
                              mu2 = NBmu[i], od2 = 0.2)
}

# The Frechet Hoeffding upper bounds  for the correlations between gamma variables
gammaShape <- seq(0.01, 15, by = 0.02)
fhub2 <- vector("numeric", length(gammaShape))
for (i in 1:length(gammaShape)) {
     fhub2[i] <- corrTG(marg1 = gm.gc(shape = gammaShape[i], rate = 1),
                        marg2 = gm.gc(shape = 0.5, rate = 1),
                        corrGauss = 1, method = "integral")
}

# Generate the plots
plot(NBmu, fhub1, type = "l", xlab = expression(E(X[2])), ylab = "FHUB")
plot(gammaShape, fhub2, type = "l", xlab = "Shape of a gamma distribution", ylab = "FHUB")

#######################################
### Section 4.1. The Atlantic fish data
#######################################

data("AtlanticFish", package = "gcKrig")

# Calculate the maximum simulated likelihood estimates with three covariates
Fitfish <- mlegc(y = AtlanticFish[, 3], x = AtlanticFish[, 4:6],
                 locs = AtlanticFish[, 1:2], longlat = TRUE,
                 marginal = negbin.gc(link = "log"),
                 corr = matern.gc(kappa = 0.5, nugget = TRUE), method = "GHK")
summary(Fitfish)

# Calculate the 95% Wald-type confidence intervals
round(Fitfish$MLE - qnorm(0.975)*sqrt(diag(solve(Fitfish$hessian))), 3)
round(Fitfish$MLE + qnorm(0.975)*sqrt(diag(solve(Fitfish$hessian))), 3)

# Calculate the 95% profile-likelihood confidence intervals
profile(Fitfish, par.index = 1, method = "GHK")
profile(Fitfish, par.index = 2, method = "GHK")
profile(Fitfish, par.index = 3, method = "GHK")
profile(Fitfish, par.index = 4, method = "GHK")
profile(Fitfish, par.index = 5, method = "GHK")

# Calculate the maximum simulated likelihood estimates with four covariates
Fitfish2 <- mlegc(y = AtlanticFish[, 3], x = AtlanticFish[, 4:7],
                  locs = AtlanticFish[, 1:2], longlat = TRUE,
                  marginal = negbin.gc(link = "log"),
                  corr = matern.gc(kappa = 0.5, nugget = TRUE), method = "GHK")
summary(Fitfish2)

####################################
### Section 4.2. The weed count data
####################################

library("parallel")

# Detect number of cores available for parallel prediction
detectCores(all.tests = FALSE, logical = TRUE)

# Load snowfall for parallel computing
library("snowfall")
data("Weed95", package = "gcKrig")

# Set the observed and prediction locations
weedobs <- subset(Weed95, dummy == 1)
weedpred <- subset(Weed95, dummy == 0)

# Parallel prediction with 4 cores using GHK method
predweed <- predgc(obs.y = weedobs$weedcount, obs.x = weedobs[, 4:5],
                   obs.locs = weedobs[, 1:2], pred.x = weedpred[, 4:5],
                   pred.locs = weedpred[, 1:2], marginal = negbin.gc(link = "log"),
                   corr = matern.gc(kappa = 0.5, nugget = TRUE), method = "GHK",
                   pred.interval = 0.95, parallel = TRUE,
                   paralleloptions = list(n.cores = 4, cluster.type = "SOCK"))
summary(predweed)
plot(predweed, plottype = "2D", xlab = "Coordinate x (m)",
     ylab = "Coordinate y (m)", col = c(3, 2), textcex = 0.8)
plot(predweed, plottype = "Predicted Variance", xlab = "Coordinate x (m)",
     ylab = "Coordinate y (m)", col = c(3, 2), textcex = 0.8)

########################################
### Section 4.3. The oil prevalence data
########################################

data("OilWell", package = "gcKrig")
gridstep <- seq(0.5, 30.5, length = 40)
locOilpred <- data.frame(Easting = expand.grid(gridstep, gridstep)[, 1],
                         Northing = expand.grid(gridstep, gridstep)[, 2])

# Predict the oil abundance with the GHK method using 4 cores in parallel
PredOil <- predgc(obs.y = OilWell[, 3], obs.locs = OilWell[, 1:2],
                  pred.locs = locOilpred, marginal = binomial.gc(link = "logit"),
                  corr = matern.gc(nugget = FALSE), obs.effort = 1,
                  pred.effort = 1, method = "GHK",
                  parallel = TRUE, paralleloptions = list(n.cores = 4))
PredMat <- summary(PredOil)

library("colorspace")

# Generate the contour plots
filled.contour(seq(0.5, 30.5, length = 40), seq(0.5, 30.5, length = 40),
               matrix(PredMat$predMean, 40), zlim = c(0, 1), col = rev(heat_hcl(12)),
               nlevels = 12, xlab = "Eastings (km)", ylab = "Northings (km)",
               plot.axes = {axis(1); axis(2); points(OilWell[, 1:2], col = 1,
                                                     cex = 0.2 + 0.4*OilWell[, 3])})

filled.contour(seq(0.5, 30.5, length = 40), seq(0.5, 30.5, length = 40),
               matrix(PredMat$predVar, 40), zlim = c(0, 0.3), col = rev(heat_hcl(12)),
               nlevels = 10, xlab = "Eastings (km)", ylab = "Northings (km)",
               plot.axes = {axis(1); axis(2); points(OilWell[, 1:2], col = 1,
                                                     cex = 0.2 + 0.4*OilWell[, 3])})

# Predict the oil abundance with the GQT method using 4 cores in parallel
# Simple change the `method' argument from GHK to GQT
PredOil2 <- predgc(obs.y = OilWell[, 3], obs.locs = OilWell[, 1:2],
                  pred.locs = locOilpred, marginal = binomial.gc(link = "logit"),
                  corr = matern.gc(nugget = FALSE), obs.effort = 1,
                  pred.effort = 1, method = "GQT",
                  parallel = TRUE, paralleloptions = list(n.cores = 4))
PredMat2 <- summary(PredOil2)

# Generate the contour plots
filled.contour(seq(0.5, 30.5, length = 40), seq(0.5, 30.5, length = 40),
               matrix(PredMat2$predMean, 40, ), zlim = c(0, 1), col = rev(heat_hcl(12)),
               nlevels = 12, xlab = "Eastings (km)", ylab = "Northings (km)",
               plot.axes = {axis(1); axis(2); points(OilWell[, 1:2], col = 1,
                                                     cex = 0.2 + 0.4*OilWell[, 3])})

filled.contour(seq(0.5, 30.5, length = 40), seq(0.5, 30.5, length = 40),
               matrix(PredMat2$predVar, 40, ), zlim = c(0, 0.3), col = rev(heat_hcl(12)),
               nlevels = 10, xlab = "Eastings (km)", ylab = "Northings (km)",
               plot.axes = {axis(1); axis(2); points(OilWell[, 1:2], col = 1,
                                                     cex = 0.2 + 0.4*OilWell[, 3])})

#############################################
### Section 5. Comparison with other packages
#############################################

# Set the sampling locations
xloc <- rep(seq(0, 1, by = 0.1), 11)
yloc <- rep(seq(0, 1, by = 0.1), each = 11)
simD <- as.matrix(dist(cbind(xloc, yloc)))

# Simulate the dataset
set.seed(321)
simData1 <- simgc(locs = cbind(xloc, yloc), sim.n = 1,
                  marginal = negbin.gc(mu = exp(1 + 0.5*xloc + yloc), od = 1),
                  corr = matern.gc(range = 0.3, nugget = 0))
simDf1 <- data.frame(xloc = xloc, yloc = yloc, data = simData1$data)

# Create a function to calculate likelihood with package 'mvtnorm' by Genz and Bretz
QMCLik <- function(v){
  mu <- exp(v[1] + v[2]*xloc + v[3]*yloc)
  lower <- qnorm(pnbinom(q = simDf1$data - 1, size = 1/v[4], mu = mu))
  upper <- qnorm(pnbinom(q = simDf1$data, size = 1/v[4], mu = mu))
  R <- exp(-simD/v[5])
  set.seed(1234)
  lik <- -log(mvtnorm::pmvnorm(lower = lower, upper = upper,
                               mean = rep(0, length(lower)), sigma = R)[1])
  return(lik)
}

# Find appropriate starting point of marginals and range
est1 <- MASS::glm.nb(data ~ xloc + yloc, data = simDf1)
start0 <- c(coef(est1), 1/est1$theta, median(simD)/2)
EPS <- .Machine$double.eps

# The GB method using package 'mvtnorm'
library("mvtnorm")
library("numDeriv")
t0 <- proc.time()
GBFit <- optim(par = start0, fn = QMCLik, gr = NULL,
               method = c("L-BFGS-B"), lower = c(-Inf, -Inf, -Inf, EPS, EPS),
               upper = c(Inf, Inf, Inf, Inf, Inf))
tGB <- proc.time() - t0
HessGB <- hessian(QMCLik, x = GBFit$par)
sdGB <- sqrt(diag(solve(HessGB)))

# The GHK method using package 'gcmr'
library("gcmr")
t0 <- proc.time()
gcmrGHK <- gcmr(data ~ xloc + yloc, data = simDf1,
                marginal = negbin.marg(link = "log"),
                cormat = matern.cormat(D = simD),
                options = gcmr.options(seed = 123, nrep = c(100, 1000)))
tgcmr <- proc.time() - t0

# The GHK method using package 'gcKrig'
library("gcKrig")
t0 <- proc.time()
gcKrigGHK <- mlegc(y = simDf1[, 3], x = simDf1[, 1:2],
                   locs = cbind(xloc,yloc), method = "GHK",
                   marginal = negbin.gc(link = "log"),
                   corr = matern.gc(nugget = FALSE),
                   ghkoptions = list(nrep = c(100, 1000), seed = 123))
tgcKrigGHK <- proc.time() - t0

# The GQT method using package 'gcKrig'
t0 <- proc.time()
gcKrigGQT <- mlegc(y = simDf1[, 3], x = simDf1[, 1:2],
                   locs = cbind(xloc,yloc), method = "GQT",
                   marginal = negbin.gc(link = "log"),
                   corr = matern.gc(nugget = FALSE))
tgcKrigGQT <- proc.time() - t0

# Table 3
est <- cbind(coef(gcmrGHK), gcKrigGHK$MLE, gcKrigGQT$MLE, GBFit$par)
sd <- cbind(sqrt(diag(vcov(gcmrGHK))), sqrt(diag(solve(gcKrigGHK$hessian))),
            sqrt(diag(solve(gcKrigGQT$hessian))), sdGB)
logLiks <- c(logLik(gcmrGHK), gcKrigGHK$log.lik, gcKrigGQT$log.lik, -GBFit$value)
times <- c(tgcmr["elapsed"], tgcKrigGHK["elapsed"], tgcKrigGQT["elapsed"], tGB["elapsed"])
Est <- matrix(paste0(format(round(est, 3), nsmall = 3),
                     " (", format(round(sd, 3), nsmall = 3), ")"), ncol = 4)
colnames(Est) <- c("gcmr", "gcKrig-GHK", "gcKrig-GQT", "mvtnorm")
rownames(Est) <- names(coef(gcmrGHK))
Table.3 <- rbind(Est,
                 logLik = format(round(logLiks, 2), nsmall = 2),
                 Time = format(round(times, 1), nsmall = 1))
print(Table.3, quote = FALSE)

##################################################
### Appendix
##################################################

# Example of specifying binomial marginal with t distribution as link
binomial_t.gc <- function(df.t = 6, size = NULL, prob = NULL)
{
  ans <- list()
  ans$discrete <- TRUE
  if(is.null(size) & is.null(prob)){
    ans$start <- function(y, x, effort) {
      mfit <- suppressWarnings(glm.fit(x, y/effort,
                                       family = binomial(link = "logit")))
      reg0 <- coef(mfit)
      glb <- qnorm(pbinom(y - 1, size = effort, prob = fitted(mfit)))
      gub <- qnorm(pbinom(y, size = effort, prob = fitted(mfit)))
      names(reg0)[1] <- "Intercept"
      return(reg0)
    }
    ans$nod <- 0
    ans$bounds <- function(y, x, pars, effort) {
      p <- pt(pars[1:ncol(x)]%*%t(x), df = df.t)
      a <- qnorm(pbinom(y - 1, size = effort, prob = p))
      b <- qnorm(pbinom(y, size = effort, prob = p))
      return(list(lower = a, upper = b))
    }
    ans$pdf <- function(y, x, pars, effort){
      p <- pt(pars[1:ncol(x)]%*%t(x), df = df.t)
      pdf <- dbinom(y, size = effort, prob = p, log = FALSE)
      return(pdf)
    }
    ans$cdf <- function(y, x, pars, effort){
      p <- pt(pars[1:ncol(x)]%*%t(x), df = df.t)
      cdf <- pbinom(y, size = effort, prob = p)
      return(cdf)
    }
  }
  if(is.numeric(size) & is.numeric(prob))
  {
    q <- function(p) qbinom(p = p, size = size, prob = prob)
    ans$margvar <- size*prob*(1-prob)
    ans$int.marg <- function (order){
      if(requireNamespace("EQL", quietly = TRUE)){
        integrate(function(x, order)
          ifelse((q(pnorm(x))==Inf | dnorm(x) < .Machine$double.eps), 0,
                 q(pnorm(x))*dnorm(x)*EQL::hermite(x, order, prob = TRUE)),
          order = order, -8, 8, subdivisions = 1500,
          rel.tol = 0.01, stop.on.error = FALSE)
      }else{
        stop("Please install {EQL} first!")
      }
    }
    ans$rsample <- function(nrep){
      rbinom(n = nrep, size = size, prob = prob)
    }
    ans$q <- q
  }
  class(ans) <- c("marginal.gc")
  return(ans)
}


# Example of specifying power exponential correlation function
powerexp.gc <- function (range = NULL, kappa = 1, nugget = TRUE)
{
  ans <- list()
  if (kappa > 2)
    stop("'Kappa' must be between 0 and 2")
  if (is.null(range)) {
    if (nugget == TRUE) {
      ans$nug <- 1
      ans$npar.cor <- 2
      ans$start <- function(D) {
        corstart <- c(median(D)/2, 0.2)
        names(corstart) <- c("range", "nugget")
        return(corstart)
      }
      ans$corr <- function(corrpar, D) {
        S <- (1 - corrpar[2]) * exp(-abs((D/corrpar[1])^(kappa))) +
          corrpar[2] * diag(NROW(D))
        return(S)
      }
    }
    else {
      ans$nug <- 0
      ans$npar.cor <- 1
      ans$start <- function(D) {
        corstart <- median(D)/2
        names(corstart) <- "range"
        return(corstart)
       }
      ans$corr <- function(corrpar, D) {
        S <- exp(-abs((D/corrpar[1])^(kappa)))
        return(S)
      }
    }
  }
  if (is.numeric(range) & is.numeric(nugget)) {
    ans$S <- function(D) (1 - nugget) * exp(-abs((D/range)^(kappa))) +
      nugget * diag(NROW(D))
  }
  class(ans) <- c("corr.gc")
  return(ans)
}
