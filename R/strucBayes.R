#' Local Bayesian modeling
#'
#' Local Bayesian modeling for the EEG data with structured spike and slab prior. Performs variable selection and prediction for local modeling using structured spike and slab prior.
#' @param y binary response \eqn{n} values
#' @param X data tensor \eqn{n x L x \tau}
#' @param dist.mat distance matrix between \eqn{L} locations
#' @param Nmcmc number of MCMC samples to generate
#' @param burnin burnin samples
#' @param thin thinning samples \eqn{>0}
#' @param X.test test data \eqn{m x p x tau}, defaults to NULL
#' @param v0 hyperparameter \eqn{v_0}, defaults to 0.005
#' @param a1 hyperparameter \eqn{a_1}, defaults to 5
#' @param a2 hyperparameter \eqn{a_2}, defaults to 50
#' @param q hyperparameter \eqn{q}, defaults to 10
#' @keywords strucBayes
#' @export
#'
#' @importFrom mvtnorm rmvnorm
#' @import MASS stats
#' @examples strucBayes()

strucBayes = function(y = NULL, # binary response n values
                      X = NULL, # data tensor n x L x tau
                      dist.mat = NULL, # distance matrix between L locations
                      Nmcmc = NULL, # number of MCMC samples to generate
                      burnin = NULL, # burnin samples
                      thin = NULL, # thinning samples
                      X.test = NULL, # test data m x p x tau
                      v0 = 0.005, # hyperparameter v0
                      a1 = 5, # hyperparameter a1
                      a2 = 50, # hyperparameter a2
                      q = 10 # hyperparameter q
){

  n = dim(X)[1]
  p = dim(X)[2]
  tau = dim(X)[3]

  #########################
  ### 1. Local modeling ###
  #########################

  c_init = 1 # initial value for c

  nusq_init = 1/rgamma(n = p, shape = a1, rate = a2) # initial value for nu^2

  # prior covariance of lambda based on distance matrix
  dmat = (dist.mat+t(dist.mat))/2
  W.mat = exp(-(dmat^2)/0.1)
  D.mat = diag(rowSums(W.mat)+0.1)
  G.mat = D.mat-W.mat
  prior.Delta = chol2inv(chol(G.mat))

  l_init = as.vector(rmvnorm(1, rnorm(p,0,0.1), prior.Delta)) # initial value for lambda

  # initial value for zeta
  zeta_init = rbinom(p, 1, pnorm(l_init))
  zeta_init[zeta_init==0] = v0

  beta_init = t(rnorm(p, 0, sqrt(zeta_init*nusq_init))) # initial value for beta

  ind = seq(burnin+1, Nmcmc, by=thin) # indices of MCMC samples to use

  # store MCMC samples of parameters
  zeta.array = array(dim = c(length(ind),p,tau))
  b.array = array(dim = c(length(ind),p,tau))
  l.array = array(dim = c(length(ind),p,tau))
  c.array = array(dim = c(length(ind),1,tau))

  alp = 0.5 # value of alpha to use in the prior mean of lambda in local models
  prior.mu = rep(0, p) # prior mean for lambda in local model at time point 1

  # sequential estimation of local models
  for(t in 1:tau){
    print(paste("Time = ", t, sep = ""))

    temp.try = tryCatch({
      sim = strucPriorT(y,
                        X[,,t],
                        c_init,
                        beta_init,
                        zeta_init,
                        nusq_init,
                        l_init,
                        prior.mu,
                        prior.Delta,
                        v0,
                        a1,
                        a2,
                        q,
                        Nmcmc,
                        ind)

      list(zeta = sim$zeta,
           b = sim$b,
           l = sim$l,
           c = sim$c,
           pr.mu = alp*colMeans(sim$l))
    }, error=function(e){
      cat("ERROR :",conditionMessage(e), "\n", "Time point failed = ", t, "\n")

      zeta = b = l = array(dim = c(length(ind),p))
      c = array(dim = c(length(ind),1))
      pr.mu = alp*prior.mu

      list(zeta = zeta,
           b = b,
           l = l,
           c = c,
           pr.mu = pr.mu)
    })

    zeta.array[,,t] = temp.try$zeta
    b.array[,,t] = temp.try$b
    l.array[,,t] = temp.try$l
    c.array[,,t] = temp.try$c
    prior.mu = temp.try$pr.mu
  }

  #############################
  ### 2. variable selection ###
  #############################

  # proportion of selections for the locations across all the time points
  loc.probs = matrix(nrow = tau, ncol = p)
  for(t in 1:tau) loc.probs[t,] = colMeans(zeta.array[,,t]==1)

  sort.probs = apply(loc.probs, 2, sort) # sort the proportions by time for each location

  # cluster sorted proportions into two clusters
  sort.clus = kmeans(t(sort.probs), centers = 2)
  ind1 = which(sort.clus$cluster==1)
  ind2 = which(sort.clus$cluster==2)
  clus1.mean = mean(loc.probs[,ind1])
  clus2.mean = mean(loc.probs[,ind2])

  # identify the locations to be selected and update the ...
  # ... coefficients for remaining locations as zero
  b.strucBayes = matrix(0, nrow = tau, ncol = p)
  if(clus1.mean==clus2.mean) b.strucBayes = apply(b.array, c(3,2), mean)
  if(clus1.mean>clus2.mean) b.strucBayes[,ind1] = apply(b.array, c(3,2), mean)[,ind1]
  if(clus2.mean>clus1.mean) b.strucBayes[,ind2] = apply(b.array, c(3,2), mean)[,ind2]

  #####################
  ### 3. Prediction ###
  #####################

  y.strucBayes.pred = NA

  if(!is.null(X.test)){
    m = dim(X.test)[1]

    # estimate c for each local model by the median of its MCMC samples
    c.est = apply(c.array, c(3,2), median)

    # local prediction probabilities for each subject
    p.strucBayes.pred = pnorm(sapply(1:tau,
                                     function(t) crossprod(t(X.test[,,t]),
                                                           b.strucBayes[t,])/c.est[t,]))

    # local response predictions
    y.matrix.pred = matrix(nrow = m, ncol = tau)
    for(i in 1:m) y.matrix.pred[i,] = as.numeric(p.strucBayes.pred[i,]>0.5)

    # weights corresponding to local predictions
    w = abs(p.strucBayes.pred-0.5)
    for(i in 1:m) w[i,] = (w[i,]^2)/sum(w[i,]^2)

    # final weighted prediciton
    y.w = sapply(1:m, function(i) sum(w[i,]*y.matrix.pred[i,]))
    y.strucBayes.pred = as.numeric(y.w>0.5)
  }

  list(b.strucBayes = b.strucBayes, # estimates of beta
       p.strucBayes.pred = p.strucBayes.pred, # local prediction probabilities
       y.strucBayes.pred = y.strucBayes.pred # final predictions
  )
}
