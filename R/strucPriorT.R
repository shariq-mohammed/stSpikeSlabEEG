#' MCMC algorithm for the estimation of local model
#'
#' Gibbs sampling, and slice sampling within Gibbs sampling for a local model
#'
#' @param y binary response n values
#' @param X data at time \eqn{t} - dimension \eqn{n x L}
#' @param c_init initial value for \eqn{c}
#' @param beta_init initial values for \eqn{\beta}
#' @param zeta_init initial values for \eqn{\zeta}
#' @param nusq_init initial values for \eqn{\nu^2}
#' @param l_init initial value for \eqn{\lambda}
#' @param prior.mu prior mean for \eqn{\lambda}
#' @param prior.Delta prior covariance for \eqn{\lambda}
#' @param v0 hyperparameter \eqn{v0}
#' @param a1 hyperparameter \eqn{a1}
#' @param a2 hyperparameter \eqn{a2}
#' @param q hyperparameter \eqn{q}
#' @param Nmcmc number of MCMC samples to generate
#' @param ind indices of MCMC samples to use
#' @keywords strucPriorT()
#' @export
#'
#' @importFrom truncnorm rtruncnorm
#' @import MASS
#' @examples strucPriorT()

########################################################
### MCMC algorithm for the estimation of local model ###
########################################################

strucPriorT = function(y, # binary response n values
                       X, # data at time t - dimension nxL
                       c_init, # initial value for c
                       beta_init, # initial values for beta
                       zeta_init, # initial values for zeta
                       nusq_init, # initial values for nu^2
                       l_init, # initial value for lambda
                       prior.mu, # prior mean for lambda
                       prior.Delta, # prior covariance for lambda
                       v0, # hyperparameter v0
                       a1, # hyperparameter a1
                       a2, # hyperparameter a2
                       q, # hyperparameter q
                       Nmcmc,  # number of MCMC samples to generate
                       ind # indices of MCMC samples to use
){

  nc = ncol(X)
  nr = nrow(X)

  c = c_init
  b = beta_init
  zeta = zeta_init
  nusq = nusq_init
  l = l_init

  uni.slice.calls = 0	# Number of calls of the slice sampling function
  uni.slice.evals = 0	# Number of density evaluations done in these calls

  # Storing all the required parameters from the MCMC samples
  c_store = matrix(nrow = Nmcmc, ncol = 1)
  zeta_store = matrix(nrow = Nmcmc, ncol = nc)
  b_store = matrix(nrow = Nmcmc, ncol = nc)
  #nusq_store = matrix(nrow = Nmcmc, ncol = nc)
  l_store = matrix(nrow = Nmcmc, ncol = nc)

  # conditional posterior mean and variance of lambdas
  mu.expr = matrix(nrow = nc, ncol = nc-1)
  cond.sigsq = matrix(nrow = nc, ncol = 1)

  for(k in 1:nc){
    mu.expr[k,] = crossprod(prior.Delta[k,-k], chol2inv(chol(prior.Delta[-k,-k])))
    cond.sigsq[k,] = prior.Delta[k,k] - crossprod(mu.expr[k,], prior.Delta[-k,k])
  }

  # posterior mean of latent parameter z
  zMean = crossprod(t(X), t(b))

  # start Gibbs sampling
  for(i in 1:Nmcmc){
    if(i %% 1000 == 0) print(i)

    # update z
    z = (y*rtruncnorm(1, 0, Inf, zMean, sd = sqrt(c))+
           (1-y)*rtruncnorm(1, -Inf, 0, zMean, sd = sqrt(c)))

    # update beta
    symX = crossprod(X, X/c)
    gm = nusq*zeta
    GMinv = diag(1/gm)
    Sig = chol2inv(chol(symX+GMinv))
    temp1 = crossprod(X, z/c)
    bMean = crossprod(t(Sig), temp1)
    bVar = Sig

    b = rmvnorm(1, mean = bMean, sigma = bVar)
    b_store[i,] = b

    # update zeta
    tempExp = exp(-(b^2)/(2*nusq))
    pr1 = (1-pnorm(l))*(tempExp^(1/v0))/sqrt(v0)
    pr2 = pnorm(l)*tempExp

    temp.probs = pr2/(pr1+pr2)
    temp.probs[is.na(temp.probs)] = 0.5
    zeta = rbinom(nc, 1, prob = temp.probs)
    zeta[zeta==0] = v0
    zeta_store[i,] = zeta

    # update lambda by slice sampling within Gibbs sampling
    for(k in 1:nc){
      cond.mu = prior.mu[k]+crossprod(mu.expr[k,], l[-k]-prior.mu[-k])

      if(zeta[k] == 1){
        l[k] = uni.slice(x0=l[k],
                         function(x){
                           log(dnorm(x, mean = cond.mu, sd = sqrt(cond.sigsq[k,]))*pnorm(x))-
                             log(pnorm(cond.mu/sqrt(1+cond.sigsq[k,])))
                         }, w = 8*sqrt(cond.sigsq[k,]))[1]
      }else{
        l[k] = uni.slice(x0=l[k],
                         function(x){
                           log(dnorm(x, mean = cond.mu, sd = sqrt(cond.sigsq[k,]))*pnorm(-x))-
                             log(pnorm(-cond.mu/sqrt(1+cond.sigsq[k,])))
                         }, w = 8*sqrt(cond.sigsq[k,]))[1]
      }
    }
    l_store[i,] = l

    # update nu^2
    nusq = 1/rgamma(nc, shape = a1+0.5, rate = a2+((b^2)/(2*zeta)))
    #nusq_store[i,] = nusq

    # update c
    zMean = crossprod(t(X), t(b))
    c = 1/rgamma(1, shape = (nr+q)/2, rate = (q+sum((z-zMean)^2))/2)
    c_store[i,] = c
  }

  # return samples of beta, zeta, lambda, and c as they will be needed for ...
  # ... variable selection and prediction
  list(b = b_store[ind,],
       zeta = zeta_store[ind,],
       l = l_store[ind,],
       #nusq = nusq_store[ind,],
       c = c_store[ind,]
  )
}
