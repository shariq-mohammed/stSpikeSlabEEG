#' Five fold cross-validation on the EEG data
#'
#' Five fold cross-validation by specifying a specific split by seed
#'
#' @param seed seed for creating splits, defaults to 113
#' @param Nmcmc number of MCMC samples to generate
#' @param burnin burnin samples
#' @param thin thinning samples \eqn{>0}
#' @param v0 hyperparameter v0, defaults to 0.005
#' @param a1 hyperparameter a1, defaults to 5
#' @param a2 hyperparameter a2, defaults to 50
#' @param q hyperparameter q, defaults to 10
#' @param nCores number of cores to use to run cross validation parallelly, defaults to 1
#' @keywords crossVal
#' @export
#'
#' @import foreach doParallel parallel
#' @examples crossVal()

crossVal = function(seed = 113, # seed for creating splits
                    Nmcmc = NULL, # number of MCMC samples to generate
                    burnin = NULL, # burnin samples
                    thin = NULL, # thinning samples
                    v0 = 0.005, # hyperparameter v0
                    a1 = 5,# hyperparameter a1
                    a2 = 50,# hyperparameter a2
                    q = 10,# hyperparameter q
                    nCores = 1 # number of cores to use to run cross validation parallelly
){

  # Create five fold split of 122 observations
  f.split = matrix(nrow = 25, ncol = 5)
  set.seed(seed)
  f.split[,1] = sample((1:122),size = 25)
  set.seed(seed)
  f.split[,2] = sample((1:122)[-f.split[,1]],size = 25)
  set.seed(seed)
  f.split[,3] = sample((1:122)[-as.vector(f.split[,1:2])],size = 25)
  set.seed(seed)
  f.split[,4] = sample((1:122)[-as.vector(f.split[,1:3])],size = 25)
  set.seed(seed)
  f.split[1:22,5] = sample((1:122)[-as.vector(f.split[,1:4])],size = 22)

  X = stSpikeSlabEEG::X_eeg # EEG data
  ind64to57 = stSpikeSlabEEG::ind64to57_eeg # location indices for which distance matrix is available
  X = X[,ind64to57,] # update EEG data for locations above
  y = t(stSpikeSlabEEG::y_eeg) # binary response
  dist.mat = stSpikeSlabEEG::dist.mat_eeg # distance matrix

  # dimensions of EEG data
  tau = dim(X)[1]
  p = dim(X)[2]
  n = dim(X)[3]

  # scale data at subject level by its Frobenius norm
  X_int = sapply(1:n, function(i) sqrt(sum(X[,,i]^2)))
  X_sc = array(NA, dim = dim(X))
  for(i in 1:n) X_sc[,,i] = X[,,i]/X_int[i]

  # begin parallel computation of five fold cross-validation
  cl = makeCluster(nCores)
  registerDoParallel(cl, cores = nCores) # register cluster
  Nreal = foreach(f.num = 1:5, .inorder=TRUE,
                  .packages = c('MASS', 'truncnorm', 'mvtnorm'),
                  .errorhandling = 'pass') %dopar% {

                    trn.ind = na.omit(as.vector(f.split[,-f.num])) # indices for training data

                    y.trn = y[trn.ind] # responses from training data
                    X.trn = aperm(X_sc[,,trn.ind], c(3,2,1)) # training EEG data
                    y.test = y[-trn.ind] # responses from test data
                    X.test = aperm(X_sc[,,-trn.ind], c(3,2,1)) # test EEG data

                    f = strucBayes(y.trn,
                                   X.trn,
                                   dist.mat,
                                   Nmcmc,
                                   burnin,
                                   thin,
                                   X.test,
                                   v0,
                                   a1,
                                   a2,
                                   q)
                    list(b.strucBayes = f$b.strucBayes,
                         p.strucBayes.pred = f$p.strucBayes.pred,
                         y.strucBayes.pred = f$y.strucBayes.pred,
                         y.test = y.test)
                  }
  stopCluster(cl)
  names(Nreal) = paste0('Fold ',1:5)

  Nreal
}
