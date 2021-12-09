
###################
# Hidden functions #
###################

# Read sample
#' Title
#'
#' @param f
#' @param data_path
#'
#' @return
#' @importFrom readr read_tsv
#' @importFrom  magrittr %>%
#' @importFrom  tibble tibble
#' @importFrom data.table fread
#' @examples
read_peaks = function(f, i, data_path){
  # read_tsv(file.path(data_path, f),
  #          col_names = F,
  #          col_types = 'ddd') %>% as.data.frame()
  if (i %% 500 == 0){
    print(sprintf('Peaks files read for %i samples', i))
  }
  data.frame(fread(file.path(data_path, f),
               colClasses=c("numeric", "numeric", "numeric"), sep="\t"))
}


# Labeler helper for plotting peptide spectra
#' Title
#'
#' @param values
#'
#' @return
#'
#' @examples
peptide_labeller = function(values){
  return(pept_labels[values])
}


#' Title
#'
#' @param lgth_g number of parameters
#' @param p peptide index
#' @param isotopes theoretical isotopic distributions
#'
#' @return Design matrix
#'
#' @examples
design = function(lgth_g, p, isotopes) {
  X = array(0,dim=c(dim(isotopes)[[2]],lgth_g))
  for (i in 1:dim(isotopes)[[2]]){
    for (j in 1:min(i,lgth_g)) {
      X[i,j] = isotopes[p,i-j+1]
    }
  }
  X
}


# Get isotopic distributions
#' Get isotopic distributions
#' @description
#' Wrapper for ballite ms_iso.
#' Calculates theoretical isotpic deistributions for a list of
#' sequences
#' @param seqs
#' @param nhyds
#'
#' @return Array of intensity distributions, without the masses
#' @export
#' @importFrom bacollite ms_iso
#'
#' @examples
get_isodists = function(seqs, nhyds){
  isotopes = array(0, dim=c(length(seqs), 5))
  for (i in 1:length(seqs)){
    isodist = ms_iso(seqs[i], nhydroxylations = nhyds[i])
    isotopes[i,] = isodist[,2]
  }
  isotopes
}

#' LS estimate of gamma
#'
#' @param X
#' @param data
#' @param weight
#'
#' @return
#'
#' @examples
LSest = function(X, data, weight) {
  # assuming diagonal weight matrix
  Xprime = diag(sqrt(weight)) %*% X
  dataprime = as.matrix(data) %*% diag(sqrt(weight))
  g <- as.vector( solve( t(Xprime)%*%Xprime ) %*% t(Xprime) %*% t(dataprime) )
  c(g,sum(( t(dataprime) - Xprime %*% g )^2))
}



## define predict sample function
#' Title
#'
#' @param Sample
#' @param Replicates
#' @param Peptides
#' @param Reliability
#' @param Logq
#'     log(q)
#' @param pars Contains the following parameters:
#' \describe{
#'     \item{data_I_s}{part of the data frame corresponding to Sample s.}
#'     \item{alpha}{estimate for fixed effect}
#'     \item{sigma2_S}{estimate for variance for random effect of sample.}
#'     \item{sigma2_R}{estimate for variance for random effect of replication.}
#'     \item{gamma}{estimate for half power on Reliability}
#'     \item{sigma2}{estimate for variance for on the different peptides. This should be a named vector.}
#'}
#' @return tibble with predicted sample PDI and standard deviation
#' @export
#' @importFrom tibble tibble
#' @examples
predict_sample <- function(Sample, Replicates, Peptides, Reliability, Logq, pars) {

  data_I_s = data.frame(
    Sample=Sample, Replicates=Replicates, Peptides=Peptides,
    Reliability=Reliability, Logq=Logq)

  alpha = pars$alpha
  sigma2_S = pars$sigma2_S
  sigma2_R = pars$sigma2_R
  gamma = pars$gamma
  sigma2 = pars$sigma2

  if (nrow(data_I_s)==0) {
    # prediction and variance without any data
    E_X = 0
    var_X = sigma2_S
  } else {
    # sanity check
    if (length(unique(data_I_s$Sample))!=1) stop("Observations must come from the same sample.")
    # build Xi
    Xi = sigma2_S*matrix(1,nrow(data_I_s),nrow(data_I_s)) +
      sigma2_R*outer(data_I_s$Replicates,data_I_s$Replicates,function(x,y){as.numeric(x==y)}) +
      diag((data_I_s$Reliability^(2*gamma))*sigma2[as.character(data_I_s$Peptides)],nrow=nrow(data_I_s))
    # prediction
    E_X = sigma2_S*sum(solve(Xi,data_I_s$Logq-alpha[as.character(data_I_s$Peptides)]))
    # variance
    var_X = sigma2_S - (sigma2_S^2)*sum(solve(Xi,rep(1,nrow(data_I_s))))
  }

  # return
  names(E_X) = names(var_X) <- NULL
  return(tibble(Prediction=E_X, sd=sqrt(var_X)))
}


#' Title
#'
#' @param m
#' @param q_data
#'
#' @return
#' @importFrom nlme fixef
#' @importFrom stats coef
#'
#' @examples
extract_estimates = function(m, q_data) {
  alpha.m = fixef(m)
  names(alpha.m) = substr(names(alpha.m), 9, 12) ## to align the code, remove Peptides from the names
  tmp = coef(m$modelStruct$reStruct, FALSE) * (m$sigma^2)
  sigma2_S.m = tmp["Sample.var((Intercept))"] ## variance parameter estimate for sample
  sigma2_R.m = tmp["Replicates.var((Intercept))"] ## variance parameter estimate for replicate
  gamma.m = coef(m$modelStruct$varStruct$A, FALSE)
  sigma2.m = m$sigma * c("Pep1"=1, coef(m$modelStruct$varStruct$B, FALSE))^2
  sigma2.m = sigma2.m[match(levels(q_data$Peptides), names(sigma2.m))]

  return(list(alpha = alpha.m, sigma2_S = sigma2_S.m, sigma2_R = sigma2_R.m,
              gamma = gamma.m, sigma2 = sigma2.m))
}

