
#' Labeler helper for plotting peptide spectra
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




#' Get isotopic distributions
#' @description
#' Wrapper for bacollite ms_iso.
#' Calculates theoretical isotopic deistributions for a list of
#' sequences
#' @param seqs
#' @param nhyds Integer, vector
#'
#' @return Array of intensity distributions, without the masses
#' @export
#' @importFrom bacollite ms_iso
#' @importFrom stringr str_count
#'
#' @examples
get_isodists = function(seqs, ndeam, nhyds, long_format=F){
  n_isopeaks = 5
  if (!long_format) {
    iso_peps = list()
    for (d in 1:ndeam) {
      isotopes = array(0, dim = c(length(seqs), n_isopeaks))
      for (i in 1:length(seqs)) {
        nQs = str_count(seqs[i], 'Q')
        if (nQs >= d) {
          isodist = ms_iso(seqs[i], ndeamidations = d, nhydroxylations = nhyds[i])
          isotopes[i, (d + 1):n_isopeaks] = isodist[1:(n_isopeaks - d), 2]
        } # else {
          # isotopes[i,] = NaN
        # }
      }
      iso_peps[[d]] = isotopes
    }
  } else {
    pep_idx = rep(1:length(seqs), each = n_isopeaks)
    mass_pos = rep(1:n_isopeaks, times = length(seqs))
    iso_peps = data.frame(
      pep_idx = pep_idx,
      mass_pos = mass_pos
    )
    for (d in 0:ndeam) {
      isotopes = array(0, dim = length(seqs)*n_isopeaks)
      for (i in 1:length(seqs)) {
        nQs = str_count(seqs[i], 'Q')
        if (nQs >= d) {
          isodist = ms_iso(seqs[i], ndeamidations = d, nhydroxylations = nhyds[i])
          isotopes[((i - 1)*5 + 1 + d):(i*5)] = isodist[1:(n_isopeaks - d),2]
        } # else {
          # isotopes[((i-1)*5+1):(i*5)] = NaN
        # }
      }
      iso_peps[[paste0('deam_', d)]] = isotopes
    }
  }
  return(iso_peps)
}



## define predict sample function
#' Title
#'
#' @param Sample
#' @param Replicates
#' @param Peptides
#' @param Reliability
#' @param resp Response variable, either q or Logq
#'
#' @param pars Contains the following parameters:
#' \describe{
#'     \item{data_I_s}{part of the data frame corresponding to Sample s.}
#'     \item{alpha}{estimate for fixed effect}
#'     \item{sigma2_S}{estimate for variance for random effect of sample.}
#'     \item{sigma2_R}{estimate for variance for random effect of replication.}
#'     \item{gamma}{estimate for half power on Reliability}
#'     \item{sigma2}{estimate for variance for on the different peptides. This should be a named vector.}
#'}
#' @return tibble with predicted sample PQI and standard deviation
#' @export
#' @importFrom tibble tibble
#' @examples
predict_sample <- function(sample, replicate, pep_number, reliability, resp, pars) {

  data_I_s = data.frame(
    sample = sample, replicate = replicate, pep_number = pep_number,
    reliability = reliability, resp = resp)

  alpha = pars$alpha
  sigma2_S = pars$sigma2_S
  sigma2_R = pars$sigma2_R
  gamma = pars$gamma
  sigma2 = pars$sigma2

  if (nrow(data_I_s) == 0) {
    # prediction and variance without any data
    E_X = 0
    var_X = sigma2_S
  } else {
    # sanity check
    if (length(unique(data_I_s$sample)) != 1) stop("Observations must come from the same sample.")
    # build Xi
    Xi = sigma2_S * matrix(1, nrow(data_I_s), nrow(data_I_s)) +
      sigma2_R * outer(data_I_s$replicate, data_I_s$replicate, function(x,y){as.numeric(x == y)}) +
      diag((data_I_s$reliability^(2*gamma)) * sigma2[as.character(data_I_s$pep_number)], nrow = nrow(data_I_s))
    # prediction
    # We don't need the estimate, as the ranef function provides it
    # E_X = sigma2_S * sum(solve(Xi, data_I_s$resp - alpha[as.character(data_I_s$pep_number)]))
    # variance
    var_X = sigma2_S - (sigma2_S^2) * sum(solve(Xi, rep(1, nrow(data_I_s))))
  }

  # return
  # names(E_X) = names(var_X) <- NULL
  # return(tibble(PQI.PredictSample = E_X, sd = sqrt(var_X)))
  return(tibble(sd = sqrt(var_X)))
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
extract_estimates = function(m) {
  q_data = m$data
  alpha.m = fixef(m)
  names(alpha.m) = substr(names(alpha.m), 9, 12) ## to align the code, remove Peptides from the names
  tmp = coef(m$modelStruct$reStruct, FALSE) * (m$sigma^2)
  sigma2_S.m = tmp["sample.var((Intercept))"] ## variance parameter estimate for sample
  sigma2_R.m = tmp["replicate.var((Intercept))"] ## variance parameter estimate for replicate
  gamma.m = coef(m$modelStruct$varStruct$A, FALSE)
  sigma2.m = m$sigma * c("pep1" = 1, coef(m$modelStruct$varStruct$B, FALSE))^2
  sigma2.m = sigma2.m[match(levels(q_data$pep_number), names(sigma2.m))]

  return(list(alpha = alpha.m, sigma2_S = sigma2_S.m, sigma2_R = sigma2_R.m,
              gamma = gamma.m, sigma2 = sigma2.m))
}

