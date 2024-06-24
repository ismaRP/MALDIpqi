

#' Fit a linear model using weights and intercept
#'
#' For use after group_by on samples and peptide index
#'
#' @param norm_int
#' @param weight
#' @param deam_0
#' @param deam_1
#' @param deam_2
#' @param return_model
#'
#' @return
#' @export
#'
#' @examples
#' q2e_vals = peaks %>% filter(n_peaks > 0) %>%
#'   group_by(sample, replicate, pep_number) %>%
#'   summarise(lm_q2e_intercept(norm_int, weight, deam_0, deam_1, deam_2))
lm_q2e_intercept = function(norm_int = NULL,
                            deam_0 = NULL, deam_1 = NULL, deam_2 = NULL,
                            data = NULL,
                            return_model = FALSE) {
  if (!is.null(data)) {
    if (all(is.na(data$norm_int))) allna = TRUE else allna = FALSE
  } else {
    if (all(is.na(norm_int))) allna = TRUE else allna = FALSE
  }
  if (allna) {
    q2e = NA
    intercept = NA
    gamma_0 = NA
    gamma_1 = NA
    gamma_2 = NA
    residual = NA
    if (return_model) return(NA)
  } else {
    lm_model = lm(
      norm_int~deam_0+deam_1+deam_2,
      na.action = na.exclude, data = data)
    if (return_model) return(lm_model)
    intercept = lm_model$coefficients[1]
    gamma_0 = lm_model$coefficients[2]
    gamma_1 = lm_model$coefficients[3]
    gamma_2 = lm_model$coefficients[4]
    if (is.na(gamma_2)) gamma_2 = 0
    q2e = 1 - (gamma_0/(gamma_0 + gamma_1 + gamma_2))
    residual = sqrt(sum(lm_model$residuals^2))/sum(!is.na(norm_int), na.rm = TRUE)
  }
  return(data.frame(q2e=q2e, intercept=intercept,
                    gamma_0=gamma_0, gamma_1=gamma_1, gamma_2=gamma_2,
                    residual=residual))
}


#' Fit a linear model using weights and intercept
#'
#' For use after group_by on samples and peptide index
#'
#' @param norm_int
#' @param weight
#' @param deam_0
#' @param deam_1
#' @param deam_2
#'
#' @return
#' @export
#'
#' @examples
#' q2e_vals = peaks %>% filter(n_peaks > 0) %>%
#'   group_by(sample, replicate, pep_number) %>%
#'   summarise(wlm_q2e_intercept(norm_int, weight, deam_0, deam_1, deam_2))
wlm_q2e_intercept = function(norm_int = NULL, weight = NULL,
                             deam_0 = NULL, deam_1 = NULL, deam_2 = NULL,
                             data = NULL,
                             return_model = FALSE) {
  if (!is.null(data)) {
    data$weight = data$weight / sum(data$weight, na.rm = TRUE)
    if (all(is.na(data$norm_int))) allna = TRUE else allna = FALSE
  } else {
    weight = weight /sum(weight, na.rm = TRUE)
    if (all(is.na(norm_int))) allna = TRUE else allna = FALSE
  }
  if (allna) {
    q2e = NA
    intercept = NA
    gamma_0 = NA
    gamma_1 = NA
    gamma_2 = NA
    residual = NA
    if (return_model) return(NA)
  } else {
    lm_model = lm(
      norm_int~deam_0+deam_1+deam_2,
      weights = weight, data=data,
      na.action = na.exclude)
    if (return_model) return(lm_model)
    intercept = lm_model$coefficients[1]
    gamma_0 = lm_model$coefficients[2]
    gamma_1 = lm_model$coefficients[3]
    gamma_2 = lm_model$coefficients[4]
    if(is.na(gamma_2)) gamma_2 = 0
    q2e = 1 - (gamma_0/(gamma_0 + gamma_1 + gamma_2))
    residual = sqrt(sum(lm_model$weights * lm_model$residuals^2))
  }
  return(data.frame(q2e=q2e, intercept=intercept,
                    gamma_0=gamma_0, gamma_1=gamma_1, gamma_2=gamma_2,
                    residual=residual))
}


#' Fit a linear model with weights and forcing a 0 intercept
#'
#' For use after group_by on samples and peptide index
#'
#' @param norm_int
#' @param weight
#' @param deam_0
#' @param deam_1
#' @param deam_2
#'
#' @return
#' @export
#'
#' @examples
#' q2e_vals = peaks %>% filter(n_peaks > 0) %>%
#'   group_by(sample, replicate, pep_number) %>%
#'   summarise(wlm_q2e(norm_int, weight, deam_0, deam_1, deam_2))
#'
wlm_q2e = function(norm_int = NULL, weight = NULL,
                   deam_0 = NULL, deam_1 = NULL, deam_2 = NULL,
                   data = NULL,
                   return_model = FALSE) {
  if (!is.null(data)) {
    data$weight = data$weight / sum(data$weight, na.rm = TRUE)
    if (all(is.na(data$norm_int))) allna = TRUE else allna = FALSE
  } else {
    weight = weight /sum(weight, na.rm = TRUE)
    if (all(is.na(norm_int))) allna = TRUE else allna = FALSE
  }
  if (allna) {
    q2e = NA
    gamma_0 = NA
    gamma_1 = NA
    gamma_2 = NA
    residual = NA
    if (return_model) return(NA)
  } else {
    lm_model = lm(
      norm_int~0+deam_0+deam_1+deam_2,
      weights = weight, data=data,
      na.action = na.exclude)
    if (return_model) return(lm_model)
    gamma_0 = lm_model$coefficients[1]
    gamma_1 = lm_model$coefficients[2]
    gamma_2 = lm_model$coefficients[3]
    if(is.na(gamma_2)) gamma_2 = 0
    q2e = 1 - (gamma_0/(gamma_0 + gamma_1 + gamma_2))
    residual = sqrt(sum(lm_model$weights * lm_model$residuals^2))
  }
  return(data.frame(q2e=q2e,
                    gamma_0=gamma_0, gamma_1=gamma_1, gamma_2=gamma_2,
                    residual=residual))
}


#' Fit a linear model without weights and forcing a 0 intercept
#'
#' For use after group_by on samples and peptide index
#'
#' @param norm_int
#' @param weight
#' @param deam_0
#' @param deam_1
#' @param deam_2
#'
#' @return
#' @export
#'
#' @examples
#' q2e_vals = peaks %>% filter(n_peaks > 0) %>%
#'    group_by(sample, replicate, pep_number) %>%
#'    summarise(lm_q2e(norm_int, weight, deam_0, deam_1, deam_2))
lm_q2e = function(norm_int = NULL,
                  deam_0 = NULL, deam_1 = NULL, deam_2 = NULL,
                  data = NULL,
                  return_model = FALSE) {
  if (!is.null(data)) {
    if (all(is.na(data$norm_int))) allna = TRUE else allna = FALSE
  } else {
    if (all(is.na(norm_int))) allna = TRUE else allna = FALSE
  }
  if (allna) {
    q2e = NA
    gamma_0 = NA
    gamma_1 = NA
    gamma_2 = NA
    residual = NA
    if (return_model) return(NA)
  } else {
    lm_model = lm(
      norm_int~0+deam_0+deam_1+deam_2,
      na.action = na.exclude, data=data)
    if (return_model) return(lm_model)
    gamma_0 = lm_model$coefficients[1]
    gamma_1 = lm_model$coefficients[2]
    gamma_2 = lm_model$coefficients[3]
    if(is.na(gamma_2)) gamma_2 = 0
    q2e = 1 - (gamma_0/(gamma_0 + gamma_1 + gamma_2))
    residual = sqrt(sum(lm_model$residuals^2))/sum(!is.na(norm_int), na.rm = TRUE)
  }
  return(data.frame(q2e=q2e,
                    gamma_0=gamma_0, gamma_1=gamma_1, gamma_2=gamma_2,
                    residual=residual))
}



#' Estimation of peptide q2e in a single linear regression with nested coefficients
#'
#' Use a single call to "lm" with nested coefficients per spectra and peptide.
#' Effectively it will calculate separate coefficient for each spectra and peptide.
#'
#' WARNING: this cannot handle too large data.frames.
#' I will try to add [biglm] at some point. Otherwise data can be grouped by.
#'
#' @param peaks_df Data.frame with peaks data and theoretical isotopic envelopes
#' @param intensity Term used as the dependent variable. It could be the unscaled
#' intensity or the normalized at the spectra and peptide level
#' @param weights Column of weights in \code{peaks_df}
#' @param intercept Logical, whether to fit an intercept or not
#' @return A data.frame with q2e estimates
#' @export
#'
#' @examples
lm_q2e_oneshot = function(peaks_df, intensity='norm_int', weights=NULL, intercept=FALSE){

  # Prepare formula
  if (intercept) {
    lm_formula = sprintf(paste0(
      "%s ~ 0 + (spectra_name : pep_idx)",
        "+(spectra_name : pep_idx : deam_0)",
        "+(spectra_name : pep_idx : deam_1)",
        "+(spectra_name : pep_idx : deam_2)"),
      intensity)
  } else {
    lm_formula = sprintf(paste0(
      "%s ~ 0",
        "+(spectra_name : pep_idx : deam_0)",
        "+(spectra_name : pep_idx : deam_1)",
        "+(spectra_name : pep_idx : deam_2)"),
      intensity)
  }
  lm_model = lm(as.formula(lm_formula), data=peaks_df)

  result = lm_model$coefficients

  q2e_data = mapply(
    function(n, x){
      spectra_name = grep("spectra_name(.*)", n, value=T)
      spectra_name = strsplit(spectra_name, 'spectra_name')[[1]][2]
      pep_idx = grep("pep_idx\\d", n, value=T)
      pep_idx = strsplit(pep_idx, 'pep_idx')[[1]][2]
      ndeam = grep("deam_\\d", n, value=T)
      if (length(ndeam) == 0) ndeam = 'intercept'
      return(data.frame(spectra_name=spectra_name, pep_idx=pep_idx,
                        ndeam=ndeam, value=x))
    },
    strsplit(names(result), ":"),
    result, SIMPLIFY = F)

  q2e_data = do.call(rbind, q2e_data)
return(q2e_data)

}

