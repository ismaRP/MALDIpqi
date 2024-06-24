
#' Calculate PQI using linear mixed effect model from q2e estimates
#'
#' @param q2e data.frame with q2e estimations per sample, replicate and peptide
#' @param logq Whether q values are log-scaled before entering the LME model
#' @param g Reliability power in the error normal distribution from linear mixed effects model
#' @param return_model
#' @param outdir Optional. directory where results tables and plots are saved
#'
#' @return list contaning the following:
#' \itemize{
#'   \item \strong{sample}: data.frame with aggregated PQI estimates per sample
#'   \item \strong{pep}: data.frame of fitted and residuals per replicate and peptide
#'   \item \strong{estimates}: list of model estimates
#' }
#' \code{sample} data.frame has the following columns:
#' \itemize{
#'   \item Sample: sample name
#'   \item Prediction: sample random effects  obtained from custom formula that. It is the log(PQI)
#'   \item sd: log(PQI) standard error
#'   \item RanefModel: sample random effects, as calculated by \code{\link[nlme]{ranef}}.
#'         It should be the same as Prediction, up to numerical precision
#'   \item PQI.Model, exponentials estimates. It is the PQI
#'
#' }
#'
#' \code{pep} data.frame contains:
#' \itemize{
#'   \item Sample, replicate and peptide IDs
#'   \item q: is the q calculated by the WLS from the isotopic distributions
#'   \item Reliability: minimum least square of q caluclated in the WLS step
#'   \item Loqq: log(q)
#'   \item resp: either q or log(q), according to whether logq was FALSE or TRUE
#'   \item Fitted: fitted q value
#'   \item Res: pearson residuals from Fitted
#'   \item Fitted0: fitted q values at the peptide level. It's equal to the peptide fixed effect
#'   \item Res0: pearson residuals from Fitted0
#' }
#'
#' \code{estimates} list contains
#' \itemize{
#'   \item alpha: peptide fixed effects
#'   \item sigma2_S: sample random effect variance
#'   \item sigma2_R: replicate random effect variance
#'   \item gamma: Reliability 2·gamma exponent
#'   \item sigma2: peptide fixed effects variance
#' }
#'
#' @importFrom nlme lme varComb varPower varIdent varFixed lmeControl
#' @importFrom nlme ranef
#' @importFrom dplyr mutate filter group_by summarise
#' @importFrom tibble as_tibble
#' @importFrom readr write_csv
#' @importFrom stats complete.cases fitted residuals
#' @export
#'
#' @examples
lme_pqi = function(q2e_vals, logq=TRUE, g=NULL, outdir=NULL,
                   return_model=F){

  if (logq){
    q2e_vals = q2e_vals %>% mutate(resp = log(q2e))
  } else {
    q2e_vals = q2e_vals %>% mutate(resp = q2e)
  }

  ## mixed effect model
  if (g == "free"){
    m = lme(
      resp~0+pep_number,
      random = ~1|sample/replicate,
      weights = varComb(
        varPower(-1/2, form = ~reliability),
        varIdent(form = ~1|pep_number)),
      control = lmeControl(maxIter = 1000, msMaxIter = 1000, msMaxEval = 1000),
      data = q2e_vals)
  } else {
    m = lme(
      resp~0+pep_number,
      random = ~1|sample/replicate,
      weight = varComb(
        varFixed(~I(1/reliability)),
        varIdent(form = ~1|pep_number)
      ),
      control=lmeControl(maxIter = 1000, msMaxIter = 1000, msMaxEval = 1000),
      data=q2e_vals)
    # m = lme(
    #   resp~0+Peptides,
    #   random=~1|Sample/Replicates,
    #   weights=varComb(varPower(fixed=g, form=~Reliability),
    #                   varIdent(form=~1|Peptides)),
    #   control=lmeControl(maxIter = 1000, msMaxIter = 1000, msMaxEval = 1000),
    #   data=q2e)
  }

  ## extract parameter estimates from lme object
  estimates_m = extract_estimates(m)
  if (g != "free"){
    estimates_m$gamma = g
  }
  prediction = predict_pqi(m, estimates_m, logq=logq)
  ## mutate into the parent dataframe fitted values and residuals for plotting
  # q2e_m = q2e %>%
  # mutate(Fitted=fitted(m),
  #        Res=residuals(m, type = "pearson"),
  #        Fitted0=fitted(m, level = 0),
  #        Res0=residuals(m, level = 0, type = "pearson")) %>%
  # as_tibble()

  # q2e_m_pred = q2e_m %>% group_by(Sample) %>%
  #   summarise(predict_sample(
  #     Sample, Replicates, Peptides, Reliability, resp,
  #     estimates_m))

  ## here I have added untransformed and transformed prediction from both the model and function
  ## Prediction & PQI.PredictSample => from the function
  ## RanefModel & PQI.Model => from the model
  # if (logq){
  #   q2e_m_pred = q2e_m_pred %>%
  #     mutate(RanefModel=ranef(m)[['Sample']][,1],
  #            PQI.Model=exp(RanefModel),
  #            PQI.PredictSample=exp(Prediction),
  #            Sample = as.character(Sample))
  # } else {
  #   q2e_m_pred = q2e_m_pred %>%
  #     mutate(RanefModel=ranef(m)[['Sample']][,1],
  #            Sample = as.character(Sample),
  #            PQI.Model=RanefModel,
  #            PQI.PredictSample=Prediction)
  # }

  if (!is.null(outdir)){
    write_csv(
      prediction$pep,
      file.path(
        outdir,
        sprintf('PQI_pep_estimates_gamma%3.3f.csv', estimates_m$gamma))
    )
    write_csv(
      prediction$sample,
      file.path(
        outdir,
        sprintf('PQI_sample_estimates_gamma%3.3f.csv', estimates_m$gamma))
    )
  }

  prediction[['estimates']] = estimates_m
  if (return_model) prediction[['model']] = m

  return(prediction)
}


#' Title
#'
#' @param model lme model object
#' @param new_q2e New q2e data for which PQI is predicted using the trained model
#' q is stored in a column called "resp"
#' @param logq whether q is in log-scale in q2e or new_q2e
#' @param estimates
#'
#' @importFrom nlme lme ranef
#' @importFrom stats predict
#' @importFrom dplyr mutate filter group_by summarise
#' @importFrom tibble as_tibble
#' @importFrom stats complete.cases fitted residuals
#' @return
#' @export
#' @examples
predict_pqi = function(model, estimates, new_q2e=NULL, logq=T){

  if (is.null(new_q2e)) {
    q2e = model$data
    pqi = q2e %>% group_by(sample) %>%
      summarise(predict_sample(
        sample, replicate, pep_number, reliability, resp, estimates)) %>%
      mutate(PQI.Model = ranef(model)[['sample']][,1])
    if (logq) {
      pqi = pqi %>% mutate(
        # PQI.PredictSample = exp(PQI.PredictSample),
        PQI.Model = exp(PQI.Model)
      )
    }
    q2e_m = q2e %>% ungroup() %>%
      mutate(Fitted = fitted(model),
             Res = residuals(model, type = "pearson"),
             Fitted0 = fitted(model, level = 0),
             Res0 = residuals(model, level = 0, type = "pearson")) %>%
      as_tibble()
    return(list('pep' = q2e_m, 'sample' = pqi))
  } else {
    new_q2e = new_q2e %>%
      mutate(sample = as.factor(sample), replicate = as.factor(replicate),
             pep_number = as.factor(pep_number))
    # new_q2e = new_q2e %>% filter(residual>0) ## filter dataset to remove 0 reliabilities that resulted in Inf when taken the reciprocal

    pqi = new_q2e %>% group_by(sample) %>%
      summarise(predict_sample(
        sample, replicate, pep_number, reliability, resp, estimates))
    q2e_m = new_q2e %>%
      mutate(
        predicted_q = predict(model, new_q2e)
      )
    if (logq){
      pqi = pqi %>% mutate(PQI.PredictSample = exp(PQI.PredictSample))
    }
    q2e_m = q2e_m %>%
      mutate(
        Pred = exp(predicted_q),
        Res = (resp - predicted_q)/sqrt(predicted_q))
    return(list('sample' = pqi))
  }

}



#' Quantile-Quantile plots of the linear mixed effect estimates per peptide
#'
#' @param pqi_m data.frame of model estimates per replicate and peptide
#' @param title Plot title
#' @param peptides_user A dataframe with peptide information. It must contain at least 3 columns,
#' peptide number or ID, name, and m/z. If NULL, default peptides are used.
#' The number or ID must have the form Pep# and be in the first column.  See \code{\link[https://github.com/ismaRP/MALDIzooMS]{getIsoPeaks}} details.
#' @param label_idx Index where to pull the labels from peptides
#' @param label_func labeller function to process labels. See \code{\link[ggplot2]{labeller}}
#' Default is label_value.
#'
#' @return
#' @importFrom ggplot2 geom_qq geom_qq_line facet_wrap
#' @importFrom ggplot2 ylab xlab ggtitle
#' @export
#'
#' @examples
pept_qqplot = function(pqi_m, title="", peptides_user=NULL, label_idx=2,
                       label_func = label_value){

  if (is.null(peptides_user)) peptides_user = peptides
  pept_labels = pull(peptides_user, label_idx)
  pept_number = pull(peptides_user, 1)
  names(pept_labels) = pept_number

  qq_plot = ggplot(pqi_m) +
    geom_qq(aes(sample=Res, color = pep_number)) +
    geom_qq_line(aes(sample = Res, color = pep_number)) +
    # facet_wrap(~Peptides, scales="free") +
    facet_wrap(
      ~pep_number,
      labeller = labeller(pep_number = as_labeller(pept_labels, label_func))) +
    theme(legend.key.size=unit(1, "cm"),
          legend.text = element_text(size = 15),
          strip.text = element_text(size = 10)) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    ylab("quantiles of standardized residuals") +
    xlab("quantiles of standard normal") +
    ggtitle(title)
  return(qq_plot)
}

#' Fitted versus residuals plot per peptide
#'
#' @param pqi_m data.frame of model estimates per replicate and peptide
#' @param title Plot title
#' @param peptides_user A dataframe with peptide information. It must contain at least 3 columns,
#' peptide number or ID, name, and m/z. If NULL, default peptides are used.
#' The number or ID must have the form Pep# and be in the first column. See \code{\link[https://github.com/ismaRP/MALDIzooMS]{getIsoPeaks}} details.
#' @param label_idx Index where to pull the labels from peptides
#' @param label_func labeller function to process labels. See \code{\link[ggplot2]{labeller}}
#' Default is label_value.
#'
#' @return
#' @importFrom ggplot2 geom_point facet_wrap
#' @importFrom ggplot2 ylab xlab ggtitle
#' @export
#'
#' @examples
fvsr = function(pqi_m, title="", peptides_user=NULL, label_idx=2,
                label_func = label_value){

  if (is.null(peptides_user)) peptides_user = peptides
  pept_labels = pull(peptides_user, label_idx)
  pept_number = pull(peptides_user, 1)
  names(pept_labels) = pept_number

  fvsr_plot = ggplot(pqi_m) +
    geom_point(aes(x = exp(Fitted), y = Res, color=pep_number),
               size = 2.5, alpha = 0.8) +
    # facet_wrap(~Peptides, scales="free") +
    facet_wrap(
      ~pep_number,
      labeller = labeller(pep_number = as_labeller(pept_labels, label_func))) +
    ylab("standardized residuals") +
    xlab("predicted q peptide") +
    theme(legend.key.size=unit(1, "cm"),
          legend.text = element_text(size = 15),
          strip.text = element_text(size = 10)) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    ggtitle(title)
  return(fvsr_plot)
}

