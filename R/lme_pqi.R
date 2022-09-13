


#' Calculate PQI using linear mixed effect model from q2e estimates
#'
#' @param q2e data.frame with q2e estimations per sample, replicate and peptide
#' @param logq Whether q values are log-scaled before entering the LME model
#' @param g Reliability power in the error normal distribution from linear mixed effects model
#' @param peptides A dataframe with peptide information. It must contain at least 3 columns,
#' peptide number or ID, name, and m/z. If NULL, default peptides are used. See \code{\link[MALDIutils]{getIsoPeaks}} details.
#' The number or ID must have the form Pep# and be in the first column.
#' @param outdir Optional. directory where results tables and plots are saved
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
#'   \item PQI.PredictSample and PQI.Model, exponentials of Prediction and RanefModel respectively.
#'         They contain the PQI
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
#' @importFrom dplyr pull mutate rename filter group_by summarise
#' @importFrom tibble as_tibble
#' @importFrom readr write_csv
#' @importFrom stats complete.cases residuals fitted
#' @export
#'
#' @examples
lme_pqi = function(q2e, logq=TRUE, g=NULL, peptides=peptides, outdir=NULL){

  n_peptides = nrow(peptides)
  n_samples = nrow(q2e)/n_peptides
  Peptides = q2e$Peptides
  Replicates = rep(c("1","2","3"), times = n_samples/3 * n_peptides)

  q2e = q2e %>%
    rename(Reliability = minLS, Sample = sample) %>%
    mutate(Sample = as.factor(Sample), Replicates = as.factor(Replicates),
           Peptides = as.factor(Peptides),
           Logq=log(q))
  q2e = q2e %>% filter(Reliability>0) ## filter dataset to remove 0 reliabilities that resulted in Inf when taken the reciprocal

  ii = c("Sample","Replicates","Peptides","q","Reliability","Logq")

  q2e = q2e[complete.cases(q2e[,ii]),ii]
  q2e$Sample = droplevels(q2e$Sample)

  if (logq){
    q2e = q2e %>% mutate(resp = Logq)
  } else {
    q2e = q2e %>% mutate(resp = q)
  }

  ## mixed effect model
  if (g == "free"){
    m = lme(
      resp~0+Peptides,
      random=~1|Sample/Replicates,
      weights=varComb(varPower(-1/2, form=~Reliability),
                      varIdent(form=~1|Peptides)),
      control=lmeControl(maxIter = 1000, msMaxIter = 1000, msMaxEval = 1000),
      data=q2e)
  } else {
    m = lme(
      Logq~0+Peptides,
      random=~1|Sample/Replicates,
      weight=varComb(
        varFixed(~I(1/Reliability)),
        varIdent(form=~1|Peptides)
      ),
      control=lmeControl(maxIter = 1000, msMaxIter = 1000, msMaxEval = 1000),
      data=q2e)
    # m = lme(
    #   resp~0+Peptides,
    #   random=~1|Sample/Replicates,
    #   weights=varComb(varPower(fixed=g, form=~Reliability),
    #                   varIdent(form=~1|Peptides)),
    #   control=lmeControl(maxIter = 1000, msMaxIter = 1000, msMaxEval = 1000),
    #   data=q2e)
  }



  ## mutate into the parent dataframe fitted values and residuals for plotting
  q2e_m = q2e %>%
    mutate(Fitted=fitted(m),
           Res=residuals(m, type = "pearson"),
           Fitted0=fitted(m, level = 0),
           Res0=residuals(m, level = 0, type = "pearson")) %>%
    as_tibble()



  ## extract parameter estimates from lme object
  estimates_m = extract_estimates(m, q2e_m)
  if (g != "free"){
    estimates_m$gamma = g
  }
  ## build a matrix containing the sample predictions and their standard errors
  # hatX_m0 = matrix(0,length(levels(q_data$Sample)),2)
  # rownames(hatX_m0) = levels(q_data$Sample)
  # colnames(hatX_m0) = c("Prediction","sd")
  # for (ii in levels(q_data$Sample)) {
  #   hatX_m0[ii,] = predict_sample(dplyr::filter(q_data,Sample==ii), alpham0, sigma2_S.m0,
  #                                              sigma2_R.m0, gamma.m0, sigma2.m0)
  # }
  # hatX_m0 = as_tibble(hatX_m0, rownames = 'Sample')

  q2e_m_pred = q2e_m %>% group_by(Sample) %>%
    summarise(predict_sample(
      Sample, Replicates, Peptides, Reliability, resp,
      estimates_m))

  ## here I have added untransformed and transformed prediction from both the model and function
  ## Prediction & PQI.PredictSample => from the function
  ## RanefModel & PQI.Model => from the model
  if (logq){
    q2e_m_pred = q2e_m_pred %>%
      mutate(RanefModel=ranef(m)[['Sample']][,1],
             PQI.Model=exp(RanefModel),
             PQI.PredictSample=exp(Prediction),
             Sample = as.character(Sample))
  } else {
    q2e_m_pred = q2e_m_pred %>%
      mutate(RanefModel=ranef(m)[['Sample']][,1],
             Sample = as.character(Sample),
             PQI.Model=RanefModel,
             PQI.PredictSample=Prediction)
  }

  if (!is.null(outdir)){
    write_csv(
      q2e_m,
      file.path(
        outdir,
        sprintf('PQI_pep_estimates_gamma%3.3f.csv', estimates_m$gamma))
    )
    write_csv(
      q2e_m_pred,
      file.path(
        outdir,
        sprintf('PQI_sample_estimates_gamma%3.3f.csv', estimates_m$gamma))
    )
  }

  return(list("pep" = q2e_m, "sample" = q2e_m_pred, "estimates" = estimates_m))
}


#' Quantile-Quantile plots of the linear mixed effect estimates per peptide
#'
#' @param pqi_m data.frame of model estimates per replicate and peptide
#' @param title Plot title
#' @param peptidesA dataframe with peptide information. It must contain at least 3 columns,
#' peptide number or ID, name, and m/z. If NULL, default peptides are used.
#' The number or ID must have the form Pep# and be in the first column.  See \code{\link[MALDIutils]{getIsoPeaks}} details.
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
pept_qqplot = function(pqi_m, title="", peptides=peptides, label_idx=2,
                       label_func = label_value){

  pept_labels = pull(peptides, label_idx)
  pept_number = pull(peptides, 1)
  names(pept_labels) = pept_number

  qq_plot = ggplot(pqi_m) +
    geom_qq(aes(sample=Res, color=Peptides)) +
    geom_qq_line(aes(sample=Res, color=Peptides)) +
    # facet_wrap(~Peptides, scales="free") +
    facet_wrap(
      ~Peptides,
      labeller = labeller(Peptides=as_labeller(pept_labels, label_func))) +
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
#' @param peptides A dataframe with peptide information. It must contain at least 3 columns,
#' peptide number or ID, name, and m/z. If NULL, default peptides are used.
#' The number or ID must have the form Pep# and be in the first column. See \code{\link[MALDIutils]{getIsoPeaks}} details.
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
fvsr = function(pqi_m, title="", peptides=peptides, label_idx=2,
                label_func = label_value){

  pept_labels = pull(peptides, label_idx)
  pept_number = pull(peptides, 1)
  names(pept_labels) = pept_number

  fvsr_plot = ggplot(pqi_m) +
    geom_point(aes(x=exp(Fitted), y=Res, color=Peptides),
               size=2.5, alpha=0.8) +
    # facet_wrap(~Peptides, scales="free") +
    facet_wrap(
      ~Peptides,
      labeller = labeller(Peptides=as_labeller(pept_labels, label_func))) +
    ylab("standardized residuals") +
    xlab("predicted q peptide") +
    theme(legend.key.size=unit(1, "cm"),
          legend.text = element_text(size = 15),
          strip.text = element_text(size = 10)) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    ggtitle(title)
  return(fvsr_plot)
}

