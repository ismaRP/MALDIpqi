


#' Title
#'
#' @param q2e
#' @param peptides
#' @param outdir
#' @param n_isopeaks
#' @param g
#'
#' @return
#' @importFrom nlme lme varComb varPower varIdent lmeControl
#' @importFrom nlme ranef
#' @importFrom dplyr pull mutate rename filter group_by summarise
#' @importFrom tibble as_tibble
#' @importFrom readr write_csv
#' @importFrom stats complete.cases residuals fitted
#' @export
#'
#' @examples
lme_pdi = function(q2e, peptides, outdir=NULL, n_isopeaks=5, g="free"){

  pept_names = pull(peptides, 1)
  pept_labels = pull(peptides, 6)
  names(pept_labels) = paste0('Pep', peptides)

  n_peptides = nrow(peptides)
  n_samples = nrow(q2e)/n_peptides
  Peptides = paste0('Pep', q2e$peptide)
  Replicates = rep(c("1","2","3"), times = n_samples/3 * n_peptides)
  PeptideMass = rep(pull(peptides, 3), each=n_samples)

  q2e = q2e %>%
    rename(Reliability = minLS, Sample = sample) %>%
    mutate(Sample = as.factor(Sample), Replicates = as.factor(Replicates),
           Peptides = as.factor(Peptides), PeptideMass = as.factor(PeptideMass),
           Logq=log(q))
  q2e = q2e %>% dplyr::filter(Reliability>0) ## filter dataset to remove 0 reliabilities that resulted in Inf when taken the reciprocal


  ii = c("Sample","Replicates","Peptides","PeptideMass","q","Reliability","Logq")

  q2e = q2e[complete.cases(q2e[,ii]),ii]
  q2e$Sample = droplevels(q2e$Sample)

  ## mixed effect model
  if (g == "free"){
    m = lme(
      Logq~0+Peptides,
      random=~1|Sample/Replicates,
      weights=varComb(varPower(-1/2, form=~Reliability),
                      varIdent(form=~1|Peptides)),
      control=lmeControl(maxIter = 1000, msMaxIter = 1000, msMaxEval = 1000),
      data=q2e)
  } else {
    # m = lme(
    #   Logq~0+Peptides,
    #   random=~1|Sample/Replicates,
    #   weight=varComb(
    #     varFixed(~I(1/Reliability)),
    #     varIdent(form=~1|Peptides)
    #   ),
    #   control=lmeControl(maxIter = 1000, msMaxIter = 1000, msMaxEval = 1000),
    #   data=q2e)
    m = lme(
      Logq~0+Peptides,
      random=~1|Sample/Replicates,
      weights=varComb(varPower(fixed=g, form=~Reliability),
                     varIdent(form=~1|Peptides)),
      control=lmeControl(maxIter = 1000, msMaxIter = 1000, msMaxEval = 1000),
      data=q2e)
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
      Sample, Replicates, Peptides, Reliability, Logq,
      estimates_m))

  ## here I have added untransformed and transformed prediction from both the model and function
  ## Prediction & PDI.PredictSample => from the function
  ## RanefModel & PDI.Model => from the model
  q2e_m_pred = q2e_m_pred %>%
    mutate(RanefModel=ranef(m)[["Sample"]][,1],
           PDI.Model=exp(RanefModel),
           PDI.PredictSample=exp(Prediction),
           Sample = as.character(Sample))

  if (!is.null(outdir)){
    write_csv(
      q2e_m,
      file.path(
        outdir,
        sprintf('PDI_pep_estimates_gamma%3.3f.csv', estimates_m$gamma))
    )
    write_csv(
      q2e_m_pred,
      file.path(
        outdir,
        sprintf('PDI_sample_estimates_gamma%3.3f.csv', estimates_m$gamma))
    )
  }

  return(list("pep" = q2e_m, "sample" = q2e_m_pred, "estimates" = estimates_m))
}


#' Title
#'
#' @param q2e
#' @param peptides
#'
#' @return
#' @importFrom dplyr pull
#' @importFrom ggplot2 ggplot geom_histogram geom_density geom_vline facet_wrap
#' @importFrom ggplot2 xlim xlab ggtitle theme
#' @export
#'
#' @examples
plot_q = function(q2e, peptides){

  pept_names = pull(peptides, 1)
  pept_labels = pull(peptides, 6)
  names(pept_labels) = paste0('Pep', peptides)

  q_hist = ggplot(q2e, aes(x=q)) +
    geom_histogram(aes(y=..density..,  fill=Peptides), colour="black")+
    geom_density(aes(fill=Peptides), color='black', alpha=0.6) +
    geom_vline(xintercept=c(0,1)) +
    xlim(-1, 2) +
    facet_wrap(
      ~Peptides,
      labeller = labeller(Peptides=as_labeller(pept_labels, label_parsed))) +
    xlab("q") +
    ggtitle("q WLS per peptide estimate") +
    theme(plot.title = element_text(size=10))

  return(q_hist)
}



#' Title
#'
#' @param pdi_m
#' @param tit
#'
#' @return
#' @importFrom ggplot2 geom_qq geom_qq_line facet_wrap
#' @importFrom ggplot2 ylab xlab ggtitle
#' @export
#'
#' @examples
pept_qqplot = function(pdi_m, tit=""){
  qq_plot = ggplot(pdi_m) +
    geom_qq(aes(sample=Res, color=Peptides)) +
    geom_qq_line(aes(sample=Res, color=Peptides)) +
    facet_wrap(~Peptides, scales="free") +
    theme(legend.key.size=unit(1, "cm"),
          legend.text = element_text(size = 15),
          strip.text = element_text(size = 15)) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    ylab("quantiles of standardized residuals") +
    xlab("quantiles of standard normal") +
    ggtitle(tit)
  return(qq_plot)
}

#' Title
#'
#' @param pdi_m
#' @param tit
#'
#' @return
#' @importFrom ggplot2 geom_point facet_wrap
#' @importFrom ggplot2 ylab xlab ggtitle
#' @export
#'
#' @examples
fvsr = function(pdi_m, tit=""){
  fvsr_plot = ggplot(pdi_m) +
    geom_point(aes(x=exp(Fitted), y=Res, color=Peptides),
               size=2.5, alpha=0.8) +
    facet_wrap(~Peptides, scales="free") +
    ylab("standardized residuals") +
    xlab("predicted q peptide") +
    theme(legend.key.size=unit(1, "cm"),
          legend.text = element_text(size = 15),
          strip.text = element_text(size = 15)) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    ggtitle(tit)
  return(fvsr_plot)
}

