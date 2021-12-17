


#' Title
#'
#' @param m
#' @param p_name
#' @param s
#' @param s_br
#' @param p
#' @param baseline
#' @param noise
#' @param s_name
#' @param n_isopeaks
#' @param min_isopeaks
#' @param tol
#'
#' @return
#' @importFrom MALDIquant match.closest as.matrix
#' @importFrom tibble tibble
#' @examples
plot_preprocessing = function(m, p_name, s, s_br, p, baseline, noise, s_name,
                              n_isopeaks, min_isopeaks, tol){
  s = as.matrix(s)
  s_br = as.matrix(s_br)
  p  = as.matrix(p)
  sel_p = rep(F, n_isopeaks)
  p_complete = rep(NA, n_isopeaks)

  sel = s[,1] > m[1]-0.5 & s[,1] < m[1]+n_isopeaks-0.5
  mz = s[sel, 1]
  int = s[sel, 2]
  b = baseline[sel]
  br = s_br[sel, 2]
  n = noise[sel]

  p = p[p[,1] > m[1]-0.5 & p[,1] < m[1]+n_isopeaks-0.5, 1]

  idx = match.closest(p, m, tolerance = tol*m)
  sel_p[idx[!is.na(idx)]] = T

  # sel_p = !is.na(match.closest(p[,1], m,
  #                              tolerance = tol*m))

  cs = cumsum(sel_p)
  first_missing = (cs[1] == 0) & (cs[min_isopeaks + 1] == min_isopeaks)
  last_missing = cs[min_isopeaks] == min_isopeaks

  if (!(first_missing | last_missing)){
    sel_p = rep(F, n_isopeaks)
  }

  p_complete[idx[!is.na(idx)]] = p[!is.na(idx)]
  p_complete[!sel_p] = NA

  detected = mz %in% p_complete

  df = tibble(mz=mz, int=int, brint=br, baseline=b, noise = n,
              peptide=p_name, spectra=s_name, peak = detected)
  return(df)
}


#' Title
#'
#'
#' @param indir Folder containing spectra.
#' @param spectra_names Character vector. Spectra names to plot
#' @param peptides
#' A dataframe with peptide information. It must contain at least 3 columns,
#' peptide number or ID, name, and m/z. IF NULL
#' @param readf
#' A string value. choose function to use to read spectra.
#' Currently restricted to one of "fread", "table" or "mzml"
#' @param smooth_method
#' Smootherning method, one of "SavitzkyGolay" or "Wavelet".
#' SavitzkyGolay uses \code{\link[MALDIquant]{smoothIntensity}} together with
#' \code{hws_smooth} parameter
#' Wavelet uses \code{\link[MALDIrppa]{wavSmoothing}} together with
#' \code{thresh.scale}
#'
#' @param thresh.scale
#' Smoothing factor for wavelet-based smoothing. Passed to \code{\link[MALDIrppa]{wavSmoothing}}
#' @param hws_smooth
#' Half-window size parameter for smoothening. Passed to \code{\link[MALDIquant]{smoothIntensity}}
#' @param iter
#' Iterations parameter for baseline detection. Passed to \code{\link[MALDIquant]{estimateBaseline}}
#' @param hws_peak
#' Half-window size parameter for local maximum detection. Passed to \code{\link[MALDIquant]{detectPeaks}}
#' @param snr
#' Signal to noise threshold for peak detection. Passed to \code{\link[MALDIquant]{detectPeaks}}
#' @param tol
#' Tolerance factor for matching detected peaks to theoretical isotopic distribution.
#' Value used is m/z of monoisotopic peak * \code{tol}.
#' @param n_isopeaks
#' Number of isotopic peaks to pick. Default is 5 and the maximum permitted.
#' @param min_isopeaks
#' If less than min_isopeaks consecutive (about 1 Da difference) isotopic peaks
#' are detected, the whole isotopic envelope is discarded. Default is 4
#' @return ggplot with raw nd preprocessed spectra
#' @export
#' @importFrom MALDIquant smoothIntensity removeBaseline detectPeaks estimateNoise estimateBaseline
#' @importFrom ggplot2 ggplot geom_line geom_point facet_wrap aes
#' @importFrom ggplot2 xlab ylab theme element_text label_parsed
#' @importFrom ggplot2 labeller as_labeller label_parsed
#' @importFrom ggpubr ggarrange
#' @importFrom MALDIutils importTsv importTable import_file.MzMl
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows
#' @importFrom readr read_lines
#'
#' @examples
plot_pept_spectra = function(indir,
                             spectra_names,
                             readf,
                             peptides=NULL,
                             smooth_method = c("SavitzkyGolay", "Wavelet"),
                             thresh.scale = 2.5,
                             hws_smooth = 8,
                             iter= 20,
                             hws_peak = 20,
                             snr = 0,
                             n_isopeaks = 5,
                             min_isopeaks = 4,
                             tol = 1.5e-4){
  switch(EXPR=readf,
         "fread" = {
           fmt = ".tab"
           read_f = importTsv
         },
         "table" = {
           fmt = ".tab"
           read_f = importTable
         },
         "mzml" = {
           fmt = ".mzML"
           read_f = import_file.MzMl

         }
  )

  switch(EXPR=smooth_method,
         "Wavelet" = {
           sm = function(thresh.scale){
             function(x) wavSmoothing(list(x), "Wavelet", thresh.scale)[[1]]
           }
           sm = sm(thresh.scale)
         },
         "SavitzkyGolay" = {
           sm = function(hws_smooth){
             function(x) smoothIntensity(x, "SavitzkyGolay", hws_smooth)
           }
           sm = sm(hws_smooth)
         }
  )
  # spectra_names = read_lines(spectra)

  spectra_paths = file.path(indir, paste0(spectra_names, fmt))

  spectra = lapply(
    spectra_paths,
    read_f
  )


  if (is.null(peptides)) {
    peptides = load("data/peptides.rda")
    # masses = c(
    #   1105.58, 2019.95, 2040.97, 2689.25,
    #   3033.50, 3093.48, 3084.42, 3116.40
    # )
    # pept_names = c(
    #   'COL1a1 508-519', 'COL1a1 270-291', 'COL1a1 375-396', 'COL1a1 934-963',
    #   'COL1a2 756-789', 'COL1a2 756-789 p', 'COL1a1 9-42', 'COL1a1 9-42 p'
    # )
    # pept_labels = c(
    #   'atop(paste("a) COL1", alpha, "1 508-519 ", "1"%*%"(Pro->Hyp)"), paste("GV", bold("Q"), "GPPGPAGEEGKR"))',
    #   'atop(paste("b) COL1", alpha, "1 270-291 ", "2"%*%"(Pro->Hyp)"), paste("GEPGPTGI", bold("Q"), "GPPGPAGEEGKR"))',
    #   'atop(paste("c) COL1", alpha, "1 375-396 ", "3"%*%"(Pro->Hyp)"), paste("TGPPGPAG", bold("Q"), "DGRPGPPGPPGAR"))',
    #   'atop(paste("d) COL1", alpha, "1 934-963 ", "2"%*%"(Pro->Hyp)"), paste("GFSGL", bold("Q"), "GPPGPPGSPGE", bold("Q"), "GPSGASGPAGPR"))',
    #   'atop(paste("e) COL1", alpha, "2 756-789 ", "5"%*%"(Pro->Hyp)"), paste("GPSGEPGTAGPPGTPGP", bold("Q"), "GLLGAPGFLGLPGSR"))',
    #   'atop(paste("f) COL1", alpha, "2 756-789 ", "5"%*%"(Pro->Hyp)"), paste("GPSGEPGTAGPPGTPGP", bold("Q"), "GFLGPPGFLGLPGSR"))',
    #   'atop(paste("g) COL1", alpha, "1 9-42 ", "5"%*%"(Pro->Hyp)"), paste("GLPGPPGAPGP", bold("Q"), "GF", bold("Q"), "GPPGEPGEPGASGPMGPR"))',
    #   'atop(paste("h) COL1", alpha, "1 9-42 ", "7"%*%"(Pro->Hyp)"), paste("GLPGPPGAPGP", bold("Q"), "GF", bold("Q"), "GPPGEPGEPGASGPMGPR"))'
    # )
    # names(pept_labels) = pept_names
  }
  pept_names = pull(peptides, 2)
  masses = pull(peptides, 3)
  pept_labels = pull(peptides, 6)
  names(pept_labels) = pept_names


  masses = matrix(masses, nrow = n_isopeaks,
                  ncol = length(masses), byrow = T)
  d = 1.00235
  masses = masses + (d * 0L:(n_isopeaks - 1L))
  masses = as.list(as.data.frame(masses))

  data = list()
  for (i in 1:length(spectra)){
    s = spectra[[i]]
    s_name = spectra_names[i]
    smoothened = sm(s)
    # smoothened = smoothIntensity(
    #   s, method="SavitzkyGolay",
    #   halfWindowSize=hws_smooth)
    baseline = estimateBaseline(smoothened, method = 'SNIP',
                                iterations = iter)
    baseline = baseline[,2]
    s_br = removeBaseline(smoothened, 'SNIP', iterations = iter)
    noise = estimateNoise(s_br, 'SuperSmoother')[,2]
    peaks = detectPeaks(s_br, halfWindowSize=hws_peak,
                        method='SuperSmoother', SNR=snr)

    spectra_df = mapply(
      plot_preprocessing,
      masses,
      pept_names,
      MoreArgs = list(
        s = s, s_br = s_br, p = peaks,
        baseline = baseline,
        noise = noise,
        s_name = s_name,
        n_isopeaks = n_isopeaks,
        min_isopeaks = min_isopeaks,
        tol=tol),
      SIMPLIFY = FALSE
    )
    spectra_df = bind_rows(spectra_df)
    data[[i]] = spectra_df
  }
  data = bind_rows(data)
  data$peptide = factor(data$peptide, levels = pept_names)
  peaks = data %>% dplyr::filter(peak == T)


  befp = ggplot(data) +
    geom_line(aes(x=mz, y=baseline, color=spectra)) +
    geom_line(aes(x=mz, y=int, color=spectra)) +
    facet_wrap(~peptide, ncol=3, scales='free',
               labeller = labeller(peptide=as_labeller(pept_labels, label_parsed))) +
    xlab('m/z') + ylab('Intensity') +
    theme(strip.text.x = element_text(size=8))
  befp

  aftp = ggplot(data) +
    geom_line(aes(x=mz, y=brint, color=spectra)) +
    geom_point(aes(x=mz, y=brint, color=spectra), data=peaks, size=3) +
    facet_wrap(~peptide, ncol=3, scales='free',
               labeller = labeller(peptide=as_labeller(pept_labels, label_parsed))) +
    xlab('m/z') + ylab('Intensity') +
    theme(strip.text.x = element_text(size=8))
  aftp

  preprocessing_plot = ggarrange(befp, aftp, nrow=2, common.legend=T,
                                 labels = "AUTO")

  return(preprocessing_plot)
}
