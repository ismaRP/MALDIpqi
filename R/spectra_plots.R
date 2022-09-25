


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
#' @param normalize
#'
#' @return
#' @importFrom MALDIquant match.closest as.matrix
#' @importFrom tibble tibble
#' @examples
plot_preprocessing = function(m, p_name, s, s_br, p, baseline, noise, s_name,
                              n_isopeaks, min_isopeaks, tol, normalize){
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

  if (normalize){
    tot_int = sum(br)
    int = int / tot_int
    b = b / tot_int

    br = br / tot_int
    n = n / tot_int
  }

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
#' @param peptides_user
#' A dataframe with peptide information. It must contain at least 3 columns,
#' peptide number or ID, name, and m/z. If NULL, default are used, see \code{\link[MALDIutils]{getIsoPeaks}} details.
#' The number or ID must have the form Pep# and be in the first column.
#' @param label_idx Column index in peptides where the label is stored
#' @param label_func labeller function to process labels. See \code{\link[ggplot2]{labeller}}
#' Default is label_value.
#' @param normalize Logical. Whether intensity hsould be normalized divinding by
#' the total intensity for each peptide peaks range. DEfault is False
#' @param readf
#' A string value. choose function to use to read spectra.
#' Currently restricted to one of "fread", "table" or "mzml"
#' @param smooth_method
#' Smootherning method, one of "SavitzkyGolay" or "Wavelet".
#' SavitzkyGolay uses \code{\link[MALDIquant]{smoothIntensity}} together with
#' \code{hws_smooth} parameter
#' Wavelet uses \code{\link[MALDIrppa]{wavSmoothing}} together with
#' \code{thresh.scale}
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
                             peptides_user=NULL,
                             label_idx=2, label_func = label_value,
                             normalize=F,
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

  smooth_method = match.arg(smooth_method)
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


  if (is.null(peptides_user)) peptides_user = peptides

  masses = pull(peptides_user, 3)
  pept_labels = pull(peptides_user, label_idx)
  pept_number = pull(peptides_user, 1)
  names(pept_labels) = pept_number

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
      pept_number,
      MoreArgs = list(
        s = s, s_br = s_br, p = peaks,
        baseline = baseline,
        noise = noise,
        s_name = s_name,
        n_isopeaks = n_isopeaks,
        min_isopeaks = min_isopeaks,
        tol=tol, normalize=normalize),
      SIMPLIFY = FALSE
    )
    spectra_df = bind_rows(spectra_df)
    data[[i]] = spectra_df
  }
  data = bind_rows(data)
  data$peptide = factor(data$peptide, levels = pept_number)
  peaks = data %>% dplyr::filter(peak == T)


  befp = ggplot(data) +
    geom_line(aes(x=mz, y=baseline, color=spectra)) +
    geom_line(aes(x=mz, y=int, color=spectra)) +
    facet_wrap(~peptide, ncol=3, scales='free',
               labeller = labeller(peptide=as_labeller(pept_labels, label_func))) +
    xlab('m/z') + ylab('Intensity') +
    theme(strip.text.x = element_text(size=8))
  befp

  aftp = ggplot(data) +
    geom_line(aes(x=mz, y=brint, color=spectra)) +
    geom_point(aes(x=mz, y=brint, color=spectra), data=peaks, size=3) +
    facet_wrap(~peptide, ncol=3, scales='free',
               labeller = labeller(peptide=as_labeller(pept_labels, label_func))) +
    xlab('m/z') + ylab('Intensity') +
    theme(strip.text.x = element_text(size=8))
  aftp

  preprocessing_plot = ggarrange(befp, aftp, nrow=2, common.legend=T,
                                 labels = "AUTO")

  return(preprocessing_plot)
}
