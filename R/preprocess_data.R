

#' preprocess_spectra
#'
#' @param s
#' \code{\link[MALDIquant]{MassSpectrum}} object
#' @param ...
#' Parameters passed to the different preprocessing functions.
#' See \code{\link{preprocess_data}}
#' @return peaks list, with mass, intensity and s2n
#' @importFrom MALDIquant smoothIntensity removeBaseline detectPeaks
#' @importFrom MALDIutils peptidePseudoClusters
#' @export
#'
#' @examples
preprocess_spectra = function(s, smoothf, iterations, halfWindowSize,
                              SNR, masses, tol, n_isopeaks, min_isopeaks){
  # prep_args = list(...)
  # s = smoothIntensity(s, method="SavitzkyGolay",
  #                     halfWindowSize=prep_args$hws_smooth)
  s = smoothf(s)
  s = removeBaseline(s, method='SNIP', iterations)
  # Estimate noise
  p = detectPeaks(s, halfWindowSize, method='SuperSmoother', SNR)
  # Extract peaks
  p = peptidePseudoClusters(p, masses, tol, n_isopeaks, min_isopeaks)
  return(p)
}




#' Preprocessing spectra for q2e estimation
#'
#' Performs smoothening, baseline removal and peak detection on MALDI samples.
#' From the peaks, isotopic peaks for a list of peptides are extracted.
#'
#' @param indir Folder containing spectra.
#' @param peptides_user
#' A dataframe with peptide information. It must contain at least 3 columns,
#' peptide number or ID, name, and m/z. IF NULL, default are used, see details.
#' The number or ID must have the form Pep# and be in the first column.
#' @param readf
#' A string value. Choose function used to read spectra.
#' Currently restricted to one of "fread", "table" or "mzml".
#' See \code{\link[MALDIutils]{preprocessData}} for more info.
#' @param outdir
#' Folder where one peaks file for each sample should be written.
#' If NULL (default), no files are written.
#' @param nchunks
#' Number of chunks all the spectra should be divided into for reading and
#' processing. If all spectra is loaded and processed at once, i.e. chunks=1,
#' it can overload RAM memory. If chunks>1 data is loaded and processed to the
#' much lighter list of peaks in chunks batches.
#' @param smooth_method
#' Smootherning method, one of "SavitzkyGolay" or "Wavelet".
#' SavitzkyGolay uses \code{\link[MALDIquant]{smoothIntensity}} together with
#' \code{hws_smooth} parameter
#' Wavelet uses \code{\link[MALDIrppa]{wavSmoothing}} together with
#' \code{thresh.scale}
#' @param hws_smooth
#' Half-window size parameter for smoothening. Passed to \code{\link[MALDIquant]{smoothIntensity}}
#' @param thresh.scale
#' Smoothing factor for wavelet-based smoothing. Passed to \code{\link[MALDIrppa]{wavSmoothing}}
#' @param iterations
#' Iterations parameter for baseline detection. Passed to \code{\link[MALDIquant]{estimateBaseline}}
#' @param halfWindowSize
#' Half-window size parameter for local maximum detection. Passed to \code{\link[MALDIquant]{detectPeaks}}
#' @param SNR
#' Signal to noise threshold for peak detection. Passed to \code{\link[MALDIquant]{detectPeaks}}
#' @param tol
#' Tolerance factor for matching detected peaks to theoretical isotopic distribution.
#' Value used is m/z of monoisotopic peak * \code{tol}.
#' @param n_isopeaks
#' Number of isotopic peaks to pick. Default is 5 and the maximum permitted.
#' @param min_isopeaks
#' If less than min_isopeaks consecutive (about 1 Da difference) isotopic peaks
#' are detected, the whole isotopic envelope is discarded. Default is 4
#' @param ncores
#' Number of cores used for the preprocessing of spectra. The cores will work in
#' parallel with the different spectra within a chunk.
#' @param iocores
#' Number of cores used for I/O operations. For some systems I/O operations involving
#' multiple cores reading or writing at the same time from disk increases time.
#' @param sep
#' Separator character for txt spectra files
#' @return A list of dataframes, 1 per sample. Each dataframe has 3 columns,
#' m/z, intensity and signal-to-noise ratio for each of the n_isopeaks from each
#' peptide. Missing peaks are NAs.
#'
#' @importFrom MALDIutils preprocessData prepFun
#' @importFrom MALDIrppa wavSmoothing
#' @importFrom parallel mcmapply detectCores mclapply
#' @importFrom readr write_tsv
#' @importFrom dplyr pull
#' @export
#' @details The default peptides are the ones from Nair et al. (2022).
#' The paper contains the details on the preprocessing procedure.
#' @references
#'
#' Nair, B. et al. (2022) ‘Parchment Glutamine Index (PQI): A novel method to estimate glutamine deamidation levels in parchment collagen obtained from low-quality MALDI-TOF data’, bioRxiv. doi:10.1101/2022.03.13.483627.
#'
#' @examples
getIsoPeaks = function(indir,
                       readf = c("fread", "table", "mzml"),
                       sep="\t",
                       peptides_user=NULL,
                       spectra=NULL,
                       outdir = NULL,
                       nchunks = 50,
                       smooth_method = c("SavitzkyGolay","Wavelet"),
                       thresh.scale = 2.5,
                       hws_smooth = 20,
                       iterations = 20,
                       halfWindowSize = 20,
                       SNR = 0,
                       tol = 1.5e-4,
                       n_isopeaks = 6,
                       min_isopeaks = 4,
                       ncores = NULL,
                       iocores = 1,
                       vch=10){

  readf = match.arg(readf)
  smooth_method = match.arg(smooth_method)
  switch(EXPR=smooth_method,
          "Wavelet" = {
            smoothf = function(thresh.scale){
              function(x) wavSmoothing(list(x), "Wavelet", thresh.scale)[[1]]
            }
            smoothf = smoothf(thresh.scale)
          },
          "SavitzkyGolay" = {
            smoothf = function(hws_smooth){
              function(x) smoothIntensity(x, "SavitzkyGolay", hws_smooth)
            }
            smoothf = smoothf(hws_smooth)
          }
  )

  if (is.null(peptides_user)) peptides_user = peptides

  masses = pull(peptides_user, 3)

  npepts = length(masses)
  spectra_f = list.files(indir)

  masses = matrix(masses, nrow = n_isopeaks, ncol = length(masses), byrow = T)
  d = 1.00235
  masses = masses + (d * 0L:(n_isopeaks - 1L))
  masses = sort(masses)


  # Load function
  prepf = prepFun(
    preprocess_spectra, smoothf = smoothf, iterations = iterations,
    halfWindowSize = halfWindowSize, SNR = SNR, masses = masses, tol = tol,
    n_isopeaks = n_isopeaks, min_isopeaks = min_isopeaks)

  iso_peaks = preprocessData(
    indir=indir, readf = readf, sep = sep, nchunks = nchunks, prepf = prepf,
    spectra = spectra, ncores = ncores, iocores = iocores, vch=vch)

  if (!is.null(outdir)){
    filenames = sub("mzML|tab|tsv", "txt", names(iso_peaks))
    invisible(mapply(
      function(x, f, outdir){
        outfile = file.path(outdir, f)
        write_tsv(x, outfile, col_names = F)
      },
      iso_peaks, filenames, MoreArgs = list(outdir=outdir)
    ))
  }
  return(iso_peaks)
}


#' Get isotopic peaks from reference dataset
#' @param list_params logical. If TRUE returns a dataframe of possible parameters
#' for which iso peaks has been precomputed for the reference dataset
#' @param params list of parameters and values for which the pre-computed iso_peaks is returned
#' Parameters supported: SNR, iterations, hws_smooth, halfWindowSize
#' @return A list of dataframes, 1 per sample. Each dataframe has 3 columns,
#' m/z, intensity and signal-to-noise ratio for each of the n_isopeaks from each
#' peptide. Missing peaks are NAs.
#' @export
#' @details It retrieves the isotopic envelopes for the default peptides for the
#' Orval Abbey dataset
#' @references
#' Bethencourt, J.H. et al. (2022) ‘Data from “A biocodicological analysis of the medieval library and archive from Orval abbey, Belgium”’, Journal of open archaeology data, 10(0). doi:10.5334/joad.89.
#' @examples
get_ref_isopeaks = function(which_params = F, params=NULL){
  iso_peaks_params = lapply(
    iso_peaks_list,
    function(x) x$params
  )
  iso_peaks_params = do.call(bind_rows, iso_peaks_params)

  if (which_params) {
    return(iso_peaks_params)
  }
  if (is.null(params)) {
    stop('Specify params')
  }
  if (is.list(params)) params = data.frame(params)
  iso_peaks_params$idx = 1:nrow(iso_peaks_params)
  params_use = merge(iso_peaks_params, params, all=F)

  return(iso_peaks_list[[params_use$idx]])
}

