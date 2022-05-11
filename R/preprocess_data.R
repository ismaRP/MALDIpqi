


#' Title
#'
#' @param x
#' @param n
#'
#' @return n list of vectos
#' @export
#'
#' @examples
chunks = function(x,n){
  split(x, cut(seq_along(x), n, labels = FALSE))
}

#' Title
#'
#' @param x
#' @param indir
#'
#' @return TRUE of FALSE for file not being empty
#'
#' @examples
detect_empty = function(x, indir){
  f = file.path(indir, x)
  t = readLines(f, 1)
  if (length(t)>0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}



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




#' getIsoPeaks
#'
#' Performs smoothening, baseline removal and peak detection on MALDI samples.
#'
#' @param indir Folder containing spectra.
#' @param peptides
#' A dataframe with peptide information. It must contain at least 3 columns,
#' peptide number or ID, name, and m/z. IF NULL
#' @param readf
#' A string value. choose function to use to read spectra.
#' Currently restricted to one of "fread", "table" or "mzml"
#' @param outdir
#' Folder where one peaks file for each sample should be written.
#' If NULL (default), no files are written.
#' @param chunks
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
#' peptide.
#'
#' @importFrom MALDIutils preprocessData prepFun
#' @importFrom MALDIrppa wavSmoothing
#' @importFrom parallel mcmapply detectCores mclapply
#' @importFrom readr write_tsv
#' @importFrom dplyr pull
#' @export
#'
#'
#' @examples
getIsoPeaks = function(indir,
                       readf = c("fread", "table", "mzml"),
                       sep="\t",
                       peptides=NULL,
                       outdir = NULL,
                       chunks = 50,
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
                       iocores = 1){

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

  if (is.null(peptides)) {
    peptides = load("data/peptides.rda")
  }
  pept_names = pull(peptides, 2)
  masses = pull(peptides, 3)

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
    indir=indir, readf = readf, sep = sep, chunks = chunks, prepf = prepf,
    ncores = ncores, iocores = iocores)

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

