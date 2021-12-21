


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



#' Title
#'
#' @param x
#' @param ...
#'
#' @return data.frame with mass, intensity and s2n values for pseudo-isotopic
#'         clusters
#' @export
#' @importFrom MALDIquant match.closest
#'
#' @examples
extractPeaks = function(x, ...){
  prep_args = list(...)
  masses = prep_args$masses
  tol = prep_args$tol
  n_isopeaks = prep_args$n_isopeaks
  min_isopeaks = prep_args$min_isopeaks
  n_pepts = length(masses)/n_isopeaks

  sel_masses = rep(NA, length(masses))
  sel_int = rep(NA, length(masses))
  sel_snr = rep(NA, length(masses))
  sel_complete = rep(F, length(masses))
  idx = match.closest(x@mass, masses, tolerance = tol*masses)

  sel_complete[idx[!is.na(idx)]] = T

  sel_complete = split(sel_complete,
                       cut(seq_along(sel_complete), n_pepts, labels = FALSE))

  sel_complete = lapply(
    sel_complete,
    function(x){
      s = cumsum(x)
      first_missing = (s[1] == 0) & (s[min_isopeaks + 1] == min_isopeaks)
      last_missing = s[min_isopeaks] == min_isopeaks
      if (!(first_missing | last_missing)){
        return(rep(F, n_isopeaks))
      } else {
        return(x)
      }
    }
  )

  sel_complete = Reduce(c, sel_complete)
  sel_masses[idx[!is.na(idx)]] = x@mass[!is.na(idx)]
  sel_masses[!sel_complete] = NA
  sel_int[idx[!is.na(idx)]] = x@intensity[!is.na(idx)]
  sel_int[!sel_complete] = NA
  sel_snr[idx[!is.na(idx)]] = x@snr[!is.na(idx)]
  sel_snr[!sel_complete] = NA
  return(data.frame(mass=sel_masses, intensity=sel_int, snr=sel_snr))
}


#' Title
#'
#' @param s
#' \code{\link[MALDIquant]{MassSpectrum}} object
#' @param ...
#' Parameters passed to the different preprocessing functions.
#' See \code{\link{preprocess_data}}
#'
#' @return peaks list, with mass, intensity and s2n
#' @importFrom MALDIquant smoothIntensity removeBaseline detectPeaks
#' @export
#'
#' @examples
preprocess_spectra = function(s, ...){
  prep_args = list(...)
  # s = smoothIntensity(s, method="SavitzkyGolay",
  #                     halfWindowSize=prep_args$hws_smooth)
  s = prep_args$sm(s)
  s = removeBaseline(s, method='SNIP', iterations=prep_args$iter)
  # Estimate noise
  p = detectPeaks(
    s, halfWindowSize=prep_args$hws_peak,
    method='SuperSmoother', SNR=prep_args$snr
  )
  # Extract peaks
  p = extractPeaks(p, ...)
  return(p)
}




#' preprocess_data
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
#' @param ncores
#' Number of cores used for the preprocessing of spectra. The cores will work in
#' parallel with the different spectra within a chunk.
#' @param iocores
#' Number of cores used for I/O operations. For some systems I/O operations involving
#' multiple cores reading or writing at the same time from disk increases time.
#'
#' @return A list of dataframes, 1 per sample. Each dataframe has 3 columns,
#' m/z, intensity and signal-to-noise ratio for each of the n_isopeaks from each
#' peptide.
#'
#' @importFrom MALDIutils importTsv importTable import_file.MzMl
#' @importFrom MALDIrppa wavSmoothing
#' @importFrom parallel mcmapply detectCores mclapply
#' @importFrom readr write_tsv
#' @importFrom dplyr pull
#' @export
#'
#'
#' @examples
preprocess_data = function(indir,
                           readf = c("fread", "table", "mzml"),
                           peptides=NULL,
                           outdir = NULL,
                           chunks = 50,
                           smooth_method = c("SavitzkyGolay","Wavelet"),
                           thresh.scale = 2.5,
                           hws_smooth = 20,
                           iter = 20,
                           hws_peak = 20,
                           snr = 0,
                           tol = 1.5e-4,
                           n_isopeaks = 6,
                           min_isopeaks = 4,
                           ncores = NULL,
                           iocores = 1){
  if (is.null(ncores)){
    ncores = detectCores() - 2
  } else {
    ncores = min(detectCores() - 2, ncores)
  }

  cat("Using ", ncores, " cores")
  readf = match.arg(readf)
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
  }
  pept_names = pull(peptides, 2)
  masses = pull(peptides, 3)

  npepts = length(masses)
  spectra_f = list.files(indir)

  # Filter out empty files
  filter_empty = lapply(
    spectra_f,
    detect_empty,
    indir
  )
  filter_empty = unlist(filter_empty)
  spectra_f = spectra_f[filter_empty]

  if (chunks > 1) {
    spectra_chunks = chunks(spectra_f, chunks)
  } else {
    spectra_chunks = list(spectra_f)
  }
  names(spectra_chunks) = NULL

  masses = matrix(masses, nrow = n_isopeaks, ncol = length(masses), byrow = T)
  d = 1.00235
  masses = masses + (d * 0L:(n_isopeaks - 1L))
  masses = sort(masses)

  iso_peaks = mapply(
    function(x, ch){
      if (ch %% 5 == 0){
        print(sprintf('Chunk %i of %i', ch, chunks))
      }
      infiles = file.path(indir, x)
      l = mclapply(
        infiles,
        read_f,
        mc.cores=iocores,
        mc.silent=T
      )
      names(l) = x
      invisible(mcmapply(
        preprocess_spectra, l,
        MoreArgs = list(
          tol=tol, n_isopeaks=n_isopeaks, min_isopeaks=min_isopeaks, snr=snr,
          sm=sm, iter=iter, hws_peak=hws_peak, masses=masses),
        mc.cores=ncores, SIMPLIFY = F
      ))
    },
    spectra_chunks,
    seq_along(spectra_chunks),
    SIMPLIFY = F
  )
  # Unlist chunks, so all spectra are in one list of depth=1
  iso_peaks = unlist(iso_peaks, recursive = F, use.names = T)
  filenames = sub("mzML|tab|tsv", "txt", names(iso_peaks))
  if (!is.null(outdir)){
    # TODO: write files
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

