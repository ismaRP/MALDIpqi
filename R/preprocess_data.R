
#' Preprocessing spectra for q2e estimation
#'
#' Performs smoothening, baseline removal and peak detection on MALDI samples.
#' From the peaks, isotopic peaks for a list of peptides are extracted.
#'
#' @param indir Folder containing spectra in mzML format.
#' @param metadata Data frame with spectra metadata with at least \code{file}
#' column. Ideally metadata has been cleaned before with [MALDIzooMS::clean_metadata]
#' @param mono_masses
#' Array with the peptides monoisotopics masses
#' @param smooth_wma_hws
#' Half-window size for WeightedMovingAverage smoothing method
#' @param smooth_sg_hws
#' Half-window size for SavitzkyGolay smoothing method
#' @param iterations
#' Iterations parameter for baseline detection.
#' @param halfWindowSize
#' Half-window size parameter for local maximum detection.
#' @param snr
#' Signal-to-noise threshold above which peaks are considered
#' @param k
#' k parameter for [MsCoreUtils::refineCentroids()]
#' @param threshold
#' threshold parameter for [MsCoreUtils::refineCentroids()]
#' @param local_bg
#' Whether to further to clean peaks of lists by modelling the local
#' background noise. See [MALDIzooMS::peaks_local_bg].
#' Ideally should work with a \code{snr} threshold of 0.
#' \code{mass_range}, \code{bg_cutoff} and \code{l_cutoff} only applied if \code{local_bg} is TRUE
#' @param mass_range
#' Mass window to both sides of a peak to be considered for backgroun modelling
#' @param bg_cutoff
#' The peaks within the mass range with intensity below the \code{bg_cutoff} quantile
#' are considered for background modelling. \code{bg_cutoff=1} keeps all peaks
#' and \code{bg_cutoff=0.5} would only keep the bottom half.
#' @param l_cutoff
#' Likelihood threshold or p-value. Peaks with a probability of being modelled as
#' background noise higher than this are filtered out.
#' @param n_isopeaks
#' Number of isotopic peaks to pick. Default is 5 and the maximum permitted.
#' @param min_isopeaks
#' If less than min_isopeaks consecutive (about 1 Da difference) isotopic peaks
#' are detected, the whole isotopic envelope is discarded. Default is 4
#' @param tolerance
#' Mass tolerance in Da between \code{mono_masses} and subsequent isotopic peaks
#' and detected peaks. See [MsCoreUtils::closest]
#' @param ppm
#' Parts-per-million added to tolerance. See [MsCoreUtils::closest]
#' @param ncores
#' Number of cores used by the [Spectra::MsBackendMzR] backend in [Spectra::peaksData]
#' @param chunk_size
#' @return A list of dataframes, 1 per sample. Each dataframe has 3 columns,
#' m/z, intensity and signal-to-noise ratio for each of the n_isopeaks from each
#' peptide. Missing peaks are NAs.
#'
#' @importFrom MALDIzooMS get_spectra_name clean_metadata
#' @importFrom MALDIzooMS smooth baseline_correction peak_detection
#' @importFrom MALDIzooMS peptide_pseudo_clusters peaks_local_bg
#' @importFrom parallel detectCores
#' @importFrom Spectra Spectra addProcessing peaksData
#' @importFrom Spectra MsBackendMzR processingChunkSize
#' @importFrom BiocParallel MulticoreParam
#' @importFrom dplyr rename bind_cols
#' @export
#' @details The default peptides are the ones from Nair et al. (2022).
#' The paper contains the details on the preprocessing procedure.
#' @references
#' Nair, B. et al. (2022) ‘Parchment Glutamine Index (PQI): A novel method to estimate glutamine deamidation levels in parchment collagen obtained from low-quality MALDI-TOF data’, bioRxiv. doi:10.1101/2022.03.13.483627.
#'
preprocess_spectra = function(
    indir, metadata,
    make_plots = FALSE,
    mono_masses = NULL,
    smooth_wma_hws = 4,
    smooth_sg_hws = 6,
    iterations = 50,
    halfWindowSize = 20,
    snr = 2, k = 0L, threshold = 0.33,
    local_bg = FALSE,
    mass_range=100, bg_cutoff=0.5, l_cutoff=1e-8,
    tolerance = 0.4, ppm=50,
    n_isopeaks = 5,
    min_isopeaks = 4,
    ncores = NULL, chunk_size=40){

  if (is.null(mono_masses)) mono_masses = peptides$mass
  if (is.null(ncores)) ncores = detectCores() - 2

  metadata = clean_metadata(metadata, indir)
  mzml_files = file.path(indir, metadata$file)

  sps_mzr = suppressMessages(
    Spectra(mzml_files, source = MsBackendMzR(), centroided = FALSE,
            BPPARAM = MulticoreParam(workers = ncores)))
  processingChunkSize(sps_mzr) = chunk_size

  sps_mzr$spectra_name = get_spectra_name(sps_mzr@backend@spectraData$dataOrigin)

  # Weighted Moving Average Smoothing
  sps_mzr = addProcessing(
    sps_mzr, MALDIzooMS::smooth, method = 'WeightedMovingAverage',
    hws = smooth_wma_hws, int_index = 'intensity', in_place = FALSE)
  # Savitzky-Golay Filter smoothing
  sps_mzr = addProcessing(
    sps_mzr, MALDIzooMS::smooth, method = 'SavitzkyGolay',
    hws = smooth_sg_hws, int_index = 'intensity', in_place = FALSE)

  if (make_plots) {
    # Baseline estimation on MA smoothed
    sps_mzr = addProcessing(
      sps_mzr, MALDIzooMS::baseline_correction, int_index = 'intensity_WeightedMovingAverage',
      keep_bl = TRUE, substract_index = 'intensity_SavitzkyGolay', in_place = FALSE,
      method = 'SNIP', iterations = iterations, decreasing = TRUE)
    # Get spectra
    sp = peaksData(sps_mzr)
    names(sp) = sps_mzr@backend@spectraData$spectra_name
    sp = as.data.frame(sp) %>% rename(spectra_name = group_name)
    sp = bind_cols(sp, separate_sample_replicate(sp$spectra_name, sep = '_'))
    # PEAKS
    sps_mzr = addProcessing(
      sps_mzr, MALDIzooMS::peak_detection, halfWindowSize = 20L,
      method = 'SuperSmoother', snr = snr, k = k, threshold = threshold,
      descending = TRUE, int_index = 'intensity_SavitzkyGolay_bl_corr_SNIP',
      add_snr=TRUE)
  } else {
    # Baseline estimation on MA smoothed
    sps_mzr = addProcessing(
      sps_mzr, MALDIzooMS::baseline_correction, int_index = 'intensity_WeightedMovingAverage',
      keep_bl = FALSE, substract_index = 'intensity_SavitzkyGolay', in_place = TRUE,
      method = 'SNIP', iterations = iterations, decreasing = TRUE)
    # PEAKS
    sps_mzr = addProcessing(
      sps_mzr, MALDIzooMS::peak_detection, halfWindowSize = 20L,
      method = 'SuperSmoother', snr = snr, k = k, threshold = threshold,
      descending = TRUE, int_index = 'intensity_SavitzkyGolay',
      add_snr=TRUE)
  }

  if (local_bg) {
    sps_mzr = addProcessing(
      sps_mzr, peaks_local_bg, mass_range = mass_range, bg_cutoff = bg_cutoff,
      l_cutoff = l_cutoff, int_index = 'intensity_SavitzkyGolay')
  }

  # GET ISOTOPIC CLUSTERS
  sps_mzr = addProcessing(
    sps_mzr, peptide_pseudo_clusters,
    mono_masses = mono_masses, n_isopeaks = n_isopeaks, min_isopeaks = min_isopeaks,
    tolerance = tolerance, ppm = ppm)


  peaks = peaksData(sps_mzr)
  names(peaks) = sps_mzr@backend@spectraData$spectra_name
  if (make_plots) {
    return(list(sp, peaks))
  } else {
    return(peaks)
  }

}

#' Normalize vector of intensities by the maximum
#'
#' @param intensity Vector of intensities
#'
#' @return
#' @export
#'
#' @examples
normalize_intensity = function(intensity) {
  if (all(is.na(intensity))) {
    norm_int = rep(NA, length(intensity))
  } else {
    norm_int = intensity/max(intensity, na.rm = TRUE)
  }
  return(norm_int)
}

#' Prepare list of peaks data into a data.frame
#' @param peaks List of peaks matrix. Names are used as spectra name
#' @param n_isopeaks Number of isotopic peaks
#' @param peptides_user
#' A dataframe with peptide information. It must contain at least 3 columns,
#' \code{pep_number} or ID, \code{mass}, \code{sequence} and \code{h_hyp} (# of hydroxyprolines).
#' IF NULL, default are used, see details.
#'
#' @return A data.frame with isotopic peaks detected from data, theoretical isotopic
#' envelopes for 1 and 2 deamidations and other associated data.
#' @details The default peptides are the ones from Nair et al. (2022).
#' The paper contains the details on the preprocessing procedure.
#'
#' @importFrom dplyr mutate rename
#' @export
#'
#' @examples
prepare_peaks = function(peaks, peptides_user, n_isopeaks) {

  peaks = as.data.frame(peaks) %>%
    rename(spectra_name = group_name)

  iso_peps = get_isodists(
    peptides_user$sequence, 2, peptides_user$n_hyp,
    long_format = T)

  n_spectra = length(unique(peaks$spectra_name))
  n_peptides = nrow(peptides_user)
  peaks = peaks %>%
    mutate(
      mass_pos = as.factor(rep(
        seq(n_isopeaks),
        times = n_spectra * n_peptides)),
      pep_idx = as.factor(rep(
        rep(seq(n_peptides), each = n_isopeaks),
        times = n_spectra)),
      pep_number = as.factor(rep(
        rep(peptides_user$pept_number, each = n_isopeaks),
        times = n_spectra)),
      ndeam = rep(
        rep(str_count(peptides_user$sequence, 'Q'), each = n_isopeaks),
        times = n_spectra),
      weight = SNR/abs(delta_mass))

  if ('pept_name' %in% colnames(peptides_user)) {
    peaks = peaks %>%
      mutate(pep_name = as.factor(rep(
        rep(peptides_user$pept_name, each = n_isopeaks),
        times = n_spectra)))
  }

  # peaks = cbind(peaks, separate_sample_replicate(peaks$spectra_name))
  peaks = bind_cols(
    peaks, separate_sample_replicate(peaks$spectra_name, sep = '_'))

  peaks = Reduce(function(x, y) merge(x, y, all = TRUE, by = c('pep_idx', 'mass_pos')),
                 list(peaks, iso_peps))
  peaks = peaks %>%
    arrange(sample, replicate, pep_idx, mass_pos) %>%
    group_by(spectra_name, pep_idx) %>%
    mutate(norm_int = normalize_intensity(intensity),
           n_peaks = sum(!is.na(intensity))) %>%
    ungroup() %>%
    mutate(intensity = replace(intensity, is.na(intensity), 0),
           norm_int = replace(norm_int, is.na(norm_int), 0)) %>%
    mutate(intensity = replace(intensity, n_peaks == 0, NA),
           norm_int = replace(norm_int, n_peaks == 0, NA))

  return(peaks)

}



#' Title
#'
#' @param x
#' @param n_isopeaks
#' @param min_isopeaks
#'
#' @return
#'
#' @examples
calc_n_frac_peaks = function(x, n_isopeaks, min_isopeaks) {
  fracs = list()
  for (n in c(0, min_isopeaks:n_isopeaks)) {
    fracs[[paste0('frac_', n)]] = sum(x == n)/length(x)
  }
  return(data.frame(fracs))
}

#' Title
#'
#' @param peaks
#' @param n_isopeaks
#' @param min_isopeaks
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_n_peaks_per_peptide = function(peaks, n_isopeaks, min_isopeaks, ...) {

  a = peaks %>%
    group_by(pep_number, ...) %>%
    summarise(calc_n_frac_peaks(n_peaks, n_isopeaks, min_isopeaks)) %>%
    ungroup() %>%
    pivot_longer(cols = starts_with('frac'), names_to='n_of_peaks',
                 names_prefix = 'frac_', values_to='fraction') %>%
    mutate(n_of_peaks = factor(as.integer(n_of_peaks), levels=c(n_isopeaks:min_isopeaks,0)),
           pep_number = as.factor(pep_number))

  ggplot(a) +
    geom_col(aes(x=pep_number, y=fraction, fill=n_of_peaks)) +
    facet_wrap(vars(...))
}



#' Title
#'
#' @param sp
#' @param peaks
#' @param peptides_user
#' @param n_isopeaks
#'
#' @return
#' @export
#'
#' @examples
plot_preprocessing = function(sp, peaks, peptides_user, n_isopeaks) {
  # Prepare peaks
  peaks = prepare_peaks(peaks, peptides_user, n_isopeaks = n_isopeaks)

  peaks_mask = list()
  sp_mask = rep(NA, nrow(sp))
  for (i in seq_along(peptides_user$mass)) {
    mono_mz = peptides_user$mass[i]
    s = (sp$mz > (mono_mz-2)) & (sp$mz < (mono_mz+6))
    sp_mask[s] = peptides_user$pept_number[i]
    p = (peaks$mz > (mono_mz-2)) & (peaks$mz < (mono_mz+6))
    peaks_mask[[i]] = p
  }
  peaks_mask = Reduce('|', peaks_mask)
  peaks = peaks[peaks_mask,]

  sp$pep_number = sp_mask
  sp = sp[!is.na(sp$pep_number),]

  sp = sp %>%
    group_by(spectra_name, pep_number) %>%
    mutate(norm_int = intensity/max(intensity, na.rm = TRUE),
           norm_int_wma = intensity_SavitzkyGolay/max(intensity, na.rm = TRUE),
           norm_int_bl = baseline_SNIP/max(intensity, na.rm = TRUE),
           norm_int_bl_corr = intensity_SavitzkyGolay_bl_corr_SNIP/max(intensity, na.rm = TRUE)) %>%
    ungroup()

  spp = ggplot(sp) +
    # Raw int
    geom_line(aes(x = mz, y = norm_int), color='black', alpha=0.4) +
    # Smooth int
    geom_line(aes(x = mz, y = norm_int_wma),
              color = 'red', alpha = 0.7) +
    # Baseline
    geom_line(aes(x = mz, y = norm_int_bl),
              color = 'darkgreen', linetype = "dashed") +
    geom_line(aes(x = mz, y = norm_int_bl_corr),
              color = 'cadetblue1', linetype = "solid") +
    # geom_line(aes(x=mz, y=b_d), color='blue') +
    # geom_text(aes(label=QCflag), x=+Inf, y=+Inf, vjust=1.3, hjust=1.2,
    #           data=sps_clusters[sele]@backend@spectraData) +
    geom_point(aes(x = mz, y = norm_int), shape = 19, size = 2,
               data = peaks) +
    facet_grid(spectra_name~pep_number, scales = 'free')
  return(spp)
}
