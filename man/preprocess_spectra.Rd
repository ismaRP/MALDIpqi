% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess_data.R
\name{preprocess_spectra}
\alias{preprocess_spectra}
\title{Preprocessing spectra for q2e estimation}
\usage{
preprocess_spectra(
  indir = NULL,
  metadata = NULL,
  mzml_files = NULL,
  spectrum_name_file = FALSE,
  sps_mzr = NULL,
  make_plots = FALSE,
  peptides_user = NULL,
  smooth_wma_hws = 4,
  smooth_sg_hws = 6,
  iterations = 50,
  halfWindowSize = 20,
  snr = 2,
  k = 0L,
  threshold = 0.33,
  local_bg = FALSE,
  mass_range = 100,
  bg_cutoff = 0.5,
  l_cutoff = 1e-08,
  tolerance = 0.4,
  ppm = 50,
  n_isopeaks = 5,
  min_isopeaks = 4,
  norm_func = NULL,
  q2e = NULL,
  ncores = NULL,
  chunk_size = 40,
  verbose = FALSE
)
}
\arguments{
\item{indir}{Folder containing spectra in mzML format.}

\item{metadata}{Data frame with spectra metadata with at least \code{file}
column. Ideally metadata has been cleaned before with \link[MALDIzooMS:clean_metadata]{MALDIzooMS::clean_metadata}}

\item{mzml_files}{Paths to mzML files}

\item{sps_mzr}{Spectra object}

\item{smooth_wma_hws}{Half-window size for WeightedMovingAverage smoothing method}

\item{smooth_sg_hws}{Half-window size for SavitzkyGolay smoothing method}

\item{iterations}{Iterations parameter for baseline detection.}

\item{halfWindowSize}{Half-window size parameter for local maximum detection.}

\item{snr}{Signal-to-noise threshold above which peaks are considered}

\item{k}{k parameter for \code{\link[MsCoreUtils:refineCentroids]{MsCoreUtils::refineCentroids()}}}

\item{threshold}{threshold parameter for \code{\link[MsCoreUtils:refineCentroids]{MsCoreUtils::refineCentroids()}}}

\item{local_bg}{Whether to further to clean peaks of lists by modelling the local
background noise. See \link[MALDIzooMS:peaks_local_bg]{MALDIzooMS::peaks_local_bg}.
Ideally should work with a \code{snr} threshold of 0.
\code{mass_range}, \code{bg_cutoff} and \code{l_cutoff} only applied if \code{local_bg} is TRUE}

\item{mass_range}{Mass window to both sides of a peak to be considered for backgroun modelling}

\item{bg_cutoff}{The peaks within the mass range with intensity below the \code{bg_cutoff} quantile
are considered for background modelling. \code{bg_cutoff=1} keeps all peaks
and \code{bg_cutoff=0.5} would only keep the bottom half.}

\item{l_cutoff}{Likelihood threshold or p-value. Peaks with a probability of being modelled as
background noise higher than this are filtered out.}

\item{tolerance}{Mass tolerance in Da between \code{mono_masses} and subsequent isotopic peaks
and detected peaks. See \link[MsCoreUtils:matching]{MsCoreUtils::closest}}

\item{ppm}{Parts-per-million added to tolerance. See \link[MsCoreUtils:matching]{MsCoreUtils::closest}}

\item{n_isopeaks}{Number of isotopic peaks to pick. Default is 5 and the maximum permitted.}

\item{min_isopeaks}{If less than min_isopeaks consecutive (about 1 Da difference) isotopic peaks
are detected, the whole isotopic envelope is discarded. Default is 4}

\item{norm_func}{Function to normalize the isotopic distribution}

\item{q2e}{If provided, it adds the theoretical isotopic distribution of peptides with
this extent of deamidation}

\item{ncores}{Number of cores used by the \link[Spectra:MsBackend]{Spectra::MsBackendMzR} backend in \link[Spectra:peaksData]{Spectra::peaksData}}

\item{spectrum_file_name}{If mzml_files are provided, whether to use file names
as spectra names. Otherwise, it is assumed the the spectra IDs are in the mzML
files' headers.}

\item{mono_masses}{Array with the peptides monoisotopics masses}
}
\value{
A list of dataframes, 1 per sample. Each dataframe has 3 columns,
m/z, intensity and signal-to-noise ratio for each of the n_isopeaks from each
peptide. Missing peaks are NAs.
}
\description{
Performs smoothening, baseline removal and peak detection on MALDI samples.
From the peaks, isotopic peaks for a list of peptides are extracted.
}
\details{
Provide the input data either using \code{metadata} and \code{indir}, or provide paths
with \code{mzml_files}. You can also provide a \code{Spectra} object directly in \code{sps_mzr}.
If data is provided using more than one of the options, the \code{sps_mzr} is used, and then the \code{mzml_files}.

The default peptides are the ones from Nair et al. (2022).
The paper contains the details on the preprocessing procedure.
}
\references{
Nair, B. et al. (2022) ‘Parchment Glutamine Index (PQI): A novel method to estimate glutamine deamidation levels in parchment collagen obtained from low-quality MALDI-TOF data’, bioRxiv. doi:10.1101/2022.03.13.483627.
}
