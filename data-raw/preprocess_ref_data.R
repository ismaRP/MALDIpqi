library(MALDIpqi)
library(MALDIutils)
library(tidyverse)
library(MALDIquant)
library(ggrepel)


orvaldata = '/home/ismael/palaeoproteomics/MALDI/orval/data'


# Run all peptides
iso_peaks_orval = getIsoPeaks(
  indir=file.path(orvaldata, 'orval_dupl_mzML'), outdir=NULL, readf="mzml",
  nchunks = 50, peptides_user = example_peptides, n_isopeaks = 5, min_isopeaks = 4, SNR = 1.5,
  smooth_method = "SavitzkyGolay", hws_smooth = 4, halfWindowSize = 20,
  ncores = 6)

orval_metadata = read_csv(file.path(orvaldata, '../tables/orval_spectra_metadata.csv'))

orval_metadata = orval_metadata %>%
  group_by(sample_name) %>% mutate(n_replicates = n()) %>%
  filter(n_replicates==3)

orval_samples = paste0(
  orval_metadata$sample_name, '_',
  orval_metadata$replicate
)


iso_peaks_orval = iso_peaks_orval[orval_samples]
