

#' Helper function to prepare data for q2e estimation
#'
#' @param peptides
#' @param outdir
#' @param data_list
#' @param indir
#' @param n_isopeaks
#'
#' @importFrom readr write_csv
#' @return
#'
#' @examples
prepare_data = function(peptides, outdir=NULL, data_list=NULL, indir=NULL,
                        n_isopeaks=5){

  #if it's a list of samples with peaks
  if (!is.null(data_list)){
    samples_repl = names(data_list) # file names as list names
    samples_repl = sort(samples_repl)
    data_list = data_list[samples_repl]
  } else if (!is.null(indir)) { # else peaks given as files in directory
    if (dir.exists(indir)){
      samples_repl = list.files(indir)
      samples_repl = sort(samples_repl)
    } else {
      stop(paste0(indir, " file does not exist."))
    }
  } else {
    stop("Either peaks or indir need to be provided")
  }

  # Remove extension and replicate number
  samples = strsplit(samples_repl, "_")
  samples = lapply(
    samples,
    function(x) paste0(x[1:length(x)-1], collapse="_")
  )
  samples = unlist(samples)
  sample_names = data.frame(sample=samples, replicate=samples_repl)
  # Set new peaks names
  # if (!is.null(data_list)) names(data_list) = samples

  n_datasets = nrow(sample_names)
  mismatch = rep(FALSE, n_datasets)

  for (i in 1:(n_datasets/3)) {
    d = sample_names$sample[(i-1)*3+1:3]

    # if (length(unlist(d))==3) d <- strsplit(txt[(i-1)*3+1:3,1],"-")
    # if (length(unlist(d))==3) d <- strsplit(txt[(i-1)*3+1:3,1],"\\.")

    if (d[[1]] != d[[2]]) mismatch[(i-1)*3+1:3] = TRUE
    if (d[[1]] != d[[3]]) mismatch[(i-1)*3+1:3] = TRUE
    if (d[[2]] != d[[3]]) mismatch[(i-1)*3+1:3] = TRUE
  }
  # Remove incomplete samples with missing replicates
  sample_names = sample_names[!mismatch, ]

  # Check if number of data sets divisble by 3
  n_datasets = nrow(sample_names)
  if (n_datasets%%3!=0) stop("UPS, number of data sets not divisble by 3")

  if (is.null(data_list)){
    # Read samples, only those that have the 3 replicates
    # data_list = lapply(
    #   sample_names$replicate,
    #   read_peaks,
    #   data_path=indir)
    data_list = mapply(
      read_peaks,
      sample_names$replicate,
      seq_along(sample_names$replicate),
      MoreArgs = list(data_path=indir),
      SIMPLIFY = F
    )
  } else{
    # Filter out replicate-incomplete samples
    data_list = data_list[sample_names$replicate]
  }
  # peptides = read_csv(pept_f, quote = "\'") %>% arrange(mass)
  n_peptides = nrow(peptides) # n of isotopics peaks to include in wls
  max_isopeaks = dim(data_list[[1]])[1]/n_peptides # n of isotopic peaks in data
  if (max_isopeaks < n_isopeaks) {
    stop(sprintf("Data contains %i isotopic peaks, but attempted to use %i",
                 max_isopeaks, n_isopeaks))
  }

  # Output data to file
  # Lines with NA will be added for the three missing replicates
  d = n_datasets*n_peptides
  data_out = data.frame(sample = rep("",d), peptide = rep(NA,d))
  data_out = cbind(data_out, array(NA, dim=c(d, 2*n_isopeaks)))

  # Add column names
  cnames <- rep("", 2*n_isopeaks)
  for (i in 1:n_isopeaks) cnames[i] = paste("peak", i, sep="")
  for (i in 1:n_isopeaks) cnames[i+n_isopeaks] = paste("noise", i, sep="")
  names(data_out)[2+1:(2*n_isopeaks)] = cnames


  # Extract data
  for (i in 1:nrow(sample_names)) {
    if (i %% 500 == 0){
      print(sprintf('Peaks extracted for %i samples', i))
    }
    data = data_list[[i]]
    indx = i
    for (j in 1:n_peptides) {
      # peak values
      data_out[(indx-1)*n_peptides+j, 2+1:n_isopeaks] = data[(j-1)*max_isopeaks+1:n_isopeaks, 2]
      # noise for each peak value
      data_out[(indx-1)*n_peptides+j, 2+(n_isopeaks+1):(2*n_isopeaks)] = data[(j-1)*max_isopeaks+1:n_isopeaks, 2] / data[(j-1)*max_isopeaks+1:n_isopeaks, 3]
      data_out[(indx-1)*n_peptides+j,2] = peptides[[j,1]]
      data_out[(indx-1)*n_peptides+j,1] = sample_names$sample[i]
    }
  }

  # data_trunc = data_out
  data_out[,-(1:2)] = round(data_out[,-c(1:2)],3)


  # Print data
  if (!is.null(outdir)){
    write_csv(
      data_out,
      file=file.path(outdir, "intensities.txt")
    )
  }

  return(list(data_out, sample_names))

}


#' Estimation of peptide q2e
#'
#' Given a list of tables, each contaning the peptide isotopic peaks per samples,
#' this function uses weighted least squares to estimate q2e
#' @param peptides A dataframe with peptide information. It must contain at least 3 columns,
#' peptide number or ID, name, and m/z. IF NULL, default are used, see details.
#' @param data_list list of isopeaks obtained from getIsoPeaks
#' @param indir Path to lsit of isopeaks file is it was saved in a file
#' @param n_isopeaks Number of isotopic peaks in isopeaks list
#' @param outdir Optional file were weighte least square estimates of q2e are saved
#'
#' @return A data.frame with gamma and q2e estimates and
#' minimum least square value of the estimate per sample, replicate and peptide
#' @importFrom readr write_csv
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_count
#' @export
#'
#' @examples
wls_q2e = function(peptides=NULL, data_list=NULL, indir=NULL,
                   n_isopeaks=5, outdir=NULL){

  if (is.null(peptides)) {
    peptides = load("data/peptides.rda")
  }

  if (is.null(data_list) & is.null(indir)) {
    stop("Either peaks or indir need to be provided")
  }

  prep_data = prepare_data(
    data_list = data_list, outdir=outdir, indir=indir,
    peptides = peptides, n_isopeaks = n_isopeaks)

  sample_names = prep_data[[2]]
  peaks_data = prep_data[[1]]

  n_peptides = nrow(peptides)
  n_datasets = nrow(peaks_data)/n_peptides

  # number of parameters to estimate (one plus number of deamidataion sites)
  lgth_gamma = str_count(pull(peptides, 4), 'Q') + 1
  ## estimates of the lgth_gammma parameters g1,...,gk
  ## ratios g1/(g1+...+gk), and min LS
  gammas = array(NA, dim=c(n_datasets, n_peptides, max(lgth_gamma)+2))

  # read data
  # for each data set (replicate within a sample), there is a line for each peptide
  # for each line, there is a number of peak values and a noise estimate for each peak value
  # peak values, columns 1:(dim(data)[2]/2)
  # noise, columns (dim(data)[2]/2+1):dim(data)[2]

  # read theoretical isotope peaks
  isodists = get_isodists(pull(peptides, 4), pull(peptides, 5))
  # use as many isotope peaks as there are peak values (excl id and peptide number)
  # nisotopes <- (dim(data)[2]-2)/2
  # check nisotopes smaller or equal to dim(isotopes)[2]
  # if (nisotopes>dim(isotopes)[2]) print("UPS, something is wrong with the number of peak values")
  max_isopeaks = (dim(peaks_data)[2]-2)/2
  if (max_isopeaks < n_isopeaks){
    stop(sprintf("Data contains %i isotopic peaks, but attempted to use %i",
                 max_isopeaks, n_isopeaks))
  }

  for (i in 1:n_datasets){
    if (i %% 500 == 0){
      print(sprintf('q calculated for %i samples', i))
    }
    for (p in 1:n_peptides) {

      d = peaks_data[(i-1)*n_peptides+p, 2+1:n_isopeaks, drop=F]
      w = peaks_data[(i-1)*n_peptides+p, 2+n_isopeaks+1:n_isopeaks, drop=F]
      n2s <- w/d

      flag <- FALSE
      # use only data sets with lgth_gamma[p] + 1 values
      if (sum(!is.na(d)==T)>=lgth_gamma[p]+1) flag = TRUE

      if (flag) {
        # design matrix
        X = design(lgth_gamma[p], p, isodists)
        # this implicitely chooses among the first 1:nisotopes rows
        aux = LSest(X[!is.na(d), , drop=FALSE], d[1, !is.na(d), drop=FALSE], n2s[1, !is.na(d), drop=FALSE])
        gammas[i, p, 1:lgth_gamma[p]] = aux[1:lgth_gamma[p]]
        gammas[i, p, max(lgth_gamma)+1] = gammas[i, p, 1]/sum(gammas[i, p, 1:lgth_gamma[p]])
        gammas[i, p, max(lgth_gamma)+2] = aux[lgth_gamma[p]+1]
      }
    }
  }

  # Split peptides estimates and print to file
  q2e_estimates = list()
  for (p in 1:n_peptides) {
    m = peptides[[p, 1]]
    df_out = as.data.frame(round(gammas[, p, ], 6))
    df_out$sample = sample_names$sample
    df_out$peptide = m
    colnames(df_out) = c('gamma1', 'gamma2', 'gamma3', 'q', 'minLS',
                         'sample', 'peptide')
    q2e_estimates[[p]] = df_out
  }
  q2e_estimates = bind_rows(q2e_estimates)
  if (!is.null(outdir)){
    write_csv(q2e_estimates, file.path(outdir, "wls_q2e.csv"))
  }
  return(q2e_estimates)
}

