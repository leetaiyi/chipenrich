# This is the main function to read files
read_bed = function(file_path) {
	if (!file.exists(file_path)) {
		stop("Can't find BED file: ",file_path);
	}

	chunk_size = 500;
	chunk = scan(file_path,what="character",nmax=chunk_size,strip.white=T,sep="\n",quiet=T);
	skip_n = suppressWarnings(min(grep("^chr(\\d+|\\w+)\\s+\\d+",chunk)) - 1);

	if (is.infinite(skip_n)) {
		stop("Error: no valid chromosomes detected within first 500 lines of BED file.");
	}

	message(sprintf("Skipping %i lines of BED header..",skip_n));

	peaks = read.table(file_path,header=F,skip=skip_n);
	peaks = peaks[,1:3];
	names(peaks) = c("chrom","start","end");

	sub_check = peaks[1:min(nrow(peaks),100),];
	if (!all(grepl("chr",sub_check$chrom))) {
		stop("First column of BED file should have chr* entries.");
	}

	if (any(sub_check$start < 0) | any(sub_check$end < 0)) {
		stop("Start/end positions of peaks should be >= 0.");
	}

	chroms = list();
	for (chr in unique(peaks$chrom)) {
		peaks_chrom = subset(peaks,chrom == chr);
		chroms[[chr]] = IRanges(start=peaks_chrom$start,end=peaks_chrom$end);
	}

	chroms = reduce_peaks(chroms)

	return(chroms);
}

# This file type is for D. melanogaster files from modENCODE
read_bedgff = function(file_path) {
  if (!file.exists(file_path)) {
    stop("Can't find BED file: ",file_path);
  }

  chunk_size = 500;
  chunk = scan(file_path,what="character",nmax=chunk_size,strip.white=T,sep="\n",quiet=T);

  skip_n = suppressWarnings(min(grep("^\\d(L|R)",chunk)) - 1);

  if (is.infinite(skip_n)) {
    stop("Error: no valid chromosomes detected within first 500 lines of BED file.");
  }

  message(sprintf("Skipping %i lines of BED header..",skip_n));

  peaks = read.table(file_path,header=F,skip=skip_n);
  peaks = peaks[,c(1,4,5)];
  names(peaks) = c("chrom","start","end");

  sub_check = peaks[1:min(nrow(peaks),100),];
  if (!all(grepl("chr",sub_check$chrom))) {
    peaks$chrom = paste('chr',peaks$chrom,sep='')
    message("Adding 'chr' to chromosome column for compatibility with locus definitions.");
  }

  if (any(sub_check$start < 0) | any(sub_check$end < 0)) {
    stop("Start/end positions of peaks should be >= 0.");
  }

  chroms = list();
  for (chr in unique(peaks$chrom)) {
    peaks_chrom = subset(peaks,chrom == chr);
    chroms[[chr]] = IRanges(start=peaks_chrom$start,end=peaks_chrom$end);
  }

  chroms = reduce_peaks(chroms)

  return(chroms);
}

# This should be deprecated because read_bed does it better
read_peaks = function(file_path) {
  d = read.table(file_path,sep="\t",header=T);

  # Check columns.
  for (col in c("chrom","start","end")) {
    if (!col %in% names(d)) {
      stop(sprintf("error reading peaks: no column named '%s' in file",col));
    }
  }

  chroms = list();
  for (chr in unique(d$chrom)) {
    peaks_chrom = subset(d,chrom == chr);
    chroms[[chr]] = IRanges(start=peaks_chrom$start,end=peaks_chrom$end);
  }

  chroms = reduce_peaks(chroms)

  return(chroms);
}

# Read peaks from a dataframe having a particular structure
load_peaks = function(dframe) {
  # Check columns.
  for (col in c("chrom","start","end")) {
    if (!col %in% names(dframe)) {
      stop(sprintf("error in peaks data frame: no column named '%s'",col));
    }
  }

  chroms = list();
  for (chr in unique(dframe$chrom)) {
    peaks_chrom = subset(dframe,chrom == chr);
    chroms[[chr]] = IRanges(start=peaks_chrom$start,end=peaks_chrom$end);
  }

  chroms = reduce_peaks(chroms)

  return(chroms);
}

reduce_peaks = function(peaks) {
	peaks = lapply(peaks, IRanges::reduce)
  peaks = lapply(peaks, unique)
	return(peaks)
}
