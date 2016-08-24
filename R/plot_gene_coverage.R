# Used in ..plot_gene_coverage(...)
avg_binned_coverage = function(gpw,bin_size=25) {
  d = gpw;
  d = d[order(d$log10_length),];
  d$group = ceiling((1:dim(d)[1])/bin_size);

  bygroup = stats::aggregate(cbind(ratio,log10_length) ~ group,d,mean);
  names(bygroup) = c("group","ratio","log_avg_length");
  bygroup;
}

plot_gene_coverage = function(peaks,locusdef="nearest_tss",genome='hg19',use_mappability=F,read_length=36,legend=T,xlim=NULL) {
	# Check genome.
  if (!genome %in% supported_genomes()) {
    stop("genome not supported: ",genome);
  }

  # Check locus definition. Should only be 1.
  if (!locusdef %in% supported_locusdefs()) {
    stop("bad locus definition requested: ",locusdef);
  }

  # Check read length.
  if (use_mappability) {
    if (!as.numeric(read_length) %in% supported_read_lengths()) {
      stop("bad read length requested: ",read_length);
    }
  }

  # Get peaks from user's file.
	if (class(peaks) == "data.frame") {
		peakobj = load_peaks(peaks);
	} else if (class(peaks) == "character") {
    if (str_sub(peaks,-4,-1) == ".gff" || str_sub(peaks,-5,-1) == '.gff3' || str_sub(peaks,-7,-1) == ".gff.gz" || str_sub(peaks,-8,-1) == '.gff3.gz') {
      message("Reading peaks file: ",peaks);
      peakobj = read_bedgff(peaks);
    } else {
      message("Reading peaks file: ",peaks);
      peakobj = read_bed(peaks);
    }
	}

  peakobj = reduce_peaks(peakobj)

	# Number of peaks in data.
  num_peaks = sum(sapply(peakobj,function(x) length(x)))

  # Load locus definitions.
  ldef_code = sprintf("locusdef.%s.%s",genome,locusdef);
  data(list=ldef_code,package = "chipenrich.data");
  ldef = get(ldef_code);

	# # Load TSS site info.
  # tss_code = sprintf("tss.%s",genome);
  # data(list=tss_code,package = "chipenrich.data");
  # tss = get(tss_code);

  # Load mappability if requested.
  if (use_mappability) {
    mappa_code = sprintf("mappa.%s.%s.%imer",genome,locusdef,read_length);
    data(list=mappa_code,package = "chipenrich.data");
    mappa = get(mappa_code);
  } else {
    mappa = NULL;
  }

  # Assign peaks to genes.
  assigned_peaks = assign_peak_segments(peakobj,ldef);
  peak_genes = unique(assigned_peaks$geneid);

  ppg = num_peaks_per_gene(assigned_peaks,ldef,mappa=NULL)

  ppg = calc_peak_gene_overlap(assigned_peaks,ppg)

	# Make plot.
	plotobj = ..plot_gene_coverage(ppg);
	return(plotobj);
}

..plot_gene_coverage = function(ppg) {

	avg_bins = avg_binned_coverage(ppg, bin_size=25)

	plotobj = xyplot(
		ratio ~ log_avg_length,
		avg_bins,
		main='Binned Locus Length versus Peak Coverage',
		xlab='log10(locus length)',
		ylab='Proportion of locus covered by peak',
		pch=20,
		col='black'
	);

	return(plotobj);
}
