# For each peak, find the nearest TSS, and distance to it.
peak_nearest_tss = function(peaks,tss,midpoint=T) {
  results = list();
  for (chrom in names(peaks)) {
    if (!chrom %in% names(tss$ranges)) {
      next;
    }

    peak_ranges = peaks[[chrom]];

    if (midpoint) {
      peak_mids = mid(peak_ranges);
      peak_ranges = IRanges(peak_mids,peak_mids);
    }

    tss_ranges = tss$ranges[[chrom]];

    nearest_tss = tss_ranges[nearest(peak_ranges,tss_ranges)];

    suppressWarnings({
      # Hide the Bioconductor 2.12 warning about distance functionality changing
      dist_to_tss = distance(peak_ranges,nearest_tss);
    })

		if (midpoint) {
			d = data.frame(
					chrom=chrom,
					peak_midpoint=start(peak_ranges),
					nearest_tss=start(nearest_tss),
					dist_to_tss=dist_to_tss,
					geneid=names(nearest_tss)
			);
		} else {
			d = data.frame(
					chrom=chrom,
					peak_start=start(peak_ranges),
					peak_end=end(peak_ranges),
					nearest_tss=start(nearest_tss),
					dist_to_tss=dist_to_tss,
					geneid=names(nearest_tss)
			);
		}

    results[[chrom]] = d;
  }

  return(Reduce(rbind,results));
}

plot_dist_to_tss = function(peaks,genome='hg19') {
	# Get peaks from user's file.
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

	# Load TSS site info.
  tss_code = sprintf("tss.%s",genome);
  data(list=tss_code,package = "chipenrich.data");
  tss = get(tss_code);

	plotobj = ..plot_dist_to_tss(peakobj,tss);
	return(plotobj);
}

..plot_dist_to_tss = function(peaks,tss) {
  # Calculate distance to each TSS.
  tss_peaks = peak_nearest_tss(peaks,tss);

  # Create distance to TSS plot.
  max_dist = max(tss_peaks$dist_to_tss);
  breaks = breaks=c(0,100,1000,5000,10000,50000,100000,max_dist);
  dist_table = table(cut(tss_peaks$dist_to_tss,breaks=breaks));
  dist_table = dist_table / sum(dist_table);

  names(dist_table) = c("< 0.1","0.1 - 1","1 - 5","5 - 10","10 - 50","50 - 100","> 100")

  pf = function(...) {
    args <- list(...);
    bar_labels = sprintf("%0.1f%%",args$y * 100);
    panel.text(seq(1,length(args$x)),args$y + 0.03,bar_labels,cex=1.5);
    panel.barchart(...);
  }

  plotobj = barchart(
    dist_table,
    panel=pf,
    horizontal=F,
    scales=list(rot=45,cex=1.6),
    col="gray",
    ylim=c(0,1),
    ylab=list(label="Proportion of Peaks",cex=1.65),
    xlab=list(label="Distance to TSS (kb)",cex=1.65),
    main=list(label="Distribution of Distance from Peaks to Nearest TSS",cex=1.45)
  );

  return(plotobj);
}
