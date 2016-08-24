..plot_expected_peaks = function(ppg) {
  # Gene/peak counts.
  num_genes = dim(ppg)[1];
  num_peaks = sum(ppg$num_peaks);

  # Compute prob of mappa for each gene.
  genome_length = sum(as.numeric(ppg$length));
  ppg$prob_mappable = ppg$length / genome_length;

  # Expected peaks (lambda)
  ppg$exp_peaks = num_peaks * ppg$prob_mappable;

  # Quantiles.
  ppg$qpois_05 = qpois(0.05,ppg$exp_peaks);
  ppg$qpois_95 = qpois(0.95,ppg$exp_peaks);
  ppg$qpois_bonf_05 = qpois(0.05 / num_genes,ppg$exp_peaks);
  ppg$qpois_bonf_95 = qpois(1 - (0.05 / num_genes),ppg$exp_peaks);

  # Max # of peaks, expected or observed.
  max_peaks = max(ppg$num_peaks,ppg$expected_peaks);
  pad = round(0.05 * max_peaks);

  # Thin out num_peaks and expected_peaks vs. length.
  d_npeak = ppg[,c("num_peaks","length")];
  d_npeak = unique(d_npeak);

  d_epeak = ppg[,c("exp_peaks","length")];
  d_epeak = unique(d_epeak);

  # Colors.
  col_npeaks = "blue";
  col_epeaks = "black";
  col_qpois_05 = "green";
  col_qpois_95 = "darkgreen";
  col_qpois_bonf_05 = "red";
  col_qpois_bonf_95 = "darkred";

  custom_key = list(
    text = list(
      c(
        'Count of Peaks',
        'Expected Peaks',
        '5th Percentile',
        '95th Percentile',
        '5th Percentile (Bonferroni)',
        '95th Percentile (Bonferroni)'
      )
    ),
    lines = list(
      pch = c(20),
      type = c("p","l","l","l","l","l"),
      lwd = 2,
      col = c(
        col_npeaks,
        col_epeaks,
        col_qpois_05,
        col_qpois_95,
        col_qpois_bonf_05,
        col_qpois_bonf_95
      )
    ),
    cex = 1.2,
    columns = 2
  );

  p = xyplot(
    num_peaks ~ length,
    d_npeak,
    par.settings = simpleTheme(pch=20,col=col_npeaks),
    ylim = c(-pad,max_peaks+pad),
    key = custom_key,
    scales = list(cex=1.4),
    xlab = list(label="Gene Locus Length",cex=1.4),
    ylab = list(label="Count of Peaks",cex=1.4)
  );

  p = p + as.layer(xyplot(
    exp_peaks ~ length,
    d_epeak,
    par.settings = simpleTheme(pch=20,col=col_epeaks,lwd=2),
    type = "l"
  ));

  for (v in c("qpois_05","qpois_95","qpois_bonf_05","qpois_bonf_95")) {
    d = ppg[,c(v,"length")];
    d = unique(round(d,2));
    d = d[order(d$length),];

    pcol = get(sprintf("col_%s",v));
    form = as.formula(sprintf("%s ~ length",v));
    p = p + as.layer(xyplot(
      form,
      d,
      par.settings = simpleTheme(pch=20,col=pcol,lwd=2),
      type = "l"
    ));
  }

  return(p);
}

plot_expected_peaks = function(peaks,locusdef="nearest_tss",genome="hg19",use_mappability=F,read_length=36,mappa_file=NULL) {
   # Check genome.
  if (!genome %in% supported_genomes()) {
    stop("genome not supported: ",genome);
  }

  # Organism - human or mouse.
  organism = genome_to_organism(genome);

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
    if (get_ext(peaks) == "bed") {
      message("Reading BED file: ",peaks);
      peakobj = read_bed(peaks);
    } else {
      message("Reading peaks file: ",peaks);
      peakobj = read_peaks(peaks);
    }
	}

	# Number of peaks in data.
  num_peaks = sum(sapply(peakobj,function(x) length(x)))

  # Load locus definitions.
  ldef_code = sprintf("locusdef.%s.%s",genome,locusdef);
  data(list=ldef_code,package = "chipenrich.data");
  ldef = get(ldef_code);

  # Load mappability if requested.
  if (!is.null(mappa_file)) {
    message("Reading user-defined gene locus mappability file: ",mappa_file);
    mappa = read_mappa(mappa_file);
  } else {
    if (use_mappability) {
      mappa_code = sprintf("mappa.%s.%s.%imer",genome,locusdef,read_length);
      data(list=mappa_code,package = "chipenrich.data");
      mappa = get(mappa_code);
      mappa = na.omit(mappa);
    } else {
      mappa = NULL;
    }
  }

  # Load TSS site info.
  tss_code = sprintf("tss.%s",genome);
  data(list=tss_code,package = "chipenrich.data");
  tss = get(tss_code);

  # Assign peaks to genes.
  message("Assigning peaks to genes..");
  assigned_peaks = assign_peaks(peakobj,ldef,tss,midpoint=T);

  # Data about the # of peaks per gene.
  ppg = num_peaks_per_gene(assigned_peaks,ldef,mappa);

  # Plot.
  p = ..plot_expected_peaks(ppg);

  return(p);
}
