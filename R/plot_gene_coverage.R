# Used in ..plot_gene_coverage(...)
avg_binned_coverage = function(gpw, bin_size = 25) {
	d = gpw
	d = d[order(d$log10_length), ]
	d$group = ceiling((1:dim(d)[1])/bin_size)

	bygroup = stats::aggregate(cbind(ratio, log10_length) ~ group, d, mean)
	names(bygroup) = c("group", "ratio", "log_avg_length")
	bygroup
}

#' Plot probability of peak being assigned to a gene vs. gene length
#'
#' Create a plot showing the probability of a gene being assigned a peak given
#' its locus length. The plot shows an empirical fit to the data using a binomial
#' smoothing spline.
#'
#' @param peaks A \code{data.frame}, or tab-delimited text file (BED, narrowPeak,
#' broadPeak, etc) with the first three columns being chrom, start, and end. The
#' data frame should have at least 3 columns: chrom, start, and end. Chrom
#' should follow UCSC convention, e.g. "chrX".
#' @param locusdef A string denoting the gene locus definition to be used, or
#' the full path to a user-defined locus definition file. A gene locus definition
#' controls how peaks are assigned to genes. See \code{\link{supported_locusdefs}}
#' for a list of supported definitions built-in. If using a user-specified file,
#' the file must have 4 columns: geneid, chrom, start, end and be tab-delimited.
#' @param genome A string indicating the genome upon which the peaks file is
#' based. Supported genomes are listed by the \code{\link{supported_genomes}} function.
#' @param use_mappability A logical variable indicating whether to adjust for
#' mappability. If enabled, this option will use our internally calculated
#' mappabilities for each gene locus given the length of reads used in the
#' experiment (see \code{read_length} option). NOTE: If providing a \code{mappa_file}
#' this parameter should be set to FALSE.
#' @param read_length If adjusting for mappability, this number specifies the
#' read length to be used. The read length given here should ideally correspond
#' to the length of reads from the original experiment.
#' @param legend If true, a legend will be drawn on the plot.
#' @param xlim Set the x-axis limit. NULL means select x-lim automatically.
#'
#' @return A trellis plot object.
#'
#' @examples
#'
#' # Spline plot for E2F4 example peak dataset.
#' data(peaks_E2F4, package = 'chipenrich.data')
#' peaks_E2F4 = subset(peaks_E2F4, peaks_E2F4$chrom == 'chr1')
#' plot_gene_coverage(peaks_E2F4)
#'
#' @export
plot_gene_coverage = function(peaks, locusdef = "nearest_tss", genome = 'hg19', use_mappability = F, read_length = 36, legend = T, xlim = NULL) {
	# Check genome.
	if (!genome %in% supported_genomes()) {
		stop("genome not supported: ",genome)
	}

	# Check locus definition. Should only be 1.
	if (!locusdef %in% supported_locusdefs()) {
		stop("bad locus definition requested: ",locusdef)
	}

	# Check read length.
	if (use_mappability) {
		if (!as.numeric(read_length) %in% supported_read_lengths()) {
			stop("bad read length requested: ",read_length)
		}
	}

	# Get peaks from user's file.
	if (class(peaks) == "data.frame") {
		peakobj = load_peaks(peaks)
	} else if (class(peaks) == "character") {
		if (str_sub(peaks,-4,-1) == ".gff" || str_sub(peaks,-5,-1) == '.gff3' || str_sub(peaks,-7,-1) == ".gff.gz" || str_sub(peaks,-8,-1) == '.gff3.gz') {
			message("Reading peaks file: ", peaks)
			peakobj = read_bedgff(peaks)
		} else {
			message("Reading peaks file: ", peaks)
			peakobj = read_bed(peaks)
		}
	}

	# Number of peaks in data.
	num_peaks = sum(sapply(peakobj,function(x) length(x)))

	# Load locus definitions.
	ldef_code = sprintf("locusdef.%s.%s", genome, locusdef)
	data(list = ldef_code, package = "chipenrich.data")
	ldef = get(ldef_code)

	# Load mappability if requested.
	if (use_mappability) {
		mappa_code = sprintf("mappa.%s.%s.%imer", genome, locusdef, read_length)
		data(list = mappa_code, package = "chipenrich.data")
		mappa = get(mappa_code)
	} else {
		mappa = NULL
	}

	# Assign peaks to genes.
	assigned_peaks = assign_peak_segments(peakobj, ldef)
	peak_genes = unique(assigned_peaks$geneid)

	ppg = num_peaks_per_gene(assigned_peaks, ldef, mappa=NULL)
	ppg = calc_peak_gene_overlap(assigned_peaks, ppg)

	# Make plot.
	plotobj = ..plot_gene_coverage(ppg)
	return(plotobj)
}

..plot_gene_coverage = function(ppg) {

	avg_bins = avg_binned_coverage(ppg, bin_size = 25)

	plotobj = xyplot(
		ratio ~ log_avg_length,
		avg_bins,
		main = 'Binned Locus Length versus Peak Coverage',
		xlab = 'log10(locus length)',
		ylab = 'Proportion of locus covered by peak',
		pch = 20,
		col = 'black'
	)

	return(plotobj)
}
