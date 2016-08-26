# Used in ..plot_spline_length(...)
add_emp_peak = function(gpw, bin_size=25) {
	d = gpw
	d = d[order(d$log10_length),]
	d$group = ceiling((1:dim(d)[1])/bin_size)

	bygroup = stats::aggregate(peak ~ group,d,mean)
	d$emp_peak = bygroup$peak[d$group]
	d
}

# Used in ..plot_spline_length(...)
avg_binned_peak = function(gpw,bin_size=25) {
	d = gpw
	d = d[order(d$log10_length),]
	d$group = ceiling((1:dim(d)[1])/bin_size)

	bygroup = stats::aggregate(cbind(peak,length) ~ group,d,mean)
	bygroup$log_avg_length = log10(bygroup$length)
	names(bygroup) = c("group","peak","avg_length","log_avg_length")
	bygroup
}

# Used in ..plot_spline_length(...)
calc_weights_gam = function(locusdef,peak_genes,mappa=NULL,...) {
	d = locusdef@dframe

	# Indicator vector for which genes have peaks.
	d$peak = as.numeric(d$geneid %in% peak_genes)

	# Compute length and log10 length for each gene.
	d$length = d$end - d$start

	# For genes that exist across multiple rows, sum their lengths.
	d = stats::aggregate(cbind(peak,length) ~ geneid,d,sum)
	d$log10_length = log10(d$length)

	# A gene could now have > 1 peak due to aggregating (above), reset it to be
	# 1 or 0.
	d$peak = as.integer(d$peak >= 1)

	# Sort by locus length.
	d = d[order(d$log10_length),]

	# If mappability was requested, add that in.
	if (!is.null(mappa)) {
		d = merge(d, mappa, by = "geneid", sort = F)
		d$orig_length = d$length
		d$length = as.numeric((d$mappa * d$length) + 1)
		d$log10_length = log10(d$length)
		d = d[order(d$log10_length), ]
	}

	# Create model.
	model = "peak ~ s(log10_length, bs = 'cr')"

	# Compute binomial spline fit.
	fit = gam(as.formula(model), data = d, family = "binomial")

	# Compute weights for each gene, based on the predicted prob(peak) for each gene.
	ppeak = fitted(fit)
	w0 = 1 / (ppeak/mean(d$peak, na.rm = T))
	w0 = w0 / mean(w0,na.rm=T)

	d$weight = w0
	d$prob_peak = ppeak
	d$resid.dev = resid(fit,type="deviance")

	cols = c("geneid", "length", "log10_length", "mappa", "orig_length", "peak", "weight", "prob_peak", "resid.dev")
	if (is.null(mappa)) {
		cols = setdiff(cols, c("mappa", "orig_length"))
	}

	d = subset(d, select = cols)
	return(d)
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
#'
#' plot_spline_length(peaks_E2F4,genome='hg19')
#'
#' # Create the plot for a different locus definition
#' # to compare the effect.
#' plot_spline_length(peaks_E2F4, locusdef = 'nearest_gene', genome = 'hg19')
#'
#' @export
plot_spline_length = function(peaks, locusdef = "nearest_tss", genome = 'hg19', use_mappability = F, read_length = 36, legend = T, xlim = NULL) {
	# Check genome.
	if (!genome %in% supported_genomes()) {
		stop("genome not supported: ",genome)
	}

	# Check locus definition. Should only be 1.
	if (!any(supported_locusdefs()$genome == genome & supported_locusdefs()$locusdef == locusdef)) {
		stop(
			sprintf("Error: invalid genome / definition combination requested: %s %s",
				genome, locusdef
			)
		)
	}

	# Check read length.
	if (use_mappability) {
		if (!any(supported_read_lengths()$genome == genome & supported_read_lengths()$locusdef == locusdef & supported_read_lengths()$read_length == as.numeric(read_length))) {
			stop(
				sprintf("Error: bad genome / locusdef / read length combination requested: %s %s %s",
					genome, locusdef, read_length)
			)
		}
	}

	# Get peaks from user's file.
	if (class(peaks) == "data.frame") {
		peakobj = load_peaks(peaks)
	} else if (class(peaks) == "character") {
		if (stringr::str_sub(peaks,-4,-1) == ".gff" || stringr::str_sub(peaks,-5,-1) == '.gff3' || stringr::str_sub(peaks,-7,-1) == ".gff.gz" || stringr::str_sub(peaks,-8,-1) == '.gff3.gz') {
			message("Reading peaks file: ", peaks)
			peakobj = read_bedgff(peaks)
		} else {
			message("Reading peaks file: ", peaks)
			peakobj = read_bed(peaks)
		}
	}

	# Number of peaks in data.
	num_peaks = sum(sapply(peakobj, function(x) length(x)))

	# Load locus definitions.
	ldef_code = sprintf("locusdef.%s.%s", genome, locusdef)
	data(list = ldef_code, package = "chipenrich.data", envir = environment())
	ldef = get(ldef_code)

	# Load TSS site info.
	tss_code = sprintf("tss.%s", genome)
	data(list = tss_code, package = "chipenrich.data", envir = environment())
	tss = get(tss_code)

	# Load mappability if requested.
	if (use_mappability) {
		mappa_code = sprintf("mappa.%s.%s.%imer", genome, locusdef, read_length)
		data(list = mappa_code, package = "chipenrich.data", envir = environment())
		mappa = get(mappa_code)
	} else {
		mappa = NULL
	}

	# Assign peaks to genes.
	assigned_peaks = assign_peaks(peakobj, ldef, tss)
	peak_genes = unique(assigned_peaks$geneid)

	# Make plot.
	plotobj = ..plot_spline_length(ldef, peak_genes, num_peaks, mappa, legend=legend, xlim=xlim)
	return(plotobj)
}

# Create diagnostic plot of Proportion of Peaks from your data against log locus length (of genes)
# along with the spline fit.
..plot_spline_length = function(locusdef, peak_genes, num_peaks, mappa = NULL, legend = T, xlim = NULL) {
	# Calculate smoothing spline fit.
	gpw = calc_weights_gam(locusdef, peak_genes, mappa = mappa) # gpw = genes, peaks, weights

	# Genome length.
	genome_length = sum(as.numeric(gpw$length))

	# Add in empirical peak by binning.
	gpw = add_emp_peak(gpw)

	# Average peak/lengths.
	avg_bins = avg_binned_peak(gpw, bin_size = 25)

	# Order by length.
	gpw = gpw[order(gpw$log10_length), ]

	# Calculate prob(all false positives).
	gpw$false_prob = 1 - (1 - (gpw$length/genome_length))^(num_peaks)

	col_rand_gene = "grey35"
	col_rand_peak = "grey74"
	col_spline = "darkorange"

	panel_func = function(x,y,...) {
		panel.xyplot(x, y, ...)
		panel.abline(h = mean(y), col = col_rand_gene, lwd = 3)
	}

	custom_key = list(
		text = list(
		c(
			"Expected Fit - Peaks Independent of Locus Length",
			"Expected Fit - Peaks Proportional to Locus Length",
			"Binomial Smoothing Spline Fit",
			"Proportion of Genes in Bin with at Least 1 Peak"
			)
		),
		lines = list(
			pch = c(20,15,17,20),
			type = c("l","l","l","p"),
			col = c(col_rand_gene, col_rand_peak, col_spline, "black"),
			lwd = 3
		),
		cex = 1.2
	)

	if (!legend) {
		custom_key = NULL
	}

	if (!is.null(mappa)) {
		xlab = expression(paste(Log[10], " Mappable Locus Length"))
	} else {
		xlab = expression(paste(Log[10], " Locus Length"))
	}

	xmin_nopad = base::ifelse(is.null(xlim[1]), floor(min(gpw$log10_length)), floor(xlim[1]))
	xmax_nopad = base::ifelse(is.null(xlim[2]), ceiling(max(gpw$log10_length)), ceiling(xlim[2]))

	scales = list(
		x = list(
			axs = 'i',
			at = seq(xmin_nopad,xmax_nopad,1)
		),
		y = list(
			axs = 'i',
			at = seq(0,1,0.2)
		)
	)

	plotobj = xyplot(
		false_prob + prob_peak ~ log10_length,
		gpw,
		xlab = list(label=xlab,cex=1.4),
		ylab = list(label = "Proportion of Peaks", cex = 1.4),
		ylim = c(-0.05,1.05),
		xlim = c(xmin_nopad - 0.5, xmax_nopad + 0.5),
		panel = panel_func,
		type = "l",
		lwd = 3,
		key = custom_key,
		scales = list(cex = 1.4),
		par.settings = simpleTheme(pch = c(15,17), col = c(col_rand_peak, col_spline))
	) + as.layer(xyplot(peak ~ log_avg_length, avg_bins, pch = 20, cex = 0.4, col = "black"))

	return(plotobj)
}
