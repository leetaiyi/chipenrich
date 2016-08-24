#!/usr/bin/env Rscript
os = Sys.info()[1]
# library(chipenrich.data)
# library(mgcv)
# library(IRanges)
# library(GenomicRanges)
# library(lattice)
# library(latticeExtra)
# library(grid)
# library(stringr)
# library(rms)
# library(parallel)

chipenrich = function(
	peaks,
	out_name = "chipenrich",
	out_path = getwd(),
	genome = "hg19",
	genesets = c(
		'GOBP',
		'GOCC',
		'GOMF'),
	locusdef = "nearest_tss",
	method = "chipenrich",
	fisher_alt = "two.sided",
	use_mappability = F,
	mappa_file = NULL,
	read_length = 36,
	qc_plots = T,
	max_geneset_size = 2000,
	num_peak_threshold = 1,
	n_cores = 1
) {

	############################################################################
	# Collect options for opts output
	l = unlist(as.list(environment()))
		opts = data.frame(
		args = names(l),
		values = sapply(unlist(l),paste,collapse=","),
		stringsAsFactors = F
	)
	rownames(opts) = 1:length(l)

	############################################################################
	# Deal with randomizations if present
	# Randomizations are accessed by appending _rndall or _rndlength to
	# the geneset names. Parsing is done here.
	rndall = all(grepl('rndall',genesets))
	rndlength = all(grepl('rndlength',genesets))
	rndloc = all(grepl('rndloc',genesets))

	if(rndall) {
		genesets = gsub('_rndall','',genesets)
	} else if (rndlength) {
		genesets = gsub('_rndlength','',genesets)
	} else if (rndloc) {
		genesets = gsub('_rndloc','',genesets)
	}

	############################################################################
	############################################################################
	# Warnings based on OS or method/ldef combinations
	############################################################################
	############################################################################

	############################################################################
	# CHECK OS and give multicore warning if Windows
	if(os == 'Windows') {
		message('Warning! Multicore enrichment is not supported on Windows.')
	}

	# Warn user if they are trying to use FET with a
	# locus definition that might lead to biased results.
	if (method == "fet") {
		if (is.character(locusdef) && !locusdef %in% c("1kb","5kb")) {
			message("Warning: Fisher's exact test should only be used with the 1kb or 5kb locus definition.")
		}

		if (user_defined_ldef) {
			message("Warning: Fisher's exact test may give biased results if the spline fit for the gene locus definitions is not flat (see QC plots.)")
		}
	}

	# Warn user if they are using the binomial test.
	if (method == "binomial") {
		message("Warning: the binomial test is provided for comparison purposes only.")
		message("This test will almost always give biased results favoring gene sets with short average locus length.")
	}

	############################################################################
	############################################################################
	# Checks and genome, locusdef, geneset, and mappa setup
	############################################################################
	############################################################################

	############################################################################
	# CHECK genome and load TSS if okay
	if (!genome %in% supported_genomes()) {
		stop("genome not supported: ",genome)
	}

	# Get the organism code once genome passes
	organism = genome_to_organism(genome)

	# Load TSS site info.
	tss_code = sprintf("tss.%s", genome)
	data(list=tss_code, package = "chipenrich.data")
	tss = get(tss_code)

	############################################################################
	# CHECK genesets and load them if okay
	# Determine if geneset codes are valid before moving on. The API for a user
	# to use their own genesets will be to put a path in the genesets argument.
	if (!check_arg(genesets, supported_genesets())) {
		bad_args = check_arg(genesets, supported_genesets(), value=T)

		# If the bad_args is a path that exists, then we know the user wants
		# to provide their own genesets
		if(file.exists(genesets)) {
			message('User-specified geneset(s)...')
		} else {
			stop("Invalid geneset(s) requested: ", paste(bad_args, collapse = ", "))
		}
	}

	# Load genesets and filter them
	if(!file.exists(genesets)) {
		# If the user does not provide a path to the geneset.
		# That is, do the normal thing
		geneset_list = list()
		for (gs in genesets) {
			geneset_code = sprintf("geneset.%s.%s", gs, organism)
			data(list = geneset_code, package = "chipenrich.data")

			geneset_list[[geneset_code]] = filter_genesets(get(geneset_code), max_geneset_size)
		}
	} else {
		# If the user provides a path to the geneset, build it and filter it
		geneset_list = list()
		geneset_code = 'user-supplied'

		geneset_list[[geneset_code]] = setup_geneset(genesets)
		geneset_list[[geneset_code]] = filter_genesets(geneset_list[[geneset_code]], max_geneset_size)
	}

	############################################################################
	# CHECK locusdefs and load if okay
	# The API for a user to use their own genesets will be to put a path in
	# the genesets argument.
	user_defined_ldef = file.exists(locusdef)
	if (user_defined_ldef) {
		# Load user-defined locus definition file.
		ldef = setup_ldef(locusdef)
	} else {
		if (!locusdef %in% supported_locusdefs()) {
			stop("Error: invalid definition requested: ",locusdef)
		}

		# Load locus definitions.
		ldef_code = sprintf("locusdef.%s.%s",genome,locusdef)
		data(list=ldef_code,package = "chipenrich.data")
		ldef = get(ldef_code)

		# Randomize locus definition if rndloc == T
		if(rndloc) {
			ldef = randomize_locusdef(ldef, 50)
		}
	}

	############################################################################
	# CHECK mappability and load if okay
	# If the user specified their own mappability, we can't use the built
	# in mappability - they must provide their own.
	user_defined_mappa = !is.null(mappa_file)

	# The user specified mappability - if they accidentally also thought
	# use_mappability had to be enabled, just disable it.
	if (user_defined_mappa & use_mappability) {
		use_mappability = FALSE
	}

	# Check for user-defined mappability with built-in mappability TRUE
	if (user_defined_ldef) {
		if (use_mappability) {
			message("Warning: built-in mappability cannot be used with a user-defined locus definition, you must calculate mappability for your definition and pass it in with the mappa_file argument.")
			use_mappability = FALSE
		}
	}

	# Check read length for using built-in mappability
	if (use_mappability) {
		if (!as.numeric(read_length) %in% supported_read_lengths()) {
			stop("Error: bad read length requested: ", read_length)
		}
	}

	# Load mappability if requested. Either user-defined or built-in.
	if (user_defined_mappa) {
		message("Reading user-specified gene locus mappability file: ", mappa_file)
		mappa = read_mappa(mappa_file)
	} else {
		if (use_mappability) {
			mappa_code = sprintf("mappa.%s.%s.%imer", genome, locusdef, read_length)
			data(list = mappa_code, package = "chipenrich.data")
			mappa = get(mappa_code)
			mappa = na.omit(mappa)
		} else {
			mappa = NULL
		}
	}

	# If they specified their own mappability, check gene names overlap
	if (user_defined_mappa && user_defined_ldef) {
		total_unique_genes = union(mappa$geneid, ldef@dframe$geneid)
		mappa_ldef_inters = intersect(mappa$geneid, ldef@dframe$geneid)
		frac_overlap = length(mappa_ldef_inters) / length(total_unique_genes)

		if (frac_overlap < 0.95) {
			stop("Error: your mappability genes and locus definition genes overlap by less than 95% (they should match almost exactly)..")
		}
	}

	############################################################################
	# CHECK method and get() it if okay
	get_test_method = function(x) {
		if (method %in% names(SUPPORTED_METHODS)) {
			return(SUPPORTED_METHODS[[method]])
		} else if (method %in% names(HIDDEN_METHODS)) {
			return(HIDDEN_METHODS[[method]])
		} else {
			stop(sprintf("Error: invalid enrichment test requested: %s, contact developer.",method))
		}
	}
	testf = get_test_method(method)
	test_func = get(testf)

	# Test name.
	method_name = METHOD_NAMES[[method]]

	############################################################################
	############################################################################
    # Start enrichment process
	############################################################################
	############################################################################

	######################################################
	# Read in and format peaks (from data.frame or file)
	if (class(peaks) == "data.frame") {
		peakobj = load_peaks(peaks)
	} else if (class(peaks) == "character") {
		if (str_sub(peaks,-4,-1) == ".gff" || str_sub(peaks,-5,-1) == '.gff3' || str_sub(peaks,-7,-1) == ".gff.gz" || str_sub(peaks,-8,-1) == '.gff3.gz') {
			message("Reading peaks file: ",peaks)
			peakobj = read_bedgff(peaks)
		} else {
			message("Reading peaks file: ",peaks)
			peakobj = read_bed(peaks)
		}
	}

	# Number of peaks in data.
	num_peaks = sum(sapply(peakobj, function(x) length(x)))

	######################################################
	# Assign peaks to genes. NOTE: If method = 'broadenrich' use
	# assign_peak_segments(), otherwise use assign_peaks().
	if(!(method == 'broadenrich' || method == 'broadenrich_splineless')) {
		message("Assigning peaks to genes with assign_peaks(...) ..")
		assigned_peaks = assign_peaks(peakobj, ldef, tss)
		message("Successfully assigned peaks..")
	} else {
		message("Assigning peaks to genes with assigned_peak_segments(...) ..")
		assigned_peaks = assign_peak_segments(peakobj, ldef)
		message("Successfully assigned peaks..")
	}

	######################################################
	# Post-process assigned peaks to add gene symbols and order the columns
	peak_genes = unique(assigned_peaks$geneid)

	# Add gene symbols to peak genes using the genes.* object
	genes_code = sprintf("genes.%s", organism)
	data(list = genes_code, package = "chipenrich.data")
	gene2symbol = get(genes_code)
	gene2symbol = change_names(gene2symbol, list(GENEID = "geneid", SYMBOL = "gene_symbol"))
	assigned_peaks = merge(assigned_peaks, gene2symbol, by="geneid", all.x=T)

	# Order the columns. NOTE: This includes the union of column names
	# when using assign_peaks() and assign_peak_segments()
	column_order = c(
		"peak_id",
		"chrom",
		"peak_start",
		"peak_end",
		"peak_midpoint",
		"geneid",
		"gene_symbol",
		"gene_locus_start",
		"gene_locus_end",
		"nearest_tss",
		"dist_to_tss",
		"nearest_tss_gene",
		"nearest_tss_gene_strand",
		"overlap_start",
		"overlap_end",
		"peak_overlap")
	column_order = intersect(column_order, names(assigned_peaks))
	assigned_peaks = assigned_peaks[, column_order]

	######################################################
	# Compute peaks per gene table
	ppg = num_peaks_per_gene(assigned_peaks, ldef, mappa)
	# This seems redundant given num_peaks_per_gene(...)
	ppg$peak = recode_peaks(ppg$num_peaks, num_peak_threshold)

	# Add relevant columns to ppg depending on the method
	if(method == 'broadenrich' || method == 'broadenrich_splineless') {
		message("Calculating peak overlaps with gene loci..")
		ppg = calc_peak_gene_overlap(assigned_peaks,ppg)
	}
	if(method == 'chipapprox') {
		message("Calculating weights for approximate method..")
		ppg = calc_approx_weights(ppg,mappa)
	}

	######################################################
	# Randomize ppg table if randomization API invoked
	# Catch randomizations if present
	if(rndall) {
		message('Randomizing across all genes.')
		ppg = randomize_ppg_all(ppg)
	} else if (rndlength) {
		message('Randomizing within length bins.')
		ppg = randomize_ppg_length(ppg)
	}

	######################################################
	# Enrichment
	# Run chipenrich method on each geneset.
	results = list()
	for (gobj in geneset_list) {
		message(sprintf("Test: %s",method_name))
		message(sprintf("Genesets: %s",gobj@type))
		message("Running tests..")
		if (testf == "test_gam") {
			rtemp = test_func(gobj,ppg,n_cores)
		}
		if (testf == "test_fisher_exact") {
			rtemp = test_func(gobj,ppg,alternative=fisher_alt)
		}
		if (testf == "test_binomial") {
			rtemp = test_func(gobj,ppg)
		}
		if (testf == "test_gam_ratio") {
			rtemp = test_func(gobj,ppg,n_cores)
		}
		if (testf == "test_approx") {
			rtemp = test_func(gobj,ppg,nwp=F,n_cores)
		}
		if (testf == "test_gam_ratio_splineless") {
			rtemp = test_func(gobj,ppg,n_cores)
		}

		# Annotate with geneset descriptions.
		rtemp$"Description" = as.character(mget(rtemp$Geneset.ID, gobj@set.name, ifnotfound=NA))
		rtemp$"Geneset.Type" = gobj@type;

		results[[gobj@type]] = rtemp;
	}
	enrich = Reduce(rbind,results)

	######################################################
	# Post-process enrichment with various orderings

	# Re-order the columns to something sensible.
	column_order = c(
		"Geneset.Type",
		"Geneset.ID",
		"Description",
		"P.value",
		"FDR",
		"Effect",
		"Odds.Ratio",
		"P.Success",
		"Status",
		"N.Geneset.Genes",
		"N.Geneset.Peak.Genes",
		"Geneset.Avg.Gene.Length",
		"Geneset.Avg.Gene.Coverage",
		"Geneset.Peak.Genes")
	column_order = intersect(column_order, names(enrich))
	enrich = enrich[, column_order]

	# Order results by p-value.
	enrich = enrich[order(enrich$P.value), ]

	# If there is a status column, re-sort so enriched terms are on top.
	if ("Status" %in% names(enrich)) {
		enrich = enrich[order(enrich$Status, decreasing=T), ]
	}

	# Pull out tests that failed.
	bad_enrich = subset(enrich,is.na(P.value))
	enrich = subset(enrich,!is.na(P.value))

	######################################################
	# Write result objects to files
	if (!is.null(out_name)) {
		filename_analysis = file.path(out_path, sprintf("%s_results.tab", out_name))
		write.table(enrich, file = filename_analysis, row.names = F, quote = F, sep = "\t")
		message("Wrote results to: ", filename_analysis)

		filename_peaks = file.path(out_path, sprintf("%s_peaks.tab", out_name))
		write.table(assigned_peaks, file = filename_peaks, row.names = F, quote = F, sep = "\t")
		message("Wrote peak-to-gene assignments to: ", filename_peaks)

		filename_opts = file.path(out_path, sprintf("%s_opts.tab", out_name))
		write.table(opts, file = filename_opts, row.names = F, quote = F, sep = "\t")
		message("Wrote run options/arguments to: ", filename_opts)

		filename_ppg = file.path(out_path, sprintf("%s_peaks-per-gene.tab", out_name))
		write.table(ppg, file = filename_ppg, row.names = F, quote = F, sep = "\t")
		message("Wrote count of peaks per gene to: ", filename_ppg)

		# If the user requested QC plots, generate those as well.
		# -- Spline fit plot
		# -- Histogram of distance from peaks to TSSs
		# -- Histogram of p-values from test
		# -- Expected # peaks vs. observed # peaks
		if (qc_plots) {
			filename_qcplots = file.path(out_path, sprintf("%s_qcplots.pdf", out_name))
			pdf(filename_qcplots)
				if (!(method=='broadenrich' || method=='broadenrich_splineless')) {
					print(..plot_spline_length(ldef, peak_genes, num_peaks, mappa=mappa))
					print(..plot_dist_to_tss(peakobj, tss))
				} else {
					print(..plot_gene_coverage(ppg))
				}
			dev.off()
			message("Wrote QC plots to: ",filename_qcplots)
		}
	}

	######################################################
	# Return objects as list
	return(list(
		peaks = assigned_peaks,
		results = enrich,
		opts = opts,
		peaks_per_gene = ppg
	))
}
