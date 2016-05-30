#!/usr/bin/env Rscript
os = Sys.info()[1]
library(chipenrich.data);
library(mgcv);
library(IRanges);
library(GenomicRanges);
library(lattice);
library(latticeExtra);
library(grid);
library(stringr);
library(rms);
library(parallel);

SUPPORTED_METHODS = list(
  'chipenrich' = "test_gam",
  'fet' = "test_fisher_exact",
  'broadenrich' = "test_gam_ratio"
);

HIDDEN_METHODS = list(
  'binomial' = "test_binomial",
  'chipapprox' = 'test_approx',
  'broadenrich_splineless' = "test_gam_ratio_splineless"
);

METHOD_NAMES = list(
	'fet' = "Fisher's Exact Test",
	'chipenrich' = "ChIP-Enrich",
	'broadenrich' = "Broad-Enrich",
	'chipapprox' = "ChIP-Enrich Approximate",
	'broadenrich_splineless' = "Broad-Enrich Splineless"
);

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

# Used in chipenrich(...) with method!='broadenrich'
assign_peaks = function(peaks,locusdef,tss,midpoint=T) {
  results = list();

  if (!midpoint) {
    warning("peak midpoints must be used currently - may be changed in future versions");
  }

  for (chrom in names(peaks)) {
    if (!chrom %in% names(locusdef@chrom2iranges)) {
      next;
    }

    peak_ranges = peaks[[chrom]];

    # Use midpoints of peak ranges.
    peak_mids = mid(peak_ranges);
    peak_mid_ranges = IRanges(peak_mids,peak_mids);

		# Find overlapping peaks.
    gene_ranges = locusdef@chrom2iranges[[chrom]];
    overlaps = findOverlaps(peak_mid_ranges,gene_ranges);
    overlap_matrix = IRanges::as.matrix(overlaps);

    if (length(overlaps) == 0) {
			# There were no peak overlaps with gene loci on this chromosome.
      next;
    }

    matched_orig = peak_ranges[overlap_matrix[,1]];
    matched_peaks = peak_mid_ranges[overlap_matrix[,1]]; # Peaks that had an overlap with a gene locus
    matched_genes = gene_ranges[overlap_matrix[,2]]; # The gene locus assigned to each peak

		# Find the nearest TSS and the distance to it from each matched peak (above).
		# tss_ranges = tss$ranges[[chrom]];
    # nearest_tss = tss_ranges[nearest(matched_peaks,tss_ranges)];
    # dist_to_tss = distance(matched_peaks,nearest_tss);

    # Find nearest TSS, distance to it (signed by up (-) or down (+) stream.
    matched_peak_grange = GRanges(chrom,matched_peaks);
    tss_hits = tss$granges[nearest(matched_peak_grange,tss$granges)];

    neg1 = start(matched_peak_grange) < start(tss_hits) & (strand(tss_hits) == "+");
    neg2 = start(matched_peak_grange) > start(tss_hits) & (strand(tss_hits) == "-");
    neg_dist = as.logical(neg1 | neg2);

	dist_obj = distanceToNearest(matched_peak_grange,tss$granges);

	dist_to_tss = dist_obj@elementMetadata$distance;

    dist_to_tss[neg_dist] = dist_to_tss[neg_dist] * -1;

    d = data.frame(
      chrom=chrom,
      peak_start=start(matched_orig),
      peak_end=end(matched_orig),
      peak_midpoint=start(matched_peaks),
      gene_locus_start=start(matched_genes),
      gene_locus_end=end(matched_genes),
      geneid=names(matched_genes),
      nearest_tss=start(tss_hits),
      nearest_tss_gene=tss_hits$geneid,
      dist_to_tss=dist_to_tss,
      nearest_tss_gene_strand=as.character(strand(tss_hits))
    );

    # A hash on chrom, peak_start, peak_end to assign unique peak_id
    d$hash = apply(cbind(as.character(d$chrom),d$peak_start,d$peak_end),1,paste,collapse='')

		# if (midpoint) {
			# d = data.frame(
				# chrom=chrom,
				# peak_midpoint=start(matched_peaks),
				# locus_start=start(matched_genes),
				# locus_end=end(matched_genes),
				# geneid=names(matched_genes),
        # nearest_tss=start(tss_hits),
        # dist_to_tss=dist_to_tss$distance,
        # strand=as.character(strand(tss_hits))
			# );
		# } else {
			# d = data.frame(
				# chrom=chrom,
				# peak_start=start(matched_peaks),
				# peak_end=end(matched_peaks),
				# locus_start=start(matched_genes),
				# locus_end=end(matched_genes),
				# dist_to_tss=dist_to_tss$distance,
				# geneid=names(matched_genes)
			# );
		# }

    results[[chrom]] = d;
  }

  if (length(results) == 0) {
    return(NULL);
  } else{
    result = Reduce(rbind,results);

    # Give each peak a unique ID based on chrom, peak_start, peak_end
	unique_peaks = apply(unique(cbind(as.character(result$chrom),result$peak_start,result$peak_end)),1,paste,collapse='')
	unique_peaks = data.frame('hash'=unique_peaks,'peak_id'=1:length(unique_peaks),stringsAsFactors=F)

	result = merge(result,unique_peaks,by='hash')

    return(result);
  }
}

# Used in chipenrich(...) with method='broadenrich'
assign_peak_segments = function(peaks,locusdef) {
  results = list();

  for (chrom in names(peaks)) {
    if (!chrom %in% names(locusdef@chrom2iranges)) {
      next;
    }

    # Pull out IRanges objects for peaks and gene loci
    peak_ranges = peaks[[chrom]];
    gene_ranges = locusdef@chrom2iranges[[chrom]];

	# Find overlapping peaks without restricting to the midpoint
	# of the peak. NOTE: Peaks may be assigned to multiple genes.
    overlaps = findOverlaps(peak_ranges,gene_ranges);
    overlap_matrix = IRanges::as.matrix(overlaps);

    if (length(overlaps) == 0) {
	  # There were no peak overlaps with gene loci on this chromosome.
      next;
    }

    # Pull out matches
    matched_peaks = peak_ranges[overlap_matrix[,1]];
    matched_genes = gene_ranges[overlap_matrix[,2]];

    # Get peak/gene start and end points for full peaks
    peak_start = start(matched_peaks)
    peak_end = end(matched_peaks)
    gene_start = start(matched_genes)
    gene_end = end(matched_genes)

    # Get overlap start and end and compute overlap in bp
    overlap_start = apply(cbind(peak_start, gene_start), 1, max)
	overlap_end = apply(cbind(peak_end, gene_end), 1, min)
	peak_overlap = overlap_end - overlap_start

    d = data.frame(
      'chrom'=chrom,
      'peak_start'=peak_start,
      'peak_end'=peak_end,
      'gene_locus_start'=gene_start,
      'gene_locus_end'=gene_end,
      'geneid'=names(matched_genes),
      'overlap_start'=overlap_start,
      'overlap_end'=overlap_end,
      'peak_overlap'=peak_overlap
    );

    # A hash on chrom, peak_start, peak_end to assign unique peak_id
    d$hash = apply(cbind(as.character(d$chrom),d$peak_start,d$peak_end),1,paste,collapse='')

    results[[chrom]] = d;
  }

  if (length(results) == 0) {
    return(NULL);
  } else{
    result = Reduce(rbind,results);

    # Give each peak a unique ID based on chrom, peak_start, peak_end
	unique_peaks = apply(unique(cbind(as.character(result$chrom),result$peak_start,result$peak_end)),1,paste,collapse='')
	unique_peaks = data.frame('hash'=unique_peaks,'peak_id'=1:length(unique_peaks),stringsAsFactors=F)

	result = merge(result,unique_peaks,by='hash')

    return(result);
  }
}

# Used in plot_expected_peaks(...) and chipenrich(...)
# Resulting columns are:
# geneid, length, log10_length, num_peaks, peak
num_peaks_per_gene = function(assigned_peaks,locusdef,mappa=NULL) {
  # Add in gene lengths to locusdef data.
  d = locusdef@dframe;
  d$length = d$end - d$start;

  # Sum up lengths for each gene.
  d = stats::aggregate(length ~ geneid,d,sum);
  d$log10_length = log10(d$length);

  # Compute the total number of peaks assigned to each gene.
  ppg = table(assigned_peaks$geneid);
  d_ppg = data.frame(geneid=names(ppg),num_peaks=as.numeric(ppg),stringsAsFactors=F);
  result = merge(d,d_ppg,by="geneid",all.x=T);
  result[is.na(result$num_peaks),]$num_peaks = 0;

  # Mappable length if requested.
	if (!is.null(mappa)) {
    message("Mappability adjustment is enabled..");
    result = merge(result,mappa,by="geneid",sort=F);
    result$length = as.numeric((result$mappa * result$length) + 1);
    result$log10_length = log10(result$length);
  }

  # Order by number of peaks in a gene.
  result = result[order(result$num_peaks,decreasing=T),];

  # Add in peak vector (0,1).
  result$peak = as.numeric(result$num_peaks >= 1);

  return(result);
}

# Used for method='broadenrich'
# Adds peak_overlap, ratio columns to peaks per gene result
calc_peak_gene_overlap = function(assigned_peaks, ppg) {
	# Sum up the lengths for each peak in a gene
	rpg = stats::aggregate(peak_overlap ~ geneid, assigned_peaks, sum)

	d_rpg = data.frame(geneid = rpg$geneid, peak_overlap = rpg$peak_overlap, stringsAsFactors=F)

	result = merge(ppg, d_rpg, by='geneid', all.x=T)
	result$peak_overlap[is.na(result$peak_overlap)] = 0

	result$ratio = result$peak_overlap / result$length

	result$ratio[result$ratio > 1] = 1

	# Order by number of peaks in a gene.
	result = result[order(result$num_peaks,decreasing=T),]

	return(result)
}

# Used for method='chipapprox'
# Identical to post-mappa part of calc_weights_gam(...)
# Adds weight, prob_peak, resid.dev columns to peaks per gene result
calc_approx_weights = function(ppg,mappa) {
  d = ppg

  # Create model.
  model = "peak ~ s(log10_length,bs='cr')";

  # Compute binomial spline fit.
  fit = gam(as.formula(model),data=d,family="binomial");

  # Compute weights for each gene, based on the predicted prob(peak) for each gene.
  ppeak = fitted(fit);
  w0 = 1 / (ppeak/mean(d$peak,na.rm=T));
  w0 = w0 / mean(w0,na.rm=T);

  d$weight = w0;
  d$prob_peak = ppeak;
  d$resid.dev = resid(fit,type="deviance");

  cols = c("geneid","length","log10_length","mappa","num_peaks","peak","weight","prob_peak","resid.dev");
  if (is.null(mappa)) {
    cols = setdiff(cols,c("mappa"));
  }
  d = subset(d,select=cols);

  return(d);
}

# Randomize the locus definition @dframe and rebuild the @granges and @chrom2iranges
randomize_locusdef = function(ldef, resolution=50) {
	message('Extracting definition data frame..')
	ldef_df = ldef@dframe

	message('Splitting data frame on chrom..')
	ldef_chrom = split(ldef_df, f=ldef_df$chrom)

	message('Working within each chrom..')
	ldef_chrom = lapply(ldef_chrom, function(lc) {
		chrom = unique(lc$chrom)

		message(paste('	On chrom',chrom))

		# Collapse consecutive rows of the ldef when the geneid is the same
			# Split by geneid
			message('		Splitting chrom on gene id..')
			lc_geneid = split(lc, lc$geneid)

			# For each geneid df, create an IRanges object to apply IRanges::reduce to
			message('		Applying IRanges::reduce to each gene id group..')
			lc_geneid_ir = lapply(lc_geneid, function(lcg) {
				geneid = unique(lcg$geneid)

				ir = IRanges(start=lcg$start, end = lcg$end, names=lcg$geneid)
				ir = IRanges::reduce(ir)
				names(ir) = rep.int(geneid, length(ir))

				return(as.data.frame(ir, stringsAsFactors=F))
			})

		# Now put all the data.framed IRange objects back together
		message('		Collapsing reduced definition..')
		lc_collapsed = Reduce(rbind, lc_geneid_ir)

		# Add chromosome column back, rename names column, and sort columns
		message('		Formatting reduced definition..')
		lc_collapsed$chrom = chrom
		colnames(lc_collapsed) = c('start','end','width','geneid','chrom')
		lc_collapsed = lc_collapsed[,c('geneid','chrom','start','end')]

		# Sort lc_collapsed by the starting position and rename rownames
		message('		Sorting reduced definition..')
		lc_collapsed = lc_collapsed[order(lc_collapsed$start),]
		rownames(lc_collapsed) = 1:nrow(lc_collapsed)

		# Form groups
		group = floor(as.numeric(rownames(lc_collapsed))+(resolution-1))/resolution
		group = floor(group)

		# Split the chromosome into parts by group
		message(paste('		Shuffling within bins of', resolution, 'genes'))
		split_lc = split(lc_collapsed, group)
		split_lc = lapply(split_lc, function(bin){
			reordering = sample(1:nrow(bin), nrow(bin))

			# Scramble geneids
			data.frame('geneid' = bin$geneid[reordering], bin[,2:ncol(bin)], stringsAsFactors=F)
		})
		lc = Reduce(rbind, split_lc)

		return(lc)
	})

	message('Done rejiggering genomic locations for data frame..')
	ldef_df = Reduce(rbind, ldef_chrom)

	message('Creating new GenomicRanges object..')
	ldef_gr = GRanges(
		seqnames = ldef_df$chrom,
		ranges = IRanges(start=ldef_df$start, end=ldef_df$end),
		names = ldef_df$geneid
	)

	message('Creating new IRanges object..')
	chroms = c(paste('chr',1:22,sep=''),'chrX','chrY')

	chr_list = chroms
	names(chr_list) = chroms

	ldef_ir = lapply(chr_list, function(chr) {
		sub_ldef_df = subset(ldef_df, chrom==chr)

		ir = IRanges(start=sub_ldef_df$start, end=sub_ldef_df$end, names=sub_ldef_df$geneid)
		return(ir)
	})

	ldef@dframe = ldef_df
	ldef@granges = ldef_gr
	ldef@chrom2iranges = ldef_ir

	return(ldef)
}

# Randomize ppg after all additions have been made across all genes
randomize_ppg_all = function(ppg) {
	ppg = ppg[order(ppg$length),]
	rownames(ppg) = 1:nrow(ppg)

	reordering = sample(1:nrow(ppg), nrow(ppg))
	ppg = data.frame('geneid'=ppg$geneid, ppg[reordering,2:ncol(ppg)], stringsAsFactors=F)

	return(ppg)
}

# Randomize ppg after all additions have been made within length bins
randomize_ppg_length = function(ppg) {
	ppg = ppg[order(ppg$length),]
	rownames(ppg) = 1:nrow(ppg)

	group = floor(as.numeric(rownames(ppg))+99)/100
	group = floor(group)

	split_ppg = split(ppg, group)
	split_ppg = lapply(split_ppg, function(bin){
		reordering = sample(1:nrow(bin), nrow(bin))

		data.frame('geneid'=bin$geneid, bin[reordering,2:ncol(bin)], stringsAsFactors=F)
	})
	ppg = Reduce(rbind, split_ppg)

	return(ppg)
}

# Used in ..plot_spline_length(...)
add_emp_peak = function(gpw,bin_size=25) {
  d = gpw;
  d = d[order(d$log10_length),];
  d$group = ceiling((1:dim(d)[1])/bin_size);

  bygroup = stats::aggregate(peak ~ group,d,mean);
  d$emp_peak = bygroup$peak[d$group];
  d;
}

# Used in plot_spline_mappa(...)
add_emp_peak_mappa = function(gpw,bin_size=25) {
  d = gpw;
  d = d[order(d$mappa),];
  d$group = ceiling((1:dim(d)[1])/bin_size);

  bygroup = stats::aggregate(peak ~ group,d,mean);
  d$emp_peak = bygroup$peak[d$group];
  d;
}

# Used in ..plot_spline_length(...)
avg_binned_peak = function(gpw,bin_size=25) {
  d = gpw;
  d = d[order(d$log10_length),];
  d$group = ceiling((1:dim(d)[1])/bin_size);

  bygroup = stats::aggregate(cbind(peak,length) ~ group,d,mean);
  bygroup$log_avg_length = log10(bygroup$length);
  names(bygroup) = c("group","peak","avg_length","log_avg_length");
  bygroup;
}

# Used in plot_spline_mappa(...)
avg_binned_peak_mappa = function(gpw,bin_size=25) {
  d = gpw;
  d = d[order(d$mappa),];
  d$group = ceiling((1:dim(d)[1])/bin_size);

  bygroup = stats::aggregate(cbind(peak,mappa) ~ group,d,mean);
  bygroup$log_avg_mappa = log10(bygroup$mappa);
  names(bygroup) = c("group","peak","avg_mappa","log_avg_mappa");
  bygroup;
}

# Used in ..plot_gene_coverage(...)
avg_binned_coverage = function(gpw,bin_size=25) {
  d = gpw;
  d = d[order(d$log10_length),];
  d$group = ceiling((1:dim(d)[1])/bin_size);

  bygroup = stats::aggregate(cbind(ratio,log10_length) ~ group,d,mean);
  names(bygroup) = c("group","ratio","log_avg_length");
  bygroup;
}

# Never used
make_gpw = function(locusdef,peak_genes) {
  d = locusdef@dframe;

  # Indicator vector for which genes have peaks.
  d$peak = as.numeric(d$geneid %in% peak_genes);

  # Compute length and log10 length for each gene.
  d$length = d$end - d$start;

  # For genes that exist across multiple rows, sum their lengths.
  d = stats::aggregate(cbind(peak,length) ~ geneid,d,sum);
  d$log10_length = log10(d$length);

  # A gene could now have > 1 peak due to aggregating (above), reset it to be
  # 1 or 0.
  d$peak = as.integer(d$peak >= 1);

  # Sort by locus length.
  d = d[order(d$log10_length),];

  return(subset(d,select=c("geneid","length","log10_length","peak")));
}

# Used in ..plot_spline_length(...)
calc_weights_gam = function(locusdef,peak_genes,mappa=NULL,...) {
  d = locusdef@dframe;

  # Indicator vector for which genes have peaks.
  d$peak = as.numeric(d$geneid %in% peak_genes);

  # Compute length and log10 length for each gene.
  d$length = d$end - d$start;

  # For genes that exist across multiple rows, sum their lengths.
  d = stats::aggregate(cbind(peak,length) ~ geneid,d,sum);
  d$log10_length = log10(d$length);

  # A gene could now have > 1 peak due to aggregating (above), reset it to be
  # 1 or 0.
  d$peak = as.integer(d$peak >= 1);

  # Sort by locus length.
  d = d[order(d$log10_length),];

  # If mappability was requested, add that in.
  if (!is.null(mappa)) {
    d = merge(d,mappa,by="geneid",sort=F);
    d$orig_length = d$length;
    d$length = as.numeric((d$mappa * d$length) + 1);
    d$log10_length = log10(d$length);
    d = d[order(d$log10_length),];
  }

  # Create model.
  model = "peak ~ s(log10_length,bs='cr')";

  # Compute binomial spline fit.
  fit = gam(as.formula(model),data=d,family="binomial");

  # Compute weights for each gene, based on the predicted prob(peak) for each gene.
  ppeak = fitted(fit);
  w0 = 1 / (ppeak/mean(d$peak,na.rm=T));
  w0 = w0 / mean(w0,na.rm=T);

  d$weight = w0;
  d$prob_peak = ppeak;
  d$resid.dev = resid(fit,type="deviance");

  cols = c("geneid","length","log10_length","mappa","orig_length","peak","weight","prob_peak","resid.dev");
  if (is.null(mappa)) {
    cols = setdiff(cols,c("mappa","orig_length"));
  }
  d = subset(d,select=cols);
  return(d);
}

calc_spline_mappa = function(locusdef,peak_genes,mappa,...) {
  d = locusdef@dframe;

  # Indicator vector for which genes have peaks.
  d$peak = as.numeric(d$geneid %in% peak_genes);

  # Compute length and log10 length for each gene.
  d$length = d$end - d$start;

  # For genes that exist across multiple rows, sum their lengths.
  d = stats::aggregate(cbind(peak,length) ~ geneid,d,sum);
  d$log10_length = log10(d$length);

  # A gene could now have > 1 peak due to aggregating (above), reset it to be
  # 1 or 0.
  d$peak = as.integer(d$peak >= 1);

  # Merge in mappability.
  d = merge(d,mappa,by="geneid",sort=F);

  # Sort by mappa.
  d = d[order(d$mappa),];

  # Create model.
  model = "peak ~ s(mappa,bs='cr')";

  # Compute binomial spline fit.
  fit = gam(as.formula(model),data=d,family="binomial");

  # Compute weights for each gene, based on the predicted prob(peak) for each gene.
  d$prob_peak = fitted(fit);
  return(subset(d,select=c("geneid","length","log10_length","mappa","peak","prob_peak")));
}

test_fisher_exact = function(geneset,gpw,alternative="two.sided") {
  # Restrict to only those genes in the genesets.
  gpw = subset(gpw,geneid %in% geneset@all.genes);

  genes = gpw$geneid;
  peaks = gpw$peak;

  results = lapply(ls(geneset@set.gene), function(go_term) {
    go_genes = geneset@set.gene[[go_term]];
    n_go_genes = length(go_genes);

    in_cat = as.numeric(genes %in% go_genes);
    n_go_peak_genes = sum((genes %in% go_genes) & (peaks == 1));
    xt = table(in_cat,peaks);

    pval = 1;
    odds_ratio = 0;
    try({
      fet_result = fisher.test(xt,conf.int=F,alternative=alternative);
      pval = fet_result$p.value;
      odds_ratio = fet_result$estimate;
    },silent=T);

		go_peak_genes = gpw[(gpw$geneid %in% go_genes) & (gpw$peak == 1),]$geneid;
		go_peak_genes_str = paste(go_peak_genes,collapse=", ");

    enr = NA;

    if (odds_ratio > 1) {
      enr = "enriched";
    }

    if (odds_ratio < 1) {
      enr = "depleted";
    }

    if (odds_ratio == 0) {
      enr = "depleted";
    }

    if (is.infinite(odds_ratio)) {
      enr = "enriched";
    }

    data.frame(
      "Geneset ID" = go_term,
      "N Geneset Genes" = n_go_genes,
      "N Geneset Peak Genes" = n_go_peak_genes,
			"Geneset Peak Genes" = go_peak_genes_str,
      "Odds.Ratio" = odds_ratio,
      "Status" = enr,
      "P-value" = pval,
      stringsAsFactors=F
    );
  });

  results = rbind.fill(results);

  results$FDR = p.adjust(results$P.value,method="BH");

  results = results[order(results$P.value),];
  return(results);
}

test_binomial = function(geneset,ppg) {
  # Restrict to genes only within the genesets.
  ppg = subset(ppg,geneid %in% geneset@all.genes);

  # Total # of peaks.
  total_peaks = sum(ppg$num_peaks);

  # If there are no total peaks, it means there were no genes w/ peaks that overlap
  # with all genes across all gene sets. We can't analyze in this case.
  if (total_peaks == 0) {
    stop("None of the peaks in this dataset were assigned to any of the genes across all pathways/genesets being tested.");
  }

  # Gene IDs must be character to be stored into an environment as names.
  ppg$geneid = as.character(ppg$geneid);

  # Index gene --> peak count.
  gene2peakc = new.env(parent=emptyenv());
  apply(ppg,1,function(x) assign(x["geneid"],x["num_peaks"],envir=gene2peakc));

  # Index gene --> log10 length.
  gene2length = new.env(parent=emptyenv());
  apply(ppg,1,function(x) assign(x["geneid"],x["length"],envir=gene2length));

  # Length of the genome.
  genome_length = sum(as.numeric(ppg$length));

  estimates = c();
  pvals = c();
  goids = c();
  gonums = c();
  n_go_peak_genes = c();
	r_go_peak_genes = c();
  for (i in 1:length(geneset@set.gene)) {
    go_term = ls(geneset@set.gene)[i];
    go_genes = as.character(geneset@set.gene[[go_term]]);

		# Eliminate GO genes that aren't in the ppg.
		#go_genes = Filter(function(x) x %in% ppg$geneid,go_genes);
    go_genes = go_genes[go_genes %in% ppg$geneid];

    go_gene_lengths = as.numeric(sapply(go_genes,function(x) get(x,gene2length)));
    go_term_length = sum(go_gene_lengths);

    go_gene_counts = as.numeric(sapply(go_genes,function(x) get(x,gene2peakc)));
    go_term_count = sum(go_gene_counts);

    binom_p = go_term_length / genome_length;
    btest = binom.test(go_term_count,total_peaks,binom_p,alternative="greater");

    goids[i] = go_term;
    pvals[i] = btest$p.value;
    estimates[i] = btest$estimate;
    gonums[i] = length(go_genes);

		go_peak_genes = ppg[(ppg$geneid %in% go_genes) & (ppg$num_peaks > 0),]$geneid;
		r_go_peak_genes[i] = paste(go_peak_genes,collapse=", ");
		n_go_peak_genes[i] = length(go_peak_genes);
  }

  fdr = p.adjust(pvals,method="BH");

  results = data.frame(
    "Geneset ID"=goids,
    "N Geneset Genes"=gonums,
    "N Geneset Peak Genes"=n_go_peak_genes,
    "P.Success"=estimates,
    "P-value"=pvals,
    "FDR"=fdr,
		"Geneset Peak Genes"=r_go_peak_genes,
    stringsAsFactors=F
  );

  results = results[order(results$P.value),];
  return(results);
}

# Extracted contents of for loop implementation in test_gam and test_gam_ratio to
# enable multicore processing of genesets. This function gets called in mclapply
# in both test_gam and test_gam_ratio, while using the correct model, and reporting
# the correct information based on the method parameter.
single_gam = function(go_id, geneset, gpw, method, model) {

	final_model = as.formula(model);

	# Genes in the geneset
	go_genes = geneset@set.gene[[go_id]];

	# Filter genes in the geneset to only those in the gpw table.
	# The gpw table will be truncated depending on which geneset type we're in.
	go_genes = go_genes[go_genes %in% gpw$geneid];

	# Background genes and the background presence of a peak
	b_genes = gpw$geneid %in% go_genes;
	sg_go = gpw$peak[b_genes];

 	# Information about the geneset
	r_go_id = go_id;
	r_go_genes_num = length(go_genes);
	r_go_genes_avg_length = mean(gpw$length[b_genes]);

    # Information about peak genes
	go_genes_peak = gpw$geneid[b_genes][sg_go==1];
	r_go_genes_peak = paste(go_genes_peak,collapse=", ");
	r_go_genes_peak_num = length(go_genes_peak);

    # Information specific to broadenrich
	if(method == 'broadenrich' || method == 'broadenrich_splineless') {
		r_go_genes_avg_coverage = mean(gpw$ratio[b_genes]);
	}

    # Small correction for case where every gene in this geneset has a peak.
	if (all(as.logical(sg_go))) {
	  cont_length = quantile(gpw$length,0.0025);

	  if(method == 'chipenrich') {
		  cont_gene = data.frame(
			geneid = "continuity_correction",
			length = cont_length,
			log10_length = log10(cont_length),
			num_peaks = 0,
			peak = 0,
			stringsAsFactors = F
		  );
	  } else if (method == 'broadenrich' || method == 'broadenrich_splineless') {
		  cont_gene = data.frame(
			geneid = "continuity_correction",
			length = cont_length,
			log10_length = log10(cont_length),
			num_peaks = 0,
			peak = 0,
			peak_overlap = 0,
			ratio = 0,
			stringsAsFactors = F
		  );
	  }

	  if ("mappa" %in% names(gpw)) {
		cont_gene$mappa = 1;
	  }
	  gpw = rbind(gpw,cont_gene);
	  b_genes = c(b_genes,1);

	  message(sprintf("Applying correction for geneset %s with %i genes...",go_id,length(go_genes)));
	}

    # Logistic regression works no matter the method because final_model is chosen above
    # and the data required from gpw will automatically be correct based on the method used.
    fit = gam(final_model,data=cbind(gpw,goterm=as.numeric(b_genes)),family="binomial");

	# Results from the logistic regression
    r_effect = coef(fit)[2];
    r_pval = summary(fit)$p.table[2,4];

	# The only difference between chipenrich and broadenrich here is
	# the Geneset Avg Gene Coverage column
	if(method == 'chipenrich') {
		out = data.frame(
			"P.value"=r_pval,
			"Geneset ID"=r_go_id,
			"N Geneset Genes"=r_go_genes_num,
			"Geneset Peak Genes"=r_go_genes_peak,
			"N Geneset Peak Genes"=r_go_genes_peak_num,
			"Effect"=r_effect,
			"Odds.Ratio"=exp(r_effect),
			"Geneset Avg Gene Length"=r_go_genes_avg_length,
			stringsAsFactors=F);
	} else if (method == 'broadenrich' || method == 'broadenrich_splineless') {
		out = data.frame(
			"P.value"=r_pval,
			"Geneset ID"=r_go_id,
			"N Geneset Genes"=r_go_genes_num,
			"Geneset Peak Genes"=r_go_genes_peak,
			"N Geneset Peak Genes"=r_go_genes_peak_num,
			"Effect"=r_effect,
			"Odds.Ratio"=exp(r_effect),
			"Geneset Avg Gene Length"=r_go_genes_avg_length,
			"Geneset Avg Gene Coverage"=r_go_genes_avg_coverage,
			stringsAsFactors=F);
	}

	return(out);
}

test_gam = function(geneset,gpw,n_cores) {
  # Restrict our genes/weights/peaks to only those genes in the genesets.
  # Here, geneset is not all combined, but GOBP, GOCC, etc.
  # i.e. A specific one.
  gpw = subset(gpw,geneid %in% geneset@all.genes);

  if (sum(gpw$peak) == 0) {
    stop("Error: no peaks in your data!");
  }

  # Construct model formula.
  model = "peak ~ goterm + s(log10_length,bs='cr')";

  # Run tests on genesets and beware of Windows!
  if(os != 'Windows' && n_cores > 1) {
	  results_list = mclapply(as.list(ls(geneset@set.gene)), function(go_id) {
		single_gam(go_id, geneset, gpw, 'chipenrich', model)
	  }, mc.cores = n_cores)
  } else {
  	  results_list = lapply(as.list(ls(geneset@set.gene)), function(go_id) {
		single_gam(go_id, geneset, gpw, 'chipenrich', model)
	  })
  }

  # Collapse results into one table
  results = Reduce(rbind,results_list)

  # Correct for multiple testing
  results$FDR = p.adjust(results$P.value, method="BH");

  # Create enriched/depleted status column
  results$Status = ifelse(results$Effect > 0, 'enriched', 'depleted')

  results = results[order(results$P.value),];

  return(results);
}

test_gam_ratio = function(geneset,gpw,n_cores) {

  # Restrict our genes/weights/peaks to only those genes in the genesets.
  gpw = subset(gpw,geneid %in% geneset@all.genes);

  if (sum(gpw$peak) == 0) {
    stop("Error: no peaks in your data!");
  }

  # Construct model formula.
  model = "goterm ~ ratio + s(log10_length,bs='cr')";

  # Run multicore tests on genesets and beware of Windows!
  if(os != 'Windows' && n_cores > 1) {
	  results_list = mclapply(as.list(ls(geneset@set.gene)), function(go_id) {
		single_gam(go_id, geneset, gpw, 'broadenrich', model)
	  }, mc.cores = n_cores)
  } else {
  	  results_list = lapply(as.list(ls(geneset@set.gene)), function(go_id) {
		single_gam(go_id, geneset, gpw, 'broadenrich', model)
	  })
  }

  # Collapse results into one table
  results = Reduce(rbind,results_list)

  # Correct for multiple testing
  results$FDR = p.adjust(results$P.value, method="BH");

  # Create enriched/depleted status column
  results$Status = ifelse(results$Effect > 0, 'enriched', 'depleted')

  results = results[order(results$P.value),];

  return(results);
}

test_gam_ratio_splineless = function(geneset,gpw,n_cores) {
  message('Using test_gam_ratio_splineless..')

  # Restrict our genes/weights/peaks to only those genes in the genesets.
  gpw = subset(gpw,geneid %in% geneset@all.genes);

  if (sum(gpw$peak) == 0) {
    stop("Error: no peaks in your data!");
  }

  # Construct model formula.
  model = "goterm ~ ratio";

  # Run multicore tests on genesets and beware of Windows!
  if(os != 'Windows' && n_cores > 1) {
	  results_list = mclapply(as.list(ls(geneset@set.gene)), function(go_id) {
		single_gam(go_id, geneset, gpw, 'broadenrich_splineless', model)
	  }, mc.cores = n_cores)
  } else {
  	  results_list = lapply(as.list(ls(geneset@set.gene)), function(go_id) {
		single_gam(go_id, geneset, gpw, 'broadenrich_splineless', model)
	  })
  }

  # Collapse results into one table
  results = Reduce(rbind,results_list)

  # Correct for multiple testing
  results$FDR = p.adjust(results$P.value, method="BH");

  # Create enriched/depleted status column
  results$Status = ifelse(results$Effect > 0, 'enriched', 'depleted')

  results = results[order(results$P.value),];

  return(results);
}

single_approx = function(go_id, geneset, gpw) {

	lrm.fast = function(x,y) {
		fit = lrm.fit(x,y);
		vv = diag(fit$var);
		cof = fit$coef;
		z = cof/sqrt(vv);
		pval = pchisq(z^2,1,lower.tail=F);
		c(cof[2],pval[2]);
	}

	# Genes in the geneset
	go_genes = geneset@set.gene[[go_id]];

	# Filter genes in the geneset to only those in the gpw table.
	# The gpw table will be truncated depending on which geneset type we're in.
	go_genes = go_genes[go_genes %in% gpw$geneid];

	# Background genes, the background presence of a peak, and the background
	# weight of peaks
	b_genes = gpw$geneid %in% go_genes;
	sg_go = gpw$peak[b_genes];
	wg_go = gpw$weight[b_genes];

 	# Information about the geneset
	r_go_id = go_id;
	r_go_genes_num = length(go_genes);
	r_go_genes_avg_length = mean(gpw$length[b_genes]);

    # Information about peak genes
	go_genes_peak = gpw$geneid[b_genes][sg_go==1];
	r_go_genes_peak = paste(go_genes_peak,collapse=", ");
	r_go_genes_peak_num = length(go_genes_peak);

	# Small correction for case where every gene in this geneset has a peak.
	if (all(as.logical(sg_go))) {
		cont_length = quantile(gpw$length,0.0025);
		cont_gene = data.frame(
			geneid = "continuity_correction",
			length = cont_length,
			log10_length = log10(cont_length),
			num_peaks = 0,
			peak = 0,
			weight = quantile(gpw$weight,0.0025),
			prob_peak = quantile(gpw$prob_peak,0.0025),
			resid.dev = quantile(gpw$resid.dev,0.0025),
			stringsAsFactors = F
		);
		if ("mappa" %in% names(gpw)) {
			cont_gene$mappa = 1;
		}
		gpw = rbind(gpw,cont_gene);
		b_genes = c(b_genes,1);

		message(sprintf("Applying correction for geneset %s with %i genes...",go_id,length(go_genes)));
	}

	# The model is still essentially peak ~ goterm + s(log10_length,bs='cr')
	# except we're using the weights that are calculated once...
	# Also, y must be binary, so as.numeric(b_genes) (goterm) must be
	# predicted against the weights...
	testm = cbind(y=as.numeric(b_genes), x=gpw$weight*gpw$peak);
	ep = lrm.fast(testm[,"x"], testm[,"y"]);

	# Results from quick regression
	r_effect = ep[1];
	r_pval = ep[2];

	out = data.frame(
		"P.value"=r_pval,
		"Geneset ID"=r_go_id,
		"N Geneset Genes"=r_go_genes_num,
		"Geneset Peak Genes"=r_go_genes_peak,
		"N Geneset Peak Genes"=r_go_genes_peak_num,
		"Effect"=r_effect,
		"Odds.Ratio"=exp(r_effect),
		"Geneset Avg Gene Length"=r_go_genes_avg_length,
		stringsAsFactors=F);

	return(out)
}

test_approx = function(geneset,gpw,nwp=F,n_cores) {

	if (!"weight" %in% names(gpw)) {
    	stop("Error: you must fit weights first using one of the calc_weights* functions.");
	}

	if (sum(gpw$peak) == 0) {
		stop("Error: no peaks in your data!");
	}

	# Restrict our genes/weights/peaks to only those genes in the genesets.
	gpw = subset(gpw,geneid %in% geneset@all.genes);

	# Re-normalize weights.
	# Not sure what nwp means, nor if we want to modify gpw like this.
	if (!nwp) {
		gpw$weight = gpw$weight / mean(gpw$weight);
	} else {
		b_haspeak = gpw$peak == 1;
		gpw$weight[b_haspeak] = gpw$weight[b_haspeak] / mean(gpw$weight[b_haspeak]);
	}

	# Run multicore tests on genesets and beware of Windows!
	if(os != 'Windows' && n_cores > 1) {
		results_list = mclapply(as.list(ls(geneset@set.gene)), function(go_id) {
			single_approx(go_id, geneset, gpw)
		}, mc.cores = n_cores)
	} else {
		results_list = lapply(as.list(ls(geneset@set.gene)), function(go_id) {
			single_approx(go_id, geneset, gpw)
		})
	}

	# Collapse results into one table
	results = Reduce(rbind,results_list)

	# Correct for multiple testing
	results$FDR = p.adjust(results$P.value, method="BH");

	# Create enriched/depleted status column
	results$Status = ifelse(results$Effect > 0, 'enriched', 'depleted')

	results = results[order(results$P.value),];

	return(results);
}

# Quick function to change column names of data frame
# (by name, not by position.)
change_names = function(dframe,name_list) {
  for (n in names(name_list)) {
    names(dframe)[names(dframe)==n] = name_list[[n]];
  }
  dframe;
}

# Remove a file extension from a file name.
strip_ext = function(file) {
  s = unlist(strsplit(file,".",fixed=T));
  s = s[1:(length(s)-1)];
  s = paste(s,collapse=".");
  return(s);
}

# Get a file's extension.
get_ext = function(file) {
  s = unlist(strsplit(file,".",fixed=T));
  return(tail(s,1));
}

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

# Do the below, but from "scratch"
plot_spline_length = function(peaks,locusdef="nearest_tss",genome='hg19',use_mappability=F,read_length=36,legend=T,xlim=NULL) {
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

	# Load TSS site info.
  tss_code = sprintf("tss.%s",genome);
  data(list=tss_code,package = "chipenrich.data");
  tss = get(tss_code);

  # Load mappability if requested.
  if (use_mappability) {
    mappa_code = sprintf("mappa.%s.%s.%imer",genome,locusdef,read_length);
    data(list=mappa_code,package = "chipenrich.data");
    mappa = get(mappa_code);
  } else {
    mappa = NULL;
  }

  # Assign peaks to genes.
  assigned_peaks = assign_peaks(peakobj,ldef,tss);
  peak_genes = unique(assigned_peaks$geneid);

	# Make plot.
	plotobj = ..plot_spline_length(ldef,peak_genes,num_peaks,mappa,legend=legend,xlim=xlim);
	return(plotobj);
}

# Create diagnostic plot of Proportion of Peaks from your data against log locus length (of genes)
# along with the spline fit.
..plot_spline_length = function(locusdef,peak_genes,num_peaks,mappa=NULL,legend=T,xlim=NULL) {
  # Calculate smoothing spline fit.
  gpw = calc_weights_gam(locusdef,peak_genes,mappa=mappa); # gpw = genes, peaks, weights

  # Genome length.
  genome_length = sum(as.numeric(gpw$length));

  # Add in empirical peak by binning.
  gpw = add_emp_peak(gpw);

  # Average peak/lengths.
  avg_bins = avg_binned_peak(gpw,bin_size=25);

  # Order by length.
  gpw = gpw[order(gpw$log10_length),];

  # Calculate prob(all false positives).
  gpw$false_prob = 1 - (1 - (gpw$length/genome_length))^(num_peaks);

  col_rand_gene = "grey35";
  col_rand_peak = "grey74";
  col_spline = "darkorange";

  panel_func = function(x,y,...) {
    panel.xyplot(x,y,...);
    panel.abline(h=mean(y),col=col_rand_gene,lwd=3);
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
      col = c(col_rand_gene,col_rand_peak,col_spline,"black"),
      lwd = 3
    ),
    cex = 1.2
  );

  if (!legend) {
    custom_key = NULL;
  }

  if (!is.null(mappa)) {
    xlab = expression(paste(Log[10]," Mappable Locus Length"))
  } else {
    xlab = expression(paste(Log[10]," Locus Length"))
  }

	xmin_nopad = base::ifelse(is.null(xlim[1]),floor(min(gpw$log10_length)),floor(xlim[1]));
	xmax_nopad = base::ifelse(is.null(xlim[2]),ceiling(max(gpw$log10_length)),ceiling(xlim[2]));

	scales = list(
		x = list(
			axs = 'i',
			at = seq(xmin_nopad,xmax_nopad,1)
		),
		y = list(
			axs = 'i',
			at = seq(0,1,0.2)
		)
	);

  plotobj = xyplot(
    false_prob + prob_peak ~ log10_length,
    gpw,
    xlab=list(label=xlab,cex=1.4),
    ylab=list(label="Proportion of Peaks",cex=1.4),
    ylim=c(-0.05,1.05),
		xlim=c(xmin_nopad - 0.5,xmax_nopad + 0.5),
    panel=panel_func,
    type="l",
    lwd=3,
    key=custom_key,
    scales=list(cex=1.4),
    par.settings=simpleTheme(pch=c(15,17),col=c(col_rand_peak,col_spline))
  ) + as.layer(xyplot(peak ~ log_avg_length,avg_bins,pch=20,cex=0.4,col="black"));

  #trellis.focus("panel",1,1,highlight=F);
  #gpw_peak = subset(gpw,peak == 1);
  #gpw_nopeak = subset(gpw,peak == 0);
  #panel.rug(gpw_peak$log10_length,regular=F,col="black");
  #panel.rug(gpw_nopeak$log10_length,regular=T,col="black");
  #trellis.unfocus();

  return(plotobj);
}

plot_spline_mappa = function(locusdef,peak_genes,mappa) {
  # Calculate spline for mappability.
  gpw = calc_spline_mappa(locusdef,peak_genes,mappa);

  # Add in empirical peak by binning.
  gpw = add_emp_peak_mappa(gpw);

  # Average peak/lengths.
  avg_bins = avg_binned_peak_mappa(gpw,bin_size=25);

  # Order by mappability.
  gpw = gpw[order(gpw$mappa),];

  panel_func = function(x,y,...) {
    panel.xyplot(x,y,...);
    panel.abline(h=mean(y),col="blue",lwd=3);
  }

  custom_key = list(
    text = list(
      c(
        "Mean(peak)",
        "Spline Fit - Mappability Only",
        "Proportion of Peaks in Bin"
      )
    ),
    lines = list(
      pch = c(20,17,20),
      type = c("l","l","p"),
      col = c("blue","green","black")
    )
  );

  plotobj = xyplot(
    prob_peak ~ mappa,
    gpw,
    xlab=list(label="Mappability",cex=1.25),
    ylab=list(label="Proportion of Peaks",cex=1.25),
    ylim=c(-0.05,1.05),
    panel=panel_func,
    type="l",
    lwd=3,
    key=custom_key,
    scales=list(cex=1.1),
    par.settings=simpleTheme(pch=c(17),col=c("green")),
  );

  trellis.focus("panel",1,1,highlight=F);
  grid.points(avg_bins$avg_mappa,avg_bins$peak,pch=20,gp=gpar(cex=0.35));
  trellis.unfocus();

   return(plotobj);
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

genome_to_organism = function(genome) {
  code = substr(as.character(genome),1,2);
  if (code == 'mm') {
    org = 'mmu';
  } else if (code == 'hg') {
    org = 'hsa';
  } else if (code == 'rn') {
    org = 'rno';
  } else if (code == 'dm') {
  	org = 'dme';
  }
  else {
    org = NULL;
  }

  if (is.null(org)) {
    stop("Error: genome requested is not supported.");
  }

  org;
}

supported_locusdefs = function() {
  piqr = data(package = "chipenrich.data");
  data_files = piqr$results[,3];

  ldefs = grep("locusdef",data_files,value=T);
  unique(str_replace(ldefs,"locusdef\\.(.+?)\\.",""));
}

supported_read_lengths = function() {
  piqr = data(package = "chipenrich.data");
  data_files = piqr$results[,3];

  mappas = grep("mappa",data_files,value=T);
  sort(unique(as.numeric(str_replace(str_extract(mappas,"(\\d+)mer"),"mer",""))));
}

supported_genesets = function() {
  piqr = data(package = "chipenrich.data");
  data_files = piqr$results[,3];

  geneset_files = grep("geneset",data_files,value=T);
  unlist(unique(Map(function(x) x[[2]],str_split(geneset_files,"\\."))));
}

supported_genomes = function() {
  piqr = data(package = "chipenrich.data");
  data_files = piqr$results[,3];

  as.character(na.omit(unique(str_match(data_files,"(locusdef|mappa)\\.(\\w+)\\.")[,3])));
}

supported_methods = function() {
  return(names(SUPPORTED_METHODS));
}

# Checks to see if all elements of a list argument are possible.
check_arg = function(arg,possible_values,value=F) {
  if (!all(arg %in% possible_values)) {
    if (value) {
      return(arg[!arg %in% possible_values]);
    } else {
      return(F);
    }
  } else {
    return(T);
  }
}

read_mappa = function(file_path) {
  if (!file.exists(file_path)) {
    stop("Error: could not find mappability file: ",file_path);
  }

	d = read.table(file_path,sep="\t",header=T,stringsAsFactors=F);

	# Check columns.
	for (col in c("geneid","mappa")) {
		if (!col %in% names(d)) {
			stop(sprintf("Error reading mappability data: no column named '%s' in file.",col));
		}
	}

	# Genes in this file should not be duplicated.
	if (sum(duplicated(d$geneid)) > 0) {
		stop("Error reading mappability data: duplicate geneids exist in file.");
	}

  # Mappability should be between 0 and 1.
	if (min(d$mappa) < 0) {
    stop("Error: mappability must be >= 0.");
  }

  if (max(d$mappa) > 1) {
    stop("Error: mappability must be <= 1.");
  }

	return(d);
}

filter_genesets = function(x,max_geneset_size = 2000) {
  x_class = class(x);
  x_class_attr = attr(x_class,"package");

  if (is.null(x_class) || is.null(x_class_attr)) {
    stop("Error: bad geneset object in filtering step..");
  }

  if (x_class != "GeneSet" || x_class_attr != "chipenrich.data") {
    stop("Error: bad geneset object in filtering step..");
  }

  g = as.list(x@set.gene);
  g = Filter(function(x) length(x) <= max_geneset_size,g);
  x@set.gene = as.environment(g);

  return(x);
}

# Recodes # of peaks to be:
# { 0, if no peaks
# { 1, if number of peaks is >= threshold argument
recode_peaks = function(num_peaks,threshold=1) {
  as.numeric(num_peaks >= threshold);
}

# Creates a LocusDefinition object given a flat file of definitions.
setup_ldef = function(filepath) {
  # Check if file exists.
  if (!file.exists(filepath)) {
    stop("Error: unable to open file for reading: ",filepath);
  }

  # Read in flat file.
  message("Reading user-specified gene locus definitions: ",filepath);
  d = read.table(filepath,sep="\t",header=T,stringsAsFactors=F);

  # Warn if missing entries, we'll ignore them though.
  num_missing = sum(!complete.cases(d));
  if (num_missing > 0) {
    message(sprintf("Warning: %i rows of user-provided locus definition were missing values, skipping these rows..",num_missing));
  }

  # Remove rows with missing data.
  d = na.omit(d);

  filename = strip_ext(filepath);

  # Create new locus definition object.
  object = new("LocusDefinition");

  # Remove duplicated rows, not duplicated genes.
  # Some genes will exist on multiple rows because they are split by
  # other transcripts, or small nuclear RNAs, etc.
  d = unique(d);

  # Create an IRanges object representing the loci for each gene on that chromosome.
  chroms = list();
  for (chr in unique(d$chrom)) {
    genes_chrom = subset(d,chrom == chr);
    chroms[[chr]] = IRanges(start=genes_chrom$start,end=genes_chrom$end,names=genes_chrom$geneid);
  }

  object@dframe = d;
  object@chrom2iranges = chroms;

  # Store as GRanges object as well, for convenience.
  object@granges = GRanges(
    seqnames=d$chrom,
    ranges=IRanges(d$start,d$end),
    names=d$geneid
  );

  # Check to make sure locus definitions are disjoint - that is, they do not overlap each other.
  # SLATED FOR DELETION PENDING CONSEQUENCES
#  if(!isDisjoint(object@granges)) {
#    stop("Error: user-provided locus definitions overlap - there should be disjoint ranges for all genes.");
#  }

  object;
}

# Creates an object that mimics the GeneSet class
# from the chipenrich.data package, as in:
#    setClass("GeneSet",representation(
#      set.gene = "environment",
#      type = "character",
#      set.name = "environment",
#      all.genes = "character",
#      organism = "character",
#      dburl = "character"
#    ),prototype(
#      set.gene = new.env(parent=emptyenv()),
#      type = "",
#      set.name = new.env(parent=emptyenv()),
#      all.genes = "",
#      organism = "",
#      dburl = ""
#    ),
#      package = "chipenrich.data"
#    )
# User supplied genesets are expected to come as two columned,
# tab-delimited text files with the first column being the geneset
# name and the second column being Entrez Gene ID
setup_geneset = function(filepath) {

    # Check if file exists.
    if (!file.exists(filepath)) {
        stop("Error: unable to open file for reading: ",filepath)
    }

    # Read in flat file.
    message("Reading user-specified gene set definitions: ",filepath)
    d = read.table(filepath,sep="\t",header=T,stringsAsFactors=F)

    # Warn if missing entries, we'll ignore them though.
    num_missing = sum(!complete.cases(d));
    if (num_missing > 0) {
        message(sprintf("Warning: %i rows of user-provided locus definition were missing values, skipping these rows..",num_missing))
    }

    # Remove rows with missing data.
    d = na.omit(d)

    filename = strip_ext(filepath)

    # Create shell for geneset list
    gs.names = unique(d[,1])
    gs = as.list(gs.names)
    names(gs) = gs.names

    gs.names.list = gs

    # Populate the shell
    gs = lapply(gs, function(g){
        return(unique(subset(d, d[,1] == g)[,2]))
    })

    # Create new GeneSet object
    object = new('GeneSet')

    # Populate the GeneSet object
    object@type = 'user-supplied'
    object@organism = 'user-supplied'
    object@dburl = 'user-supplied'

    object@set.gene = as.environment(gs)
    object@all.genes = as.character(sort(Reduce(function(x,y) union(x,y),object@set.gene)))

    object@set.name = as.environment(gs.names.list)

    return(object)
}

# d_fwf = function(d) {
  # d = colwise(as.character)(d);
  # widths = sapply(d,function(x) max(nchar(x)));
  # widths + 2;

  # final = NULL;
  # for (i in 1:dim(d)[2]) {
    # n = names(d)[i];
    # sf = sprintf("%%-0%is",widths[i]);
    # v = sprintf(sf,d[,n]);
    # final = cbind(final,v);
  # }

  # final = data.frame(final,stringsAsFactors=F);
  # names(final) = names(d);

  # final;
# }

chipenrich = function(
  peaks,
  out_name = "chipenrich",
  out_path = getwd(),
  genome = "hg19",
  genesets = c(
    'GOBP',
    'GOCC',
    'GOMF'
  ),
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

  if(os == 'Windows') {
  	message('Warning! Multicore enrichment is not supported on Windows.')
  }

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

  # Assign all chipenrich opts to a special environment
  # This prevents huge opts outputs for particular modes of running chipenrich
  opts_env = new.env()
  assign('os' = os, envir = opts_env)
  assign('peaks' = peaks, envir = opts_env)
  assign('out_name' = out_name, envir = opts_env)
  assign('out_path' = out_path, envir = opts_env)
  assign('genome' = genome, envir = opts_env)
  assign('genesets' = genesets, envir = opts_env)
  assign('locusdef' = locusdef, envir = opts_env)
  assign('method' = method, envir = opts_env)
  assign('fisher_alt' = fisher_alt, envir = opts_env)
  assign('use_mappability' = use_mappability, envir = opts_env)
  assign('mappa_file' = mappa_file, envir = opts_env)
  assign('read_length' = read_length, envir = opts_env)
  assign('qc_plots' = qc_plots, envir = opts_env)
  assign('max_geneset_size' = max_geneset_size, envir = opts_env)
  assign('num_peak_threshold' = num_peak_threshold, envir = opts_env)
  assign('n_cores' = n_cores, envir = opts_env)

  opts_l = unlist(as.list(opts_env))

  opts = data.frame(
  	args = names(opts_l),
  	values = sapply(unlist(opts_l), paste, collapse=","),
  	stringsAsFactors = F)
  rownames(opts) = 1:length(opts_l)

  # Check genome.
  if (!genome %in% supported_genomes()) {
    stop("genome not supported: ",genome);
  }

  # Organism - human or mouse.
  organism = genome_to_organism(genome);

  # Check genesets. The API for a user to use their own genesets
  # will be to put a path in the genesets argument.
  if (!check_arg(genesets,supported_genesets())) {
    bad_args = check_arg(genesets,supported_genesets(),value=T);

    # If the bad_args is a path that exists, then we know the user wants
    # to provide their own genesets
    if(file.exists(genesets)) {
        message('Will use genesets specific by the user...')
    } else {
        stop("bad geneset(s) requested: ",paste(bad_args,collapse=", "))
    }
  }

  # Check and load locus definition.
  user_defined_ldef = file.exists(locusdef);
  if (user_defined_ldef) {
    # Load user-defined locus definition file.
    ldef = setup_ldef(locusdef);
  } else {
    if (!locusdef %in% supported_locusdefs()) {
      stop("Error: invalid definition requested: ",locusdef);
    }

    # Load locus definitions.
    ldef_code = sprintf("locusdef.%s.%s",genome,locusdef);
    data(list=ldef_code,package = "chipenrich.data");
    ldef = get(ldef_code);

	# Randomize locus definition if rndloc == T
	if(rndloc) {
		ldef = randomize_locusdef(ldef, 50)
	}
  }

  # If the user specified their own mappability, we can't use the built in mappability - they must provide their own.
  user_defined_mappa = !is.null(mappa_file);

  # The user specified mappability - if they accidentally also thought use_mappability had to be enabled, just disable it.
  if (user_defined_mappa & use_mappability) {
    use_mappability = FALSE;
  }

  if (user_defined_ldef) {
    if (use_mappability) {
      message("Warning: built-in mappability cannot be used with a user-defined locus definition, you must calculate mappability for your definition and pass it in with the mappa_file argument.");
      use_mappability = FALSE;
    }
  }

  # Check read length.
  if (use_mappability) {
    if (!as.numeric(read_length) %in% supported_read_lengths()) {
      stop("Error: bad read length requested: ",read_length);
    }
  }

  # Load mappability if requested.
  if (user_defined_mappa) {
    message("Reading user-specified gene locus mappability file: ",mappa_file);
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

  # If they specified their own mappability, we should check that the gene names overlap.
  if (user_defined_mappa && user_defined_ldef) {
    total_unique_genes = union(mappa$geneid,ldef@dframe$geneid);
    mappa_ldef_inters = intersect(mappa$geneid,ldef@dframe$geneid);
    frac_overlap = length(mappa_ldef_inters) / length(total_unique_genes);

    if (frac_overlap < 0.95) {
      stop("Error: your mappability genes and locus definition genes overlap by less than 95% (they should match almost exactly)..");
    }
  }

  # Testing function.
  get_test_method = function(x) {
    if (method %in% names(SUPPORTED_METHODS)) {
      return(SUPPORTED_METHODS[[method]]);
    } else if (method %in% names(HIDDEN_METHODS)) {
      return(HIDDEN_METHODS[[method]]);
    } else {
      stop(sprintf("Error: invalid enrichment test requested: %s, contact developer.",method));
    }
  }
  testf = get_test_method(method);
  test_func = get(testf);

	# Test name.
	method_name = METHOD_NAMES[[method]];

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

  # Warn user if they are trying to use FET with a
  # locus definition that might lead to biased results.
  #if ((method == "fet") & (!locusdef %in% c("1kb","5kb"))) {
  if (method == "fet") {
    if (is.character(locusdef) && !locusdef %in% c("1kb","5kb")) {
      message("");
      message("Warning: Fisher's exact test should only be used with the 1kb or 5kb locus definition.");
      message("");
    }

    if (user_defined_ldef) {
      message("");
      message("Warning: Fisher's exact test may give biased results if the spline fit for the gene locus definitions is not flat (see QC plots.)");
      message("");
    }
  }

  # Warn user if they are using the binomial test.
  if (method == "binomial") {
    message("");
    message("Warning: the binomial test is provided for comparison purposes only.");
    message("This test will almost always give biased results favoring gene sets with short average locus length.");
    message("");
  }

	# Number of peaks in data.
  num_peaks = sum(sapply(peakobj,function(x) length(x)))

  # Load genesets and filter them
  if(!file.exists(genesets)) {
      # If the user does not provide a path to the geneset.
      # That is, do the normal thing
      geneset_list = list();
      for (gs in genesets) {
        geneset_code = sprintf("geneset.%s.%s",gs,organism);
        data(list=geneset_code,package = "chipenrich.data");

        geneset_list[[geneset_code]] = filter_genesets(get(geneset_code),max_geneset_size);
      }
  } else {
      # If the user provides a path to the geneset, build it and filter it
      geneset_list = list()
      geneset_code = 'user-supplied'

      geneset_list[[geneset_code]] = setup_geneset(genesets)
      geneset_list[[geneset_code]] = filter_genesets(geneset_list[[geneset_code]],max_geneset_size)
  }

  # Load TSS site info.
  tss_code = sprintf("tss.%s",genome);
  data(list=tss_code,package = "chipenrich.data");
  tss = get(tss_code);

  # Assign peaks to genes. NOTE: Depending on method,
  # peaks are assigned using assign_peaks(...) or
  # assign_peak_segments(...).
  if(!(method == 'broadenrich' || method == 'broadenrich_splineless')) {
	message("Assigning peaks to genes with assign_peaks(...) ..");
	assigned_peaks = assign_peaks(peakobj,ldef,tss,midpoint=T);
	message("Successfully assigned peaks..")
  } else {
  	message("Assigning peaks to genes with assigned_peak_segments(...) ..");
  	assigned_peaks = assign_peak_segments(peakobj,ldef)
  	message("Successfully assigned peaks..")
  }

  peak_genes = unique(assigned_peaks$geneid);

	# Add gene symbols to peak genes.
	genes_code = sprintf("genes.%s",organism);
	data(list=genes_code,package = "chipenrich.data");
	gene2symbol = get(genes_code);
	gene2symbol = change_names(gene2symbol,list(GENEID="geneid",SYMBOL="gene_symbol"));
	assigned_peaks = merge(assigned_peaks,gene2symbol,by="geneid",all.x=T);

  if(!(method == 'broadenrich' || method == 'broadenrich_splineless')) {

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
    "nearest_tss_gene_strand"
	);
  } else {
  	column_order = c(
  		"peak_id",
  		"chrom",
      "peak_start",
      "peak_end",
      "geneid",
      "gene_symbol",
    	"gene_locus_start",
    	"gene_locus_end",
    	"overlap_start",
    	"overlap_end",
    	"peak_overlap"
  	);
  }
  assigned_peaks = assigned_peaks[,column_order];

  ppg = num_peaks_per_gene(assigned_peaks,ldef,mappa);
  # This seems redundant given num_peaks_per_gene(...)
  ppg$peak = recode_peaks(ppg$num_peaks,num_peak_threshold);

  # Add relevant columns to ppg depending on the method
  if(method == 'broadenrich' || method == 'broadenrich_splineless') {
	message("Calculating peak overlaps with gene loci..")
	ppg = calc_peak_gene_overlap(assigned_peaks,ppg);
  }
  if(method == 'chipapprox') {
  	message("Calculating weights for approximate method..")
  	ppg = calc_approx_weights(ppg,mappa);
  }

  # Catch randomizations if present
  if(rndall) {
  	message('Randomizing across all genes.')
  	ppg = randomize_ppg_all(ppg)
  } else if (rndlength) {
  	message('Randomizing within length bins.')
  	ppg = randomize_ppg_length(ppg)
  }

  # Run chipenrich method on each geneset.
  results = list();
  for (gobj in geneset_list) {
		message(sprintf("Test: %s",method_name));
		message(sprintf("Genesets: %s",gobj@type));
		message("Running tests..");
    if (testf == "test_gam") {
      rtemp = test_func(gobj,ppg,n_cores);
    }
    if (testf == "test_fisher_exact") {
      rtemp = test_func(gobj,ppg,alternative=fisher_alt);
    }
    if (testf == "test_binomial") {
      rtemp = test_func(gobj,ppg);
    }
    if (testf == "test_gam_ratio") {
      rtemp = test_func(gobj,ppg,n_cores);
    }
    if (testf == "test_approx") {
      rtemp = test_func(gobj,ppg,nwp=F,n_cores);
    }
    if (testf == "test_gam_ratio_splineless") {
      rtemp = test_func(gobj,ppg,n_cores);
    }

    # Annotate with geneset descriptions.
    rtemp$"Description" = as.character(mget(rtemp$Geneset.ID,gobj@set.name,ifnotfound=NA));
    rtemp$"Geneset.Type" = gobj@type;

    results[[gobj@type]] = rtemp;
  }
  enrich = Reduce(rbind,results);

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
    "Geneset.Peak.Genes"
  );
  column_order = intersect(column_order,names(enrich));
  enrich = enrich[,column_order];

  # Order results by p-value.
  enrich = enrich[order(enrich$P.value),];

  # If there is a status column, re-sort so enriched terms are on top.
  if ("Status" %in% names(enrich)) {
    enrich = enrich[order(enrich$Status,decreasing=T),];
  }

  # Pull out tests that failed.
  bad_enrich = subset(enrich,is.na(P.value));
  enrich = subset(enrich,!is.na(P.value));

  # Write results to file.
  if (!is.null(out_name)) {
    filename_analysis = file.path(out_path,sprintf("%s_results.tab",out_name));
    write.table(enrich,file=filename_analysis,row.names=F,quote=F,sep="\t");
    message("Wrote results to: ",filename_analysis);

    filename_peaks = file.path(out_path,sprintf("%s_peaks.tab",out_name));
		write.table(assigned_peaks,file=filename_peaks,row.names=F,quote=F,sep="\t");
		message("Wrote peak-to-gene assignments to: ",filename_peaks);

    filename_opts = file.path(out_path,sprintf("%s_opts.tab",out_name));
		write.table(opts,file=filename_opts,row.names=F,quote=F,sep="\t");
		message("Wrote run options/arguments to: ",filename_opts);

    filename_ppg = file.path(out_path,sprintf("%s_peaks-per-gene.tab",out_name));
		write.table(ppg,file=filename_ppg,row.names=F,quote=F,sep="\t");
		message("Wrote count of peaks per gene to: ",filename_ppg);

    # If the user requested QC plots, generate those as well.
    # -- Spline fit plot
    # -- Histogram of distance from peaks to TSSs
    # -- Histogram of p-values from test
    # -- Expected # peaks vs. observed # peaks
    if (qc_plots) {
      filename_qcplots = file.path(out_path,sprintf("%s_qcplots.pdf",out_name));
      pdf(filename_qcplots);
			if (!(method=='broadenrich' || method=='broadenrich_splineless')) {
				print(..plot_spline_length(ldef,peak_genes,num_peaks,mappa=mappa));
				print(..plot_dist_to_tss(peakobj,tss));
			} else {
				print(..plot_gene_coverage(ppg));
			}
#     print(..plot_expected_peaks(ppg));
      dev.off();

      message("Wrote QC plots to: ",filename_qcplots);
    }
  }

  return(list(
    peaks = assigned_peaks,
    results = enrich,
    opts = opts,
    peaks_per_gene = ppg
  ));
}
