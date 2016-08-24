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

# Used in plot_spline_mappa(...)
# Which is not called
add_emp_peak_mappa = function(gpw,bin_size=25) {
  d = gpw;
  d = d[order(d$mappa),];
  d$group = ceiling((1:dim(d)[1])/bin_size);

  bygroup = stats::aggregate(peak ~ group,d,mean);
  d$emp_peak = bygroup$peak[d$group];
  d;
}

# Used in plot_spline_mappa(...)
# Which is not called
avg_binned_peak_mappa = function(gpw,bin_size=25) {
  d = gpw;
  d = d[order(d$mappa),];
  d$group = ceiling((1:dim(d)[1])/bin_size);

  bygroup = stats::aggregate(cbind(peak,mappa) ~ group,d,mean);
  bygroup$log_avg_mappa = log10(bygroup$mappa);
  names(bygroup) = c("group","peak","avg_mappa","log_avg_mappa");
  bygroup;
}

# Used in plot_spline_mappa()
# Which is not called
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
	if(method == 'chipenrich') {
		if (all(as.logical(sg_go))) {
			cont_length = quantile(gpw$length,0.0025);

			cont_gene = data.frame(
				geneid = "continuity_correction",
				length = cont_length,
				log10_length = log10(cont_length),
				num_peaks = 0,
				peak = 0,
				stringsAsFactors = F)

			if ("mappa" %in% names(gpw)) {
				cont_gene$mappa = 1;
			}
			gpw = rbind(gpw,cont_gene);
			b_genes = c(b_genes,1);

			message(sprintf("Applying correction for geneset %s with %i genes...",go_id,length(go_genes)));
		}
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

# Not used or exported
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

  l = unlist(as.list(environment()));
    opts = data.frame(
    args = names(l),
    values = sapply(unlist(l),paste,collapse=","),
    stringsAsFactors = F
  );
  rownames(opts) = 1:length(l);

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
	assigned_peaks = assign_peaks(peakobj,ldef,tss);
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
