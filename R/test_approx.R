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
