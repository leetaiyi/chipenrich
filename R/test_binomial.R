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
