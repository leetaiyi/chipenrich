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
