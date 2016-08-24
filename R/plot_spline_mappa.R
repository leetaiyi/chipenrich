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
