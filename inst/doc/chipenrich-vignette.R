## ------------------------------------------------------------------------
library(chipenrich)

## ---- eval=FALSE---------------------------------------------------------
#  chipenrich(peaks, out_name = "chipenrich", out_path = getwd(),
#    genome = "hg19", genesets = c("GOBP", "GOCC", "GOMF"),
#    locusdef = "nearest_tss", method = "chipenrich",
#    fisher_alt = "two.sided", use_mappability = F, mappa_file = NULL,
#    read_length = 36, qc_plots = T, max_geneset_size = 2000,
#    num_peak_threshold = 1, n_cores = 1)

## ------------------------------------------------------------------------
data(peaks_E2F4, package = 'chipenrich.data')
data(peaks_H3K4me3_GM12878, package = 'chipenrich.data')

head(peaks_E2F4)

## ---- echo=FALSE---------------------------------------------------------
peaks_E2F4 = subset(peaks_E2F4, peaks_E2F4$chrom == 'chr1')
peaks_H3K4me3_GM12878 = subset(peaks_H3K4me3_GM12878, peaks_H3K4me3_GM12878$chrom == 'chr1')

## ------------------------------------------------------------------------
supported_genomes()

## ------------------------------------------------------------------------
supported_locusdefs()

## ------------------------------------------------------------------------
supported_genesets()

## ------------------------------------------------------------------------
supported_methods()

## ------------------------------------------------------------------------
supported_read_lengths()

## ---- fig.align='center', fig.cap='E2F4 peak distances to TSS', fig.height=6, fig.width=6, fig.show='hold'----
plot_dist_to_tss(peaks = peaks_E2F4, genome = 'hg19')

## ---- fig.align='center', fig.cap='E2F4 spline without mappability', fig.height=6, fig.width=6, fig.show='hold'----
plot_spline_length(peaks = peaks_E2F4, locusdef = 'nearest_tss', genome = 'hg19')

## ---- fig.align='center', fig.cap='E2F4 spline with mappability', fig.height=6, fig.width=6, fig.show='hold'----
plot_spline_length(peaks = peaks_E2F4, locusdef = 'nearest_tss',  genome = 'hg19', use_mappability = T, read_length = 24)

## ---- fig.align='center', fig.cap='H3K4me3 gene coverage', fig.height=6, fig.width=6, fig.show='hold'----
plot_gene_coverage(peaks = peaks_H3K4me3_GM12878, locusdef = 'nearest_tss',  genome = 'hg19')

## ------------------------------------------------------------------------
# Without mappability
gs_path = system.file('extdata','vignette_genesets.txt', package='chipenrich')
results = chipenrich(peaks = peaks_E2F4, genesets = gs_path,
	locusdef = "nearest_tss", qc_plots = F, out_name = NULL, n_cores = 1)
results.ce = results$results
print(results.ce[1:5,1:5])

## ------------------------------------------------------------------------
# With mappability
results = chipenrich(peaks = peaks_E2F4, genesets = gs_path,
	locusdef = "nearest_tss", use_mappability=T, read_length=24, qc_plots = F,
	out_name = NULL,n_cores=1)
results.cem = results$results
print(results.cem[1:5,1:5])

## ------------------------------------------------------------------------
results = chipenrich(peaks = peaks_H3K4me3_GM12878, genesets = gs_path,
	method='broadenrich', locusdef = "nearest_tss", qc_plots = F,
	out_name = NULL, n_cores=1)
results.be = results$results
print(results.be[1:5,1:5])

## ------------------------------------------------------------------------
results = chipenrich(peaks = peaks_E2F4, genesets = gs_path, locusdef = "5kb",
	method = "fet", fisher_alt = "two.sided", qc_plots = F, out_name = NULL)
results.fet = results$results
print(results.fet[1:5,1:5])

## ------------------------------------------------------------------------
head(results$peaks)

## ------------------------------------------------------------------------
head(results$peaks_per_gene)

## ------------------------------------------------------------------------
head(results$results)

