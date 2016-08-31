context('Test peaks per gene functions')

################################################################################
# Test human locus definitions

# Load the required data for the test
hg19_ldefs = grep('locusdef.hg19', data(package = 'chipenrich.data')$results[,3], value=T)
data(list = hg19_ldefs, package = 'chipenrich.data')
tss = get(data('tss.hg19', package = 'chipenrich.data'))

file = system.file('extdata', 'test_assign.bed', package = 'chipenrich')
peaks = read_bed(file)

# Expected values for peak assignments
expected_numpeaks = list(
	'locusdef.hg19.10kb' = c(3),
	'locusdef.hg19.10kb_and_more_upstream' = c(1),
	'locusdef.hg19.1kb' = c(2),
	'locusdef.hg19.5kb' = c(3),
	'locusdef.hg19.exon' = c(1,1),
	'locusdef.hg19.intron' = c(1,1),
	'locusdef.hg19.nearest_gene' = c(3,1),
	'locusdef.hg19.nearest_tss' = c(3,1)
)

# Do the tests in a loop for each locus definition of hg19
for(ldef in hg19_ldefs) {
	test_assign = assign_peaks(peaks = peaks, locusdef = get(ldef), tss = tss)

	test_ppg = num_peaks_per_gene(assigned_peaks = test_assign, locusdef = get(ldef), mappa = NULL)

	expect_true(all(subset(test_ppg, num_peaks != 0)$num_peaks == expected_numpeaks[[ldef]]),
		info = sprintf('Test num_peaks_per_gene(): %s', ldef))
}

expected_overlaps = list(
	'locusdef.hg19.10kb' = c(700,81),
	'locusdef.hg19.10kb_and_more_upstream' = c(918),
	'locusdef.hg19.1kb' = c(327),
	'locusdef.hg19.5kb' = c(700),
	'locusdef.hg19.exon' = c(100,100),
	'locusdef.hg19.intron' = c(100,100),
	'locusdef.hg19.nearest_gene' = c(700,1000),
	'locusdef.hg19.nearest_tss' = c(700,1000)
)

# Do the tests in a loop for each locus definition of hg19
for(ldef in hg19_ldefs) {
	test_assign = assign_peak_segments(peaks = peaks, locusdef = get(ldef))

	test_ppg = num_peaks_per_gene(assigned_peaks = test_assign, locusdef = get(ldef), mappa = NULL)
	test_ppg = calc_peak_gene_overlap(assigned_peaks = test_assign, ppg = test_ppg)

	expect_true(all(subset(test_ppg, peak_overlap != 0)$peak_overlap == expected_overlaps[[ldef]]),
		info = sprintf('Test num_peaks_per_gene(): %s', ldef))
}
