SUPPORTED_METHODS = list(
  'chipenrich' = "test_gam",
  'fet' = "test_fisher_exact",
  'broadenrich' = "test_gam_ratio"
)

HIDDEN_METHODS = list(
  'binomial' = "test_binomial",
  'chipapprox' = 'test_approx',
  'broadenrich_splineless' = "test_gam_ratio_splineless"
)

METHOD_NAMES = list(
	'fet' = "Fisher's Exact Test",
	'chipenrich' = "ChIP-Enrich",
	'broadenrich' = "Broad-Enrich",
	'chipapprox' = "ChIP-Enrich Approximate",
	'broadenrich_splineless' = "Broad-Enrich Splineless"
)
