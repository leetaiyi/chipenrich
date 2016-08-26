context('Test User Supplied Genesets')

# Get some peaks to do enrichment on
data(peaks_H3K4me3_GM12878,package='chipenrich.data')

# Test building the geneset itself
geneset_file = system.file('extdata', 'test_geneset.txt', package='chipenrich')
test_geneset = setup_geneset(geneset_file)

test_that('Object of GeneSet class is returned',{
    expect_equal(class(test_geneset)[1], 'GeneSet')
})

test_that('Three genesets are created',{
    expect_equal(length(test_geneset@set.gene), 3)
})

test_that('Test names of genesets',{
    expect_equal(all(ls(test_geneset@set.gene) == c('GO:0035909','GO:0045822','GO:0045823')), TRUE)
})

# # Test running chipenrich
# enrich_test = chipenrich(
#     peaks = peaks_H3K4me3_GM12878,
#     out_name = 'user_geneset_test',
#     genome = 'hg19',
#     genesets = geneset_file,
#     locusdef = 'nearest_tss',
#     method = 'chipenrich',
#     qc_plots = F,
#     n_cores = 2)
#
# test_that('Test that there are three results',{
#     expect_equal(nrow(enrich_test$results),3)
# })
#
# test_that('Genesets in enrichment results match those given by user',{
#     expect_equal(all(enrich_test$result$Geneset.ID == c('GO:0035909', 'GO:0045823', 'GO:0045822')), TRUE)
# })
