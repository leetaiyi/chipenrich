#' Display supported locus definitions
#'
#' @return A vector indicating the supported locus definitions.
#'
#' @examples
#'
#' supported_locusdefs()
#'
#' @export
supported_locusdefs = function() {
	piqr = data(package = "chipenrich.data")
	data_files = piqr$results[,3]

	ldefs = grep("locusdef", data_files, value = T)
	combos = Reduce(rbind,sapply(ldefs, strsplit, '[.]'))[,c(2,3)]
	df = data.frame(
		'genome' = combos[,1],
		'locusdef' = combos[,2],
		stringsAsFactors = F)
	return(df)
}

#' Display supported read lengths for mappability
#'
#' @return A vector indicating supported mappability read lengths.
#'
#' @examples
#'
#' supported_read_lengths()
#'
#' @export
supported_read_lengths = function() {
	piqr = data(package = "chipenrich.data")
	data_files = piqr$results[,3]

	mappas = grep("mappa", data_files, value = T)
	combos = Reduce(rbind,sapply(mappas, strsplit, '[.]'))[,c(2,3,4)]
	df = data.frame(
		'genome' = combos[,1],
		'locusdef' = combos[,2],
		'read_length' = combos[,3],
		stringsAsFactors = F)
	return(df)
}

#' Display supported genesets for gene set enrichment.
#'
#' @return A vector indicating supported gene sets for enrichment.
#'
#' @examples
#'
#' supported_genesets()
#'
#' @export
supported_genesets = function() {
	piqr = data(package = "chipenrich.data")
	data_files = piqr$results[,3]

	geneset_files = grep("geneset", data_files, value = T)
	combos = Reduce(rbind,sapply(geneset_files, strsplit, '[.]'))[,c(2,3)]
	df = data.frame(
		'geneset' = combos[,1],
		'organism' = combos[,2],
		stringsAsFactors = F)
	df = df[order(df$organism), ]
	return(df)
}

#' Display supported genomes.
#'
#' @return A vector indicating supported genomes.
#'
#' @examples
#'
#' supported_genomes()
#'
#' @export
supported_genomes = function() {
	piqr = data(package = "chipenrich.data")
	data_files = piqr$results[,3]

	as.character(stats::na.omit(unique(stringr::str_match(data_files,"(locusdef|mappa)\\.(\\w+)\\.")[,3])))
}

#' Display supported gene set enrichment methods.
#'
#' @return A vector indicating supported methods for gene set enrichment.
#'
#' @examples
#'
#' supported_methods()
#'
#' @export
supported_methods = function() {
	return(names(SUPPORTED_METHODS))
}
