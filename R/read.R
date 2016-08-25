#' Read BEDX+Y files and convert into GRanges
#'
#' Given a \code{file_path}, read in a delimited file, assuming it is BEDX+Y,
#' keep only the first three columns: chrom, start, end, and output \code{GRanges}.
#'
#' Typically, this function will not be used alone, but inside \code{chipenrich()}.
#'
#' @param file_path A path to a valid BEDX+Y file.
#'
#' @return A \code{GRanges} that is unstranded, and contains only chrom, start, and end.
#'
#' @examples
#'
#' # Example of BED3 with no header
#' file = system.file('extdata', 'test_assign.bed', package = 'chipenrich')
#' peaks = read_bed(file)
#'
#' # Example of BED3 with header
#' file = system.file('extdata', 'test_header.bed', package = 'chipenrich')
#' peaks = read_bed(file)
#'
#' # Example of narrowPeak with header
#' file = system.file('extdata', 'test.narrowPeak', package = 'chipenrich')
#' peaks = read_bed(file)
#'
#' @export
read_bed = function(file_path) {
	if (!file.exists(file_path)) {
		stop("Can't find BED file: ",file_path)
	}

	# Determine how many lines to skip
	chunk_size = 500;
	chunk = scan(file_path,what="character",nmax=chunk_size,strip.white=T,sep="\n",quiet=T)
	skip_n = suppressWarnings(min(grep("^chr(\\d+|\\w+)\\s+\\d+",chunk)) - 1)

	if (is.infinite(skip_n)) {
		stop("Error: no valid chromosomes detected within first 500 lines of BED file.")
	}

	message(sprintf("Skipping %i lines of BED header..",skip_n))

	# Read the file, subset to first 3 columns, and rename columns
	peaks = read.table(file_path,header=F,skip=skip_n)
	peaks = peaks[,1:3]
	names(peaks) = c("chrom","start","end")

	sub_check = peaks[1:min(nrow(peaks),100),]
	if (!all(grepl("chr",sub_check$chrom))) {
		stop("First column of BED file should have chr* entries.")
	}

	if (any(sub_check$start < 0) | any(sub_check$end < 0)) {
		stop("Start/end positions of peaks should be >= 0.")
	}

	# Convert to list of GenomicRanges
	gr = GenomicRanges::GRanges(
		seqnames = peaks$chrom,
		ranges = IRanges::IRanges(start = peaks$start, end = peaks$end))

	# Reduce peaks
	gr = GenomicRanges::reduce(gr)
	gr$name = paste('peak:', 1:length(gr), sep='')

	return(gr)
}

#' Read BEDGFF files and convert into GRanges
#'
#' Given a \code{file_path}, read in a delimited file, assuming it is BEDGFF (as
#' is output by modENCODE for D. Melanogaster TF ChIP-seq experiments), keep
#' only chrom, start, and end columns, and output \code{GRanges}.
#'
#' Typically, this function will not be used alone, but inside \code{chipenrich()}.
#'
#' @param file_path A path to a valid BEDGFF file (as from modENCODE).
#'
#' @return A \code{GRanges} that is unstranded, and contains only chrom, start, and end.
#'
#' @examples
#'
#' # Example of GFF3
#' file = system.file('extdata', 'test.gff3', package = 'chipenrich')
#' peaks = read_bedgff(file)
#'
#' # Example of gzipped GFF3
#' file = system.file('extdata', 'test.gff3.gz', package = 'chipenrich')
#' peaks = read_bedgff(file)
#'
#' @export
read_bedgff = function(file_path) {
	if (!file.exists(file_path)) {
		stop("Can't find BED file: ",file_path)
	}

	chunk_size = 500;
	chunk = scan(file_path,what="character",nmax=chunk_size,strip.white=T,sep="\n",quiet=T)

	skip_n = suppressWarnings(min(grep("^\\d(L|R)",chunk)) - 1)

	if (is.infinite(skip_n)) {
		stop("Error: no valid chromosomes detected within first 500 lines of BED file.")
	}

	message(sprintf("Skipping %i lines of BED header..",skip_n))

	peaks = read.table(file_path,header=F,skip=skip_n)
	peaks = peaks[,c(1,4,5)]
	names(peaks) = c("chrom","start","end")

	sub_check = peaks[1:min(nrow(peaks),100),]
	if (!all(grepl("chr",sub_check$chrom))) {
		peaks$chrom = paste('chr',peaks$chrom,sep='')
		message("Adding 'chr' to chromosome column for compatibility with locus definitions.")
	}

	if (any(sub_check$start < 0) | any(sub_check$end < 0)) {
		stop("Start/end positions of peaks should be >= 0.")
	}

	# Convert to list of GenomicRanges
	gr = GenomicRanges::GRanges(
		seqnames = peaks$chrom,
		ranges = IRanges::IRanges(start = peaks$start, end = peaks$end))

	# Reduce peaks
	gr = GenomicRanges::reduce(gr)
	gr$name = paste('peak:', 1:length(gr), sep='')

	return(gr)
}

#' Convert BEDX+Y data.frames and into GRanges
#'
#' Given a \code{data.frame} in BEDX+Y format, keep only the first three
#' columns: chrom, start, end, and output \code{GRanges}.
#'
#' Typically, this function will not be used alone, but inside \code{chipenrich()}.
#'
#' @param dframe A BEDX+Y style \code{data.frame}.
#'
#' @return A \code{GRanges} that is unstranded, and contains only chrom, start, and end.
#'
#' @examples
#'
#' # Example of BED3 with no header
#' data(peaks_H3K4me3_GM12878, package='chipenrich.data')
#' peaks = load_peaks(peaks_H3K4me3_GM12878)
#'
#' @export
load_peaks = function(dframe) {
	# Check columns.
	for (col in c("chrom","start","end")) {
		if (!col %in% names(dframe)) {
			stop(sprintf("error in peaks data frame: no column named '%s'",col))
		}
	}

	# Convert to list of GenomicRanges
	gr = GenomicRanges::GRanges(
		seqnames = dframe$chrom,
		ranges = IRanges::IRanges(start = dframe$start, end = dframe$end))

	# Reduce peaks
	gr = GenomicRanges::reduce(gr)
	gr$name = paste('peak:', 1:length(gr), sep='')

	return(gr)
}
