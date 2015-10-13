\name{plot_dist_to_tss}
\alias{plot_dist_to_tss}

\title{
Plot histogram of distance from peak to nearest TSS
}

\description{
Create a histogram of the distance from each peak to the nearest transcription start site (TSS) of any gene. 
}

\usage{
plot_dist_to_tss(peaks, genome = "hg19")
}

\arguments{
  \item{peaks}{
Either a data frame, a BED file, a .broadPeak file, or a .narrowPeak file. The data frame should have at least 3 columns: chrom, start, and end. Chrom should follow UCSC convention, e.g. "chrX". The file should be tab-delimited.
}
  \item{genome}{
A string indicating the genome upon which the peaks file is based. Supported genomes are listed by the \code{\link{supported_genomes}} function. 
}
}

\value{
A trellis plot object.
}

\author{
Ryan Welch \email{welchr@umich.edu}
}

\seealso{
\code{\link{chipenrich}}
}

\examples{
library(chipenrich.data)
library(chipenrich)

# Create histogram of distance from peaks to nearest TSS.
data(peaks_E2F4)
plot_dist_to_tss(peaks_E2F4)
}
