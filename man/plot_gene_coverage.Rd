\name{plot_gene_coverage}
\alias{plot_gene_coverage}

\title{
Plot gene locus length vs. gene coverage by a peak
}

\description{
Create a plot showing the relationship between the length of a gene locus and the proportion of the locus covered by a peak.
}

\usage{
plot_gene_coverage(peaks,locusdef="nearest_tss",
genome='hg19',use_mappability=F,read_length=36,legend=T,xlim=NULL)
}

\arguments{
  \item{peaks}{
Either a data frame, a BED file, a .broadPeak file, or a .narrowPeak file. The data frame should have at least 3 columns: chrom, start, and end. Chrom should follow UCSC convention, e.g. "chrX". The file should be tab-delimited.
}
  \item{locusdef}{
A string denoting the locus definition to be used. A locus definition controls how peaks are assigned to genes. See \code{\link{supported_locusdefs}} for a list of supported definitions. }
  \item{genome}{
A string indicating the genome upon which the peaks file is based. Supported genomes are listed by the \code{\link{supported_genomes}} function. 
}
  \item{use_mappability}{
If true, each gene's locus length is corrected for by mappability. 
}
  \item{read_length}{
If using mappability (see above), the read length should match the length of reads used in the original experiment.  
}
  \item{legend}{
If true, a legend will be drawn on the plot. 
}
  \item{xlim}{
Set the x-axis limit. NULL means select x-lim automatically. 
}
}

\details{
The x-axis is the log10 of the gene locus length (defined by the locus definition). The y-axis is the proportion of the gene locus covered by a peak. Each plotted black dot is the mean (over both values) of 25 genes (ordered by log10 of the locus length).
}

\value{
A trellis plot object.
}

\author{
Raymond Cavalcante \email{rcavalca@umich.edu}
}

\seealso{
\code{\link{chipenrich}}
}

\examples{
library(chipenrich.data)
library(chipenrich)

# Create histogram of distance from peaks to nearest TSS.
data(peaks_E2F4)
plot_gene_coverage(peaks_E2F4)
}
