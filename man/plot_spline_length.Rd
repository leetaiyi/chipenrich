\name{plot_spline_length}
\alias{plot_spline_length}

\title{
Plot probability of peak being assigned to a gene vs. gene length
}

\description{
Create a plot showing the probability of a gene being assigned a peak given its locus length. The plot shows an empirical fit to the data using a binomial smoothing spline. 
}

\usage{
plot_spline_length(peaks,locusdef="nearest_tss",
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
The x-axis represents the log10 of the gene locus length (defined by the locus definition.) The y-axis is the probability of a gene being assigned a peak. Each black dot on the plot is a bin of 25 genes, the y-axis value is then the proportion of peaks that were assigned to genes in that bin, and the x-axis value is the average of the locus lengths of genes in that bin. 

The spline curve is fit using a binomial smoothing spline model, see \code{\link{chipenrich}} for more information. This curve is created by modeling presence of peak (a 0/1 binary variable denoting whether the gene was assigned a peak) given the gene locus length. 

The random genes curve represents the model where each gene has the same probability of being assigned a peak given the total number of peaks in the dataset. 

The "random peaks across genome" curve represents the model where the probability of a gene being assigned at least 1 peak is dependent on its length: \deqn{p(peak\textbar L_{locus}) = 1 - (1 - \frac{L_{locus}}{L_{genome}})^{n_{peaks}}} \eqn{L_{locus}} is the length of the gene locus, \eqn{L_{genome}} is the length of the sampleable genome, and \eqn{n_{peaks}} is the number of peaks in the dataset. 
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

# Spline plot for E2F4 example peak dataset. 
data(peaks_E2F4)
plot_spline_length(peaks_E2F4,genome='hg19')

# Create the plot for a different locus definition
# to compare the effect. 
plot_spline_length(peaks_E2F4,locusdef='nearest_gene',genome='hg19');

\donttest{
# Create the plot, but write the result to a PDF instead of displaying it interactively. 
pdf("spline_plot.pdf");
p = plot_spline_length(peaks_E2F4,genome='hg19');
dev.off();
}
}
