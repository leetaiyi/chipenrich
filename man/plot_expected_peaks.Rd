\name{plot_expected_peaks}
\alias{plot_expected_peaks}

\title{
Plot the expected number of peaks given locus length for each gene. 
}

\description{
Create a plot showing the expected count of peaks per gene given their locus lengths.

Also plotted are confidence intervals for the expected count, and the actual observed number of peaks per gene. 
}

\usage{
plot_expected_peaks(peaks, 
  locusdef = "nearest_tss", 
  genome = "hg19", 
  use_mappability = F, 
  read_length = 36, 
  mappa_file = NULL
)
}

\arguments{
  \item{peaks}{
    Either a data frame, a BED file, a .broadPeak file, or a .narrowPeak file. The data frame should have at least 3 columns: chrom, start, and end. Chrom should follow UCSC convention, e.g. "chrX". The file should be tab-delimited.
  }

  \item{locusdef}{
    A string denoting the locus definition to be used. A locus definition controls how peaks are assigned to genes. See \code{\link{supported_locusdefs}} for a list of supported definitions.
  }

  \item{genome}{
    A string indicating the genome upon which the peaks file is based. Supported genomes are listed by the \code{\link{supported_genomes}} function. 
  }

  \item{use_mappability}{
    If true, each gene's locus length is corrected for by mappability. 
  }

  \item{read_length}{
    If using mappability (see above), the read length should match the length of reads used in the original experiment.  
  }

  \item{mappa_file}{
    Path to a file containing user-specified gene locus mappability. The file should contain two columns: geneid and mappa. Gene IDs should be Entrez gene IDs. Mappability values should be between 0 and 1.
  }
}

\details{
The x-axis is gene locus length (for the defined \code{locusdef}.) The y-axis is count of peaks. Each blue dot represents the observed count of peaks assigned to a gene. 

The black line represents the expected number of peaks given locus length. 

Also drawn are the 5 and 95\% percentiles of a Poisson distribution for the expected number of peaks, and the 5 and 95\% percentiles adjusted for the number of genes (Bonferroni adjustment - e.g. 0.05 / # of genes.)
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

# Expected peak count plot for the E2F4 dataset. 
data(peaks_E2F4)
plot_expected_peaks(peaks_E2F4,genome='hg19')

# Create the plot for a different locus definition
# to compare the effect. 
plot_expected_peaks(peaks_E2F4,locusdef='nearest_gene',genome='hg19');

\donttest{
# Create the plot, but write the result to a PDF 
# instead of displaying it interactively. 
pdf("expected_peak_plot.pdf");
p = plot_expected_peaks(peaks_E2F4,genome='hg19');
dev.off();
}
}
