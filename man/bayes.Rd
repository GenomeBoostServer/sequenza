\name{baf.bayes}
\alias{baf.bayes}
\alias{mufreq.bayes}
\title{Model Allele-Specific Copy Numbers with Given Cellularity and Ploidy}

\description{
  For a given pair of cellularity and ploidy parameters, returns the most likely 
  allele-specific copy numbers with corresponding log-posterior probability of fit,
  based on B-allele frequencies and depth ratios.
}

\usage{
  baf.bayes(Bf, depth.ratio, cellularity, ploidy, avg.depth.ratio,
            sd.Bf = 0.1, sd.ratio = 0.5, weight.Bf = 1, weight.ratio = 1,
            CNt.min = 0, CNt.max = 7, CNn = 2,
            priors.table = data.frame(CN = CNt.min:CNt.max, value = 1),
            ratio.priority = FALSE)
            
  mufreq.bayes(mufreq, depth.ratio, cellularity, ploidy, avg.depth.ratio,
            weight.mufreq = 100, weight.ratio = 100, CNt.min = 1, CNt.max = 7,
            CNn = 2, priors.table = data.frame(CN = CNt.min:CNt.max, value = 1))
}

\arguments{
  \item{Bf}{Vector of B-allele frequencies (values from 0 to 0.5).}
  \item{mufreq}{Vector of mutation frequencies (values from 0 to 1).}
  \item{depth.ratio}{Vector of depth ratios.}
  \item{sd.ratio}{Standard deviation observed in depth ratio measures per segment.}
  \item{sd.Bf}{Standard deviation observed in B-allele frequency measures per segment.}
  \item{weight.Bf}{Vector of weights for B-allele frequency values.}
  \item{weight.mufreq}{Vector of weights for mutation frequency values.}
  \item{weight.ratio}{Vector of weights for depth ratio values.}
  \item{cellularity}{Fraction of tumor cells in the sample.}
  \item{ploidy}{2 * ratio between tumor and normal cell DNA content.}
  \item{avg.depth.ratio}{Average normalized depth ratio.}
  \item{CNt.min}{Minimum copy number to consider.}
  \item{CNt.max}{Maximum copy number to consider.}
  \item{CNn}{Copy number in normal genome (typically 2).}
  \item{priors.table}{Data frame with columns \code{CN} and \code{value}, specifying
    copy numbers and their prior weights. Default weight is 1.}
  \item{ratio.priority}{Logical; if TRUE, use only depth ratio for copy number state,
    and B-allele frequency for determining number of B-alleles.}
}

\value{
  \item{CNt}{Copy number in tumor cells at tested position.}
  \item{A}{Number of A alleles at tested position.}
  \item{B}{Number of B alleles at tested position.}
  \item{CNn}{Copy number in normal cells (equals input CNn).}
  \item{Mt}{Number of mutated alleles at tested position.}
  \item{LPP}{Log-posterior probability of model fit.}
}

\details{
  Both functions use a naive Bayesian approach to calculate posterior probabilities
  of the data fitting the model given the cellularity and DNA-content parameters.
}

\seealso{
   \code{\link{baf.model.fit}}, \code{\link{mufreq.model.fit}}
}

\examples{
  \dontrun{
    # Load and process example data
    data.file <- system.file("extdata", "example.seqz.txt.gz",
                            package = "sequenza")
    seqz.data <- read.seqz(data.file, chr_name = "1")
    gc.stats <- gc.sample.stats(data.file)
    gc.vect <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
    
    # Normalize coverage
    seqz.data$adjusted.ratio <- seqz.data$depth.ratio /
                               gc.vect[as.character(seqz.data$GC.percent)]
    
    # Process segments
    breaks <- find.breaks(seqz.data[seqz.data$zygosity.normal == "het",],
                         gamma = 80, kmin = 10)
    seg.data <- segment.breaks(seqz.data, breaks = breaks)
    
    # Detect copy numbers
    cn.alleles <- baf.bayes(
        Bf = seg.data$Bf,
        depth.ratio = seg.data$depth.ratio,
        cellularity = 0.8,
        ploidy = 2,
        avg.depth.ratio = 1
    )
    
    # Process mutations
    mutations <- mutation.table(seqz.data,
        mufreq.threshold = 0.15,
        min.reads = 40,
        segments = seg.data
    )
    
    # Detect mutated alleles
    mut.alleles <- mufreq.bayes(
        mufreq = mutations$F,
        depth.ratio = mutations$adjusted.ratio,
        cellularity = 0.8,
        ploidy = 2,
        avg.depth.ratio = 1
    )
  }
}