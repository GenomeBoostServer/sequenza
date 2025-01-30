\name{baf.model.fit}
\alias{baf.model.fit}
\alias{mufreq.model.fit}
\title{Model Fitting Using Maximum a Posteriori Inference}

\description{
  Computes the log-posterior probability distribution for the specified range
  of cellularity and ploidy parameters.
}

\usage{
    baf.model.fit(cellularity = seq(0.3, 1, by = 0.01),
        ploidy = seq(1, 7, by = 0.1), mc.cores = getOption("mc.cores", 2L),
        ...)
    mufreq.model.fit(cellularity = seq(0.3, 1, by = 0.01),
        ploidy = seq(1, 7, by = 0.1), mc.cores = getOption("mc.cores", 2L),
        ...)
}

\arguments{
    \item{cellularity}{Vector of cellularity values to test.}
    \item{ploidy}{Vector of ploidy values to test.}
    \item{mc.cores}{Number of cores to use for parallel processing.}
    \item{...}{Additional arguments passed to \code{\link{mufreq.bayes}} or
        \code{\link{baf.bayes}}.}
}

\value{
    A list containing three items:
    \item{ploidy}{Tested values of the ploidy parameter.}
    \item{cellularity}{Tested values of the cellularity parameter.}
    \item{lpp}{Log-posterior probability matrix for each cellularity/ploidy pair.}
}

\details{
    \code{baf.model.fit} uses \code{\link{baf.bayes}} to infer
    the log-posterior probability of model fit using all possible
    combinations of cellularity and ploidy values provided.
    
    \code{mufreq.model.fit} fits the mutation/depth ratio model using
    \code{\link{mufreq.bayes}}. Note that mutation data typically has more
    background noise compared to segmented B-allele frequency, potentially
    leading to less accurate results.
    
    \code{baf.model.fit} is the default method used to infer cellularity and
    ploidy on segmented chromosomes.
}

\seealso{
  \code{\link{cp.plot}} for visualization of results,
  \code{\link{get.ci}} for extracting confidence intervals.
}

\examples{
    \dontrun{
    # Load example data
    data.file <- system.file("extdata", "example.seqz.txt.gz",
                            package = "sequenza")
    seqz.data <- read.seqz(data.file, chr_name = "1")
    
    # Process and segment data
    gc.stats <- gc.sample.stats(data.file)
    gc.vect <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
    seqz.data$adjusted.ratio <- seqz.data$depth.ratio / 
                               gc.vect[as.character(seqz.data$GC.percent)]
    
    # Detect and process segments
    breaks <- find.breaks(seqz.data[seqz.data$zygosity.normal == "het",],
                         gamma = 80, kmin = 10)
    seg.data <- segment.breaks(seqz.data, breaks = breaks)
    
    # Filter segments and fit model
    seg.filt <- seg.data[seg.data$end.pos - seg.data$start.pos > 3e6,]
    weights <- (seg.filt$end.pos - seg.filt$start.pos) / 1e6
    
    CP <- baf.model.fit(
        Bf = seg.filt$Bf,
        depth.ratio = seg.filt$depth.ratio, 
        weight.ratio = weights,
        weight.Bf = weights,
        sd.ratio = seg.filt$sd.ratio,
        sd.Bf = seg.filt$sd.BAF,
        avg.depth.ratio = 1,
        cellularity = seq(0.1, 1, 0.01),
        ploidy = seq(0.5, 3, 0.05)
    )

    # Get confidence intervals
    confint <- get.ci(CP)
    }
}
