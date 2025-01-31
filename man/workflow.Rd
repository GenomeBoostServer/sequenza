\name{sequenza}
\alias{sequenza.extract}
\alias{sequenza.fit}
\alias{sequenza.results}
\title{Sequenza workflow functions for copy number and mutation analysis}

\description{
  Three main interface functions that implement the standard Sequenza analysis workflow:
  1. Data extraction and preprocessing
  2. Model fitting for cellularity and ploidy
  3. Results generation and visualization
}

\usage{
  sequenza.extract(file, window = 1e6, overlap = 1,
    slide_win = 100, peak_wins = seq(from = 50, to = 300, by = 25),
    normalization.method = "mean", ignore.normal = FALSE,
    verbose = TRUE, chromosome.list = NULL, breaks = NULL,
    min.mut.freq = 0.1, min.reads = 40, min.reads.normal = 10,
    min.reads.baf = 1, max.mut.types = 1, min.type.freq = 0.9,
    min.fw.freq = 0, assembly = "hg38", gc.stats = NULL,
    do_raster = FALSE, smooth_gc = FALSE, min_times_gc = 5,
    gc_grid = 250, parallel = 1, weighted.mean = TRUE, ...)

  sequenza.fit(sequenza.extract, female = TRUE, N.ratio.filter = 10,
    N.BAF.filter = 1, segment.filter = 3e6, mufreq.threshold = 0.1,
    XY = c(X = "X", Y = "Y"), cellularity = seq(0.1,1,0.01),
    ploidy = seq(1, 7, 0.1), ratio.priority = FALSE, method = "baf",
    priors.table = data.frame(CN = 2, value = 2), chromosome.list = 1:24,
    mc.cores = getOption("mc.cores", 2L), verbose = TRUE)

  sequenza.results(sequenza.extract, cp.table = NULL, sample.id,
    out.dir = getwd(), cellularity = NULL, ploidy = NULL,
    female = TRUE, CNt.max = 20, ratio.priority = FALSE,
    XY = c(X = "X", Y = "Y"), chromosome.list = 1:24)
}

\arguments{
  \item{file}{the name of the seqz file to read.}
  \item{window}{size of windows used when plotting mean and quartile ranges of depth ratios and B-allele frequencies. Smaller windows will take more time to compute.}
  \item{overlap}{integer specifying the number of overlapping windows.}
  \item{slide_win}{size of the sliding window for depth ratio calculation.}
  \item{peak_wins}{vector of window sizes for peak detection in depth ratio.}
  \item{normalization.method}{string defining the operation to perform during the GC-normalization process. Possible values are \code{mean} (default) and \code{median}. A \code{median} normalization is preferable with noisy data.}
  \item{ignore.normal}{boolean, when set to TRUE the process will ignore the normal coverage and perform the analysis by using the normalized tumor coverage.}
  \item{verbose}{logical, indicating whether to print progress information during processing.}
  \item{chromosome.list}{vector containing the index or the names of the chromosome to include in the model fitting.}
  \item{breaks}{Optional data.frame in the format chrom, start.pos, end.pos, defining a pre-existing segmentation. When the argument is set the built-in segmentation will be skipped in favor of the suggested breaks.}
  \item{min.mut.freq}{Minimum mutation frequency threshold.}
  \item{min.reads}{Minimum number of reads above the quality threshold to accept mutation calls.}
  \item{min.reads.normal}{Minimum number of reads required to determine genotype in the normal sample.}
  \item{min.reads.baf}{Threshold for depth of positions included in average BAF calculation per segment.}
  \item{max.mut.types}{Maximum number of different base substitutions per position (1-3).}
  \item{min.type.freq}{minimum frequency of aberrant types.}
  \item{min.fw.freq}{Minimum frequency of variant reads in the forward strand. When set to 0, excludes variants with strand frequency outside [0,1].}
  \item{assembly}{Genome assembly version (see \code{\link{aspcf}} or \code{\link{pcf}}).}
  \item{gc.stats}{GC statistics object from \code{\link{gc.sample.stats}}. If \code{NULL}, computed from input file.}
  \item{do_raster}{Logical; whether to generate raster plots.}
  \item{smooth_gc}{Logical; whether to smooth GC content calculations.}
  \item{min_times_gc}{minimum times for GC normalization (default: 5).}
  \item{gc_grid}{grid size for GC binning (default: 250).}
  \item{parallel}{Number of parallel processes to use (pbapply backend).}
  \item{weighted.mean}{Logical; whether to use read depth as weights for depth ratio and BAF calculations.}
  \item{sequenza.extract}{a list of objects as output from the \code{sequenza.extract} function.}
  \item{method}{method to use to fit the data; possible values are \code{baf} to use \code{\link{baf.model.fit}} or \code{mufreq} to use the \code{\link{mufreq.model.fit}} function to fit the data.}
  \item{cp.table}{a list of objects as output from the \code{sequenza.fit} function.}
  \item{female}{logical, indicating whether the sample is male or female, to properly handle the X and Y chromosomes. Implementation only works for the human normal karyotype.}
  \item{CNt.max}{maximum copy number to consider in the model.}
  \item{N.ratio.filter}{threshold of minimum number of observation of depth ratio in a segment.}
  \item{N.BAF.filter}{threshold of minimum number of observation of B-allele frequency in a segment.}
  \item{segment.filter}{threshold segment length (in base pairs) to filter out short segments, that can cause noise when fitting the cellularity and ploidy parameters. The threshold will not affect the allele-specific segmentation.}
  \item{XY}{character vector of length 2 specifying the labels used for the X and Y chromosomes.}
  \item{cellularity}{vector of candidate cellularity parameters.}
  \item{ploidy}{vector candidate ploidy parameters.}
  \item{priors.table}{data frame with the columns \code{CN} and \code{value}, containing the copy numbers and the corresponding weights. To every copy number is assigned the value 1 as default, so every values different then 1 will change the corresponding weight.}
  \item{ratio.priority}{logical, if TRUE only the depth ratio will be used to determine the copy number state, while the Bf value will be used to determine the number of B-alleles.}
  \item{sample.id}{identifier of the sample, to be used as a prefix for saved objects.}
  \item{out.dir}{output directory where the files and objects will be saved.}
  \item{mc.cores}{legacy argument to set the number of cores, but it refers to the \code{cl} of \code{\link{pblapply}}. It uses \code{\link{mclapply}} when set to an integer.}
  \item{segment_weights}{List of weights for segment scoring:
    \itemize{
      \item{fit}: Weight for raw fit quality of segments (default: 0.8). Higher values favor segments that better match the data.
      \item{penalty}: Weight for consecutive outliers penalty (default: 0.4). Higher values more strongly penalize long runs of outliers.
      \item{elbow}: Weight for proximity to the elbow point in the fit curve (default: 0.1). Higher values favor selecting window sizes near where marginal improvements in fit begin to diminish.
      \item{segments}: Weight for segment count before the elbow point (default: 0.1). Higher values favor solutions with more segments up to the elbow point, promoting finer segmentation where supported by the data.
      \item{window}: Window size penalty weight (default: 0.05). Higher values more strongly prefer smaller window sizes when fit scores are similar.
    }
    
    The weights control the tradeoff between different aspects of segmentation:
    \itemize{
      \item Higher \code{fit} and lower \code{penalty} weights will favor solutions that match the data closely
      \item Higher \code{elbow} weights help avoid oversegmentation by identifying diminishing returns
      \item Higher \code{segments} weights allow more segments when supported by clear evidence
      \item Higher \code{window} weights prefer more granular analysis when possible
    }
  }
}

\details{
  The workflow consists of three main steps:

  1. \code{sequenza.extract}: Processes raw sequencing data to:
     * Normalize depth ratios for GC-content bias
     * Perform allele-specific segmentation
     * Filter mutations
     * Calculate binned data for visualization
     Uses pbapply for parallel processing when parallel > 1.

  2. \code{sequenza.fit}: Estimates sample cellularity and ploidy by:
     * Analyzing depth ratios and B-allele frequencies
     * Supporting both BAF and mutation frequency methods
     * Providing progress tracking for long computations
     * Handling memory efficiently for large datasets

  3. \code{sequenza.results}: Generates comprehensive results including:
     * Segment and mutation tables
     * Genome-wide visualizations
     * Chromosome-level plots
     * Copy number statistics
     * Model fit diagnostics
     With progress tracking and efficient memory management.
}

\seealso{
  \code{\link{genome.view}}, \code{\link{baf.bayes}}, \code{\link{cp.plot}}, \code{\link{get.ci}}.
}

\examples{
  \dontrun{

data.file <-  system.file("extdata", "example.seqz.txt.gz",
              package = "sequenza")
test <- sequenza.extract(data.file)
test.CP   <- sequenza.fit(test)
sequenza.results(test, test.CP, out.dir = "example",
                 sample.id = "example")

   }
}