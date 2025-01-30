\name{read.seqz}
\alias{read.seqz}
\title{Read and Process Sequenza Files}

\description{
    Efficiently reads and processes sequenza (seqz) files, providing support for both 
    chromosome-specific and coordinate-specific data extraction. The function implements 
    optimized reading strategies for both standard and tabix-indexed files.
}

\usage{
    read.seqz(file, n_lines = NULL, col_types = "ciciidddcddccc", chr_name = NULL,
        buffer = 33554432, parallel = 1,
        col_names = c("chromosome", "position", "base.ref", "depth.normal",
            "depth.tumor", "depth.ratio", "Af", "Bf", "zygosity.normal",
            "GC.percent", "good.reads", "AB.normal", "AB.tumor",
            "tumor.strand"), ...)
}

\arguments{
    \item{file}{Path to the input seqz file (can be gzipped)}
    \item{col_types}{Column type specification string for data parsing. Defaults to 
        "ciciidddcddccc" matching standard seqz file structure}
    \item{chr_name}{Optional chromosome selector. For standard files, extracts data 
        for the specified chromosome. For tabix-indexed files, supports genomic 
        coordinates (e.g., "chr5:1-1000000")}
    \item{n_lines}{vector of length 2 specifying the first and last line to
        read from the file. If specified, only the selected portion of the
        file will be used.}
    \item{buffer}{maximal size of each chunk in bytes(see
        \code{\link{chunk.apply}}).}
    \item{parallel}{integer, number of threads used to process a seqz file
        (see \code{\link{chunk.apply}}).}
    \item{col_names}{names of the columns of the seqz file. The default
        corresponds to the column names of a seqz file.}
    \item{...}{any arguments accepted by \code{read_tsv}.}
}

\details{
    This function provides efficient access to seqz files through:
    \itemize{
        \item Optimized chunk-based reading for memory efficiency
        \item Parallel processing support for large files
        \item Direct coordinate access for tabix-indexed files
        \item Automatic format detection and handling
    }
    
    The function automatically detects file type (standard or tabix-indexed) and 
    applies the most efficient reading strategy.
}

\format{
    A seqz file is a tab-separated text file with 14 columns and a header row.
    The first 3 columns are derived from the original \code{pileup}
    file and contain:
    \describe{
        \item{chromosome}{the chromosome name}
        \item{position}{the base position}
        \item{base.ref}{the base in the reference genome.
            Note that this is NOT necessarily the same base as in
            the normal specimen.}
    The remaining 10 columns contain the following information:
        \item{depth.normal}{read depth observed in the normal sample}
        \item{depth.tumor}{read depth observed in the tumor sample}
        \item{depth.ratio}{ratio of \code{depth.tumor} and \code{depth.normal}}
        \item{Af}{A-allele frequency observed in the tumor sample}
        \item{Bf}{B-allele frequency observed in the tumor sample in
            heterozygous positions}
        \item{zygosity.normal}{zygosity of the reference sample. "hom"
            corresponds to AA or BB, whereas "het" corresponds to AB or BA}
        \item{GC.percent}{GC-content (percent), calculated from the reference
            genome in fixed nucleotide windows }
        \item{good.reads}{number of reads that passed the quality threshold
            (threshold specified in the pre-processing software), in the
            tumor specimen}
        \item{AB.normal}{base(s) found in the germline sample; for heterozygous
            positions AB are sorted using the values of Af and Bf respectively}
        \item{AB.tumor}{base(s) found in the tumor sample not present in the
            normal specimen. The field include all the variants found in the
            tumor alignment, separated by a colon. Each variant contains the
            base and the observed frequency}
        \item{tumor.strand}{frequency of the variant nucleotides detected on
            the forward orientation. The field have a consistent structure with
            \code{AB.tumor}, indicating the fraction, relative to the total
            number of reads presenting the specific variant, orientated
            in the forward direction}
    }
}

\seealso{
  \code{read_delim}.
}

\examples{
   \dontrun{
    data_file <- system.file("extdata", "example.seqz.txt.gz", package = "sequenza")

    # Read chromosome 1 from a seqz file
    seqz_data <- read.seqz(data_file, chr_name = "1")

    # Read a specific region from a tabix-indexed file
    seqz_data <- read.seqz(data_file, chr_name = "1:1000000-2000000")

    # Process large files with parallel processing
    seqz_data <- read.seqz(data_file, parallel = 4)

    # Compare different reading methods
    system.time(seqz_data <- read.seqz(data_file, chr_name = "X"))
    system.time(seqz_data <- read.seqz(data_file, 
               n_lines = c(gc.stats$file.metrics$start[1], 
                         gc.stats$file.metrics$end[1])))
   }
}
