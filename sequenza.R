#!/usr/bin/env Rscript

# Accept R_LIBS path as command line argument (optional)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && nchar(args[1]) > 0) {
    R_LIB_PATH <- args[1]
    .libPaths(c(R_LIB_PATH, .libPaths()))
    Sys.setenv(R_LIBS = paste(R_LIB_PATH, Sys.getenv("R_LIBS"), sep = ":"))
    cat("==> Using custom R library path:", R_LIB_PATH, "\n")
}

library(sequenza)
library(parallel)

# SEQZ_FILE is a required input; there is no meaningful default.
seqz_file <- Sys.getenv("SEQZ_FILE")
if (!nzchar(seqz_file)) {
    stop("SEQZ_FILE environment variable is not set; cannot run analysis.")
}

# SAMPLE_ID defaults to the seqz file name with the seqz extension stripped.
sample_id <- Sys.getenv("SAMPLE_ID")
if (!nzchar(sample_id)) {
    sample_id <- sub("\\.(small\\.)?seqz\\.gz$", "", basename(seqz_file))
    cat("==> SAMPLE_ID not set; derived from seqz file:", sample_id, "\n")
}

# OUTPUT_DIR defaults to the current working directory.
out_dir <- Sys.getenv("OUTPUT_DIR")
if (!nzchar(out_dir)) {
    out_dir <- "."
    cat("==> OUTPUT_DIR not set; defaulting to current directory\n")
}

n_cores <- as.integer(Sys.getenv("N_CORES"))
n_cores_fit <- as.integer(Sys.getenv("N_CORES_FIT"))

# Optional: local cytoBand file path for offline use. When empty, sequenza
# falls back to UCSC download (hgdownload.cse.ucsc.edu/.../cytoBand.txt.gz).
cytoband_file <- Sys.getenv("CYTOBAND_FILE")
if (!nzchar(cytoband_file)) {
    cytoband_file <- NULL
} else {
    cat("==> Using local cytoBand file:", cytoband_file, "\n")
}

# Reference assembly name ("hg38" or "hg19"), passed to sequenza.extract
# as the `assembly` argument. Common aliases are normalized automatically.
ref_type <- switch(tolower(trimws(Sys.getenv("ASSEMBLY"))),
    "hg38" = , "grch38" = , "assembly38" = "hg38",
    "hg19" = , "grch37" = , "assembly19" = , "b37" = "hg19",
    NA_character_)

if (is.na(ref_type)) {
    raw <- Sys.getenv("ASSEMBLY")
    if (nzchar(raw)) {
        stop("ASSEMBLY='", raw, "' is not recognized. ",
             "Use hg38 (GRCh38) or hg19 (GRCh37/b37).")
    }
    ref_type <- "hg38"
    cat("==> WARNING: ASSEMBLY not set; defaulting to", ref_type, "\n")
} else {
    cat("==> Reference type:", ref_type, "\n")
}

# Sample sex from the GENDER env var. sequenza.extract expects a logical
# `female`, so the raw value is mapped: female/f/xx -> TRUE, male/m/xy -> FALSE.
gender <- Sys.getenv("GENDER")
female <- switch(tolower(trimws(gender)),
    "female" = , "f" = , "xx" = TRUE,
    "male" = , "m" = , "xy" = FALSE,
    NA)
if (is.na(female)) {
    female <- FALSE
    cat("==> WARNING: GENDER not set or unrecognized ('", gender,
        "'); defaulting female to FALSE\n", sep = "")
} else {
    cat("==> Sample sex (GENDER='", gender, "') -> female = ", female,
        "\n", sep = "")
}

if (is.na(n_cores) || n_cores < 1) {
    n_cores <- max(1L, detectCores(logical = FALSE) - 1L, na.rm = TRUE)
}
if (is.na(n_cores_fit) || n_cores_fit < 1) {
    n_cores_fit <- max(1L, detectCores(logical = FALSE) - 1L, na.rm = TRUE)
}

cat("==> Using", n_cores, "cores for extract,", n_cores_fit, "cores for fit\n")
cat("==> Loading seqz file:", seqz_file, "\n")

seqz_data <- sequenza.extract(seqz_file, verbose = TRUE, parallel = n_cores, gamma = 80, kmin = 10,
                              cytoband_file = cytoband_file,
                              assembly = ref_type, female = female)

# Edge case guard: chr-subset / shallow data can yield zero segments passing
# sequenza.fit's BAF-method filters (N.ratio > 10, N.BAF > 1, length >= 3 Mb).
# Calling sequenza.fit on empty input crashes baf.bayes workers with
# "argument is of length zero". Pre-check and emit a header-only
# _segments.txt placeholder so downstream consumers can parse it as a
# 0-row table.
all_segs <- do.call(rbind, seqz_data$segments)
seg_len <- all_segs$end.pos - all_segs$start.pos
n_eligible <- sum(all_segs$N.ratio > 10 & all_segs$N.BAF > 1 & seg_len >= 3e6,
                  na.rm = TRUE)
if (n_eligible == 0) {
    cat("==> WARNING: 0 segments pass BAF filters",
        "(N.ratio > 10, N.BAF > 1, length >= 3 Mb).",
        "Skipping sequenza.fit and emitting header-only _segments.txt.\n")
    # Mirror sequenza.results column layout: seg.tab columns + baf.bayes columns.
    header_cols <- c(colnames(all_segs), "CNt", "A", "B", "LPP")
    writeLines(paste(header_cols, collapse = "\t"),
               con = file.path(out_dir, paste0(sample_id, "_segments.txt")))
    quit(save = "no", status = 0)
}

cat("==> Fitting cellularity and ploidy...\n")
CP <- sequenza.fit(seqz_data, mc.cores = n_cores_fit)

cat("==> Generating results...\n")
sequenza.results(
    sequenza.extract = seqz_data,
    cp.table = CP,
    sample.id = sample_id,
    out.dir = out_dir
)

cat("==> Sequenza analysis completed!\n")
cat("Results saved to:", out_dir, "\n")
