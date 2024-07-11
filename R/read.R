read.seqz <- function(file, n_lines = NULL, col_types = "ciciidddcddccc",
    chr_name = NULL, buffer = 33554432, parallel = 1,
    col_names = c("chromosome", "position", "base.ref", "depth.normal",
                  "depth.tumor", "depth.ratio", "Af", "Bf", "zygosity.normal",
                  "GC.percent", "good.reads", "AB.normal", "AB.tumor",
                  "tumor.strand"), ...) {

    if (is.null(n_lines)) {
        skip <- 1
        n_max <- Inf
    } else {
        n_lines <- round(sort(n_lines), 0)
        skip <- n_lines[1]

        n_max <- n_lines[2] - skip + 1
    }
    chr_name <- as.character(chr_name)
    tbi <- file.exists(paste(file, "tbi", sep = "."))
    if (tbi) {
        read.seqz.tbi(file, split_chr_coord(chr_name), col_names)
    } else {
        read.seqz.chr(file, chr_name = chr_name, col_types = col_types,
            col_names = col_names, skip = skip, buffer = buffer, parallel = parallel)
    }
}

read.seqz.chr <- function(file, chr_name, col_names,
    col_types, skip, buffer, parallel) {
    #con <- gzfile(file, "rb")
    #suppressWarnings(skip_line <- readLines(con, n = 1))
    #remove(skip_line)
    #parse_chunck <- function(x, chr_name, col_names, col_types) {
    #    x <- read_tsv(file = paste(mstrsplit(x), collapse = "\n"),
    #        col_types = col_types, skip = 0, n_max = Inf,
    #        col_names = col_names, progress = FALSE)
    #    x[x$chromosome == chr_name, ]
    #}
    #res <- chunk.apply(input = con, FUN = parse_chunck, chr_name = chr_name,
    #    col_names = col_names, col_types = col_types, CH.MAX.SIZE = buffer,
    #    parallel = parallel)
    #close(con)
    #res
    if (!is.null(chr_name)) {
        f <- function(x, pos) {
            subset(x, chromosome == chr_name)
        }
    } else {
        f <- function(x, pos) {
            x
        }
    }
    read_tsv_chunked(
        file, DataFrameCallback$new(f), col_types = col_types, skip = skip,
        col_names = col_names)
}

read.seqz.tbi <- function(file, chr_name, col_names) {
    #res <- tabix.read(file, chr_name)
    #res <- read_tsv(file = paste(mstrsplit(res), collapse = "\n"),
    #    col_types = col_types, skip = 0, n_max = Inf,
    #    col_names = col_names, progress = FALSE)
    res <- tabix.read.table(file, chr_name, col.names = TRUE, stringsAsFactors = FALSE)
    colnames(res) <- col_names
    as_tibble(
        res
    )
}