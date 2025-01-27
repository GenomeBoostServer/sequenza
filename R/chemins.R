
select_chromosomes_with_centromere <- function(assembly, all_sequences) {
    golden_path <- paste("http://hgdownload.cse.ucsc.edu", "goldenPath", assembly,
        "database", "cytoBand.txt.gz", sep = "/")
    arms <- get_assembly(golden_path, prefix = "chr")
    chromosomes <- unique(arms[arms$arm != "n", ]$chromosome)
    select_sequences <- all_sequences %in% chromosomes
    if (sum(select_sequences) == 0) {
        select_sequences <- paste0("chr", all_sequences) %in% chromosomes
    }
    all_sequences[select_sequences]
}

get_arms <- function(chromosome) {
    chrom <- chromosome$chromosome[1]
    start <- min(chromosome$start)
    end <- max(chromosome$end)
    acen <- chromosome$type %in% "acen"
    if (sum(acen) == 0) {
        res <- tibble(chromosome = chrom, start = start, end = end, arm = "n")
    } else {
        p_end <- min(chromosome$start[acen])
        q_start <- max(chromosome$end[acen])
        res <- tibble(chromosome = chrom, start = start, end = p_end, arm = "p")
        res <- add_row(res, chromosome = chrom, start = q_start, end = end, arm = "q")
    }
    res
}

get_assembly <- function(name, url = NULL, prefix = "chr") {
    if (is.null(url)) {
        # fetch assembly by name
        file_path <- name
    } else {
        file_path <- url
    }
    cytobands <- read_tsv(file_path, col_types = "ciicc", col_names = c("chromosome",
        "start", "end", "cytoband", "type"))
    arms <- do.call(rbind, lapply(split(cytobands, cytobands$chromosome), get_arms))
    arms$chromosome <- str_replace(arms$chromosome, "^chr", prefix)
    return(arms)
}

slide_matrix <- function(x, position = NULL, w = 100, smooth = TRUE, method = c("kstest",
    "meandiff", "both"), verbose = TRUE) {

    if (is.null(position)) {
        position <- seq_len(length(x))
    }
    if (length(method) > 1) {
        method <- method[1]
    }
    if (method == "kstest") {
        method <- 1
    } else if (method == "meandiff") {
        method <- 2
    } else if (method == "both") {
        method <- 3
    } else {
        method <- 1
    }
    .Call(`_sequenza_slide_matrix`, x, position, w, smooth, method, verbose)
}


get_peaks <- function(x, w = 100, position = NULL) {
    if (is.null(position)) {
        position <- seq_len(length(x))
    }
    .Call(`_sequenza_get_peaks`, x, position, w)
}

get_gaps_peaks <- function(x, w = 100, position = NULL, arms) {
    apply(arms, 1, FUN = function(arm) {
        index <- between(position, as.numeric(arm["start"]), as.numeric(arm["end"]))
        if (sum(index) > 0) {
            coords <- get_peaks(x = x[index], position = position[index], w = w)
            if (length(coords) > 0) {
                pos_start <- coords[-length(coords)]
                pos_end <- coords[-1]
                breaks <- tibble(start.pos = pos_start, end.pos = pos_end)

                breaks <- add_row(breaks, start.pos = min(position[index]), end.pos = coords[1] -
                  1, .before = 1)
                breaks <- add_row(breaks, start.pos = coords[length(coords)] + 1,
                  end.pos = max(position[index]))
            } else (breaks <- tibble(start.pos = min(position[index]), end.pos = max(position[index])))
            breaks
        }
    })
}
