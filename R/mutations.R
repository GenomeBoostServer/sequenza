mut.fractions <- function(AB.tumor, Af, tumor.strand) {
    # Input validation
    if (length(AB.tumor) != length(Af) || length(AB.tumor) != length(tumor.strand)) {
        stop("Input vectors must have equal lengths")
    }
    
    # Calculate complement frequency
    F <- 1 - Af
    
    # Helper function to split strings
    split_string <- function(x) {
        unlist(strsplit(as.character(x), split = "[:]"))
    }
    
    # Parse base mutations
    base.mut <- lapply(AB.tumor, split_string)
    base.fw <- lapply(tumor.strand, split_string)
    
    # Function to extract frequency from split data
    frequencify <- function(x) {
        # Split the data and handle potential errors
        split_data <- unlist(strsplit(unlist(x), split = ".", fixed = TRUE))
        
        # Validate split data
        if (length(split_data) < 2) {
            warning(paste("Invalid data format:", x))
            return(NULL)
        }
        
        # Extract base name and value
        base.name <- substr(split_data[1], 1, nchar(split_data[1]) - 1)
        base.val <- as.numeric(paste("0", split_data[2], sep = "."))
        
        setNames(base.val, base.name)
    }
    
    # Calculate base frequencies
    base.freqs <- lapply(base.mut, frequencify)
    fw.freqs <- lapply(base.fw, frequencify)
    
    # Count base mutations
    n.base.mut <- sapply(base.mut, length)
    
    # Find maximum frequency
    max.fq <- function(x) {
        # Handle potential NA or NULL values
        if (is.null(base.freqs[[x]]) || is.null(fw.freqs[[x]])) {
            return(rep(NA, 4))
        }
        
        freq.rel <- base.freqs[[x]] / F[x]
        f.max <- which.max(freq.rel)
        
        c(
            freq.rel[f.max], 
            names(base.freqs[[x]])[f.max],
            base.freqs[[x]][f.max], 
            fw.freqs[[x]][f.max]
        )
    }
    
    # Apply max frequency function
    max.freqs <- do.call(rbind, lapply(seq_along(F), max.fq))
    
    # Create result dataframe with error checking
    result <- data.frame(
        base.count = as.integer(n.base.mut),
        maj.base.freq = as.numeric(max.freqs[, 1]),
        base = as.character(max.freqs[, 2]),
        freq = as.numeric(max.freqs[, 3]),
        fw.freq = as.numeric(max.freqs[, 4])
    )
    
    # Validate result dataframe
    if (nrow(result) != length(AB.tumor)) {
        warning(paste(
            "Mismatch in result rows. Expected:", 
            length(AB.tumor), 
            "Actual:", 
            nrow(result)
        ))
    }
    
    return(result)
}


mutation.table <- function(
    seqz.tab, mufreq.treshold = 0.15, min.reads = 40, min.reads.normal = 10,
    max.mut.types = 3, min.type.freq = 0.9, min.fw.freq = 0,
    segments = NULL
) {
    chroms <- unique(seqz.tab$chromosome)
    hom.filt <- seqz.tab$zygosity.normal == "hom" & seqz.tab$AB.tumor !=
        "."
    seqz.tab <- seqz.tab[hom.filt, ]
    reads.filt <- seqz.tab$good.reads >= min.reads & seqz.tab$depth.normal >=
        min.reads.normal
    seqz.tab <- seqz.tab[reads.filt, ]
    mufreq.filt <- seqz.tab$Af <= (1 - mufreq.treshold)
    seqz.tab <- seqz.tab[mufreq.filt, ]
    if (!is.null(segments)) {
        for (i in 1:nrow(segments)) {
            pos.filt <- seqz.tab$chromosome == segments$chromosome[i] &
                seqz.tab$position >= segments$start.pos[i] &
                seqz.tab$position <= segments$end.pos[i]
            seqz.tab$adjusted.ratio[pos.filt] <- segments$depth.ratio[i]
        }
    }
    seqz.dummy <- data.frame(
        chromosome = chroms, position = 1, GC.percent = NA, good.reads = NA,
        adjusted.ratio = NA, F = 0, mutation = "NA", stringsAsFactors = FALSE
    )
    if (nrow(seqz.tab) >=
        1) {
        mu.fracts <- mut.fractions(
            AB.tumor = seqz.tab$AB.tumor, Af = seqz.tab$Af, tumor.strand = seqz.tab$tumor.strand
        )
        mufreq.filt <- mu.fracts$freq >= mufreq.treshold
        type.filt <- mu.fracts$base.count <= max.mut.types
        prop.filt <- mu.fracts$maj.base.freq >= min.type.freq
        if (!is.na(min.fw.freq)) {
            fw.2 <- 1 - min.fw.freq
            fw.2 <- sort(c(fw.2, min.fw.freq))
            fw.filt <- mu.fracts$fw.freq >= fw.2[1] & mu.fracts$fw.freq <=
                fw.2[2]
            mufreq.filt <- mufreq.filt & type.filt & prop.filt &
                fw.filt
        } else {
            mufreq.filt <- mufreq.filt & type.filt & prop.filt
        }
        mut.type <- paste(seqz.tab$AB.normal, mu.fracts$base, sep = ">")
        seqz.tab <- seqz.tab[, c(
            "chromosome", "position", "GC.percent", "good.reads",
            "adjusted.ratio"
        )]
        seqz.tab <- cbind(seqz.tab, F = mu.fracts$freq, mutation = mut.type)
        rbind(seqz.tab[mufreq.filt, ], seqz.dummy)
    } else {
        seqz.dummy
    }
}
