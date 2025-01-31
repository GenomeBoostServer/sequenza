# Cached environment for repeated reads
.seqz_cache <- new.env(parent = emptyenv())

read.seqz <- function(
  file, n_lines = NULL, col_types = "ciciidddcddccc", chr_name = NULL,
  buffer = 33554432, parallel = 1, col_names = c(
    "chromosome", "position", "base.ref",
    "depth.normal", "depth.tumor", "depth.ratio", "Af", "Bf", "zygosity.normal",
    "GC.percent", "good.reads", "AB.normal", "AB.tumor", "tumor.strand"
  ), cache = TRUE,
  ...
) {
  # Validate and normalize inputs
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  if (parallel < 1) {
    stop("parallel must be >= 1")
  }
  if (buffer < 1024) {
    warning("Very small buffer size may impact performance")
  }

  # Check cache for repeated reads
  if (cache) {
    cache_key <- digest::digest(list(file, n_lines, col_types, chr_name, buffer))
    if (exists(cache_key, envir = .seqz_cache)) {
      return(get(cache_key, envir = .seqz_cache))
    }
  }

  # Process line ranges more efficiently
  range_info <- if (!is.null(n_lines)) {
    if (length(n_lines) != 2) {
      stop("n_lines must be NULL or length 2")
    }
    list(skip = n_lines[1], n_max = diff(round(sort(n_lines))) + 1)
  } else {
    list(skip = 1, n_max = Inf)
  }

  # Determine reading strategy
  chr_name <- as.character(chr_name)
  result <- if (file.exists(paste(file, "tbi", sep = "."))) {
    read.seqz.tbi(file, split_chr_coord(chr_name), col_names)
  } else {
    read.seqz.chr(file,
      chr_name = chr_name, col_types = col_types, col_names = col_names,
      skip = range_info$skip, buffer = buffer, parallel = parallel
    )
  }

  # Cache result if enabled
  if (cache) {
    assign(cache_key, result, envir = .seqz_cache)
  }

  result
}

read.seqz.chr <- function(file, chr_name, col_names, col_types, skip, buffer, parallel) {
  # Initialize connection with cleanup
  con <- gzfile(file, "rb")
  on.exit(close(con))

  # Skip header efficiently
  if (skip > 0) {
    suppressWarnings(readLines(con, n = skip))
  }

  # Pre-allocate and optimize chunk processing
  parse_chunk <- function(x, chr_name, col_names, col_types) {
    # Process chunk data as a single string first
    chunk_text <- paste(mstrsplit(x), collapse = "\n")

    # Process chunk with optimized settings
    chunk_data <- read_tsv(
      file = chunk_text, col_types = col_types, col_names = col_names,
      skip = 0, n_max = Inf, progress = FALSE, show_col_types = FALSE, lazy = TRUE # Enable lazy reading
    )

    # Efficient chromosome filtering using data.table-style optimization
    if (!is.null(chr_name)) {
      idx <- chunk_data$chromosome == chr_name
      chunk_data[idx, , drop = FALSE]
    } else {
      chunk_data
    }
  }

  # Process chunks with correct parameters
  results <- chunk.apply(
    input = con, FUN = parse_chunk, chr_name = chr_name, col_names = col_names,
    col_types = col_types, CH.MAX.SIZE = buffer, CH.PARALLEL = parallel
  )

  # Convert to tibble efficiently
  as_tibble(results)
}

read.seqz.tbi <- function(file, chr_name, col_names) {
  # Handle different coordinate formats correctly
  tabix_range <- if (!is.null(chr_name)) {
    if (grepl(":", chr_name)) {
      # Parse coordinates in format 'chr:start-end'
      parts <- strsplit(chr_name, "[:-]")[[1]]
      if (length(parts) != 3) {
        stop(
          "Invalid coordinate format. Expected 'chr:start-end', got: ",
          chr_name
        )
      }
      chr_name # Keep original format if it includes coordinates
    } else {
      # For chromosome-only queries, add full range specification
      sprintf("%s:0-", chr_name)
    }
  } else {
    NULL
  }

  # Read tabix data efficiently
  res <- tabix.read.table(file, tabix_range, col.names = TRUE, stringsAsFactors = FALSE)

  # Set column names efficiently
  setNames(as_tibble(res), col_names)
}

# Helper function for coordinate splitting
split_chr_coord <- function(chr_name) {
  if (is.null(chr_name) || !grepl(":", chr_name)) {
    return(chr_name)
  }
  parts <- strsplit(chr_name, "[:-]")[[1]]
  if (length(parts) == 3) {
    # Validate numeric parts
    if (!all(grepl("^[0-9]+$", parts[2:3]))) {
      stop("Invalid coordinate format. Start and end positions must be numeric.")
    }
    return(chr_name) # Return original format if valid
  } else {
    stop("Invalid coordinate format. Expected 'chr:start-end', got: ", chr_name)
  }
}
