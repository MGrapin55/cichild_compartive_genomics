# Query NCBI SRA for genomic, paired-end runs from a list of species.
# Install once if needed: install.packages("rentrez")

library(rentrez)

# One species name per line in species.txt
species <- c("Serranochromis robustus")
#species <- trimws(readLines("path/to/species.txt", warn = FALSE))
species <- species[nzchar(species)]

# Retrieve RunInfo in batches to avoid oversized requests.
fetch_runinfo <- function(ids, batch_size = 200L) {
  batches <- split(ids, ceiling(seq_along(ids) / batch_size))
  
  results <- lapply(batches, function(batch) {
    csv_text <- entrez_fetch(
      db = "sra",
      id = batch,
      rettype = "runinfo",
      retmode = "text"
    )
    read.csv(text = csv_text, stringsAsFactors = FALSE, check.names = FALSE)
  })
  
  do.call(rbind, results)
}

# Return the first matching column name, ignoring capitalization.
find_column <- function(data, candidates) {
  matches <- match(tolower(candidates), tolower(names(data)))
  matches <- matches[!is.na(matches)]
  if (length(matches)) names(data)[matches[1]] else NULL
}

all_runs <- list()

for (sp in species) {
  query <- sprintf(
    '"%s"[Organism] AND genomic[Source] AND paired[Layout]',
    gsub('"', '\\"', sp)
  )
  
  search <- entrez_search(
    db = "sra",
    term = query,
    #retmax = 100000
  )
  
  if (search$count > length(search$ids)) {
    stop(
      sprintf(
        "The result limit was reached for %s (%s of %s records retrieved).",
        sp, length(search$ids), search$count
      )
    )
  }
  
  if (length(search$ids) == 0L) next
  
  runinfo <- fetch_runinfo(search$ids)
  
  # Confirm the requested criteria using the returned run metadata.
  source_col <- find_column(runinfo, "LibrarySource")
  layout_col <- find_column(runinfo, "LibraryLayout")
  
  if (!is.null(source_col) && !is.null(layout_col)) {
    keep <- toupper(runinfo[[source_col]]) == "GENOMIC" &
      toupper(runinfo[[layout_col]]) == "PAIRED"
    runinfo <- runinfo[!is.na(keep) & keep, , drop = FALSE]
  }
  
  if (nrow(runinfo) == 0L) next
  
  run_col <- find_column(runinfo, "Run")
  sample_col <- find_column(runinfo, "BioSample")
  project_col <- find_column(runinfo, "BioProject")
  study_col <- find_column(runinfo, c("SRAStudy", "Study"))
  pubmed_col <- find_column(runinfo, c("PubmedId", "PubMedID"))
  authors_col <- find_column(runinfo, "Authors")
  year_col <- find_column(runinfo, "Year")
  organism_col <- find_column(runinfo, c("ScientificName", "Organism"))
  
  get_value <- function(column, default = NA_character_) {
    if (is.null(column)) rep(default, nrow(runinfo)) else runinfo[[column]]
  }
  
  all_runs[[sp]] <- data.frame(
    PUBMED_ID = suppressWarnings(as.integer(get_value(pubmed_col))),
    AUTHORS = as.character(get_value(authors_col)),
    YEAR = suppressWarnings(as.integer(get_value(year_col))),
    # SRA study accessions (for example, SRP...) are alphanumeric.
    STUDY = as.character(get_value(study_col)),
    SPECIES = as.character(get_value(organism_col, sp)),
    GENOMIC_DNA = TRUE,
    PAIRED_END = TRUE,
    SRR_ACCESSION = as.character(get_value(run_col)),
    BIOSAMPLE_ACCESSION = as.character(get_value(sample_col)),
    BIOPROJECT_ACCESSION = as.character(get_value(project_col)),
    stringsAsFactors = FALSE
  )
}

if (length(all_runs)) {
  sra_results <- do.call(rbind, all_runs)
  rownames(sra_results) <- NULL
} else {
  sra_results <- data.frame(
    PUBMED_ID = integer(),
    AUTHORS = character(),
    YEAR = integer(),
    STUDY = character(),
    SPECIES = character(),
    GENOMIC_DNA = logical(),
    PAIRED_END = logical(),
    SRR_ACCESSION = character(),
    BIOSAMPLE_ACCESSION = character(),
    BIOPROJECT_ACCESSION = character()
  )
}

write.csv(sra_results, "ncbi_sra_results.csv", row.names = FALSE)