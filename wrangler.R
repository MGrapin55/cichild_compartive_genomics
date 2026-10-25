# ==============================================================================
# 1. SETUP & INITIALIZATION
# ==============================================================================

# Query NCBI SRA for genomic, paired-end runs from a list of species.
# Install once if needed: install.packages("rentrez")

library(rentrez)

# Define the target species names (one name per element)
species <- c("Serranochromis robustus", "Metriaclima mbenji")

# Alternative: Read species names dynamically from an external text file
# species <- trimws(readLines("path/to/species.txt", warn = FALSE))

# Clean the vector by removing any empty or whitespace-only strings
species <- species[nzchar(species)]


# ==============================================================================
# 2. FUNCTIONS SECTION
# ==============================================================================

# Purpose: NCBI limits individual search returns. This function implements 
# pagination using `retstart` and `retmax` inside a `while` loop to guarantee 
# all matching SRA run IDs are fetched, even for large species datasets.
search_sra_ids <- function(term, page_size = 500L) {
  first_page <- entrez_search(
    db = "sra",
    term = term,
    retmax = page_size,
    retstart = 0L
  )
  ids <- first_page$ids
  
  while (length(ids) < first_page$count) {
    next_page <- entrez_search(
      db = "sra",
      term = term,
      retmax = page_size,
      retstart = length(ids)
    )
    
    if (length(next_page$ids) == 0L) {
      stop(
        sprintf(
          "SRA search stopped after %s of %s results.",
          length(ids), first_page$count
        )
      )
    }
    
    ids <- c(ids, next_page$ids)
  }
  
  ids
}

# Purpose: Prevents HTTP query errors by breaking large ID lists into smaller 
# chunks (batches of 200), downloading their `runinfo` metadata in CSV format, 
# and row-binding them together.
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

# Purpose: NCBI schema headers can vary in capitalization across different queries 
# (e.g., `BioProject` vs `bioproject`). This helper safely finds the correct column 
# name regardless of casing.
find_column <- function(data, candidates) {
  matches <- match(tolower(candidates), tolower(names(data)))
  matches <- matches[!is.na(matches)]
  if (length(matches)) names(data)[matches[1]] else NULL
}

# Purpose: Connects sequencing projects to published literature. It tries a direct 
# database link first (`entrez_link`), and if that fails, executes a fallback text 
# search on PubMed using the BioProject accession.
linked_pubmed_ids <- function(project_accessions) {
  projects <- unique(project_accessions[!is.na(project_accessions) &
                                        nzchar(project_accessions)])
  
  setNames(lapply(projects, function(project) {
    project_search <- entrez_search(
      db = "bioproject",
      term = project,
      retmax = 1L
    )
    if (length(project_search$ids) == 0L) return(character())
    
    # Strategy 1: Try the direct database-to-database link (elink)
    links <- entrez_link(
      dbfrom = "bioproject",
      id = project_search$ids[1],
      db = "pubmed"
    )
    linked_ids <- as.character(unlist(links$links, use.names = FALSE))
    
    # Strategy 2: Fallback to text searching PubMed if elink is empty
    if (length(linked_ids) == 0L) {
      fallback_search <- entrez_search(
        db = "pubmed",
        term = paste0(project, "[ALL]"),
        retmax = 1L
      )
      linked_ids <- fallback_search$ids
    }
    
    linked_ids
  }), projects)
}

# Purpose: Downloads article metadata in raw XML format, parses out authors and 
# publication years using `xml2`, and structures them into a clean tabular format.
fetch_pubmed_metadata <- function(ids, batch_size = 100L) {
  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  if (length(ids) == 0L) {
    return(data.frame(
      PUBMED_ID = character(),
      AUTHORS = character(),
      YEAR = integer(),
      stringsAsFactors = FALSE
    ))
  }
  
  batches <- split(ids, ceiling(seq_along(ids) / batch_size))
  articles <- unlist(lapply(batches, function(batch) {
    xml_text <- entrez_fetch(
      db = "pubmed",
      id = batch,
      rettype = "xml",
      parsed = FALSE
    )
    document <- xml2::read_xml(xml_text)
    xml2::xml_find_all(document, ".//PubmedArticle")
  }), recursive = FALSE)
  
  if (length(articles) == 0L) {
    return(data.frame(
      PUBMED_ID = character(),
      AUTHORS = character(),
      YEAR = integer(),
      stringsAsFactors = FALSE
    ))
  }
  
  xml_value <- function(node, xpath) {
    values <- xml2::xml_find_all(node, xpath)
    if (length(values) == 0L) NA_character_ else xml2::xml_text(values[[1]])
  }
  
  get_author <- function(author) {
    collective <- xml_value(author, "./CollectiveName")
    if (!is.na(collective) && nzchar(collective)) return(collective)
    
    last_name <- xml_value(author, "./LastName")
    initials <- xml_value(author, "./Initials")
    if (is.na(initials)) initials <- xml_value(author, "./ForeName")
    if (is.na(last_name)) return(NA_character_)
    if (is.na(initials)) last_name else paste(last_name, initials)
  }
  
  rows <- lapply(articles, function(article) {
    author_nodes <- xml2::xml_find_all(
      article,
      ".//MedlineCitation/Article/AuthorList/Author"
    )
    authors <- vapply(author_nodes, get_author, character(1))
    authors <- authors[!is.na(authors) & nzchar(authors)]
    
    year <- xml_value(article, ".//JournalIssue/PubDate/Year")
    if (is.na(year)) year <- xml_value(article, ".//ArticleDate/Year")
    if (is.na(year)) {
      medline_date <- xml_value(article, ".//JournalIssue/PubDate/MedlineDate")
      year_match <- regmatches(medline_date, regexpr("[0-9]{4}", medline_date))
      year <- if (length(year_match)) year_match else NA_character_
    }
    
    data.frame(
      PUBMED_ID = xml_value(article, ".//MedlineCitation/PMID"),
      AUTHORS = if (length(authors)) paste(authors, collapse = "; ") else NA_character_,
      YEAR = suppressWarnings(as.integer(year)),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, rows)
}


# ==============================================================================
# 3. MAIN PROCESSING SECTION
# ==============================================================================

all_runs <- list()

for (sp in species) {
  query <- sprintf(
    '"%s"[Organism] AND genomic[Source] AND paired[Layout]',
    gsub('"', '\\"', sp)
  )
  
  sra_ids <- search_sra_ids(query)
  if (length(sra_ids) == 0L) next
  
  runinfo <- fetch_runinfo(sra_ids)
  
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
  
  projects <- as.character(get_value(project_col))
  pubmed_ids <- as.character(get_value(pubmed_col))
  pubmed_ids[is.na(pubmed_ids) | !nzchar(trimws(pubmed_ids))] <- NA_character_
  
  project_links <- linked_pubmed_ids(projects)
  missing_pubmed <- which(is.na(pubmed_ids))
  
  # Keep one citation per run to match the existing PUBMED_ID output column.
  for (row in missing_pubmed) {
    project <- projects[row]
    if (!is.na(project) && project %in% names(project_links)) {
      linked_ids <- project_links[[project]]
      if (length(linked_ids)) pubmed_ids[row] <- linked_ids[1]
    }
  }
  
  publication_metadata <- fetch_pubmed_metadata(pubmed_ids)
  metadata_rows <- match(pubmed_ids, publication_metadata$PUBMED_ID)
  metadata_authors <- publication_metadata$AUTHORS[metadata_rows]
  metadata_years <- publication_metadata$YEAR[metadata_rows]
  
  authors <- as.character(get_value(authors_col))
  missing_authors <- is.na(authors) | !nzchar(trimws(authors))
  authors[missing_authors] <- metadata_authors[missing_authors]
  
  years <- suppressWarnings(as.integer(get_value(year_col)))
  missing_years <- is.na(years)
  years[missing_years] <- metadata_years[missing_years]
  
  all_runs[[sp]] <- data.frame(
    PUBMED_ID = suppressWarnings(as.integer(pubmed_ids)),
    AUTHORS = authors,
    YEAR = years,
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