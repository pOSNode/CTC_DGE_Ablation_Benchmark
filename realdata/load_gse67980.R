# Load GSE67980 (Miyamoto et al. 2015) prostate CTC read counts.
# Returns a gene x sample integer matrix plus per-sample annotation.
GSE67980_FILES <- c(
  counts = "GSE67980_readCounts.txt.gz",
  props  = "GSE67980_sampleProperties.txt.gz")
GSE67980_URL <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE67nnn/GSE67980/suppl/"

#' Download GSE67980 from GEO if it is not already present.
fetch_gse67980 <- function(dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  for (f in GSE67980_FILES) {
    dest <- file.path(dir, f)
    if (!file.exists(dest)) {
      message("downloading ", f, " from GEO ...")
      utils::download.file(paste0(GSE67980_URL, f), dest, mode = "wb",
                           quiet = TRUE)
    }
  }
  invisible(dir)
}

load_gse67980 <- function(dir = file.path(Sys.getenv("BENCH_ROOT",
                            "/Users/owaissiddiqi/Documents/methods_paper"),
                            "realdata")) {
  fetch_gse67980(dir)
  cnt <- read.delim(gzfile(file.path(dir, "GSE67980_readCounts.txt.gz")),
                    check.names = FALSE, stringsAsFactors = FALSE)
  ann_cols <- c("ID", "Entrez GeneID", "uniGene", "symbol", "name")
  genes <- cnt[, intersect(ann_cols, names(cnt)), drop = FALSE]
  m <- as.matrix(cnt[, setdiff(names(cnt), ann_cols), drop = FALSE])
  storage.mode(m) <- "integer"
  rownames(m) <- make.unique(ifelse(is.na(genes$symbol) | genes$symbol == "",
                                    genes$ID, genes$symbol))

  props <- read.delim(gzfile(file.path(dir, "GSE67980_sampleProperties.txt.gz")),
                      check.names = FALSE, stringsAsFactors = FALSE)
  names(props)[1:2] <- c("title", "source")
  props <- props[match(colnames(m), props$title), , drop = FALSE]
  props$donor <- props[["characteristics: donor"]]
  list(counts = m, samples = props, genes = genes)
}
