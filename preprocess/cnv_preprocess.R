## ========== WORKING DIRECTORY ==========
setwd("~/Desktop/RA/BrCa/Code/PhylogicNDT/phylogic_inputs")

cnv_dir      <- "~/Desktop/RA/BrCa/Data/cnv_segment"      # Directory of raw CNV files (.tsv / .tsv.gz)
summary_csv  <- "~/Desktop/RA/BrCa/Results/brca-garvan-wgts_annotated_summary.csv"  # Summary table containing purity
out_dir_allelic <- "seg_alleliccapseg"   # Output folder for alleliccapseg-format files

dir.create(out_dir_allelic, showWarnings = FALSE, recursive = TRUE)

## ========== READ PURITY SUMMARY ==========
summary_df <- read.csv(summary_csv, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
if (!"Purity_value" %in% names(summary_df)) {
  summary_df$Purity_value <- suppressWarnings(as.numeric(sub(" .*", "", summary_df$Purity)))
}
id_col     <- "WGS-tumor-lib"   # Sample ID column (ensure it matches your summary table)
purity_col <- "Purity_value"
stopifnot(id_col %in% names(summary_df), purity_col %in% names(summary_df))

## ========== CONVERSION: STANDARDIZE TO alleliccapseg FORMAT ==========
# Output columns:
#   Chromosome, Start.bp, End.bp, mu.minor, sigma.minor, mu.major, sigma.major
normalize_to_allelic <- function(df, default_sigma = 0.20) {
  idx <- function(name) which(tolower(names(df)) == tolower(name))
  has <- function(name) length(idx(name)) == 1
  
  # chr/start/end
  chr_col <- NULL; for (c in c("chromosome", "chr")) if (has(c)) { chr_col <- idx(c); break }
  if (is.null(chr_col)) stop("Missing Chromosome column (chromosome/chr)")
  start_col <- NULL; for (k in c("start", "start.bp")) if (has(k)) { start_col <- idx(k); break }
  if (is.null(start_col)) stop("Missing start position column (start/start.bp)")
  end_col <- NULL; for (k in c("end", "end.bp")) if (has(k)) { end_col <- idx(k); break }
  if (is.null(end_col)) stop("Missing end position column (end/end.bp)")
  
  # A1/A2 from either modal.a1/modal.a2 or purple (major/minorAlleleCopyNumber)
  if (has("modal.a1") && has("modal.a2")) {
    a1_col <- idx("modal.a1"); a2_col <- idx("modal.a2")
  } else if (has("majorallelecopynumber") && has("minorallelecopynumber")) {
    a1_col <- idx("majorallelecopynumber"); a2_col <- idx("minorallelecopynumber")
  } else {
    stop("Missing allelic copy number columns (modal.a1/modal.a2 or major/minorAlleleCopyNumber)")
  }
  
  # deviation (optional)
  maj_dev <- if (has("majorallelecopynumberdeviation")) idx("majorallelecopynumberdeviation") else NA
  min_dev <- if (has("minorallelecopynumberdeviation")) idx("minorallelecopynumberdeviation") else NA
  
  out <- data.frame(
    Chromosome = df[[chr_col]],
    Start.bp   = as.integer(df[[start_col]]),
    End.bp     = as.integer(df[[end_col]]),
    mu.major   = as.numeric(df[[a1_col]]),
    mu.minor   = as.numeric(df[[a2_col]]),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  
  # Remove "chr" prefix and keep standard chromosomes
  out$Chromosome <- sub("^chr", "", out$Chromosome, ignore.case = TRUE)
  keep_chr <- c(as.character(1:22), "X", "Y")
  out <- out[out$Chromosome %in% keep_chr, , drop = FALSE]
  
  # Ensure major >= minor
  swap <- out$mu.major < out$mu.minor
  if (any(swap)) {
    tmp <- out$mu.major[swap]; out$mu.major[swap] <- out$mu.minor[swap]; out$mu.minor[swap] <- tmp
  }
  
  # sigma: use deviation if available, otherwise assign default
  if (!is.na(maj_dev) && !is.na(min_dev)) {
    out$sigma.major <- suppressWarnings(as.numeric(df[[maj_dev]]))
    out$sigma.minor <- suppressWarnings(as.numeric(df[[min_dev]]))
    out$sigma.major[!is.finite(out$sigma.major)] <- default_sigma
    out$sigma.minor[!is.finite(out$sigma.minor)] <- default_sigma
  } else {
    out$sigma.major <- default_sigma
    out$sigma.minor <- default_sigma
  }
  
  # Sort
  chr_order <- function(v) as.integer(factor(v, levels = c(as.character(1:22), "X", "Y")))
  out <- out[order(chr_order(out$Chromosome), out$Start.bp, out$End.bp), ]
  
  # Final column order for PhylogicNDT
  out <- out[, c("Chromosome","Start.bp","End.bp","mu.minor","sigma.minor","mu.major","sigma.major")]
  out
}

## ========== BATCH PROCESSING ==========
cnv_files <- list.files(cnv_dir, pattern="\\.tsv(\\.gz)?$", full.names=TRUE)
if (length(cnv_files) == 0L) stop(sprintf("No .tsv/.tsv.gz files found in directory: %s", cnv_dir))

for (fp in cnv_files) {
  sid <- sub("\\..*$", "", basename(fp))  # Use filename prefix as sample ID
  if (!(sid %in% summary_df[[id_col]])) {
    message(sprintf("[SKIP] %s — Sample ID not found in summary: %s", basename(fp), sid)); next
  }
  
  seg_raw <- tryCatch({
    if (grepl("\\.gz$", fp, ignore.case = TRUE)) {
      read.delim(gzfile(fp), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      read.delim(fp, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    }
  },
  error = function(e) { message("[READ ERROR] ", fp, " : ", e$message); return(NULL) }
  )
  if (is.null(seg_raw)) next
  
  seg_allelic <- tryCatch(normalize_to_allelic(seg_raw, default_sigma = 0.20),
                          error = function(e){ message("[FORMAT ERROR] ", fp, " : ", e$message); return(NULL) })
  if (is.null(seg_allelic)) next
  
  out1 <- file.path(out_dir_allelic, paste0(sid, ".seg.tsv"))
  write.table(seg_allelic, out1, sep="\t", quote=FALSE, row.names=FALSE)
  
  message(sprintf("[OK] %-18s | allelic=%s", sid, basename(out1)))
}

message("All processing completed.")