## ====== BASIC PATHS (edit as needed) ======
library(data.table)
library(GenomicRanges)

maf_dir  <- "/Users/jiajunzhan/Desktop/RA/BrCa/Data/sample_maf_files"
seg_dir  <- "/Users/jiajunzhan/Desktop/RA/BrCa/Code/Evolution/data/cnv_with_purity"  # Raw PURPLE segments (may include purity)
## Output root visible to container (you mapped this to /work)
out_root <- "/Users/jiajunzhan/Desktop/RA/BrCa/Code/PhylogicNDT/phylogic_inputs"
out_maf  <- file.path(out_root, "maf_with_cn")
out_seg  <- file.path(out_root, "seg_timing")        # Canonical 5-column seg (Cluster expects Start.bp/End.bp)
sif_file <- file.path(out_root, "Patient.sif")

dir.create(out_maf, showWarnings = FALSE, recursive = TRUE)
dir.create(out_seg, showWarnings = FALSE, recursive = TRUE)

## ====== Helper: normalize segment headers -> 5 columns (Chromosome/Start/End/A1.Seg.CN/A2.Seg.CN) ======
to_seg5 <- function(df){
  nm <- tolower(names(df)); names(df) <- nm
  # Support PURPLE or already 5-column input
  if (all(c("chromosome","start","end","majorallelecopynumber","minorallelecopynumber") %in% nm)) {
    x <- data.frame(
      Chromosome  = df$chromosome,
      Start       = as.integer(df$start),
      End         = as.integer(df$end),
      `A1.Seg.CN` = as.numeric(df$minorallelecopynumber),
      `A2.Seg.CN` = as.numeric(df$majorallelecopynumber),
      check.names = FALSE
    )
  } else if (all(c("chromosome","start","end","a1.seg.cn","a2.seg.cn") %in% nm)) {
    x <- data.frame(
      Chromosome  = df$chromosome, Start=as.integer(df$start), End=as.integer(df$end),
      `A1.Seg.CN` = as.numeric(df$`a1.seg.cn`), `A2.Seg.CN` = as.numeric(df$`a2.seg.cn`),
      check.names = FALSE
    )
  } else if (all(c("chromosome","start.bp","end.bp","a1.seg.cn","a2.seg.cn") %in% nm)) {
    x <- data.frame(
      Chromosome  = df$chromosome, Start=as.integer(df$`start.bp`), End=as.integer(df$`end.bp`),
      `A1.Seg.CN` = as.numeric(df$`a1.seg.cn`), `A2.Seg.CN` = as.numeric(df$`a2.seg.cn`),
      check.names = FALSE
    )
  } else {
    stop("Segment file is missing required columns: Chromosome/Start/End/A1.Seg.CN/A2.Seg.CN or PURPLE columns")
  }
  # Standardize
  x$Chromosome <- sub("^chr","",x$Chromosome,ignore.case=TRUE)
  keep_chr <- c(as.character(1:22),"X","Y")
  x <- x[x$Chromosome %in% keep_chr, , drop=FALSE]
  x[["A1.Seg.CN"]] <- pmax(0, as.numeric(x[["A1.Seg.CN"]]))
  x[["A2.Seg.CN"]] <- pmax(0, as.numeric(x[["A2.Seg.CN"]]))
  ord <- as.integer(factor(x$Chromosome, levels=c(as.character(1:22),"X","Y")))
  x <- x[order(ord, x$Start, x$End), ]
  x
}

## ====== Helper: read seg and return GRanges + purity (if available) ======
read_seg_as_gr <- function(seg_fp){
  seg_raw <- fread(seg_fp, sep="\t", header=TRUE, check.names=FALSE, showProgress=FALSE)
  seg5    <- to_seg5(seg_raw)
  # purity column if present
  nm <- tolower(names(seg_raw))
  purity <- NA_real_
  pur_col <- c("purity","sample_purity","tumor_purity")
  hit <- pur_col[pur_col %in% nm]
  if (length(hit)) {
    v <- as.numeric(seg_raw[[ hit[1] ]])
    purity <- suppressWarnings(as.numeric(v[!is.na(v)][1]))
  }
  seg5$`A1.Seg.CN` <- pmax(0, seg5$`A1.Seg.CN`)
  seg5$`A2.Seg.CN` <- pmax(0, seg5$`A2.Seg.CN`)
  gr <- GRanges(
    seqnames = seg5$Chromosome,
    ranges   = IRanges(seg5$Start, seg5$End),
    a1       = seg5$`A1.Seg.CN`,
    a2       = seg5$`A2.Seg.CN`
  )
  list(gr=gr, seg5=seg5, purity=purity)
}

## ====== Helper: annotate MAF with local_cn_a1/a2 (map sites to segments) ======
annotate_maf_with_cn <- function(maf_dt, gr_seg){
  # Harmonize common headers
  cn <- tolower(names(maf_dt))
  ren <- c(Chromosome="chromosome",
           Start_position="start_position",
           Reference_Allele="Reference_Allele",
           Tumor_Seq_Allele2="Tumor_Seq_Allele2",
           t_ref_count="t_ref_count",
           t_alt_count="t_alt_count")
  for (std in names(ren)) {
    if (!(tolower(std) %in% cn) && (ren[std] %in% cn)) {
      setnames(maf_dt, ren[std], std)
    }
  }
  need <- c("Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2")
  if (!all(need %in% names(maf_dt))) {
    stop(paste("MAF missing required columns:", paste(setdiff(need, names(maf_dt)), collapse=", ")))
  }
  maf_dt[, Chromosome := sub("^chr","",Chromosome, ignore.case=TRUE)]
  keep_chr <- c(as.character(1:22),"X","Y")
  maf_dt   <- maf_dt[Chromosome %in% keep_chr]
  gr_maf <- GRanges(seqnames = maf_dt$Chromosome,
                    ranges   = IRanges(as.integer(maf_dt$Start_Position),
                                       as.integer(maf_dt$Start_Position)))
  hit <- findOverlaps(gr_maf, gr_seg, select="first")
  a1 <- rep(NA_integer_, nrow(maf_dt))
  a2 <- rep(NA_integer_, nrow(maf_dt))
  ok <- !is.na(hit)
  a1[ok] <- mcols(gr_seg)$a1[hit[ok]]
  a2[ok] <- mcols(gr_seg)$a2[hit[ok]]
  # Default to 1/1 if not overlapping any segment
  a1[is.na(a1)] <- 1L
  a2[is.na(a2)] <- 1L
  maf_dt[, local_cn_a1 := a1]
  maf_dt[, local_cn_a2 := a2]
  maf_dt
}

## ====== Iterate over all MAFs: produce withCN MAFs, canonical segs, and SIF ======
maf_files <- list.files(maf_dir, pattern="\\.maf(\\.gz)?$", full.names=TRUE)
stopifnot(length(maf_files) > 0)

sif_rows <- list()

for (mf in maf_files) {
  sid <- sub("-somatic-PASS.*$","", basename(mf))
  sid <- sub("\\.maf(\\.gz)?$","", sid)
  
  ## Locate matching seg file (according to your naming pattern)
  cand <- c(
    file.path(seg_dir, paste0(sid, "_with_purity.tsv")),
    file.path(seg_dir, paste0(sid, ".tsv")),
    file.path(seg_dir, paste0(sid, ".seg.tsv"))
  )
  seg_fp <- cand[file.exists(cand)][1]
  if (is.na(seg_fp)) {
    message(sprintf("[SKIP] No segment file found for sample: %s", sid))
    next
  }
  
  # Read seg -> GRanges + purity
  seg_obj <- tryCatch(read_seg_as_gr(seg_fp),
                      error=function(e){ message("[SEG READ ERROR] ", seg_fp, " : ", e$message); NULL })
  if (is.null(seg_obj)) next
  gr_seg <- seg_obj$gr
  seg5   <- seg_obj$seg5
  purity <- seg_obj$purity
  
  # Write canonical seg (Start.bp/End.bp) for clustering
  seg5_cluster <- seg5
  names(seg5_cluster)[names(seg5_cluster)=="Start"] <- "Start.bp"
  names(seg5_cluster)[names(seg5_cluster)=="End"]   <- "End.bp"
  seg_out <- file.path(out_seg, paste0(sid, ".seg.tsv"))
  fwrite(seg5_cluster, seg_out, sep="\t", quote=FALSE)
  message(sprintf("[SEG] Wrote: %s", seg_out))
  
  # Read MAF, standardize headers, add local_cn_*, and write out
  maf_dt <- fread(mf, sep="\t", header=TRUE, quote="", check.names=FALSE, showProgress=FALSE)
  maf_dt <- annotate_maf_with_cn(maf_dt, gr_seg)
  maf_out <- file.path(out_maf, paste0(sid, "-somatic-PASS.withCN.maf"))
  fwrite(maf_dt, maf_out, sep="\t", quote=FALSE)
  message(sprintf("[MAF] Wrote: %s (rows=%d)", maf_out, nrow(maf_dt)))
  
  # Record SIF row (container-visible paths)
  maf_ct <- file.path("./myinput/maf_with_cn", basename(maf_out))
  seg_ct <- file.path("./myinput/maf_with_cn/seg_timing", basename(seg_out))
  sif_rows[[length(sif_rows)+1]] <- data.frame(
    sample_id = sid,
    maf_fn    = maf_ct,
    seg_fn    = seg_ct,
    purity    = purity,
    timepoint = 0,
    check.names = FALSE
  )
}

## Write SIF
if (length(sif_rows)) {
  sif <- data.table::rbindlist(sif_rows, use.names=TRUE, fill=TRUE)
  fwrite(sif, sif_file, sep="\t", quote=FALSE)
  message(sprintf("[SIF] Wrote: %s (%d samples)", sif_file, nrow(sif)))
} else {
  message("[SIF] Nothing to write (no paired maf/seg found)")
}



library(data.table)

out_maf <- "/Users/jiajunzhan/Desktop/RA/BrCa/Code/PhylogicNDT/phylogic_inputs/maf_with_cn"

maf_files <- list.files(out_maf, pattern="\\.maf(\\.gz)?$", full.names=TRUE)
stopifnot(length(maf_files) > 0)

for (mf in maf_files) {
  dt <- fread(mf, sep = "\t", header = TRUE, quote = "", check.names = FALSE, showProgress = FALSE)
  
  # Normalize headers to exactly what PhylogicNDT expects
  # Start_Position -> Start_position (lowercase p)
  if ("Start_Position" %in% names(dt) && !("Start_position" %in% names(dt))) {
    setnames(dt, "Start_Position", "Start_position")
  }
  if ("End_Position" %in% names(dt) && !("End_position" %in% names(dt))) {
    setnames(dt, "End_Position", "End_position")
  }
  if ("HGVSp_Short" %in% names(dt) && !("Protein_Change" %in% names(dt))) {
    setnames(dt, "HGVSp_Short", "Protein_Change")
  }
  # Map common alternates to standard column names (only if the standard is missing)
  ren_map <- list(
    Chromosome        = c("chromosome"),
    Reference_Allele  = c("reference_allele"),
    Tumor_Seq_Allele2 = c("tumor_seq_allele2"),
    t_ref_count       = c("t.ref.count","t_ref_cnt"),
    t_alt_count       = c("t.alt.count","t_alt_cnt")
  )
  for (std in names(ren_map)) {
    if (!(std %in% names(dt))) {
      for (alt in ren_map[[std]]) {
        if (alt %in% names(dt)) { setnames(dt, alt, std); break }
      }
    }
  }
  
  # Ensure copy-number columns exist
  if (!("local_cn_a1" %in% names(dt))) dt[, local_cn_a1 := NA_integer_]
  if (!("local_cn_a2" %in% names(dt))) dt[, local_cn_a2 := NA_integer_]
  
  # Overwrite the original file
  fwrite(dt, mf, sep = "\t", quote = FALSE)
  message("Fixed header: ", basename(mf))
}

## Patch SIF 'maf_fn' paths for timing analysis
library(data.table)

sif <- fread("/Users/jiajunzhan/Desktop/RA/BrCa/Code/PhylogicNDT/phylogic_inputs/Patient.sif")

# Assume maf_fn should point to ./myoutput/<sample_id>/<sample_id>.mut_ccfs.txt
sif[, maf_fn := file.path("./myoutput", sample_id, paste0(sample_id, ".mut_ccfs.txt"))]

# Write a new SIF (avoid overwriting the original)
fwrite(sif,
       "/Users/jiajunzhan/Desktop/RA/BrCa/Code/PhylogicNDT/phylogic_inputs/Patient_timing.sif",
       sep = "\t", quote = FALSE)

