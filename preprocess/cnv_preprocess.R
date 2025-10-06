## ========== 路径 ==========
setwd("~/Desktop/RA/BrCa/Code/PhylogicNDT/phylogic_inputs")

cnv_dir      <- "~/Desktop/RA/BrCa/Data/cnv_segment"      # 原始 CNV 段文件目录（.tsv / .tsv.gz）
summary_csv  <- "~/Desktop/RA/BrCa/Results/brca-garvan-wgts_annotated_summary.csv"  # 含 purity 的汇总表
out_dir_allelic <- "seg_alleliccapseg"   # 只导出 alleliccapseg

dir.create(out_dir_allelic, showWarnings = FALSE, recursive = TRUE)

## ========== 读取 purity 汇总 ==========
summary_df <- read.csv(summary_csv, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
if (!"Purity_value" %in% names(summary_df)) {
  summary_df$Purity_value <- suppressWarnings(as.numeric(sub(" .*", "", summary_df$Purity)))
}
id_col     <- "WGS-tumor-lib"   # 你的样本ID列（注意和表头一致）
purity_col <- "Purity_value"
stopifnot(id_col %in% names(summary_df), purity_col %in% names(summary_df))

## ========== 工具：标准化到 alleliccapseg ==========
# 输出列：
#   Chromosome, Start.bp, End.bp, mu.minor, sigma.minor, mu.major, sigma.major
normalize_to_allelic <- function(df, default_sigma = 0.20) {
  idx <- function(name) which(tolower(names(df)) == tolower(name))
  has <- function(name) length(idx(name)) == 1
  
  # chr/start/end
  chr_col <- NULL; for (c in c("chromosome", "chr")) if (has(c)) { chr_col <- idx(c); break }
  if (is.null(chr_col)) stop("找不到 Chromosome 列（chromosome/chr）")
  start_col <- NULL; for (k in c("start", "start.bp")) if (has(k)) { start_col <- idx(k); break }
  if (is.null(start_col)) stop("找不到起始位点列（start/start.bp）")
  end_col <- NULL; for (k in c("end", "end.bp")) if (has(k)) { end_col <- idx(k); break }
  if (is.null(end_col)) stop("找不到结束位点列（end/end.bp）")
  
  # A1/A2：modal.a1/modal.a2 或 PURPLE: major/minorAlleleCopyNumber
  if (has("modal.a1") && has("modal.a2")) {
    a1_col <- idx("modal.a1"); a2_col <- idx("modal.a2")
  } else if (has("majorallelecopynumber") && has("minorallelecopynumber")) {
    a1_col <- idx("majorallelecopynumber"); a2_col <- idx("minorallelecopynumber")
  } else {
    stop("找不到等位拷贝数（modal.a1/modal.a2 或 major/minorAlleleCopyNumber）")
  }
  
  # deviation（可选）
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
  
  # 标准化 chr & 过滤
  out$Chromosome <- sub("^chr", "", out$Chromosome, ignore.case = TRUE)
  keep_chr <- c(as.character(1:22), "X", "Y")
  out <- out[out$Chromosome %in% keep_chr, , drop = FALSE]
  
  # 确保 major >= minor
  swap <- out$mu.major < out$mu.minor
  if (any(swap)) {
    tmp <- out$mu.major[swap]; out$mu.major[swap] <- out$mu.minor[swap]; out$mu.minor[swap] <- tmp
  }
  
  # sigma: 有 deviation 就用；否则给默认常数
  if (!is.na(maj_dev) && !is.na(min_dev)) {
    out$sigma.major <- suppressWarnings(as.numeric(df[[maj_dev]]))
    out$sigma.minor <- suppressWarnings(as.numeric(df[[min_dev]]))
    out$sigma.major[!is.finite(out$sigma.major)] <- default_sigma
    out$sigma.minor[!is.finite(out$sigma.minor)] <- default_sigma
  } else {
    out$sigma.major <- default_sigma
    out$sigma.minor <- default_sigma
  }
  
  # 排序
  chr_order <- function(v) as.integer(factor(v, levels = c(as.character(1:22), "X", "Y")))
  out <- out[order(chr_order(out$Chromosome), out$Start.bp, out$End.bp), ]
  
  # 列顺序符合 PhylogicNDT 预期
  out <- out[, c("Chromosome","Start.bp","End.bp","mu.minor","sigma.minor","mu.major","sigma.major")]
  out
}

## ========== 批量处理 ==========
cnv_files <- list.files(cnv_dir, pattern="\\.tsv(\\.gz)?$", full.names=TRUE)
if (length(cnv_files) == 0L) stop(sprintf("目录里没有 .tsv/.tsv.gz：%s", cnv_dir))

for (fp in cnv_files) {
  sid <- sub("\\..*$", "", basename(fp))  # 去后缀当样本ID
  if (!(sid %in% summary_df[[id_col]])) {
    message(sprintf("[跳过] %s —— 汇总表中找不到样本ID：%s", basename(fp), sid)); next
  }
  
  seg_raw <- tryCatch({
    if (grepl("\\.gz$", fp, ignore.case = TRUE)) {
      read.delim(gzfile(fp), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      read.delim(fp, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    }
  },
  error = function(e) { message("[读取失败] ", fp, " : ", e$message); return(NULL) }
  )
  if (is.null(seg_raw)) next
  
  seg_allelic <- tryCatch(normalize_to_allelic(seg_raw, default_sigma = 0.20),
                          error = function(e){ message("[格式错误] ", fp, " : ", e$message); return(NULL) })
  if (is.null(seg_allelic)) next
  
  out1 <- file.path(out_dir_allelic, paste0(sid, ".seg.tsv"))
  write.table(seg_allelic, out1, sep="\t", quote=FALSE, row.names=FALSE)
  
  message(sprintf("[OK] %-18s | allelic=%s", sid, basename(out1)))
}

message("全部完成。")