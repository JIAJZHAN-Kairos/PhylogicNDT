## ====== 基本路径（按需改） ======
library(data.table)
library(GenomicRanges)

maf_dir  <- "/Users/jiajunzhan/Desktop/RA/BrCa/Data/sample_maf_files"
seg_dir  <- "/Users/jiajunzhan/Desktop/RA/BrCa/Code/Evolution/data/cnv_with_purity"  # 原始 PURPLE 段（含 purity）
## 容器可见的输出根目录（你已经用 -v 映射到了 /work）
out_root <- "/Users/jiajunzhan/Desktop/RA/BrCa/Code/PhylogicNDT/phylogic_inputs"
out_maf  <- file.path(out_root, "maf_with_cn")
out_seg  <- file.path(out_root, "seg_timing")        # 规范后的 5 列 seg（Cluster 用 Start.bp/End.bp）
sif_file <- file.path(out_root, "Patient.sif")

dir.create(out_maf, showWarnings = FALSE, recursive = TRUE)
dir.create(out_seg, showWarnings = FALSE, recursive = TRUE)

## ====== 工具：把各种 seg 表头 -> 5 列（Chromosome/Start/End/A1.Seg.CN/A2.Seg.CN）======
to_seg5 <- function(df){
  nm <- tolower(names(df)); names(df) <- nm
  # 兼容 PURPLE 或已是 5 列的情况
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
    stop("seg 缺少必须列：Chromosome/Start/End/A1.Seg.CN/A2.Seg.CN 或 PURPLE 列")
  }
  # 标准化
  x$Chromosome <- sub("^chr","",x$Chromosome,ignore.case=TRUE)
  keep_chr <- c(as.character(1:22),"X","Y")
  x <- x[x$Chromosome %in% keep_chr, , drop=FALSE]
  x[["A1.Seg.CN"]] <- pmax(0, as.numeric(x[["A1.Seg.CN"]]))
  x[["A2.Seg.CN"]] <- pmax(0, as.numeric(x[["A2.Seg.CN"]]))
  ord <- as.integer(factor(x$Chromosome, levels=c(as.character(1:22),"X","Y")))
  x <- x[order(ord, x$Start, x$End), ]
  x
}

## ====== 工具：读 seg 并返回 GRanges + purity 值（若有）======
read_seg_as_gr <- function(seg_fp){
  seg_raw <- fread(seg_fp, sep="\t", header=TRUE, check.names=FALSE, showProgress=FALSE)
  seg5    <- to_seg5(seg_raw)
  # purity（若原表有一列 purity/ Purity）
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

## ====== 工具：把 maf 加上 local_cn_a1/a2（按位点落在哪个段）======
annotate_maf_with_cn <- function(maf_dt, gr_seg){
  # 兼容不同表头
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
    stop(paste("MAF 缺少列：", paste(setdiff(need, names(maf_dt)), collapse=", ")))
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
  # 默认落空位点用 1/1
  a1[is.na(a1)] <- 1L
  a2[is.na(a2)] <- 1L
  maf_dt[, local_cn_a1 := a1]
  maf_dt[, local_cn_a2 := a2]
  maf_dt
}

## ====== 扫描所有 MAF，逐个处理，写出 withCN + 规范 seg + SIF ======
maf_files <- list.files(maf_dir, pattern="\\.maf(\\.gz)?$", full.names=TRUE)
stopifnot(length(maf_files) > 0)

sif_rows <- list()

for (mf in maf_files) {
  sid <- sub("-somatic-PASS.*$","", basename(mf))
  sid <- sub("\\.maf(\\.gz)?$","", sid)
  
  ## 找对应 seg 文件（按你当前目录规则）
  cand <- c(
    file.path(seg_dir, paste0(sid, "_with_purity.tsv")),
    file.path(seg_dir, paste0(sid, ".tsv")),
    file.path(seg_dir, paste0(sid, ".seg.tsv"))
  )
  seg_fp <- cand[file.exists(cand)][1]
  if (is.na(seg_fp)) {
    message(sprintf("[跳过] 未找到 seg：%s", sid))
    next
  }
  
  # 读 seg -> GRanges + purity
  seg_obj <- tryCatch(read_seg_as_gr(seg_fp),
                      error=function(e){ message("[seg 读取失败] ", seg_fp, " : ", e$message); NULL })
  if (is.null(seg_obj)) next
  gr_seg <- seg_obj$gr
  seg5   <- seg_obj$seg5
  purity <- seg_obj$purity
  
  # 写出规范 seg（Start.bp/End.bp，供 Cluster 用）
  seg5_cluster <- seg5
  names(seg5_cluster)[names(seg5_cluster)=="Start"] <- "Start.bp"
  names(seg5_cluster)[names(seg5_cluster)=="End"]   <- "End.bp"
  seg_out <- file.path(out_seg, paste0(sid, ".seg.tsv"))
  fwrite(seg5_cluster, seg_out, sep="\t", quote=FALSE)
  message(sprintf("[SEG] 写出：%s", seg_out))
  
  # 读 MAF + 标准化列名 + 加 local_cn_* 并写出
  maf_dt <- fread(mf, sep="\t", header=TRUE, quote="", check.names=FALSE, showProgress=FALSE)
  maf_dt <- annotate_maf_with_cn(maf_dt, gr_seg)
  maf_out <- file.path(out_maf, paste0(sid, "-somatic-PASS.withCN.maf"))
  fwrite(maf_dt, maf_out, sep="\t", quote=FALSE)
  message(sprintf("[MAF] 写出：%s（行数=%d）", maf_out, nrow(maf_dt)))
  
  # 记录 SIF 行（容器内路径）
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

## 汇总写 SIF
if (length(sif_rows)) {
  sif <- data.table::rbindlist(sif_rows, use.names=TRUE, fill=TRUE)
  fwrite(sif, sif_file, sep="\t", quote=FALSE)
  message(sprintf("[SIF] 写出：%s（%d 个样本）", sif_file, nrow(sif)))
} else {
  message("[SIF] 没有样本可写（可能没匹配到 seg 或 maf）")
}



library(data.table)

out_maf <- "/Users/jiajunzhan/Desktop/RA/BrCa/Code/PhylogicNDT/phylogic_inputs/maf_with_cn"

maf_files <- list.files(out_maf, pattern="\\.maf(\\.gz)?$", full.names=TRUE)
stopifnot(length(maf_files) > 0)

for (mf in maf_files) {
  dt <- fread(mf, sep = "\t", header = TRUE, quote = "", check.names = FALSE, showProgress = FALSE)
  
  # 统一表头为 PhylogicNDT 期望的精确拼写
  # Start_Position -> Start_position（小写 p）
  if ("Start_Position" %in% names(dt) && !("Start_position" %in% names(dt))) {
    setnames(dt, "Start_Position", "Start_position")
  }
  if ("End_Position" %in% names(dt) && !("End_position" %in% names(dt))) {
    setnames(dt, "End_Position", "End_position")
  }
  if ("HGVSp_Short" %in% names(dt) && !("Protein_Change" %in% names(dt))) {
    setnames(dt, "HGVSp_Short", "Protein_Change")
  }
  # 其他关键列也校正一次（防止大小写/下划线差异）
  # 这些判断是“有等价列但没有标准列名”才改名，不会破坏已经正确的文件
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
  
  # 确保拷贝数列存在
  if (!("local_cn_a1" %in% names(dt))) dt[, local_cn_a1 := NA_integer_]
  if (!("local_cn_a2" %in% names(dt))) dt[, local_cn_a2 := NA_integer_]
  
  # 写回原文件（覆盖）
  fwrite(dt, mf, sep = "\t", quote = FALSE)
  message("Fixed header: ", basename(mf))
}


## 修改 SIF 文件中的 maf_fn 路径 for timing analysis
library(data.table)

sif <- fread("/Users/jiajunzhan/Desktop/RA/BrCa/Code/PhylogicNDT/phylogic_inputs/Patient.sif")

# 假设 maf_fn 要改成 ./myoutput/<sample_id>/<sample_id>.mut_ccfs.txt
sif[, maf_fn := file.path("./myoutput", sample_id, paste0(sample_id, ".mut_ccfs.txt"))]

# 写回一个新文件（避免覆盖原始）
fwrite(sif,
       "/Users/jiajunzhan/Desktop/RA/BrCa/Code/PhylogicNDT/phylogic_inputs/Patient_timing.sif",
       sep = "\t", quote = FALSE)
