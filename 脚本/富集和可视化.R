suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(jsonlite)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# =========================================================
# 0) 只需要改这里：工作目录（你已给出）
# =========================================================
work_dir <- "F:/h何昊硕士毕业论文开展/data/KEGG-GO_db"
setwd(work_dir)

# 输入文件（都在 work_dir 下）
emapper_file  <- "ROC22.egg.emapper.annotations"
json_file     <- "osa00001.json"     # <-- 关键：水稻 KEGG BRITE JSON
gene_file     <- "DEG_all_union.in_universe.txt"
universe_file <- "universe_ids.txt"

# 展示多少条
kegg_show_n <- 12
go_show_n   <- 12

# =========================================================
# 1) 输出目录：每次运行独立文件夹，防止覆盖
# =========================================================
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(work_dir, paste0("OUT_osa_", run_id))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("work_dir:", work_dir, "\n")
cat("out_dir :", out_dir, "\n")

# 检查文件存在
stopifnot(file.exists(emapper_file),
          file.exists(json_file),
          file.exists(gene_file),
          file.exists(universe_file))

# =========================================================
# 2) 读取 gene list 与 universe
# =========================================================
genes_raw <- read.table(gene_file, header=FALSE, sep="", quote="", comment.char="",
                        stringsAsFactors=FALSE)[,1] %>%
  as.character() %>% trimws() %>% unique()
genes_raw <- genes_raw[!is.na(genes_raw) & genes_raw != ""]
cat("输入基因数(genes_raw):", length(genes_raw), "\n")
if (length(genes_raw) == 0) stop("targets_intersect_DEG.txt 为空或读入失败")

universe_ids <- read.table(universe_file, header=FALSE, sep="", quote="", comment.char="",
                           stringsAsFactors=FALSE)[,1] %>%
  as.character() %>% trimws() %>% unique()
universe_ids <- universe_ids[!is.na(universe_ids) & universe_ids != ""]
cat("背景基因数(universe_ids):", length(universe_ids), "\n")
if (length(universe_ids) == 0) stop("universe_ids.txt 为空或读入失败")

# =========================================================
# 3) 稳健读取 eggNOG emapper（兼容 # 表头）
# =========================================================
read_emapper_annotations <- function(path) {
  x <- readLines(path, n = 5000, warn = FALSE)
  hdr_i <- which(grepl("^#?query\\b", x))[1]
  if (is.na(hdr_i)) stop("未在文件前5000行找到表头(query)。请检查 emapper 文件格式。")
  
  header_line <- sub("^#+", "", x[hdr_i])
  cols <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  
  df <- read.delim(path, sep="\t", header=FALSE, quote="",
                   comment.char="", stringsAsFactors=FALSE,
                   check.names=FALSE, skip=hdr_i)
  colnames(df) <- cols
  df
}

cat("读取 EggNOG:", emapper_file, "\n")
emapper <- read_emapper_annotations(emapper_file)

# 清理空值
emapper[emapper == ""] <- NA
emapper[emapper == "-"] <- NA

need_cols <- c("query", "Preferred_name", "GOs", "KEGG_ko")
if (!all(need_cols %in% colnames(emapper))) {
  stop("缺少关键列。当前列名为：\n", paste(colnames(emapper), collapse=", "))
}

# gene_info / gene2go / gene2ko
gene_info_all <- emapper %>%
  dplyr::select(GID = query, GENENAME = Preferred_name) %>%
  dplyr::mutate(GENENAME = ifelse(is.na(GENENAME), GID, GENENAME)) %>%
  dplyr::distinct()

gene2go_all <- emapper %>%
  dplyr::select(GID = query, GO = GOs) %>%
  dplyr::filter(!is.na(GO)) %>%
  tidyr::separate_rows(GO, sep = ",") %>%
  dplyr::mutate(GO = stringr::str_extract(GO, "GO:[0-9]+"),
                EVIDENCE = "IEA") %>%
  dplyr::filter(!is.na(GO)) %>%
  dplyr::distinct()

gene2ko_all <- emapper %>%
  dplyr::select(GID = query, Ko = KEGG_ko) %>%
  dplyr::filter(!is.na(Ko)) %>%
  tidyr::separate_rows(Ko, sep = ",") %>%
  dplyr::mutate(Ko = stringr::str_remove(Ko, "^ko:"),
                Ko = stringr::str_extract(Ko, "^K[0-9]+")) %>%
  dplyr::filter(!is.na(Ko)) %>%
  dplyr::distinct()

cat("注释概况：基因数", nrow(gene_info_all),
    "；带GO基因", dplyr::n_distinct(gene2go_all$GID),
    "；带KO基因", dplyr::n_distinct(gene2ko_all$GID), "\n")

# universe 与 genes 必须能映射到注释
universe_ids <- intersect(universe_ids, gene_info_all$GID)
cat("与注释匹配后的 universe 基因数:", length(universe_ids), "\n")
if (length(universe_ids) == 0) stop("universe_ids 与 emapper$query 交集为0，请检查ID体系")

genes_in_anno <- intersect(genes_raw, gene_info_all$GID)
cat("genes_raw 中能匹配注释的基因数:", length(genes_in_anno), "\n")
if (length(genes_in_anno) == 0) stop("targets_intersect_DEG 与 emapper$query 完全不匹配")

# =========================================================
# 4) 解析本地 osa00001.json -> pathway2name / ko2pathway（不做白名单）
# =========================================================
get_kegg_info_local <- function(json_file) {
  kegg <- jsonlite::fromJSON(json_file, simplifyVector = FALSE)
  if (is.null(kegg[["children"]])) stop("JSON 结构异常：无 children")
  
  path_list <- list()
  ko_list <- list()
  p_count <- 0
  
  # 适配 BRITE 树：A -> B -> C -> pathway 节点
  for (a in seq_along(kegg[["children"]])) {
    B <- kegg[["children"]][[a]][["children"]]
    if (is.null(B)) next
    for (b in seq_along(B)) {
      C <- B[[b]][["children"]]
      if (is.null(C)) next
      for (c in seq_along(C)) {
        node <- C[[c]]
        info <- node[["name"]]
        # 关键：支持 osa00010 / ko00010 / taes00010 等
        pid <- stringr::str_extract(info, "[a-z]{2,4}[0-9]{5}")
        pname <- info %>%
          stringr::str_remove(" \\[PATH:[a-z]{2,4}[0-9]{5}\\]") %>%
          stringr::str_remove("^[0-9]{5} ")
        
        if (!is.na(pid)) {
          p_count <- p_count + 1
          path_list[[p_count]] <- data.frame(Pathway=pid, Name=pname, stringsAsFactors=FALSE)
          
          kids <- node[["children"]]
          if (!is.null(kids)) {
            raw <- vapply(kids, function(x) x[["name"]], FUN.VALUE=character(1))
            kos <- stringr::str_extract(raw, "K[0-9]+")
            kos <- kos[!is.na(kos)]
            if (length(kos) > 0) {
              ko_list[[length(ko_list)+1]] <- data.frame(Ko=kos, Pathway=pid, stringsAsFactors=FALSE)
            }
          }
        }
      }
    }
  }
  
  p2n <- dplyr::bind_rows(path_list) %>% dplyr::distinct()
  k2p <- if (length(ko_list) == 0) {
    data.frame(Ko=character(), Pathway=character(), stringsAsFactors=FALSE)
  } else {
    dplyr::bind_rows(ko_list) %>% dplyr::distinct()
  }
  
  if (nrow(p2n) == 0) stop("从 JSON 中未解析到任何 Pathway。请确认 json_file 是 osa00001.json 且文件未损坏。")
  list(p2n=p2n, k2p=k2p)
}

cat("解析 KEGG JSON:", json_file, "\n")
kegg_data <- get_kegg_info_local(json_file)

pathway2name <- kegg_data$p2n %>%
  dplyr::mutate(Name = gsub(" \\[BR:[a-z]{2,4}[0-9]{5}\\]", "", Name)) %>%
  dplyr::distinct()

ko2pathway <- kegg_data$k2p %>% dplyr::distinct()

cat("解析结果：通路数", nrow(pathway2name), "；KO-Pathway数", nrow(ko2pathway), "\n")
if (nrow(ko2pathway) == 0) stop("KO-Pathway 映射为空：请检查 JSON 是否包含 KO 列表（children）。")

# TERM2GENE（全量：不白名单）
gene2pathway_all <- gene2ko_all %>%
  dplyr::inner_join(ko2pathway, by="Ko") %>%
  dplyr::select(GID, Pathway) %>%
  dplyr::distinct()

pathway2gene_all <- gene2pathway_all %>%
  dplyr::transmute(Pathway, GID) %>%
  dplyr::distinct()

write.table(pathway2name, file = file.path(out_dir, "ALL_pathway2name_osa.txt"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(pathway2gene_all, file = file.path(out_dir, "ALL_pathway2gene_osa.txt"),
            sep="\t", quote=FALSE, row.names=FALSE)

# GO TERM2GENE（全量）
go_term2gene_all <- gene2go_all %>% dplyr::transmute(Term=GO, Gene=GID) %>% dplyr::distinct()
write.table(go_term2gene_all, file = file.path(out_dir, "ALL_GO_term2gene.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# =========================================================
# 5) 富集（KEGG/GO）与作图
# =========================================================
# universe 必须与 TERM2GENE 匹配
universe_kegg <- intersect(universe_ids, unique(pathway2gene_all$GID))
genes_kegg <- intersect(genes_in_anno, universe_kegg)
cat("KEGG：universe=", length(universe_kegg), "；genes=", length(genes_kegg), "\n")
if (length(genes_kegg) == 0) stop("KEGG：genes 与 universe 交集为0（可能是 KO 注释不足或ID不一致）")

ekk <- clusterProfiler::enricher(
  gene = genes_kegg,
  universe = universe_kegg,
  TERM2GENE = pathway2gene_all,
  TERM2NAME = pathway2name,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 1,
  maxGSSize = 5000
)

ekk_df <- as.data.frame(ekk)
write.table(ekk_df, file = file.path(out_dir, "KEGG_enrich_osa.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

if (!is.null(ekk) && nrow(ekk_df) > 0) {
  p1 <- enrichplot::dotplot(ekk, showCategory=kegg_show_n, label_format=45) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size=9))
  ggplot2::ggsave(file.path(out_dir, "KEGG_dotplot_osa.pdf"), p1, width=10, height=8)
  
  dfp <- ekk_df %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = kegg_show_n) %>%
    dplyr::mutate(Description = stringr::str_wrap(Description, width = 45),
                  Description = factor(Description, levels = rev(unique(Description))),
                  neglog10_padj = -log10(p.adjust))
  
  p2 <- ggplot2::ggplot(dfp, ggplot2::aes(x = neglog10_padj, y = Description)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = "-log10(p.adjust)", y = NULL)
  
  ggplot2::ggsave(file.path(out_dir, "KEGG_barplot_osa.pdf"), p2, width=10, height=8)
} else {
  cat("KEGG：无显著富集条目，跳过绘图。\n")
}

# GO：背景用 universe_ids 与 GO_term2gene 交集
universe_go <- intersect(universe_ids, unique(go_term2gene_all$Gene))
genes_go <- intersect(genes_in_anno, universe_go)
cat("GO：universe=", length(universe_go), "；genes=", length(genes_go), "\n")
if (length(genes_go) == 0) stop("GO：genes 与 universe 交集为0，请检查ID")

# 给 GO 补 TERM 名称（避免 select 冲突：不 attach AnnotationDbi）
go_term2name <- NULL
if (requireNamespace("GO.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
  go_ids <- unique(go_term2gene_all$Term)
  go_map <- tryCatch({
    AnnotationDbi::select(GO.db::GO.db, keys = go_ids, columns = "TERM", keytype = "GOID")
  }, error = function(e) NULL)
  
  if (!is.null(go_map) && nrow(go_map) > 0) {
    go_term2name <- go_map %>%
      dplyr::select(Term = GOID, Name = TERM) %>%
      dplyr::distinct()
  }
}

ego <- clusterProfiler::enricher(
  gene = genes_go,
  universe = universe_go,
  TERM2GENE = go_term2gene_all,
  TERM2NAME = go_term2name,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

ego_df <- as.data.frame(ego)
write.table(ego_df, file = file.path(out_dir, "GO_enrich.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

if (!is.null(ego) && nrow(ego_df) > 0) {
  p3 <- enrichplot::dotplot(ego, showCategory=go_show_n, label_format=45) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size=9))
  ggplot2::ggsave(file.path(out_dir, "GO_dotplot.pdf"), p3, width=10, height=9)
  
  dfp <- ego_df %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = go_show_n) %>%
    dplyr::mutate(Description = stringr::str_wrap(Description, width = 45),
                  Description = factor(Description, levels = rev(unique(Description))),
                  neglog10_padj = -log10(p.adjust))
  
  p4 <- ggplot2::ggplot(dfp, ggplot2::aes(x = neglog10_padj, y = Description)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = "-log10(p.adjust)", y = NULL)
  
  ggplot2::ggsave(file.path(out_dir, "GO_barplot.pdf"), p4, width=10, height=9)
} else {
  cat("GO：无显著富集条目，跳过绘图。\n")
}

# =========================================================
# 6) 输出子集明细与对象（便于复现）
# =========================================================
write.table(gene_info_all %>% dplyr::filter(GID %in% genes_in_anno),
            file = file.path(out_dir, "SUBSET_gene_info.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(gene2ko_all %>% dplyr::filter(GID %in% genes_in_anno),
            file = file.path(out_dir, "SUBSET_gene2ko.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(gene2go_all %>% dplyr::filter(GID %in% genes_in_anno),
            file = file.path(out_dir, "SUBSET_gene2go.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(gene2pathway_all %>% dplyr::filter(GID %in% genes_in_anno),
            file = file.path(out_dir, "SUBSET_gene2pathway_osa.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

saveRDS(list(
  pathway2name = pathway2name,
  pathway2gene = pathway2gene_all,
  ekk = ekk,
  ego = ego
), file = file.path(out_dir, "objects_osa.rds"))

cat("\n================= 完成 =================\n")
cat("输出目录:", out_dir, "\n")
cat("KEGG 背景: ALL_pathway2gene_osa.txt / ALL_pathway2name_osa.txt\n")
cat("KEGG 富集: KEGG_enrich_osa.tsv + PDF\n")
cat("GO   富集: GO_enrich.tsv + PDF\n")
cat("对象: objects_osa.rds\n")


