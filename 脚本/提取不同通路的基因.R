suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(tidyr)
})

# =========================
# 参数窗口（你主要改这里）
# =========================

# 工作目录（KEGG-GO_db 总目录）
work_dir <- "F:/h何昊硕士毕业论文开展/data/KEGG-GO_db"

# 自动选择 OUT_osa_* 的策略： "latest" = 最近修改时间；"manual" = 手动指定 out_dir
out_dir_mode <- "latest"

# 若 out_dir_mode == "manual"，在此指定具体 OUT_osa_* 目录
out_dir_manual <- "F:/h何昊硕士毕业论文开展/data/KEGG-GO_db/OUT_osa_xxx"

# 输入文件名（位于 out_dir 下）
input_file <- "GO_enrich.tsv"

# 选择条目的方式：
# 1) "exact"：按 Description 精准匹配（忽略大小写）
# 2) "regex"：按 Description 正则匹配（忽略大小写）
match_mode <- "exact"

# 你要提取的通路/条目（对应 GO_enrich.tsv 的 Description 列）
# - exact 模式：写成向量，多条就多写几个
# - regex 模式：写成一个正则字符串（例如 "lignin.*process|phenylpropanoid"）
target_terms <- c("lignin biosynthetic process", "lignin metabolic process")

# 输出文件名前缀（便于区分不同通路）
output_prefix <- "LIGNIN"

# geneID 分隔符（clusterProfiler 常见是 "/"）
gene_sep <- "/"

# 每次运行创建独立输出文件夹（TRUE 强烈建议保持）
create_unique_run_dir <- TRUE

# =========================
# 主程序（一般不需要改）
# =========================

setwd(work_dir)

# 1) 确定 out_dir
if (out_dir_mode == "latest") {
  out_dirs <- list.dirs(work_dir, full.names = TRUE, recursive = FALSE)
  out_dirs <- out_dirs[grepl("^OUT_osa_", basename(out_dirs))]
  if (length(out_dirs) == 0) stop("未找到 OUT_osa_* 输出目录，请确认你已运行富集脚本。")
  out_dir <- out_dirs[which.max(file.info(out_dirs)$mtime)]
} else if (out_dir_mode == "manual") {
  out_dir <- out_dir_manual
  if (!dir.exists(out_dir)) stop("手动指定的 out_dir 不存在：", out_dir)
} else {
  stop("out_dir_mode 只能是 'latest' 或 'manual'")
}
cat("使用输出目录：", out_dir, "\n")

# 2) 读入 GO_enrich.tsv
in_path <- file.path(out_dir, input_file)
if (!file.exists(in_path)) stop("未找到输入文件：", in_path)

go <- read.delim(in_path, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
need_cols <- c("ID","Description","geneID")
if (!all(need_cols %in% colnames(go))) {
  stop("输入文件缺少必要列：", paste(setdiff(need_cols, colnames(go)), collapse = ", "))
}

# 3) 选择目标条目
if (match_mode == "exact") {
  target_lower <- tolower(target_terms)
  sel <- go %>% filter(tolower(Description) %in% target_lower)
} else if (match_mode == "regex") {
  # regex 模式下，target_terms 允许写多个，但会合并成 OR
  pattern <- paste(target_terms, collapse = "|")
  sel <- go %>% filter(str_detect(tolower(Description), tolower(pattern)))
} else {
  stop("match_mode 只能是 'exact' 或 'regex'")
}

if (nrow(sel) == 0) {
  stop("没有匹配到任何目标条目。请检查 target_terms / match_mode 以及 Description 字段内容。")
}

# 4) 每次运行建立独立输出目录，避免覆盖
run_dir <- out_dir
if (create_unique_run_dir) {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  run_dir <- file.path(out_dir, paste0("EXTRACT_", output_prefix, "_", ts))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
}
cat("本次输出目录：", run_dir, "\n")

# 5) 提取基因（geneID 按 gene_sep 拆分）
by_term <- sel %>%
  mutate(gene = str_split(geneID, fixed(gene_sep))) %>%
  unnest(gene) %>%
  mutate(gene = str_trim(gene)) %>%
  filter(gene != "") %>%
  transmute(Term = Description, Gene = gene) %>%
  distinct()

genes_union <- by_term %>%
  pull(Gene) %>%
  unique() %>%
  sort()

cat("匹配到的条目数：", nrow(sel), "\n")
print(sel[, intersect(c("ID","Description","p.adjust","Count"), colnames(sel)), drop = FALSE])
cat("合并去重后的基因数：", length(genes_union), "\n")

# 6) 输出文件
writeLines(genes_union, file.path(run_dir, paste0(output_prefix, "_genes_union.txt")))
write.table(by_term, file.path(run_dir, paste0(output_prefix, "_genes_byTerm.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sel, file.path(run_dir, paste0(output_prefix, "_enrich_rows.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("完成。输出文件已写入：", run_dir, "\n")