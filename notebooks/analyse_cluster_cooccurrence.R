
#
# analyse_cluster_cooccurrence.R
#
# Quantifies whether the largest clusters in a genome x cluster presence/
# absence matrix tend to co-occur on the same genomes more than expected by
# chance. Complements analyse_cluster_phylogeny.py's tree-based view with a
# purely presence/absence-based statistical view.
#
# Outputs:
#   <prefix>_jaccard_matrix.tsv   pairwise Jaccard similarity, top N clusters
#   <prefix>_jaccard_heatmap.png  heatmap of the above
#   <prefix>_pairwise_fisher.tsv  pairwise Fisher's exact test, BH-corrected
#   <prefix>_summary.txt          observed vs. permuted mean Jaccard, p-value
#
# Example usage:
#   Rscript analyse_cluster_cooccurrence.R \
#     --presence_matrix results/cluster_spread/presence_absence_matrix.tsv \
#     --top_n 20 \
#     --permutations 1000 \
#     --output_prefix results/cluster_cooccurrence/top20

suppressMessages({
  library(optparse)
  library(pheatmap)
})

option_list <- list(
  make_option("--presence_matrix", type = "character",
              help = "presence_absence_matrix.tsv from analyse_cluster_spread.py (clusters x genomes)"),
  make_option("--top_n", type = "integer", default = 20,
              help = "Number of largest clusters (by genome count) to analyse [default %default]"),
  make_option("--permutations", type = "integer", default = 1000,
              help = "Permutations for the null co-occurrence test. Runtime scales with permutations x n_genomes, so increase gradually [default %default]"),
  make_option("--output_prefix", type = "character",
              help = "Prefix for output files, e.g. results/cluster_cooccurrence/top20")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$presence_matrix) || is.null(opt$output_prefix)) {
  stop("--presence_matrix and --output_prefix are required. Run with --help for usage.")
}

# --- Load matrix: rows = clusters, columns = genomes, values 0/1 ---
mat_raw <- read.delim(opt$presence_matrix, row.names = 1, check.names = FALSE)
M_full <- as.matrix(mat_raw)
storage.mode(M_full) <- "numeric"

cat(sprintf("Loaded matrix: %d clusters x %d genomes\n", nrow(M_full), ncol(M_full)))

# --- Select top N clusters by genome count ---
genome_counts <- rowSums(M_full)
top_clusters <- names(sort(genome_counts, decreasing = TRUE))[1:min(opt$top_n, nrow(M_full))]
M <- M_full[top_clusters, , drop = FALSE]
n <- nrow(M)

cat(sprintf("Analysing top %d clusters (genome counts: %d-%d)\n",
            n, min(rowSums(M)), max(rowSums(M))))

# --- Vectorized pairwise Jaccard similarity via matrix multiplication ---
# inter[i,j] = number of genomes where both cluster i and cluster j are present
# union[i,j] = |A| + |B| - inter, by inclusion-exclusion
compute_jaccard_matrix <- function(mat) {
  inter <- mat %*% t(mat)
  row_sums <- rowSums(mat)
  union_mat <- outer(row_sums, row_sums, "+") - inter
  jmat <- inter / union_mat
  diag(jmat) <- NA
  jmat
}

jmat <- compute_jaccard_matrix(M)
rownames(jmat) <- top_clusters
colnames(jmat) <- top_clusters

write.table(round(jmat, 3), file = paste0(opt$output_prefix, "_jaccard_matrix.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

# --- Pairwise Fisher's exact test for co-occurrence ---
pair_idx <- t(combn(n, 2))
fisher_results <- data.frame(
  cluster_a = character(nrow(pair_idx)), cluster_b = character(nrow(pair_idx)),
  odds_ratio = numeric(nrow(pair_idx)), p_value = numeric(nrow(pair_idx)),
  jaccard = numeric(nrow(pair_idx)), stringsAsFactors = FALSE
)

for (k in seq_len(nrow(pair_idx))) {
  i <- pair_idx[k, 1]; j <- pair_idx[k, 2]
  a <- factor(M[i, ] == 1, levels = c(FALSE, TRUE))
  b <- factor(M[j, ] == 1, levels = c(FALSE, TRUE))
  tab <- table(a, b)
  ft <- fisher.test(tab)
  fisher_results$cluster_a[k] <- top_clusters[i]
  fisher_results$cluster_b[k] <- top_clusters[j]
  fisher_results$odds_ratio[k] <- unname(ft$estimate)
  fisher_results$p_value[k] <- ft$p.value
  fisher_results$jaccard[k] <- jmat[i, j]
}

fisher_results$p_adj <- p.adjust(fisher_results$p_value, method = "BH")
fisher_results <- fisher_results[order(fisher_results$p_adj), ]

write.table(fisher_results, file = paste0(opt$output_prefix, "_pairwise_fisher.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Global permutation test ---
# Null model: each cluster's presence is randomly reshuffled across genomes,
# preserving that cluster's own prevalence (row sum), but destroying any
# genome-level association between clusters. Tests whether the *observed*
# mean pairwise Jaccard among the top clusters is higher than this null.
observed_mean_jaccard <- mean(jmat[upper.tri(jmat)], na.rm = TRUE)

permute_rows <- function(mat) {
  t(apply(mat, 1, sample))
}

set.seed(1)
null_means <- numeric(opt$permutations)
for (p in seq_len(opt$permutations)) {
  M_perm <- permute_rows(M)
  jmat_perm <- compute_jaccard_matrix(M_perm)
  null_means[p] <- mean(jmat_perm[upper.tri(jmat_perm)], na.rm = TRUE)
}

perm_p <- (sum(null_means >= observed_mean_jaccard) + 1) / (opt$permutations + 1)

summary_lines <- c(
  sprintf("Clusters analysed: %d (top %d by genome count)", n, opt$top_n),
  sprintf("Genomes in matrix: %d", ncol(M)),
  sprintf("Cluster genome-count range: %d-%d", min(rowSums(M)), max(rowSums(M))),
  "",
  sprintf("Observed mean pairwise Jaccard similarity: %.4f", observed_mean_jaccard),
  sprintf("Null mean Jaccard across %d permutations: %.4f (sd %.4f)",
          opt$permutations, mean(null_means), sd(null_means)),
  sprintf("Permutation p-value (observed >= null): %.4g", perm_p),
  "",
  "Interpretation: a low p-value means these top clusters co-occur on the same",
  "genomes more than expected if each cluster were independently scattered at",
  "its own observed prevalence - i.e. the visual pattern in the iTOL tree is",
  "unlikely to be a coincidence of shared prevalence alone. Check",
  "*_pairwise_fisher.tsv for which specific cluster pairs drive this."
)
writeLines(summary_lines, con = paste0(opt$output_prefix, "_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

# --- Heatmap ---
pheatmap(jmat,
         main = sprintf("Pairwise Jaccard similarity, top %d clusters", n),
         filename = paste0(opt$output_prefix, "_jaccard_heatmap.png"),
         display_numbers = TRUE, number_format = "%.2f",
         na_col = "grey90")

cat(sprintf("\nDone. Wrote %s_jaccard_matrix.tsv, %s_pairwise_fisher.tsv, %s_summary.txt, %s_jaccard_heatmap.png\n",
            opt$output_prefix, opt$output_prefix, opt$output_prefix, opt$output_prefix))
