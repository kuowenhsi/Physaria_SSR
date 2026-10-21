library(ape)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg[1])
  script_path <- gsub("~\\+~", " ", script_path, fixed = FALSE)
  dirname(normalizePath(script_path, mustWork = FALSE))
} else {
  getwd()
}

repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
tree_dir <- file.path(repo_root, "figure", "Physaria_Neis_tree")
input_tree <- file.path(tree_dir, "Physaria_Neis_neighbor_joining_tree_relabeled.nwk")
output_plot <- file.path(
  tree_dir,
  "Physaria_Neis_neighbor_joining_tree_relabeled_without_EX_01.png"
)

nj_tree <- read.tree(input_tree)

if (!"EX_01" %in% nj_tree$tip.label) {
  stop("The relabeled tree does not contain an EX_01 tip.")
}

nj_tree <- drop.tip(nj_tree, "EX_01")

png(
  filename = output_plot,
  width = 2,
  height = 3.5,
  units = "in",
  res = 300
)
par(mar = c(0.6, 0.6, 1.2, 0.6))
plot(
  nj_tree,
  type = "phylogram",
  cex = 0.4,
  font = 1,
  no.margin = FALSE,
  direction = "downwards",
  label.offset = 0.06
)
add.scale.bar(length = 0.2, lwd = 1)
dev.off()
