##############################################################
# PCA Visualization in R
# ------------------------------------------------------
# Description:
# This script  plots the first two PCs, coloring samples by group (e.g., SR vs ST).
##############################################################

setwd("")

#Load Covariance Matrix
cov <- as.matrix(read.table("angsd_output.cov"))
e <- eigen(cov)

#Load Sample Names
# Each row in "names.txt" should correspond to one individual in the covariance matrix.
names <- read.table("names.txt")
rownames(cov) <- names$V1  # Assign sample names to the covariance matrix rows

# Define Sample Group Colors 
#samples starting with "SR" = red, "ST" = blue.
group_colors <- c("SR" = "red", "ST" = "blue")

#scale axes slightly beyond the data range
x_limits <- range(e$vectors[, 1]) * 1.3
y_limits <- range(e$vectors[, 2]) * 1.3

# This step classifies each sample into SR or ST based on its name.
sample_groups <- ifelse(grepl("^SR", rownames(cov)), "SR", "ST")

#PCA Scatter Plot
plot(
  e$vectors[, 1], e$vectors[, 2],
  col = group_colors[sample_groups],
  pch = 16,
  cex = 0.7,
  xlab = paste0("PC1 (", round(100 * e$values[1] / sum(e$values), 1), "%)"),
  ylab = paste0("PC2 (", round(100 * e$values[2] / sum(e$values), 1), "%)"),
  main = "PCA (ANGSD Output)",
  xlim = x_limits,
  ylim = y_limits
)

# Add Legend 
legend(
  "topright",
  legend = names(group_colors),
  col = group_colors,
  pch = 16,
  bty = "n",         
  pt.cex = 1.2,
  cex = 0.9
)
