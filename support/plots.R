# All PCA calculations and result plotting are done in R.
read_snapshot <- function(folder, stage) {
  prefix <- file.path(folder, sprintf("stage_%02d", stage))
  individuals <- read.csv(paste0(prefix, "_individuals.csv"))
  con <- gzfile(paste0(prefix, "_dosage.csv.gz"), "rt")
  on.exit(close(con))
  dosage <- as.matrix(read.csv(con, header=FALSE, check.names=FALSE))
  stopifnot(nrow(dosage) == nrow(individuals))
  # Unscaled, centered genotype PCA, matching the original lab.
  # A sample-by-sample Gram matrix avoids an expensive full SNP SVD.
  centered <- sweep(dosage, 2, colMeans(dosage), "-")
  eig <- eigen(tcrossprod(centered), symmetric=TRUE)
  pcs <- matrix(0, nrow(dosage), 3)
  k <- min(3L, nrow(dosage), ncol(dosage))
  pcs[, seq_len(k)] <- sweep(eig$vectors[, seq_len(k), drop=FALSE], 2,
                            sqrt(pmax(0, eig$values[seq_len(k)])), "*")
  colors <- apply(pcs, 2, function(x) {
    span <- diff(range(x))
    if (span < 1e-10) rep(0.5, length(x)) else (x-min(x))/span
  })
  list(individuals=individuals, pcs=pcs,
       colors=rgb(colors[,1], colors[,2], colors[,3]),
       variance=if (sum(pmax(eig$values, 0)) > 0)
         100*head(c(pmax(eig$values, 0), 0, 0), 2)/sum(pmax(eig$values, 0)) else c(0,0))
}

plot_snapshot <- function(snapshot, layers, steps) {
  old <- par(mfrow=c(2,2), mar=c(3,3,3,1), oma=c(0,0,2,0))
  on.exit(par(old))
  dat <- snapshot$individuals
  plot(snapshot$pcs[,1], snapshot$pcs[,2], pch=21, bg=snapshot$colors,
       xlab=sprintf("PC1 (%.1f%%)", snapshot$variance[1]),
       ylab=sprintf("PC2 (%.1f%%)", snapshot$variance[2]), main="Genetic PCA")
  palette <- hcl.colors(100, "Viridis")
  paint_map <- function(layer, colors, title, background=gray.colors(100, start=0.97, end=0.55)) {
    nr <- nrow(layer); nc <- ncol(layer)
    plot(NA, xlim=c(0,nc), ylim=c(nr,0), asp=1, axes=FALSE,
         xlab="", ylab="", main=title, xaxs="i", yaxs="i")
    image(seq_len(nc)-0.5, seq_len(nr)-0.5, t(layer),
          col=background, zlim=c(0,1), add=TRUE)
    points(dat$x, dat$y, pch=21, bg=colors, cex=0.8)
    box()
  }
  paint_map(layers$capacity, snapshot$colors, "Genetic PCA on landscape")
  het <- pmin(1, pmax(0, dat$heterozygosity))
  paint_map(layers$connectivity, palette[1+floor(99*het)], "Individual heterozygosity")
  legend("bottomleft", legend=c("0","0.5","1"), fill=palette[c(1,50,100)],
         horiz=TRUE, bty="n", cex=0.7)
  phenotype <- pmin(1, pmax(0, dat$phenotype))
  paint_map(layers$environment, palette[1+floor(99*phenotype)], "Phenotype", background=palette)
  legend("bottomleft", legend=c("0","0.5","1"), fill=palette[c(1,50,100)],
         horiz=TRUE, bty="n", cex=0.7)
  mtext(sprintf("Step %d | %d individuals", steps, nrow(dat)), outer=TRUE)
}
