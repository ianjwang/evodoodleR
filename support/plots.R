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

# Match the editor's linear RGB interpolation; values always use the full 0–1 scale.
landscape_colors <- function(values, anchors) {
  channels <- grDevices::colorRamp(anchors, space="rgb")(pmin(1,pmax(0,values)))
  grDevices::rgb(round(channels), maxColorValue=255)
}

plot_snapshot <- function(snapshot, layers, steps) {
  old <- par(mfrow=c(2,2), mar=c(5,3,3,1), oma=c(0,0,2,0))
  on.exit(par(old))
  dat <- snapshot$individuals
  plot(snapshot$pcs[,1], snapshot$pcs[,2], pch=21, bg=snapshot$colors,
       xlab=sprintf("PC1 (%.1f%%)", snapshot$variance[1]),
       ylab=sprintf("PC2 (%.1f%%)", snapshot$variance[2]), main="Genetic PCA")
  environment_anchors <- c("#542788","#998ec3","#f7f7f7","#f1a340","#b35806")
  phenotype_palette <- landscape_colors(seq(0,1,length.out=1025),environment_anchors)
  het_anchors <- hcl.colors(256,"Viridis")
  het_palette <- landscape_colors(seq(0,1,length.out=1025),het_anchors)
  paint_map <- function(layer, colors, title, background=gray.colors(100, start=0.97, end=0.55)) {
    nr <- nrow(layer); nc <- ncol(layer)
    plot(NA, xlim=c(0,nc), ylim=c(nr,0), asp=1, axes=FALSE,
         xlab="", ylab="", main=title, xaxs="i", yaxs="i")
    image(seq_len(nc)-0.5, seq_len(nr)-0.5, t(layer),
          col=background, zlim=c(0,1), add=TRUE)
    points(dat$x, dat$y, pch=21, bg=colors, cex=0.8)
    box()
  }
  # A continuous legend below the map, with absolute values (never rescaled per run).
  colorbar <- function(palette, label) {
    edges <- grconvertX(seq(0.15,0.85,length.out=length(palette)+1),"npc","user")
    ys <- grconvertY(c(-0.12,-0.06),"npc","user")
    rect(head(edges,-1),ys[1],tail(edges,-1),ys[2],col=palette,border=NA,xpd=NA)
    rect(edges[1],ys[1],tail(edges,1),ys[2],border="#555555",xpd=NA)
    text(grconvertX(c(0.15,0.5,0.85),"npc","user"),
         grconvertY(-0.19,"npc","user"),labels=c("0","0.5","1"),cex=0.7,xpd=NA)
    text(grconvertX(0.5,"npc","user"),grconvertY(-0.28,"npc","user"),
         labels=label,cex=0.7,xpd=NA)
  }
  paint_map(layers$capacity, snapshot$colors, "Genetic PCA on landscape")
  paint_map(layers$connectivity,landscape_colors(dat$heterozygosity,het_anchors),
            "Individual heterozygosity")
  colorbar(het_palette,"Observed heterozygosity")
  paint_map(layers$environment,landscape_colors(dat$phenotype,environment_anchors),
            "Phenotype and environment",background=phenotype_palette)
  colorbar(phenotype_palette,"Phenotype / environmental optimum")
  mtext(sprintf("Step %d | %d individuals", steps, nrow(dat)), outer=TRUE)
}
