DIFplot = function(d, 
  what=c('distances','statistics','pvalues','significance'),
  pam=c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
     "fdr"), cluster = TRUE, alpha=0.05, ...) {
  what = match.arg(what)
  pam = match.arg(pam)
  alpha = alpha/2
  dist = abs(d$Delta_R)
  o = hclust(as.dist(dist))$order
  lbl = d$items$item_id
  stat = abs(d$DIF_pair)
  outl = lbl[o]
  if (cluster) {stat=stat[o,o]; lbl=lbl[o]}
  pval = 1 - pnorm(stat)
  u0 = pval[lower.tri(pval)]
  u1 = p.adjust(u0, method=pam)
  pval[lower.tri(pval)] = u1

  if(what=='distances') {
    if (cluster) {dist = dist[o,o]}
    rownames(dist) = colnames(dist) = d$items$item_id
    diag(dist) = NA
    pheatmap::pheatmap(dist, main='PDIF: raw differences', cluster_rows=FALSE, cluster_cols=FALSE)
  }
  if(what=='statistics') {
    rownames(stat) = colnames(stat) = lbl
    diag(stat) = NA
    pheatmap::pheatmap(stat, main='PDIF: standardized differences', cluster_rows=FALSE, cluster_cols=FALSE)
  } 
  if(what=='pvalues') {
    rownames(pval) = colnames(pval) = lbl
    diag(pval) = NA
    ttl = 'PDIF: p-values'
    if (pam != 'none') ttl=paste0(ttl, ' (below diagonal adjusted by ',pam,')')
    pheatmap::pheatmap(pval,
      main = ttl,
      cluster_rows=FALSE, cluster_cols=FALSE,
      color = colorRampPalette(RColorBrewer::brewer.pal(n=7, name="RdYlBu"))(100))
  }
  if(what=='significance') {
    ttl = paste0('PDIF: significance at alpha=',2*alpha)
    if (pam != 'none') ttl=paste0(ttl, ' (b/d adjusted by ',pam,')')
    v = (0 + (pval<alpha))
    rownames(v) = colnames(v) = lbl
    diag(v) = NA
    pheatmap::pheatmap(v, main=ttl, cluster_rows=FALSE, cluster_cols=FALSE, legend=FALSE)
  }
  return(outl)
}



