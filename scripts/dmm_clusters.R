# DMM refers to https://bioconductor.org/packages/release/bioc/vignettes/DirichletMultinomial/inst/doc/DirichletMultinomial.pdf

library(DirichletMultinomial)
library(lattice)
library(xtable)
library(parallel)
requireNamespace("ComplexHeatmap")


# load metagenome abundance profiling, row:sample, col:MAG
tab <- read.csv(
    '../data/relative_abundance_genomes.csv',
    header=T,
    row.names=1,
    check.names=F
)
gtdb <- read.csv(
    '../data/genome_taxonomy_annotation.csv',
    header=T,
    row.names=1,
    check.names=F
)
meta <- read.csv(
    '../data/metadata.csv',
    header=T,
    row.names=1,
    check.names=F
) %>% .[rownames(tab), ]

# convert relative abundance to count data due to DMN only works on integer
# multiply by 1M, although it doesn't match to inStrain quick profile, which calculating rel.abund = (genome length * genome coverage) / total metagenome read base pairs
count_data <- mapply(tab, FUN=function(x) round(x*1e+6))

#================
# 1. check distribution of reads from each MAG, on a log scale
#================
cnts <- log10(colSums(count_data))
pdf("dmm_MAGs_counts.pdf")
densityplot(
    cnts,
    xlim=range(cnts),
    xlab="MAG representation (log 10 count)"
)
dev.off()

#================
# 2. fits for dirichlet multinomial model
#================
fit <- mclapply(
    1:15,  # cluster range from 1 to 15
    dmn,   # function of dirichlet model, followed by its kwargs
    count=count_data,
    verbose=TRUE    
)
save(fit, file='../data/dmm_fit_models.rds')

# optimize number of clusters by measuring Laplace (also works for goodness of fit: AIC, BIC) of fitted model
# choose the k-clusters from which resulting the minimum laplace measurement
lplc <- sapply(fit, laplace)  # used in bayesian settings to approx marginal likelihood for model comparison, while accounting for uncerteinty in params est.
bic <- sapply(fit, BIC)  # best for prediction tasks and larger datasets
aic <- sapply(fit, AIC)  # favored in identifying the most plausible model given the data, esp. for smaller datasets
pdf("../data/dmm_min_laplace.pdf")
plot(
    lplc,
    type="b",
    xlab="Number of Dirichlet Components",
    ylab="Model Fit (Laplace)"
)
plot(
    bic,
    type="b",
    xlab="Number of Dirichlet Components",
    ylab="Model Fit (BIC)"
)
plot(
    aic,
    type="b",
    xlab="Number of Dirichlet Components",
    ylab="Model Fit (AIC)"
)
dev.off()

best <- fit[[which.min(lplc)]]
# report the weight:pi and homogeneity: theta (larger -> more homogeneous) of the fitted model
print(mixturewt(best))

#================
# 3. visualize clustering
#================
# assign community types (estimated Dirichlet components)
community_types <- as.factor(mixture(best, assign=TRUE))
print(table(community_types))

# posterior mean difference between the best and the single-component dirichlet multinomial model
# measures how each component differs from the population average
pop.comp <- fitted(fit[[1]], scale=TRUE)  # scale by theta
best.comp <- fitted(best, scale=TRUE)
colnames(best.comp) <- paste("Community_type", 1:ncol(best.comp), sep="_")

# summarize taxonomic contributions to each dirichlet components
diff <- rowSums(
    abs(best.comp - as.vector(pop.comp))
)
ro <- order(diff, decreasing=TRUE)
# cumulative contribution
cdiff <- cumsum(diff[ro]) / sum(diff)
df_contribs <- cbind(Mean=pop.comp[ro], best.comp, diff=diff[ro], cdiff=cdiff) %>% as.data.frame

# select K features with contribution > 0.01
cutoff <- 0.01
# feats <- rownames(df_contribs)[which(df_contribs$diff > cutoff)]
# feats_splits <- apply(best.comp[feats, ], 1, which.max)
# alternatively, select K features dominant in each community type
sp_per_type <- cbind(tab, community_types) %>% group_by(community_types) %>% summarise_all(mean) %>% select(-community_types)
feats <- which(colSums(sp_per_type) > cutoff) %>% names
feats_splits <- apply(sp_per_type[, feats], 2, FUN=which.max)
print(table(feats_splits))

write.csv(df_contribs, '../data/dmm_tax_contribs.csv', quote=FALSE)
pdf("../data/dmm_tax_contribs.pdf")
plot(df_contribs$diff, type="b", xlab="N largest features", ylab="Difference from single-component model")
abline(v=length(feats), h=cutoff, col="grey")
plot(df_contribs$cdiff, type="b", xlab="N largest features", ylab="Cumulative difference from single-component model")
abline(v=length(feats), h=df_contribs[dplyr::last(feats), "cdiff"], col="grey")
dev.off()


# mat <- as.matrix(t(apply(tab[, feats], 2, scale)))
mat <- as.matrix(t(tab[, feats]))
col_fun <- circlize::colorRamp2(c(0, 0.1, 1), c("#4393C3", "#ffffff", "#D6604D"))
ht <- ComplexHeatmap::Heatmap(
    mat,
    col=col_fun,
    clustering_method_rows='ward.D2',
    clustering_method_columns='ward.D2',
    row_split=feats_splits,
    column_split=community_types,
    column_labels=rep('', nrow(tab)),
    row_labels=sapply(gtdb[feats, 'GTDB_tax'], function(x){ifelse(endsWith(x, 's__'), strsplit(x, 'g__')[[1]][2], strsplit(x, 's__')[[1]][2])}),
    # heatmap_legend_param=list(title='relative abundance (scaled)'),
    heatmap_legend_param=list(title='relative abundance'),
    top_annotation=ComplexHeatmap::HeatmapAnnotation(
        df=subset(cbind(meta, community_types), select=-SubjectID)
    )
)

pdf("../data/dmm_heatmap.pdf", width=15, height=10)
ht
dev.off()

# ballon plot for age bins & community types
age_bins <- cut(meta$age_months, breaks=c(0, 6, 12, 24, 36), include.lowest=T, right=F)
dstats <- cbind(meta,community_types, age_bins) %>% 
    group_by(community_types, age_bins, lifestyle) %>% 
    dplyr::summarise(size=dplyr::n()) %>%
    group_by(age_bins, lifestyle) %>%
    mutate(total_per_group=sum(size), freq=size/total_per_group)
dstats$community_types <- factor(dstats$community_types, levels=c(3,1,2))
dstats$lifestyle <- factor(dstats$lifestyle, levels=c("Industrialized", "Transitional", "Non-industrialized"))

pdf('../data/dmm_vs_meta.pdf', width=10, height=6)
ggpubr::ggballoonplot(
    # as.data.frame(table(community_types, age_bins)),
    dstats,
    x="community_types",
    y="age_bins",
    facet.by="lifestyle",
    fill="size",
    size="freq",
    size.range=c(1,25),
    ylab="Age months",
    xlab="Community types"
) + ggplot2::scale_fill_distiller(direction=1) + ggplot2::theme_linedraw() + ggplot2::labs(fill="Actual sample size", size="Freq in age groups")
dev.off()
