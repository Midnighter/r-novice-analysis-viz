## ---- include=FALSE------------------------------------------------------
source("../bin/chunk-options.R")
knitr_fig_path("03-explore-gene-expression-")

## ---- eval=FALSE---------------------------------------------------------
## source("https://bioconductor.org/biocLite.R")
## biocLite(c("edgeR", "org.EcK12.eg.db"))
## install.packages("locfit")

## ---- results="hide"-----------------------------------------------------
library(edgeR)
library(ggplot2)
library(org.EcK12.eg.db)

## ---- eval=FALSE---------------------------------------------------------
## edgeRUsersGuide()

## ------------------------------------------------------------------------
wulffenTable <- read.table("data/GSE71562.csv", header = TRUE, row.names = 1, sep = ",")
head(wulffenTable)

## ------------------------------------------------------------------------
samples <- read.table("data/pheno.csv", header = TRUE, row.names = 1, sep = ",")
samples

## ------------------------------------------------------------------------
wulffen <- DGEList(counts=wulffenTable, genes=rownames(wulffenTable),
                   samples=samples)
wulffen <- calcNormFactors(wulffen)

## ------------------------------------------------------------------------
wulffenCpm <- cpm(wulffen)

## ------------------------------------------------------------------------
scores <- prcomp(log2(t(wulffenCpm) + 0.25))$x

## ------------------------------------------------------------------------
pcaDf <- merge(scores, samples, by=0)

## ------------------------------------------------------------------------
ggplot(pcaDf, aes(PC1, PC2, label=time, color=replicate)) +
    geom_text()

## ------------------------------------------------------------------------
wulffenShort <- wulffen[, wulffen$samples$time %in% c("t0", "t10")]
design <- model.matrix(~as.character(time), data=wulffenShort$samples)
colnames(design) <- c("(Intercept)", "t10")
design

## ------------------------------------------------------------------------
wulffenShort <- estimateDisp(wulffenShort, design)
fit <- glmFit(wulffenShort, design)
lrt <- glmLRT(fit)
topTags(lrt)

## ---- echo=FALSE, eval=FALSE---------------------------------------------
## df <- merge(topTags(lrt, n=Inf), read.csv("data/ecoli.csv"))
## write.csv(df, file="data/deg.csv", quote=FALSE, row.names=FALSE)

## ------------------------------------------------------------------------
library(org.EcK12.eg.db)
symbol2entrez <- mapIds(org.EcK12.eg.db, rownames(lrt), "ENTREZID", keytype="SYMBOL")
universe <- na.omit(symbol2entrez)

## ------------------------------------------------------------------------
fdr <- p.adjust(lrt$table[,"PValue"], "fdr")
allSymbols <- rownames(lrt$table)
deSymbols <- allSymbols[fdr < 0.05 & !is.na(symbol2entrez)]
deEntrez <- symbol2entrez[deSymbols]

## ------------------------------------------------------------------------
goTable <- goana(deEntrez, universe=universe, species="EcK12")
head(goTable[order(goTable$P.DE),])

