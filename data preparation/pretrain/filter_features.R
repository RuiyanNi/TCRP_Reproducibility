setwd('~/Desktop/MBP1413Project')

library(Biobase)
library(SummarizedExperiment)
library(S4Vectors)
library(PharmacoGx)
library(simpIntLists)
library(gprofiler2)


drugs <- c('Cetuximab', 'Erlotinib', 'Paclitaxel', 'Tamoxifen', 'Trametinib')
data("HumanBioGRIDInteractionOfficial")
EGFR <- c()
ESR1 <- c()
MAP2K1 <- c()
MAP2K2 <- c()
CDK1 <- c()
CDK2 <- c()
for (item in HumanBioGRIDInteractionOfficial) {
  if (item$name == "EGFR") {EGFR <- c(EGFR, item)}
  if (item$name == "ESR1") {ESR1 <- c(ESR1, item)}
  if (item$name == "MAP2K1") {MAP2K1 <- c(MAP2K1, item)}
  if (item$name == "MAP2K2") {MAP2K2 <- c(MAP2K2, item)}
  if (item$name == "CDK1") {CDK1 <- c(CDK1, item)} # mitosis
  if (item$name == "CDK2") {CDK2 <- c(CDK2, item)} # mitosis
}
interactions <- c(EGFR$interactors, ESR1$interactors, 
                  MAP2K1$interactors, MAP2K2$interactors, 
                  CDK1$interactors, CDK2$interactors)
pp_feats <- read.csv("data/gene_symbol_ppi.csv")
pp_feats <- colnames(pp_feats)

feats <- c(interactions, pp_feats)
conversion <- gconvert(query = feats, organism = "hsapiens",target="ENSG", mthreshold = Inf, filter_na = TRUE)
write.csv(conversion, file="data/processed/conversion.20220326.csv",row.names = FALSE)


GDSC1<- readRDS('./data/GDSC1.rds')
GDSCexpression1 <- summarizeMolecularProfiles(
  GDSC1,
  cellNames(GDSC1),
  mDataType="rna",
  features=fNames(GDSC1,'rna'),
  verbose=FALSE
)
expr1 <- assay(GDSCexpression1,1)
expr1 <- expr1[rownames(expr1) %in% conversion$target,]
convert1 <- gconvert(query = rownames(expr1), organism = "hsapiens",target="HGNC", mthreshold = Inf, filter_na = TRUE)
expr1 <- expr1[rownames(expr1) %in% convert1$input,]
rownames(expr1) <- convert1$target

exprSds <- rowSds(expr1,na.rm=T)
exprSdsDecile <- quantile(exprSds, probs = seq(0, 1, 1/10))
exprToKeep <- exprSds > exprSdsDecile[2]
expr1 <- expr1[exprToKeep,]

saveRDS(expr1,"data/processed/GDSC1.expr.RDS")

GDSCmutation1 <- summarizeMolecularProfiles(
  GDSC1,
  cellNames(GDSC1),
  mDataType="mutation",
  features=fNames(GDSC1,'mutation'),
  summary.stat = "and",
  verbose=FALSE
) # mutation data
mutations1 <- assay(GDSCmutation1,1)
rn <- rownames(mutations1)
mutations1 <- as.data.frame(apply(as.matrix(mutations1),2,as.numeric))
rownames(mutations1) <- rn

GDSCmutation_exome1 <- summarizeMolecularProfiles(
  GDSC1,
  cellNames(GDSC1),
  mDataType="mutation_exome",
  features=fNames(GDSC1,'mutation_exome'),
  summary.stat = "and",
  verbose=FALSE
) # exome mutation data
mutation_exome1 <- assay(GDSCmutation_exome1,1)
rn <- rownames(mutation_exome1)
mutation_exome1 <- as.data.frame(apply(as.matrix(mutation_exome1),2,as.numeric))
rownames(mutation_exome1) <- rn

mutations1 <- rbind(mutations1, mutation_exome1)
mutations1 <- mutations1[rownames(mutations1) %in% conversion$input,]

mutationCount <- rowSums(mutations1,na.rm=T)
mutationToKeep <- mutationCount > 10
mutations1 <- mutations1[mutationToKeep,]

mutations1[is.na(mutations1)] <- 0

saveRDS(mutations1,"data/processed/GDSC1.mutations.RDS")

GDSC1.feats <- rbind(expr1, mutations1)
saveRDS(GDSC1.feats,"data/processed/final/GDSC1.combined.RDS")


GDSC2<- readRDS('./data/GDSC2.rds')
GDSCexpression <- summarizeMolecularProfiles(
  GDSC2,
  cellNames(GDSC2),
  mDataType="rna",
  features=fNames(GDSC2,'rna'),
  verbose=FALSE
)
expr <- assay(GDSCexpression,1)
expr <- expr[rownames(expr) %in% conversion$target,]
convert <- gconvert(query = rownames(expr),organism = "hsapiens",target="HGNC", mthreshold = Inf, filter_na = TRUE)
expr <- expr[rownames(expr) %in% convert$input,]
rownames(expr) <- convert$target

exprSds <- rowSds(expr,na.rm=T)
exprSdsDecile <- quantile(exprSds, probs = seq(0, 1, 1/10))
exprToKeep <- exprSds > exprSdsDecile[2]
expr <- expr[exprToKeep,]

saveRDS(expr,"data/processed/GDSC2.expr.RDS")

GDSCmutation <- summarizeMolecularProfiles(
  GDSC2,
  cellNames(GDSC2),
  mDataType="mutation",
  features=fNames(GDSC2,'mutation'),
  summary.stat = "and",
  verbose=FALSE
) # mutation data
mutations <- assay(GDSCmutation,1)
rn <- rownames(mutations)
mutations <- as.data.frame(apply(as.matrix(mutations),2,as.numeric))
rownames(mutations) <- rn

GDSCmutation_exome <- summarizeMolecularProfiles(
  GDSC2,
  cellNames(GDSC2),
  mDataType="mutation_exome",
  features=fNames(GDSC2,'mutation_exome'),
  summary.stat = "and",
  verbose=FALSE
) # exome mutation data
mutation_exome <- assay(GDSCmutation_exome,1)
rn <- rownames(mutation_exome)
mutation_exome <- as.data.frame(apply(as.matrix(mutation_exome),2,as.numeric))
rownames(mutation_exome) <- rn

mutations <- rbind(mutations, mutation_exome)
mutations <- mutations[rownames(mutations) %in% conversion$input,]

mutationCount <- rowSums(mutations,na.rm=T)
mutationToKeep <- mutationCount > 10
mutations <- mutations[mutationToKeep,]

mutations[is.na(mutations)] <- 0

saveRDS(mutations,"data/processed/GDSC2.mutations.RDS")

GDSC2.feats <- rbind(expr, mutations)
saveRDS(GDSC2.feats,"data/processed/final/GDSC2.combined.RDS")
