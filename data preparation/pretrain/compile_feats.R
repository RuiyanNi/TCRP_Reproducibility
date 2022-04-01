library(Biobase)
library(SummarizedExperiment)
library(S4Vectors)
library(PharmacoGx)
library(tidyverse)

setwd('~/Desktop/MBP1413Project')

resp <- readRDS('data/processed/final/drugResp_auc.RDS')
GDSC1.feats <- readRDS('data/processed/final/GDSC1.combined.RDS')
GDSC2.feats <- readRDS('data/processed/final/GDSC2.combined.RDS')


for (drug in c("Erlotinib",  "Paclitaxel", "Tamoxifen",  "Trametinib")) { 
  cur_resp <- resp[drug,]
  cur_resp <- cur_resp[,!is.na(cur_resp)]
  celllines <- colnames(cur_resp) # cell lines with auc data
  allFeats <- GDSC2.feats[,colnames(GDSC2.feats) %in% celllines]
  saveRDS(allFeats, file=paste0("data/processed/final/", drug, "_feats.RDS"))
}


drug <- 'Cetuximab'
cur_resp <- resp[drug,]
cur_resp <- cur_resp[,!is.na(cur_resp)]
celllines <- colnames(cur_resp) # cell lines with auc data
allFeats <- GDSC1.feats[,colnames(GDSC1.feats) %in% celllines]
saveRDS(allFeats, file=paste0("data/processed/final/", drug, "_feats.RDS"))
