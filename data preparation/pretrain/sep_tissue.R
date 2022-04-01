setwd('~/Desktop/MBP1413Project/data/processed/final')

resp <- read.csv("drugResp_auc.csv")
resp <- tibble::column_to_rownames(resp,var = "cell.line")
drugs <- c('Cetuximab','Paclitaxel','Tamoxifen','Trametinib','Erlotinib')

for (drug_name in drugs) {
  feats <- readRDS(paste0(drug_name,"_feats.RDS"))
  drug.metadata <- read.csv(paste0(drug_name,"_meta.csv"))
  drug.metadata <- drug.metadata[order(drug.metadata$tissue),]
  
  cur_resp <- resp[,drug_name]
  names(cur_resp) <- rownames(resp)
  cur_resp <- cur_resp[!is.na(cur_resp)]
  
  avail_tissue <- unique(drug.metadata$tissue)
  for (tissue in avail_tissue) {
    cellline <- drug.metadata[drug.metadata$tissue==tissue,'cellline_name']
    tissue.resp <- as.data.frame(cur_resp[(names(cur_resp) %in% cellline)])
    colnames(tissue.resp) <- "AUC"
    
    tissue.feats <- feats[,colnames(feats) %in% cellline]
    tissue.feats.desc <- apply(tissue.feats, MARGIN=1, FUN=function(g) {
      ind <- length(levels(as.factor(g))) <= 2
      ret <- ifelse(ind, "_mutation", "_expression")
      ret
    })
    tissue.desc <- paste0(rownames(tissue.feats), tissue.feats.desc)
    filename.feats <- paste(tissue, drug_name, "feature.csv", sep="_")
    filename.desc <- paste(tissue, drug_name, "description.csv", sep="_")
    outdir.feats <- paste0('by_tissue/',filename.feats)
    #outdir.desc <- paste0('feature_by_tissue/',filename.desc)
    tissue.feats <- t(tissue.feats)
    colnames(tissue.feats) <- tissue.desc
    tissue.feats <- tissue.feats[,order(colnames(tissue.feats))]
    
    tissue.all <- cbind(tissue.resp, tissue.feats)
    tissue.all <- tissue.all[complete.cases(tissue.all),]

    write.csv(tissue.all, file = outdir.feats, row.names = TRUE)
    #write.csv(tissue.desc, file = outdir.desc, row.names = FALSE)
  }
}
