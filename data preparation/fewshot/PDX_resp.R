library(Xeva)
library(SummarizedExperiment)
setwd('~/Desktop/MBP1413Project')

pdxe <- readRDS('data/Xeva_PDXE.rds')
fullmod <- modelInfo(pdxe)
length(unique(fullmod$model.id))  # 4706
length(unique(fullmod$patient.id)) # 277

drugs <- c('cetuximab', 'erlotinib', 'paclitaxel', 'tamoxifen', 'trametinib')
sPdxe <- subsetXeva(pdxe, ids=drugs, id.name="drug")

mod <- modelInfo(sPdxe)
length(unique(mod$model.id)) # 392
length(unique(mod$patient.id)) # 149

compute_minChangeVol <- function(dV, days) {
  over10days <- which(days >= 10)
  if (length(over10days) == 0) { return(NA) }  # no data are over 10 days
  if (length(over10days) == length(days)) { over10days <- over10days[2:length(over10days)] } # all data are over 10 days
  dVt <- dV[over10days]
  #percentageChange <- (abs((Vt-Vinit))/Vinit) * 100
  minChange <- min(abs(dVt))
  return(minChange)
}

res <- sapply(rownames(mod), FUN=function(model) {
  model.data <- getExperiment(sPdxe, model.id = model)
  compute_minChangeVol(model.data$volume.normal, model.data$time)
})

df <- data.frame(cbind(names(res), res))
colnames(df) <- c('model', 'delta_vol')
write.csv(df, file="data/processed/PDXE_dVol.csv", row.names = FALSE)

NAmods <- names(res[which(is.na(res))]) # none of these have volume time > 10
model.data <- getExperiment(sPdxe, model.id = NAmods[7])
model.data

#library(readxl)
#sheet <- read_xlsx("data.xlsx", sheet="PCT curve metrics")
#sSheet <- sheet[sheet$Treatment %in% drugs,]
