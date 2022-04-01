library(survival)
library(survminer)
library(openxlsx)
#library(Xeva)
library(stringr)
library(mosaic)

setwd('~/Desktop/MBP1413Project/data/output/')
files <- list.files('prediction_old')
corr.out <- grep("PDX_corr_",files,value=TRUE)
pred.out <- grep("PDX_pred_",files,value=TRUE)
drugs <- c('cetuximab', 'erlotinib', 'paclitaxel', 'tamoxifen', 'trametinib')
drugs.acroynm <- c('ceab','erib','pael','taen','trib')
#pdxe <- readRDS('~/Desktop/MBP1413Project/data/Xeva_PDXE.rds')
pdx.ids <- readRDS('PDX.ids.RDS')
meta <- read.xlsx('41591_2015_BFnm3954_MOESM10_ESM.xlsx', sheet = 'PCT curve metrics')
meta$Model <- str_replace( meta$Model, '-', '\\.')

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
# use k=8, num_trial=20

k <- 10
i <- 3
drug_name <- drugs[i]
cur.samples <- unlist(pdx.ids[drug_name])
pred <- read.csv(paste0('prediction_old/PDX_pred_', firstup(drug_name),'_k',k,'.csv'))
rownames(pred) <- cur.samples

sMeta <- meta[meta$Treatment == drug_name,]
sMeta$Model <- paste(sMeta$Model, drugs.acroynm[i], sep=".")

#sPdxe <- subsetXeva(pdxe, ids=drug_name, id.name="drug")
cur.drug.outcome <- as.data.frame(matrix(NA, nrow = length(cur.samples), ncol = 6))
for (i in 1:length(cur.samples)) {
  #cur.experiment <- getExperiment(sPdxe, model.id = cur.samples[i])
  surv.time <- sMeta[sMeta$Model == cur.samples[i], 'Day_Last']
  double.before <- sMeta[sMeta$Model == cur.samples[i], 'TimeToDouble'] > sMeta[sMeta$Model == cur.samples[i], 'Day_Last']
  status <- ifelse(double.before,0,1)
  vol <- pred[rownames(pred) == cur.samples[i], 20]
  vol.class <- ifelse( vol > 0.3, 'PD', 'SD' )
  cur.drug.outcome[i,1] <- cur.samples[i]
  cur.drug.outcome[i,2] <- as.numeric(surv.time)
  cur.drug.outcome[i,3] <- as.numeric(status)
  cur.drug.outcome[i,4] <- as.numeric(vol)
  cur.drug.outcome[i,5] <- vol.class
  actual.outcome <- ifelse(sMeta[sMeta$Model == cur.samples[i], 'BestAvgResponse'] > 0.3 * 100, 'Progressive', 'Stable')
  cur.drug.outcome[i,6] <- actual.outcome
}
colnames(cur.drug.outcome) <- c('id','time','status','deltaV','class','actual.outcome')
write.csv(cur.drug.outcome, file=paste0("survival/",drug_name,".survival.csv"))
fit <- survfit(Surv(time, status) ~ class, data = cur.drug.outcome)

#### SURVIVAL PLOT

pdf(paste0('res/', drug_name, '.survplot.pdf'), width = 16, height = 9)
print(ggsurvplot(fit, pval = T,
           legend="bottom",
           legend.title="",
           legend.labs = c("predicted drug response >= 30 (28 PDX models)", 
                           "predicted drug response < 30 (8 PDX models)"),
           title = '2-shot trametinib log(rank P) = 0.2',
           data = cur.drug.outcome) + guides(colour = guide_legend(nrow = 2)))
dev.off()

#### WATER FALL PLOT
cur.drug.outcome <- cur.drug.outcome[order(cur.drug.outcome$deltaV, decreasing=TRUE),]
x <- 1:nrow(cur.drug.outcome)
wilcox.test(deltaV ~ actual.outcome, data = cur.drug.outcome, exact = FALSE, alternative = "greater")
p <- ggplot(cur.drug.outcome, aes(x=x, y=deltaV * 100, fill=actual.outcome, color=actual.outcome)) +
  geom_bar(stat="identity", width=0.7, position = position_dodge(width=0.4)) + 
  scale_fill_discrete(name="Actual outcome") + 
  scale_color_discrete(guide="none") +
  ylab(expression(atop("Predicted drug response", paste('(%',Delta, Vol,')')))) + 
  labs(title='1 shot paclitaxel Rank test P = 0.8169') + 
  theme_classic() %+replace%
  theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

pdf(paste0('plots/',drug_name, '_waterfall.pdf'), width=16, heigh=9)
print(p)
dev.off()

#### ODDS RATIO #SD+ / SD- * PD+ / PD-
contingency_tab <- function(df) {
  t(table(df[,c('class','actual.outcome')]))
}
for (drug_name in drugs) {
  cat(drug_name)
  survival.info <- read.csv(paste0("survival/",drug_name,".survival.csv"))
  #oddsRatio(contingency_tab(survival.info), verbose=TRUE)
  cat(drug_name)
  print(contingency_tab(survival.info))
}

OR.summary <- as.data.frame(matrix(NA, nrow=5, ncol=3))
OR.summary[1,] <- c(0.5357, 0.1545, 1.857)
OR.summary[2,] <- c(0.8571, 0.06527, 11.26)
OR.summary[3,] <- c(0.8409, 0.072, 9.821)
OR.summary[4,] <- c(1.182, 0.09735, 14.35)
OR.summary[5,] <- c(0.4286, 0.04458, 4.12)
colnames(OR.summary) <- c('OR','lb','ub')
rownames(OR.summary) <- drugs
write.csv(OR.summary, file="res/or_summary.csv")

or.p <- ggplot(OR.summary, aes(x = OR, y = firstup(rownames(OR.summary)))) + 
    geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
    geom_errorbarh(aes(xmax = ub, xmin = lb), size = .5, height = .2, color = "black") +
    geom_point(size = 1.5, color = "black") +
    theme_classic()+
    ylab("") +
    xlab("Odds ratio (95% CI)")

pdf(file='plots/oddsRatio.pdf')
print(or.p)
dev.off()

  