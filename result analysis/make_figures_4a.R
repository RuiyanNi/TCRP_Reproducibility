library(reshape2)
library(ggplot2)
library(boot)
library(stringr)
library(tibble)
library(RColorBrewer)
library(openxlsx)

setwd('~/Desktop/MBP1413Project/data/output/')
files <- list.files('prediction_old')
corr.out <- grep("PDX_corr_",files,value=TRUE)
pred.out <- grep("PDX_pred_",files,value=TRUE)
drug <- c("Cetuximab", "Erlotinib", "Paclitaxel", "Tamoxifen", "Trametinib")

###### MAKE FIG4A.
comp_ci <- function(corr.mat, R=1000) {
  corr.mat.reshape <- melt(corr.mat)
  rs <- corr.mat.reshape$value
  
  r <- function(data, indices) {
    d <- data[indices] # allows boot to select sample
    return(mean(d))
  }
  
  boot.res <- boot(rs,r,R=1000)
  ci <- boot.ci(boot.res,type='basic')$basic[4:5]

  return( c(boot.res$t0, ci) )
}

kshot.pear.info <- list()
kshot.spear.info <- list()

for (k in 1:10) {
  fs <- grep(paste0("k",k,".csv"),corr.out,value=TRUE)
  pear.corr.mat.k <- matrix(NA, nrow=20, ncol=length(fs))
  spear.corr.mat.k <- pear.corr.mat.k
  dn <- rep(NA,length(fs)) 
  for (i in 1:length(fs)) {
    df <- read.csv(paste0('prediction_old/',fs[i]))
    pear.corr.mat.k[,i] <- df$Pearson
    spear.corr.mat.k[,i] <- df$Spearman
    dn[i] <- substr(fs[i], 10, gregexpr("\\.",fs[i])[[1]][1]-1)
  }
  colnames(pear.corr.mat.k) <- dn
  colnames(spear.corr.mat.k) <- dn
  pear.corr.mat.k <- as.data.frame(pear.corr.mat.k)
  spear.corr.mat.k <- as.data.frame(spear.corr.mat.k)
  
  write.csv(pear.corr.mat.k, file=paste0('res/pear_corr_k',k,'.csv'), row.names = FALSE)
  write.csv(spear.corr.mat.k, file=paste0('res/spear_corr_k',k,'.csv'), row.names = FALSE)

  pear.ci <- as.data.frame(t(sapply( pear.corr.mat.k, FUN=comp_ci )))
  colnames(pear.ci) <- c('mu', 'lb', 'ub')
  pear.ci$k <- k
  rownames(pear.ci) <- str_replace(rownames(pear.ci), paste0("_k",k), "")
  pear.ci <- rownames_to_column(pear.ci,'drug')
  kshot.pear.info[[k]] <- pear.ci

  spear.ci <- as.data.frame(t(sapply( spear.corr.mat.k, FUN=comp_ci )))
  colnames(spear.ci) <- c('mu', 'lb', 'ub')
  spear.ci$k <- k
  rownames(spear.ci) <- str_replace(rownames(spear.ci), paste0("_k",k), "")
  spear.ci <- rownames_to_column(spear.ci,'drug')
  kshot.spear.info[[k]] <- spear.ci

}
kshot.pear.info <- do.call(rbind.data.frame, kshot.pear.info)
kshot.spear.info <- do.call(rbind.data.frame, kshot.spear.info)


####
zero.shot <- read.xlsx('res/zero_shot.xlsx')
zero.shot <- cbind(zero.shot, cbind(rep(0,5), rep(0,5)))
zero.shot.pear <- zero.shot[,c('drug','Pearson','1','2')]
zero.shot.spear <- zero.shot[,c('drug','Spearman','1','2')]

zero.shot.pear$k=0
zero.shot.spear$k=0

colnames(zero.shot.pear) <- c('drug','mu','lb','ub','k') 
colnames(zero.shot.spear) <- c('drug','mu','lb','ub','k') 

kshot.pear.info <- read.csv('res/pear_corr_ci.csv')
kshot.spear.info <- read.csv('res/spear_corr_ci.csv')

kshot.pear.info <- rbind(zero.shot.pear, kshot.pear.info)
kshot.spear.info <- rbind(zero.shot.spear, kshot.spear.info)

kshot.pear.info[1:5,c('lb','ub')] <- kshot.pear.info[1:5,'mu']
kshot.spear.info[1:5,c('lb','ub')] <- kshot.spear.info[1:5,'mu']

#write.csv(kshot.pear.info, file = "res/pear_corr_ci.csv", row.names=FALSE)
#write.csv(kshot.spear.info, file = "res/spear_corr_ci.csv", row.names=FALSE)


kshot.pear.info.reshape <- melt(kshot.pear.info,id.vars = c('mu','lb','ub','k'),value.name = 'drug')
kshot.spear.info.reshape <- melt(kshot.spear.info,id.vars = c('mu','lb','ub','k'),value.name = 'drug')

pear.plot <- ggplot(kshot.pear.info.reshape, aes(x=k, y=mu, colour=drug)) + 
  geom_line() +
  geom_point(colour='black') +
  geom_errorbar(aes(ymin=lb, ymax=ub), width=.2, position=position_dodge(0.1)) +
  theme_classic() + theme(legend.position = "bottom")+
  ylab("correlation (predicted, actual)") + xlab("Number of samples for drug response") +
  scale_x_continuous(labels = 0:10, breaks = 0:10) +
  scale_color_brewer(palette="Set1")

spear.plot <- ggplot(kshot.spear.info.reshape, aes(x=k, y=mu, colour=drug)) + 
  geom_line() +
  geom_point(colour='black') +
  geom_errorbar(aes(ymin=lb, ymax=ub), width=.2, position=position_dodge(0.1)) +
  theme_classic() + theme(legend.position = "bottom")+
  ylab("correlation (predicted, actual)") + xlab("Number of samples for drug response") +
  scale_x_continuous(labels = 1:10, breaks = 1:10) +
  scale_color_brewer(palette="Set1")

pdf(file="plots/corr-plots.pdf",width = 16, height = 9,onefile = TRUE)
print(pear.plot)
print(spear.plot)
dev.off()
