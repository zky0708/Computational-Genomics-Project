### Systematic analysis of tumor purity
### Author: Kaiyi Zhu
### Date: 9 May 2016
### PART2 -- Comparing purity estimates from each method

load("purity.ABSOLUTE.rda")
load("purity.slides.rda")
load("purity.ESTIMATE.rda")
load("purity.DPC.rda")

cancerTypes <- c("BLCA","BRCA","COAD","HNSC","KIRP",
                 "LGG","LIHC","LUAD","LUSC","PCPG","PRAD",
                 "SKCM","STAD","THCA","UCEC")

# Filter out matched normal samples in tumor purity estimates based on ESTIMATE and biospecimen slides
purity.ESTIMATE.tumor = sapply(purity.ESTIMATE, function(x) x[which(substr(names(x), 14, 14) == "0")])
purity.slides.tumor = sapply(purity.slides, function(x) x[which(substr(names(x), 14, 14) == "0")])


## create a table containing tumor purity estimates of all samples by all methods
totalPurities <- c()

for(cancer.type in cancerTypes) {
  a <- purity.ESTIMATE.tumor[[cancer.type]]
  b <- purity.ABSOLUTE[[cancer.type]]
  c <- purity.DPC[[cancer.type]]; 
  d <- purity.slides.tumor[[cancer.type]]
  
  # unify the barcode formats of each sample
  if(length(a)>0) names(a) <- substr(gsub("[.]", "-", names(a)),1,12)
  if(length(c)>0) names(c) <- substr(names(c), 1, 12)
  if(length(d)>0) names(d) <- substr(names(d), 1, 12)
  samples <- Reduce(union, x <- c(names(a), names(b), names(c), names(d)))
  
  tmp <- matrix(NA, nrow = length(samples), ncol = 5)
  rownames(tmp) <- samples
  colnames(tmp) <- c("CancerType", "ESTIMATE", "ABSOLUTE", "DPC", "Pathology")
  tmp[, 1] <- cancer.type
  tmp[names(a),2] <- round(a, digits = 4)
  tmp[names(b),3] <- round(b, digits = 4)
  tmp[names(c),4] <- round(c, digits = 4)
  tmp[names(d),5] <- round(d, digits = 4)
  
  totalPurities <- rbind(totalPurities, tmp)
}

total <- apply(totalPurities[,-1], 2, as.numeric)
rownames(total) <- rownames(totalPurities)
# colMeans(total, na.rm = TRUE)
# apply(total, 2, function(x) sd(x, na.rm = TRUE))

# Use the median as the consensus purity
totalPurities <- cbind(totalPurities, apply(total, 1, function(x) round(median(x, na.rm = TRUE), digits = 4)))
colnames(totalPurities)[6] <- "Consensus"
# save(totalPurities, file = "purity_all_samples.rda")
write.csv(totalPurities, file = "purity_all_samples.csv", quote = FALSE)

## Get the sample size available for each method and each cancer type
numMatrix <- matrix(NA, ncol = 6, nrow = length(cancerTypes))
rownames(numMatrix) <- cancerTypes
colnames(numMatrix) <- c("CancerTypeName","ESTIMATE(geneExpr)","ABSOLUTE(CNA)","DPC(mutation)","Pathology","Consensus")
for(cancer.type in cancerTypes){
  numMatrix[cancer.type, "CancerTypeName"] <- unique(TSS$Study.Name[which(TSS$Cancer.Types == cancer.type)])
  id <- which(totalPurities[,"CancerType"] == cancer.type)
  numMatrix[cancer.type,"ESTIMATE(geneExpr)"] <- sum(!is.na(totalPurities[id,"ESTIMATE"]))
  numMatrix[cancer.type,"ABSOLUTE(CNA)"] <- sum(!is.na(totalPurities[id,"ABSOLUTE"]))
  numMatrix[cancer.type,"DPC(mutation)"] <- sum(!is.na(totalPurities[id,"DPC"]))
  numMatrix[cancer.type,"Pathology"] <- sum(!is.na(totalPurities[id,"Pathology"]))
  numMatrix[cancer.type,"Consensus"] <- sum(!is.na(totalPurities[id,"Consensus"]))
}
write.csv(numMatrix, file = "samples_summary.csv", quote = FALSE)


## Compute pairwise correlation
minSampleSize = 100

corMatrix <- matrix(NA, nrow = 6, ncol = length(cancerTypes)+1)
colnames(corMatrix) <- c(cancerTypes, "Total")
rownames(corMatrix) <- c("ESTIMATE_ABSOLUTE","ESTIMATE_DPC","ABSOLUTE_DPC",
                         "ESTIMATE_Pathology","ABSOLUTE_Pathology", "DPC_Pathology")

corMatrix["ESTIMATE_ABSOLUTE","Total"] <- cor(total[,"ESTIMATE"], total[,"ABSOLUTE"], use = "na.or.complete")
corMatrix["ESTIMATE_DPC","Total"] <- cor(total[,"ESTIMATE"], total[,"DPC"], use = "na.or.complete")
corMatrix["ABSOLUTE_DPC","Total"] <- cor(total[,"DPC"], total[,"ABSOLUTE"], use = "na.or.complete")
corMatrix["ESTIMATE_Pathology","Total"] <- cor(total[,"ESTIMATE"], total[,"Pathology"], use = "na.or.complete")
corMatrix["ABSOLUTE_Pathology","Total"] <- cor(total[,"ABSOLUTE"], total[,"Pathology"], use = "na.or.complete")
corMatrix["DPC_Pathology","Total"] <- cor(total[,"DPC"], total[,"Pathology"], use = "na.or.complete")

for(cancer.type in cancerTypes) {
  a <- purity.ESTIMATE.tumor[[cancer.type]]
  b <- purity.ABSOLUTE[[cancer.type]]
  c <- purity.DPC[[cancer.type]]
  d <- purity.slides.tumor[[cancer.type]]
  
  if(length(a)>0) names(a) <- substr(gsub("[.]", "-", names(a)),1,12)
  if(length(c)>0) names(c) <- substr(names(c), 1, 12)
  if(length(d)>0) names(d) <- substr(names(d), 1, 12)
  
  common_ab <- intersect(names(a), names(b))
  common_ac <- intersect(names(a), names(c))
  common_ad <- intersect(names(a), names(d))
  common_bc <- intersect(names(b), names(c))
  common_bd <- intersect(names(b), names(d))
  common_cd <- intersect(names(c), names(d))
  
  if(length(common_ab) >= minSampleSize){
    corMatrix["ESTIMATE_ABSOLUTE",cancer.type] <- cor(a[common_ab], b[common_ab], use = "na.or.complete")
  }
  if(length(common_ac) >= minSampleSize){
    corMatrix["ESTIMATE_DPC",cancer.type] <- cor(a[common_ac], c[common_ac], use = "na.or.complete")
  }
  if(length(common_ad) >= minSampleSize){
    corMatrix["ESTIMATE_Pathology",cancer.type] <- cor(a[common_ad], d[common_ad], use = "na.or.complete")
  }
  if(length(common_bc) >= minSampleSize){
    corMatrix["ABSOLUTE_DPC",cancer.type] <- cor(b[common_bc], c[common_bc], use = "na.or.complete")
  }
  if(length(common_bd) >= minSampleSize){
    corMatrix["ABSOLUTE_Pathology",cancer.type] <- cor(b[common_bd], d[common_bd], use = "na.or.complete")
  }
  if(length(common_cd) >= minSampleSize){
    corMatrix["DPC_Pathology",cancer.type] <- cor(c[common_cd], d[common_cd], use = "na.or.complete")
  }
  
}


## Draw heatmap of Pearson correlation
library("fields")
par(oma = c(0,0,0,0))
m <- nrow(corMatrix)
n <- ncol(corMatrix)
# image.plot(1:n, 1:m, t(corMatrix[m:1,]), axes = FALSE, xlab = "", ylab = "")
# text(1:n, -0.2, labels = colnames(corMatrix), srt = 90, xpd = TRUE, cex = 1)
# text(-1, 1:m, labels = rev(rownames(corMatrix)), xpd = TRUE, cex = 0.8)
image.plot(1:m, 1:n, corMatrix[,n:1], axes = FALSE, xlab = "", ylab = "")
text(-0.3, 1:n, labels = rev(colnames(corMatrix)), xpd = TRUE, cex = 1)
text(1:m, -3.1, labels = rownames(corMatrix), srt = 90, xpd = TRUE, cex = 1)


## Draw pairwise scatterplots
groupScatterplot("ESTIMATE","ABSOLUTE")
groupScatterplot("ESTIMATE","DPC")
groupScatterplot("ESTIMATE","Pathology")
groupScatterplot("ABSOLUTE","DPC")
groupScatterplot("ABSOLUTE","Pathology")
groupScatterplot("DPC","Pathology")

# function for creating group scatterplots 
groupScatterplot <- function(method1, method2) {
  par(mfrow = c(4,4), oma = c(0,3,3,0),mai = c(0.3,0.25,0.25,0.15))
  for(cancer.type in cancerTypes) {
    id <- which(totalPurities[,"CancerType"] == cancer.type)
    x <- as.numeric(totalPurities[id, method1])
    y <- as.numeric(totalPurities[id, method2])
    corr <- cor(x, y, use = "na.or.complete")
    plot(x, y, pch = 20, col = "blue", cex = 0.5,
         xlim = c(0,1), ylim = c(0,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "",
         main = paste(cancer.type, ", r = ", round(corr, digits = 4), sep = ""), line = 0.2)
    axis(1, at = c(0,0.5,1))
    axis(2, at = c(0,0.5,1))
    abline(a = 0, b = 1, lty = 2)
  }
  mtext(method1, side = 3, outer = TRUE)
  mtext(method2, side = 2, outer = TRUE, line = 1)
}

