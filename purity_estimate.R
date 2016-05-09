### Systematic analysis of tumor purity
### Author: Kaiyi Zhu
### Date: 9 May 2016
### PART1 -- Get tumor purity estimates

load("TCGA_barcode_TissueSourceSite.rda")
cancerTypes <- c("BLCA","BRCA","COAD","HNSC","KIRP",
                 "LGG","LIHC","LUAD","LUSC","PCPG","PRAD",
                 "SKCM","STAD","THCA","UCEC")

########################## Gene expression -- ESTIMATE ##########################
library("data.table")
library("estimate")
library("preprocessCore")
library("CePa")

allGeneExpr <- fread("unc.edu_PANCAN_IlluminaHiSeq_RNASeqV2.geneExp_whitelisted.tsv", 
                     header = TRUE, sep = "\t", data.table = FALSE)
geneSymbol <- do.call(rbind, strsplit(allGeneExpr[,"gene_id"], split = "\\|"))[,1]
geneSymbol[16301] <- "SLC35E2B"   # There are duplicate names of "SLC35E2"
entrezID <- do.call(rbind, strsplit(allGeneExpr[,"gene_id"], split = "\\|"))[,2]
geneSymbol[which(geneSymbol == "?")] <- entrezID[which(geneSymbol == "?")]
names(entrezID) <- geneSymbol

purity.ESTIMATE <- list()
for(cancer.type in cancerTypes) {
  tss.code <- TSS[which(TSS[,"Cancer.Types"] == cancer.type),"TSS.Code"]
  ge <- allGeneExpr[, which(substr(colnames(allGeneExpr), 6,7) %in% tss.code)]
  rownames(ge) <- geneSymbol
  # write.table(ge, file = paste(cancer.type, "expression_raw.txt", sep = "_"), sep = "\t", quote = FALSE)
  
  ge1 <- log2(ge + 1)
  ge1.log.qn <- normalize.quantiles(data.matrix(ge1)) 
  rownames(ge1.log.qn) <- rownames(ge1)
  colnames(ge1.log.qn) <- colnames(ge1)
  # write.table(ge1.log.qn, file = paste(cancer.type, "expression_log_qn.txt", sep = "_"), sep = "\t", quote = FALSE)
  
  filterCommonGenes(input.f = paste(cancer.type, "expression_log_qn.txt", sep = "_"), 
                    output.f = paste(cancer.type, "10412genes.gct", sep = "_"))
  estimateScore(paste(cancer.type, "10412genes.gct", sep = "_"), paste(cancer.type,"estimate_score.gct", sep = "_"), platform = "illumina")
  scores <- read.gct(paste(cancer.type,"estimate_score.gct", sep = "_"))
  
  # The formula calculating purity estimate from ESTIMATE score comes from ESTIMATE paper.
  purity.ESTIMATE[[cancer.type]] <- cos(0.6049872018 + 0.0001467884*scores["ESTIMATEScore",])
}
save(purity.ESTIMATE, file = "purity.ESTIMATE.rda")


########################## Copy Number Abbreviation -- ABSOLUTE ##########################
ABSOLUTE.data <- read.table(file = "Purity_Ploidy_All_Samples_4-17-15.txt", 
                            header = TRUE, sep="\t", quote = "", stringsAsFactors = FALSE)
ABSOLUTE.data <- ABSOLUTE.data[setdiff(1:nrow(ABSOLUTE.data), which(is.na(ABSOLUTE.data[,"purity"]))),]

purity.ABSOLUTE <- list()
for(cancer.type in cancerTypes) {
  id <- which(substr(ABSOLUTE.data[,"individual_id"], 6,7) %in% TSS[which(TSS[,"Cancer.Types"] == cancer.type),"TSS.Code"])
  tmp <- ABSOLUTE.data[id, "purity"]
  names(tmp) <- ABSOLUTE.data[id, "individual_id"]
  purity.ABSOLUTE[[cancer.type]] <- tmp
}
save(purity.ABSOLUTE, file = "purity.ABSOLUTE.rda")


########################## Somatic Mutations -- DPC ##########################
library("data.table")
maf = fread("tcga_pancancer_082115.vep.filter_whitelisted.maf", header= TRUE, sep = "\t")
tumor.types = unique(maf$tumor_type)
# Separate data into each cancer type
for(i in tumor.types){
  data = maf[which(maf$tumor_type == i),]
  save(data, file = paste("maf_per_cancer_type/",i,".mutations.rda", sep = ""))
}

library("matrixStats")
# 'dpc4project.R' a simplified version of the original 'dpc.R' so that it only outputs tumor purity estimate
source("dpc4project.R")

purity.DPC <- list()
# No read counts in GBM, KICH, PAAD, OV; no file of MESO
for(cancer.type in setdiff(cancerTypes,c("GBM","KICH","MESO","OV","PAAD"))) {
  load(paste("maf_per_cancer_type/",tolower(cancer.type),".mutations.rda", sep = ""))
  samples = unique(data$Tumor_Sample_Barcode)
  purity <- rep(NA, length(samples))
  names(purity) <- samples
  
  for(i in samples) {
    id <- which(data$Tumor_Sample_Barcode == i)
    a <- data$t_alt_count[id]   # variant counts
    d <- data$t_depth[id]       # total counts (= variant counts + reference counts)
    # filter out mutations with a=NA, d=0
    id <- id[setdiff(seq(1,length(id)),union(which(is.na(a)), which(d==0)))]
    if(length(id) == 0)
      next
    else
      purity[i] <- dpc4project(data$t_alt_count[id], data$t_depth[id], paste(cancer.type,".density.txt", sep = ""))
  }
  purity.DPC[[cancer.type]] <- purity[!is.na(purity)]
}
save(purity.DPC, file = "purity.DPC.rda")


########################## Pathological Review -- biospecimen slides ##########################
purity.slides <- list()
for(cancer.type in cancerTypes) {
  slide <- read.table(file = paste("biospecimen_slides/nationwidechildrens.org_biospecimen_slide_",cancer.type,".txt", sep = ""), 
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  slide <- slide[-1,]
  samples <- unique(slide$bcr_sample_barcode)
  purity <- rep(NaN, length(samples))
  names(purity) <- samples
  for(j in samples) {
    id <- which(slide$bcr_sample_barcode == j)
    purity[j] <- mean(as.numeric(slide$percent_tumor_cells[id]), na.rm = TRUE)/100
  }
  purity.slides[[cancer.type]] <- purity[!is.nan(purity)]
}
save(purity.slides, file = "purity.slides.rda")
