### Systematic analysis of tumor purity
### PART3 -- Clinical Association Examination
### Author: Kaiyi Zhu
### Date: 9 May 2016

allClinical <- read.csv(file = "clinical_PANCAN_patient_with_followup.tsv", sep = "\t", stringsAsFactors = FALSE)
# Separate data to each cancer type
types <- unique(allClinical$acronym)
clinical_data <- list()
for(i in types) {
  id <- which(allClinical$acronym == i)
  clinical_data[[i]] <- allClinical[id, ]
}
save(clinical_data, file = "clinical_per_cancer_type.rda")

load("clinical_per_cancer_type.rda")
load("purity_all_samples.rda")

nonsense <- c("","[Not Available]","[Not Applicable]","[Unknown]","[Not Evaluated]","[Discrepancy]", NA)
clinical_features <- colnames(clinical_data[[1]])
survival_features <- c("vital_status","days_to_death","days_to_last_followup",
                       "days_to_initial_pathologic_diagnosis","days_to_birth")
id_features <- c("bcr_patient_uuid","bcr_patient_barcode","patient_id","acronym","tissue_source_site","form_completion_date",
                 "year_of_initial_pathologic_diagnosis","system_version","surgical_procedure_purpose_other_text",
                 "state_province_country_of_procurement", "city_of_procurement")
continuous_features <- c("age_at_initial_pathologic_diagnosis","lymph_node_examined_count","number_of_lymphnodes_positive_by_he","mitoses_count",
                         "weiss_score","height","weight","age_began_smoking_in_years","number_pack_years_smoked","karnofsky_performance_score", "total_aor_lnr",
                         "days_to_new_tumor_event_after_initial_treatment", "fluorescence_in_situ_hybridization_diagnostic_procedure_chromosome_17_signal_result_range",
                         "her2_neu_and_centromere_17_copy_number_analysis_input_total_number_count","days_to_performance_status_assessment")

## Repeat tests of clinical associations with the five methods
pValue.clinicals <- vector("list", length = 5)
names(pValue.clinicals) <- c("ESTIMATE", "ABSOLUTE", "DPC", "Pathology","Consensus")

for(i in names(pValue.clinicals)){
  pValue.clinical <- matrix(NA, nrow = length(clinical_features), ncol = length(cancerTypes))
  colnames(pValue.clinical) <- cancerTypes
  rownames(pValue.clinical) <- clinical_features
  
  for(cancer.type in cancerTypes){
    clinical <- clinical_data[[cancer.type]]
    rownames(clinical) <- clinical$bcr_patient_barcode
    # Remove features with less than 100 samples
    tmp = apply(clinical, 2, function(x) sum(x %in% nonsense))
    features <- names(which(tmp < (nrow(clinical)-100)))
    # Exclude clinical related to survival and profile information
    features <- setdiff(features, union(survival_features, id_features))
    # Separate into continuous and discrete clinical features
    features.cont <- intersect(features, continuous_features)
    features.disc <- setdiff(features, features.cont)
    
    id <- names(which(totalPurities[,"CancerType"] == cancer.type))
    purity <- as.numeric(totalPurities[id, i])
    names(purity) <- id
    purity <- purity[!is.na(purity)]
    
    # For continuous features, test using Spearman coefficient
    for(j in features.cont){
      samples <- intersect(names(purity), clinical$bcr_patient_barcode[which(!(clinical[,j] %in% nonsense))])
      if(length(samples) > minSampleSize)
        pValue.clinical[j,cancer.type] <- cor.test(purity[samples], as.numeric(clinical[samples, j]), method = "spearman", exact = TRUE, use = "na.or.complete")$p.value
    }
    # For binary or categorical features, test using analysis of variance (ANOVA)
    for(k in features.disc){
      samples <- intersect(names(purity), clinical$bcr_patient_barcode[which(!(clinical[,k] %in% nonsense))])
      if(length(samples) > minSampleSize) {
        if(length(table(clinical[samples,k])) < 2)
          next
        else{
          tmp <- data.frame(Purity = purity[samples], Clinical = clinical[samples, k])
          pValue.clinical[k,cancer.type] <- (summary(aov(Purity~Clinical, data=tmp))[[1]]$`Pr(>F)`)[1]
          # plot(Purity~Clinical, data = tmp)
        }
      }
    }
  }
  
  # Remove features with NAs for all cancer types
  pValue.clinical <- pValue.clinical[which(apply(pValue.clinical,1,function(x) sum(is.na(x))) < length(cancerTypes)), ]
  # Apply multiple test correction
  pValue.clinicals[[i]] <- apply(pValue.clinical, 2, function(x) p.adjust(x, method = "fdr"))
  write.csv(pValue.clinicals[[i]], file = paste("clinical_",i,".csv", sep = ""))
}

threshold <- 0.01
pValue.consensus <- pValue.clinicals$Consensus
apply(pValue.consensus, 2, function(x) sum(x < threshold, na.rm = TRUE))


## Create boxplots of purity levels in significant features
# The following is an example for drawing boxplots of histological subtypes
feature <- "histological_type"
types <- names(which(pValue.consensus[feature,] < threshold))
for(i in types) {
  clinical <- clinical_data[[i]]
  rownames(clinical) <- clinical$bcr_patient_barcode
  id <- names(which(totalPurities[,"CancerType"] == i))
  purity <- as.numeric(totalPurities[id, "Consensus"])
  names(purity) <- id
  purity <- purity[!is.na(purity)]
  samples <- intersect(names(purity), clinical$bcr_patient_barcode[which(!(clinical[,feature] %in% nonsense))])

  # Don't include categories that have less than 10 samples
  tmp0 = samples[which(clinical[samples, feature] %in% names(which(table(clinical[samples,feature]) >= 10)))]
  tmp <- data.frame(Purity = purity[tmp0], Clinical = clinical[tmp0, feature])
  
  boxplot(Purity~Clinical, data = tmp, main = paste(i, ", p = ", signif(pValue.consensus[feature, i], digits = 2), sep = ""), 
          xlab =feature, ylab = "Consensus Purity", ylim = c(0,1), pch = 20, cex = 0.5)
}


## Survival analysis
library("survival")

pValue.survival <- matrix(NA, nrow=length(cancerTypes), ncol = 5)
rownames(pValue.survival) <- cancerTypes
colnames(pValue.survival) <- c("ESTIMATE", "ABSOLUTE", "DPC", "Pathology","Consensus")

for(cancer.type in cancerTypes){
  id <- names(which(totalPurities[,"CancerType"] == cancer.type))
  clinical <- clinical_data[[cancer.type]]
  
  SurvTime <- as.numeric(clinical$days_to_death)
  SurvTime[is.na(SurvTime)] <- as.numeric(clinical$days_to_last_followup[is.na(SurvTime)])
  names(SurvTime) <- clinical$bcr_patient_barcode
  Status <- clinical$vital_status == "Dead"
  names(Status) <- clinical$bcr_patient_barcode
  
  for(i in colnames(pValue.survival)){
    purity <- as.numeric(totalPurities[id, i])
    names(purity) <- id
    purity <- purity[!is.na(purity)]
    samples <- intersect(names(purity), clinical$bcr_patient_barcode)
    if(length(samples) > minSampleSize) {
      res.cox <- summary(coxph(Surv(SurvTime[samples], Status[samples]) ~ purity[samples]))
      pValue.survival[cancer.type, i] = res.cox[["waldtest"]]["pvalue"]
    }
  }
}
# write.csv(pValue.survival, file = "purity_prognosis.csv", quote = FALSE)

## plot Kaplan-Meier curve
cancer.type <- names(which(pValue.survival[,"Consensus"] < threshold))  # only SKCM is significant
id <- names(which(totalPurities[,"CancerType"] == cancer.type))
clinical <- clinical_data[[cancer.type]]

SurvTime <- as.numeric(clinical$days_to_death)
SurvTime[is.na(SurvTime)] <- as.numeric(clinical$days_to_last_followup[is.na(SurvTime)])
names(SurvTime) <- clinical$bcr_patient_barcode
Status <- clinical$vital_status == "Dead"
names(Status) <- clinical$bcr_patient_barcode

purity <- as.numeric(totalPurities[id, "Consensus"])
names(purity) <- id
purity <- purity[!is.na(purity)]
# Use 1st and 4th quantile to plot 
high <- names(which(purity >= quantile(purity, probs = 0.75)))
low <- names(which(purity <= quantile(purity, probs = 0.25)))

samples <- intersect(union(high,low), clinical$bcr_patient_barcode)
purity.binary <- rep("low", length(samples))
names(purity.binary) <- samples
purity.binary[intersect(high, samples)] <- "high"

surv.purity <- survfit(Surv(SurvTime[samples], Status[samples]) ~ strata(purity.binary))
par(mfrow = c(1,1))
plot(surv.purity, lty = c(1,3), xlab = "Time (days)", ylab="Survival Probability", main = paste(cancer.type,", p =", round(pValue.survival[cancer.type,"Consensus"], digits = 3)))
legend(5000, 1, c("high purity","low purity"), lty = c(1,3), bty = "n")


## Find most appropriate method for each cancer type
scoreMatrix <- matrix(0, nrow = length(cancerTypes), ncol = 4)
colnames(scoreMatrix) <- c("ESTIMATE","ABSOLUTE","DPC","Pathology")
rownames(scoreMatrix) <- cancerTypes

# View survival as one clinical feature
types <- cancerTypes[which(pValue.survival[,"Consensus"] < threshold)]
for(t in types) {
  methods <- which(pValue.survival[t, 1:4] <= pValue.survival[t,"Consensus"])
  if(length(methods) > 0){
    for(m in methods)
      scoreMatrix[t, m] <- scoreMatrix[t, m]+1
  }
}

# Other clinical features
pValue.clinical <- pValue.clinicals$Consensus
features <- rownames(pValue.clinical)
for(cancer.type in cancerTypes) {
  features.signif <- features[which(pValue.clinical[, cancer.type] < threshold)]
  for(m in c("ESTIMATE","ABSOLUTE","DPC","Pathology")) {
    tmp <- pValue.clinicals[[m]][,cancer.type]
    tmp <- tmp[!is.na(tmp)]
    features.signif2 <- intersect(features.signif, names(tmp))
    if(length(features.signif2) > 0){
      for(f in features.signif2) {
        if(tmp[f] <= pValue.clinical[f, cancer.type])
          scoreMatrix[cancer.type, m] <- scoreMatrix[cancer.type, m] + 1
      }
    }
  }
}

