library(data.table)
library(plyr)
library(stringr)

lambda = 0.5
maf.file.name <- 'acc_tcga.maf.txt'
cancer.name <- toupper(str_split(maf.file.name, pattern = '_')[[1]][1]) 
gene.feature.score <- as.data.frame(fread(paste0(cancer.name, '_geneFeatureScore.txt')))

# calculate the pValue using the ks.test
gene.feature.score$pValue <- NA

for (g in 1:nrow(gene.feature.score)){
  
  if(gene.feature.score$seriRate[g] == 0){
    gene.feature.score$pValue[g] <- 1
  }else{
    if(gene.feature.score$mutationScore[g] >= 0){
      weighted_FIscore <- gene.feature.score$mutationScore[g]*gene.feature.score$seriRate[g]*gene.feature.score$seriNumPatient[g]
      final_score <- lambda*gene.feature.score$estimateFI[g] + (1-lambda)*gene.feature.score$w_nei_score[g]
      pValue_norm <- pnorm(weighted_FIscore,final_score,gene.feature.score$sigma[1])
      gene.feature.score$pValue[g] <- 1 - pValue_norm
    }else if(gene.feature.score$mutationScore[g] < 0){
      gene.feature.score$pValue[g] <- 1
    }
  }
}


gene.feature.score$qValue <- stats::p.adjust(gene.feature.score$pValue, method="fdr",length(gene.feature.score$pValue))
gene.feature.score <- gene.feature.score[order(gene.feature.score$qValue,decreasing = F),]
#gene.feature.score$gene <- rownames(gene.feature.score)
gene.feature.score <- subset(gene.feature.score, select = c('gene', 'pValue', 'qValue'))

write.csv(gene.feature.score,file = paste(cancer.name,
                                          ".results",lambda*10,(1-lambda)*10,'.csv',sep = ""))