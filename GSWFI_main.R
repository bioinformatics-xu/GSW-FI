library(data.table)
library(plyr)
library(ggplot2)
library(stringr)
library(MASS)
library(gamlss)

lambda <- 0.5

maf.file.name <- 'brac.maf.txt'

print(maf.file.name)

M <- as.data.frame(fread(file = maf.file.name)) 
cancer.name <- toupper(str_split(maf.file.name, pattern = '_')[[1]][1]) 
dict <- as.data.frame(fread("mutation_type_dictionary_file.txt"))

#GENE
if (("gene" %in% colnames(M))&("Hugo_Symbol" %in% colnames(M))){
  cat("NOTE: Both 'gene' and 'Hugo_Symbol' are present in mutation_file. Using 'gene'. \n")
}else if("Hugo_Symbol" %in% colnames(M)){
  M$gene <- M$Hugo_Symbol
}else if("gene" %in% colnames(M)){
  
}else{
  stop("mutation_file lacks 'gene' or 'Hugo_Symbol' column.")
}


if("gene" %in% colnames(M)){
  if((str_sub(M$gene[1],1,4) == "ENSG")&("Hugo_Symbol" %in% colnames(M))){
    M$gene <- M$Hugo_Symbol
  }
}

#PATIENT
if (("patient" %in% colnames(M))&("Tumor_Sample_Barcode" %in% colnames(M))){
  cat("NOTE: Both 'patient' and 'Tumor_Sample_Barcode' are present in mutation_file. Using 'patient'. \n")
}else if("patient" %in% colnames(M)){
  #ok
}else if("Tumor_Sample_Barcode" %in% colnames(M)){
  M$patient <- M$Tumor_Sample_Barcode
}else{
  stop("mutation_file lacks 'patient' or 'Tumor_Sample_Barcode' column.")
}

if (length(unique(M$patient)) < 2){
  stop("identification is not applicable to single patients. \n")
}

#Variant_Classification
if(!("Variant_Classification" %in% colnames(M))&("type" %in% colnames(M))){
  cat('NOTE:  Both Variant_Classification and type are present in mutation_file. Using Variant_Classification \n')
  M$Variant_Classification <- M$type
}else if("Variant_Classification" %in% colnames(M)){
  #ok
}else if("type" %in% colnames(M)){
  M$Variant_Classification <- M$type
}else{
  stop("mutation.file is missing Variant_Classification")
}

#Chromosome
if(!("chr" %in% tolower(colnames(M))) & "chromosome" %in% tolower(colnames(M))){
  colnum <- which(tolower(colnames(M)) %in% "chromosome")
  M$chr <- M[[colnum]]
}
if('start' %in% tolower(colnames(M)) & 'end' %in% tolower(colnames(M))){
  colnum <- which(tolower(colnames(M)) %in% "start")
  M$start <- M[[colnum]]
  
  colnum <- which(tolower(colnames(M)) %in% "end")
  M$end <- M[[colnum]]
}else if("start_position" %in% tolower(colnames(M)) & "end_position" %in% tolower(colnames(M))){
  colnum <- which(tolower(colnames(M)) %in% "start_position")
  M$start <- M[[colnum]]
  
  colnum <- which(tolower(colnames(M)) %in% "end_position")
  M$end <- M[[colnum]]
}else{
  if(!("start" %in% tolower(colnames(M))) & "start_position" %in% tolower(colnames(M))){
    colnum <- which(tolower(colnames(M)) %in% "start_position")
    M$start <- M[[colnum]]
  }
  if(!("end" %in% tolower(colnames(M))) & "end_position" %in% tolower(colnames(M))){
    colnum <- which(tolower(colnames(M)) %in% "end_position")
    M$end <- M[[colnum]]
  }
}

if(!("ref_allele" %in% colnames(M)) & "reference_allele" %in% tolower(colnames(M))){
  colnum <- which(tolower(colnames(M)) %in% "reference_allele")
  M$ref_allele <- M[[colnum]]
}
if(!("newbase" %in% tolower(colnames(M)))){
  if("tumor_seq_allele1" %in% tolower(colnames(M))){
    colnum <- which(tolower(colnames(M)) %in% "tumor_seq_allele1")
    M$newbase <- M[[colnum]]
    if("tumor_seq_allele2" %in% tolower(colnames(M))){
      colnum.allele1 <- which(tolower(colnames(M)) %in% "tumor_seq_allele1")
      idx <- which(M$ref_allele == M[[colnum.allele1]])
      colnum.allele2 <- which(tolower(colnames(M)) %in% "tumor_seq_allele2")
      M$newbase[idx] <- M[[colnum.allele2]][idx]
    }
  }
}

if(length(which(!(M$Variant_Classification %in% unique(dict$Variant_Classification))))>0){
  M <- M[-which(!(M$Variant_Classification %in% unique(dict$Variant_Classification))),]
}


if ("effect" %in% colnames(M)){
  # cat("Will use the pre-existing effect column.\n")
  M$effect <- gsub("^flank.*","noncoding",M$effect)
  M_effectnames <- c("noncoding","silent","nonsilent","null")
  if (!any(unique(M$effect) %in% M_effectnames)){
    stop("in mutation_file, 'effect' must be one of
                 noncoding/silent/nonsilent/null")
  }
}else{
  if(!("Variant_Classification" %in% colnames(M))&("type" %in% colnames(M))){
    # cat("after if variant_classification in colnames ")
    M$Variant_Classification <- M$type
  }
  if(!("Variant_Classification" %in% colnames(M))){
    stop("mutation_file is missing Variant_Classification")
  }
  flag_num <- match(toupper(M$Variant_Classification),
                    toupper(dict$Variant_Classification),nomatch =nrow(dict)+1)
  dict <- rbind(dict,dict[nrow(dict),])
  dict[nrow(dict)+1,] <- "unknown"
  
  M$effect <- dict$effect[flag_num]
  bad <- which(M$effect == "unknown")
  if(length(bad)>0){
    cat(sprintf("WARNING: %d/%d mutations could not be mapped to
                            effect using mutation_type_dictionary_file:\n",
                length(bad),length(M$effect)))
    #table(bad_variant <- M$Variant_Classification[bad])
    cat("They will be removed from the analysis. \n")
    M <- M[-bad]
  }
  if(nrow(M) == 0){
    stop("No mutations left!")
  }
}

# take the number of non-synonymous mutations for genes
M$nonsyeffect <- 0
M$nonsyeffect[which(M$effect %in% c("null","nonsilent"))] <- 1
nonsyeffect <- as.matrix(tapply(M$nonsyeffect,M$gene,sum)) 

#take the mutation number of every gene as an explaination variable
M$alleffect <- 1
alleffect <- as.matrix(tapply(M$alleffect,M$gene,sum)) 

# take the ratio of null and nonsilent mutation as an explaination variable
nonsyeffect_rate <- nonsyeffect/alleffect

#tabke the null and nonsilent mutation ratio and mutation number of genes
gene.effect <- data.frame(gene=rownames(nonsyeffect_rate), totalMutNum=alleffect,seriMutNum=nonsyeffect)

mutationdata <- subset(M,select = c("gene","chr","start","end","Variant_Classification","Variant_Type","ref_allele","newbase","patient","effect"))
mutationdata$mutation <- paste("hg19",mutationdata$chr,mutationdata$start,mutationdata$ref_allele,mutationdata$new,sep = ",")
# take the mutation of maf into the form of mutation accessor

mutationdata <- mutationdata[which(mutationdata$gene != "TTN"),]
# remove the TNN gene

# list all the files of mutation accessor
MA.scores.direct <- "./MA_scores_rel3_hg19_full"
MA.scores.filename <- list.files(MA.scores.direct)

# set the gene effect silent noncoding null nonsilent
mutation_score_dictionary <- fread(file = "mutation_type_dictionary_file.txt")

# build a data.frame mutationdata.score to store the mutaiton score
mutationdata.score <- as.data.frame(matrix(data = NA,ncol = ncol(mutationdata)+2))
colnames(mutationdata.score) <- c(colnames(mutationdata),"mutation","score")
mutationdata.score <- mutationdata.score[-1,]

# build a data.frame MA.score.effect to store the mutation score of MutationAccessor
MA.score.effect <- as.data.frame(matrix(data = NA,ncol = ncol(mutationdata)+2))
MA.score.effect <- MA.score.effect[-1,]

for(i in c(1:22,"X","Y","M")){
  # take the chr i to deal with the mutation score for every mutation
  print(i)
  # take the mutation score of chromosome i
  mutationdata.i <- mutationdata[which(mutationdata$chr == i),]
  
  if(nrow(mutationdata.i) > 0){
    # take the mutation data of mutationdata.i into the form of MutationAccessor
    mutationdata.i$mutation <- paste("hg19",i,mutationdata.i$start,mutationdata.i$ref_allele,mutationdata.i$new,sep = ",")
    
    MA.scores.filename.i <- paste(MA.scores.direct,"/","MA_scores_rel3_hg19_chr",i,"_full.csv",sep = "")
    MA.scores <- fread(file = MA.scores.filename.i)
    # read the mutation data of MutationAccessor
    
    # assign mutation score for missense mutation using the MutationAccessor FI score
    mutation.intersect <- intersect(mutationdata.i$mutation,MA.scores$Mutation)
    mutationdata.filter <- mutationdata.i[mutationdata.i$mutation %in% mutation.intersect,]
    MA.scores.filter <- MA.scores[MA.scores$Mutation %in% mutation.intersect,]
    flag.match <- match(mutationdata.filter$mutation,MA.scores.filter$Mutation)
    mutationdata.filter$score <- MA.scores.filter$`FI score`[flag.match]
    
    MA.score.effect <- rbind(MA.score.effect,mutationdata.filter)
    # take the effective  mutation score of MutationAccessor, and store in MA.score.effect
  }   
}


if(nrow(MA.score.effect) < nrow(M)*0.5){
  bad <- nrow(M) - nrow(MA.score.effect)
  cat(sprintf("WARNING: %d/%d mutations could not be mapped to
                    mutation functional impact score using mutation_accessor_file:\n",
              bad,nrow(M)))
  #stop("The maf file does not match with the mutation impact file of hg19")
}

MAscoreIndex$Mrow[maf] <- nrow(M)
MAscoreIndex$MArow[maf] <- nrow(MA.score.effect)

# remove the mutation data of NA in MA.score.effect
effect.score <- MA.score.effect[which(!is.na(MA.score.effect$score)),]

# match the mutation of mutationdata with the mutation of effect.score
dict <- as.data.frame(fread("mutation_type_dictionary_file.txt"))
flag_num <- match(toupper(mutationdata$mutation),
                  toupper(effect.score$mutation),nomatch = nrow(effect.score)+1)

# assign mutation FI score to mutationdata
MA.score <- c(effect.score$score,NA)
mutationdata$MAscore <- MA.score[flag_num]


# calculate the mean mutation FI score of every mutation effect
mutationdata_ex_NA <- mutationdata[which(!is.na(mutationdata$MAscore)),]
mutEffectMean <- tapply(mutationdata_ex_NA[,"MAscore"], mutationdata_ex_NA[,"effect"], mean)

mutEffectScore <- data.frame(effect=c("noncoding","nonsilent","null","silent"),
                             score=c(1,2,3,0)) 
flag_effect <- mutEffectScore$effect %in% names(mutEffectMean)

# assign the mutation FI score of effect which is not in the MutationAccessor file
if(sum(flag_effect) < 4){
  mutEffectMean <- data.frame(effect=c(names(mutEffectMean),
                                       as.vector(mutEffectScore$effect)[!flag_effect]),
                              score=c(mutEffectMean,as.vector(mutEffectScore$score)[!flag_effect]))
}else{
  mutEffectMean <- data.frame(effect=names(mutEffectMean),score=mutEffectMean)
}

# assign the mutation FI score of NA in mutationdata by the mean of mutation effect
mutationdata_NA <- mutationdata[which(is.na(mutationdata$MAscore)),]
flag_match_effect <- match(mutationdata_NA$effect,mutEffectMean$effect)
mutationdata_NA$MAscore <- mutEffectMean$score[flag_match_effect]

# rbind the mutation score with NA and without NA
mutationdata_new <- rbind(mutationdata_ex_NA,mutationdata_NA)

allGenes <- unique(mutationdata_new$gene)
allPatients <- unique(mutationdata_new$patient)
gePatMat <- as.data.frame(matrix(data = 0,nrow = length(allGenes),
                                 ncol = length(allPatients)))
dimnames(gePatMat) <- list(allGenes,allPatients)
for (j in 1:length(allPatients)) {
  patName <- colnames(gePatMat)[j]
  patFI <- mutationdata_new[which(mutationdata_new$patient %in% patName),]
  
  geneFIscore <- as.data.frame(tapply(patFI$MAscore, patFI$gene, sum))
  existFlag <- match(rownames(gePatMat),rownames(geneFIscore))
  existFlag1 <- which(!is.na(existFlag))
  gePatMat[existFlag1,j] <- geneFIscore[existFlag[existFlag1],1]
}


sdFI <- as.data.frame(apply(gePatMat,1,sd)) 
sdFI$gene <- rownames(sdFI)
sdFI <- sdFI[order(sdFI$`apply(gePatMat, 1, sd)`,decreasing = T),]
colnames(sdFI) <- c("sd","gene")

meanFI <- as.data.frame(apply(gePatMat,1,mean)) 
meanFI$gene <- rownames(meanFI)
meanFI <- meanFI[order(meanFI$`apply(gePatMat, 1, mean)`,decreasing = T),]
colnames(meanFI) <- c("mean","gene")

meanSdFI <- merge(sdFI,meanFI,by="gene")
meanSdFI <- meanSdFI[order(meanSdFI$sd,decreasing = T),]
meanSdFI <- meanSdFI[order(meanSdFI$mean,decreasing = T),]


gene.score <- ddply(mutationdata_new, "gene", summarise, sum(MAscore))
colnames(gene.score) <- c("gene","FI")

# take the exp of mutation FI score
mutationdata_new$expScore <- exp(mutationdata_new$MAscore)
# scale the exp mutation FI score 
mutationdata_new$scaleExpScore <- mutationdata_new$expScore/max(mutationdata_new$expScore)

# calculate the total FI score of every gene
gene.score <- ddply(mutationdata_new, "gene", summarise, sum(MAscore))
# order the FI score of all genes
gene.score.order <- gene.score[order(gene.score[,2],decreasing = T),]

colnames(gene.score) <- c("gene","mutationScore")
quantile(gene.score$mutationScore)

# take all the features of GLM model
gene.feature <- as.data.frame(fread("geneCominedFeature.txt"))
colnames(gene.feature)[1] <- "gene"
gene.feature.score <- merge(gene.score,gene.feature,by="gene")
gene.feature.score <- merge(gene.feature.score,gene.effect,by="gene")
gene.feature.score <- merge(gene.feature.score,sdFI,by="gene")

rownames(gene.feature.score) <- gene.feature.score$gene
gene.feature.score <- gene.feature.score[,-1]

# build a GLM model
glm.model <- gamlss(formula = mutationScore ~.,sigma.formula = ~1,
                    family = NO(mu.link = "identity",sigma.link = "identity"),
                    data = gene.feature.score)
summary(glm.model)

# take the explaination data
explainationData <- matrix(data = rep(1,times=nrow(gene.feature.score)))
explainationData <- cbind(explainationData,gene.feature.score[,2:ncol(gene.feature.score)])

# calculate the mean mutation FI score of every gene as background FIS
mu.coeff <- matrix(data = glm.model$mu.coefficients,ncol = 1)
mu.g <- as.matrix(explainationData) %*% mu.coeff

sigma.g <- glm.model$sigma.coefficients

gene.feature.score$estimateFI <- mu.g[,1]


# Building the FIS circle
V <- subset(gene.feature.score,select = c("AvgCodingLen","expr","reptime","hic",
                                          "constraint_score","totalMutNum",
                                          "seriMutNum","sd","Expression_hub","Regulator",
                                          "Genomic_abnormalities_CNV","meEthylation"))
Z <- scale(V)
ng <- nrow(gene.feature.score)
nv <- ncol(V)
neighbour_size <- 100
qual_min <- 0.1
gene.feature.score$finalFI <- NA
gene.feature.score$w_nei_score <- NA
gene.feature.score$neighborNum <- NA

for (g in 1:nrow(gene.feature.score)){
  if (g %% 1000 == 0) print(g)
  df2= (Z-matrix(rep(Z[g,],ng),nrow=ng,ncol=nv,byrow = TRUE))^2
  dfz  = df2
  dfz[is.nan(dfz)] <-0
  dist2 <- apply(dfz,1,sum)/apply(!is.nan(df2),1,sum)
  dist2 <- round(dist2,15)
  ord <- order(dist2,decreasing = FALSE,na.last = TRUE)
  ord <- c(g,ord[ord!=g])
  
  targx <- ord[1]
  
  neighbour <- vector()
  
  for (ni in 1:(nrow(gene.feature.score)-1))
  {
    #if(ni %% 1000 == 0) print(ni)
    gidx = ord[ni+1]
    p_norm <- pnorm(abs(gene.feature.score$mutationScore[targx]-gene.feature.score$mutationScore[gidx]),0,sigma.g)
    qual_left <- min(p_norm,1-p_norm)
    qual <- 2*qual_left
    #print(qual)
    
    if(qual < qual_min){
      next
    }else{
      neighbour <- c(neighbour,ord[ni+1])
    }
    
    if(length(neighbour) >= neighbour_size){
      break
    }
    
  }
  
  if(length(neighbour) == 0){
    cat(sprintf("%s has 0 neighbour ",rownames(gene.feature.score)[g]))
    gene.feature.score$finalFI[g] <- gene.feature.score$estimateFI[g]
    
  }else{
    Z_neighbour <- Z[neighbour,]
    distance_nei <- dist2[neighbour]
    FI_neighbour <- gene.feature.score$mutationScore[neighbour]
    
    #w_nei_score <- sum(FI_neighbour/sqrt(distance_nei))/sum(1/sqrt(distance_nei))
    if(length(which(distance_nei == 0)) > 0){
      flag0 <- which(distance_nei == 0)
      distance_nei <- distance_nei[-flag0]
      FI_neighbour <- FI_neighbour[-flag0]
    }            
    w_nei_score <- sum(FI_neighbour/distance_nei)/sum(1/distance_nei)
    
    final_score <- lambda*gene.feature.score$estimateFI[g] + (1-lambda)*w_nei_score
    gene.feature.score$finalFI[g] <- final_score
  }
  
  gene.feature.score$w_nei_score[g] <- w_nei_score
  gene.feature.score$neighborNum[g] <- length(neighbour)
  
} # of for (gi in 1:ng)

# calculate the correlation coefficient of predict FI score and real FI score
MAscoreIndex$cor[maf] <- cor.test(gene.feature.score$mutationScore,gene.feature.score$finalFI)$estimate[1]

# calculate the pValue using the ks.test
gene.feature.score$pValue <- NA
seriRate <- (gene.feature.score$seriMutNum)/(gene.feature.score$totalMutNum)

seriNumPatient <- exp(gene.feature.score$seriMutNum/length(unique(mutationdata$patient)))
quantile(seriRate)
quantile(seriNumPatient)

for (g in 1:nrow(gene.feature.score)){
  
  if(gene.feature.score$seriMutNum[g] == 0){
    gene.feature.score$pValue[g] <- 1
  }else{
    if(gene.feature.score$mutationScore[g] >= 0){
      weighted_FIscore <- gene.feature.score$mutationScore[g]*seriRate[g]*seriNumPatient[g]
      final_score <- lambda*gene.feature.score$estimateFI[g] + (1-lambda)*gene.feature.score$w_nei_score[g]
      gene.feature.score$finalFI[g] <- final_score
      pValue_norm <- pnorm(weighted_FIscore,final_score,sigma.g)
      gene.feature.score$pValue[g] <- 1 - pValue_norm
    }else if(gene.feature.score$mutationScore[g] < 0){
      gene.feature.score$pValue[g] <- 1
    }
  }
}


gene.feature.score$qValue <- stats::p.adjust(gene.feature.score$pValue, method="fdr",length(gene.feature.score$pValue))
gene.feature.score <- gene.feature.score[order(gene.feature.score$qValue,decreasing = F),]
gene.feature.score$gene  <- rownames(gene.feature.score)

gene.feature.score <- subset(gene.feature.score, select = c('gene', 'pValue', 'qValue'))
write.csv(gene.feature.score,file = paste(cancer.name,
                                          ".results",lambda*10,(1-lambda)*10,'.csv',sep = ""),row.names = F)



