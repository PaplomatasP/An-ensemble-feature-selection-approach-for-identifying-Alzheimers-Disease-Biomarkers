
PvalueCalc=function(data){

object_BPSC <- as.matrix(SummarizedExperiment::assay(data, "normcounts"))

controllds <- which(data$label == levels(factor(data$label)[1]))
design <- model.matrix(~ data$label)
resbp <- BPSC::BPglm(data = object_BPSC, controlIds = controllds,
                     design = design, coef = 2, estIntPar = FALSE, useParallel =TRUE )
FDR <- p.adjust(resbp$PVAL, method = "BH")
result_BPSC <- list(gene_names = names(resbp$PVAL), pvalue = resbp$PVAL, FDR = FDR)
return(result_BPSC)

}

StatisticalPvalueFilter = function(data,Labels, threshold) {
  options(scipen = 999)
  
  #Labels=data[,ncol(data)]
  data = as.data.frame(t(data[,-ncol(data)])  )
  data1 = data_fitting(data, Labels)
  
  
  
  PvalueData = PvalueCalc(data1)
  
PvalueData  <- as.data.frame(PvalueData)

  rownames(PvalueData)=PvalueData[,1]
  count = 0
  PvalueTreshold = list()
  for (i in seq(1, nrow(PvalueData))) {
    if (PvalueData[, 2][i] <= threshold) {
      count = count + 1
      
      PvalueTreshold[[count]] = rownames(PvalueData)[i]
    }
  }
  PvalueTreshold = as.data.frame(cbind(PvalueTreshold))
  row.names(PvalueTreshold) = PvalueTreshold[, 1]
  
  PvalueData1 <- data.frame(PvalueData[order(PvalueData$pvalue, decreasing = FALSE), drop = FALSE, ])
  iG1 = subset(PvalueData1 , rownames(PvalueData) %in% rownames(PvalueTreshold))
  iG = cbind((iG1[,-c(1,3)]))
  rownames(iG) = rownames(iG1)
  
  iG<-as.data.frame(round(iG,digits = 6) )
  
  # neudata = RDS_file1[,colnames(RDS_file1) %in% row.names(PvalueTreshold) ]
  # 
  # neudata = as.data.frame(neudata)
  # 
  # neudata$Labels = as.factor(Labels)


return(iG)
}
StatistData=StatisticalPvalueFilter(dataGSE103334,Labels,0.05)
saveRDS(StatistData,"StatisticalGenes.rds")
