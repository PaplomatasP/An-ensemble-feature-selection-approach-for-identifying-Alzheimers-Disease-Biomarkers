library(caret)

dataGSE103334=dataGSE103334[,-ncol(dataGSE103334)]
ML_model=function(data,partition,Labels) {
  data$labels=as.factor(Labels)
  
  partionData=createDataPartition(data$labels,p=partition,list=FALSE)
  
  
  testData<<-data[-partionData,]
  trainData=data[partionData,]
  
  
  Model=train(labels ~ .,data=trainData ,method = "rf",metric = "Accuracy",
              preProc = c("center","scale"),trControl = trainControl(method = "cv",
                                                                     number = 10, classProbs = TRUE),na.action = na.omit )
  

 
  return(Model)
}
Labels=metadata1[,2]
dataGSE103334=dataGSE103334[,-ncol(dataGSE103334)]
initialModel=ML_model(dataGSE103334 ,0.8,Labels)
VarimpRforest=varImp(initialModel)
df_imps1 <- VarimpRforest[["importance"]][1]
df_imps1 <- subset(df_imps1, df_imps1[,1]>0)
df_imps1 = data.frame(df_imps1[order(df_imps1, decreasing = TRUE), drop = FALSE, ])

saveRDS(initialModel,"initialModel.rds")
saveRDS(df_imps1,"RF variable importance.rds")
