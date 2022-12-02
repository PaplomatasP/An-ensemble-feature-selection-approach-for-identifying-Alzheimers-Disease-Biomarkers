library(ggvenn)
library(votesys)

# ImportanceGenes=ImportanceGenes[order(ImportanceGenes$Overall,decreasing = TRUE),,drop=FALSE]
# ImportanceGenes=as.data.frame(ImportanceGenes)

#Ranking 
raw <- c(
  rownames(`RF variable importance`)[1:200], 
  rownames(StatisticalGenes)[1:200], 
  rownames(VariableGenes)[1:200]
) 
raw1=unique.default(sapply(raw, unique))
raw <- matrix(raw, ncol = 3 , byrow = TRUE)

vote <- create_vote(raw, xtype = 2, candidate = c(raw1))
y <- borda_method(vote,modified = TRUE)

RankingResults=as.data.frame(sort(y[["other_info"]][["count_max"]],decreasing = TRUE ) ) 
RankingResultsDF=data.frame(RankingResults[1:200,])
rownames(RankingResultsDF)=rownames(RankingResults)[1:200]


#Ranking Barplot
ColorFun <- colorRampPalette( c( "#CCCCCC" , "#104E8B" ) )
ColorPaleta <- ColorFun( n = nrow( x = RankingResults ) )
RankingResults$Color <-
  as.character(
    x = cut(
      x = rank( x = RankingResults[,1])  # used to assign order in the event of ties
      , breaks = nrow( x = RankingResults)  # same as the 'n' supplied in ColorFun
      , labels = ColorPaleta  # label the groups with the color in ColorPaleta
    )
  )

barplot( height = RankingResults[,1]
         , names.arg = rownames(RankingResults)
         , las = 2
         , col = RankingResults$Color
         , border = NA  # eliminates borders around the bars
         , main = "Ranking of  Genes"
         , ylab = "Count of Votes",ylim = c(0,8)
         , xlab = "Genes ID " , cex.names=0.3
)


#Ven Circulis
set.seed(20190708)
x <- list(
  VariableGenes=rownames(VariableGenes)[1:200], 
  RF_Variable_Importance=rownames(`RF variable importance`)[1:200], 
  StatisticalGenes=rownames(StatisticalGenes)[1:200]
)


ggvenn(
  x, 
  fill_color = c("#0073C2FF","#CD534CFF" , "#EFC000FF"),
  stroke_size = 0.2, set_name_size = 1
)
library("ggVennDiagram")
p=ggVennDiagram(x, label_alpha = 0,set_size = 4,category.names = c("DUBStepR",
                                                                     "RF variable importance",
                                                                     "BPSC "))
p + #scale_fill_distiller(palette = "Reds", direction = 1)
  scale_fill_distiller(palette = "RdBu")
