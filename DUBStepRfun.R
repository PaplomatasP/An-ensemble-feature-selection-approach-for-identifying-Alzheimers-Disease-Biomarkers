data_fitting = function(data, group, normal.Method) {
  if (missing(normal.Method)) {
    normcounts <- data
    
  }
  else{
    normcounts <- normaliz_method(data, normal.Method)
  }
  SingleCell <-
    SingleCellExperiment::SingleCellExperiment(assays = list(counts = data, normcounts = normcounts))
  gene_df <- data.frame(Gene = rownames(SingleCell))
  cell_df <- data.frame(label = group, cell = colnames(SingleCell))
  rownames(gene_df) <- gene_df$Gene
  rownames(cell_df) <- cell_df$cell
  SingleCell <-
    SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = data, normcounts = normcounts),
      colData = cell_df,
      rowData = gene_df
    )
}
DUBStepRfun=function(data,Labels){
  obj = as.data.frame(t(data[,-ncol(data)])  )
  scale.factor <- mean(colSums(obj))
  
  obj <- Seurat::LogNormalize(obj,
                              scale.factor = scale.factor)
  
  obj=as.data.frame(obj)
  data1 = data_fitting(obj, Labels)
  object_Seurat <- as.matrix(SummarizedExperiment::assay(data1, "counts"))
  names(Labels) <- colnames(obj)
  MetaData <- data.frame(groups = Labels)
  Input <- Seurat::CreateSeuratObject(counts = object_Seurat, meta.data = MetaData)
  dubstepR.out <- DUBStepR::DUBStepR(input.data = Input@assays$RNA@data, min.cells = 0.05*ncol(Input), optimise.features = T, k =10, num.pcs = 10, error = 0)
  Input@assays$RNA@var.features <- dubstepR.out$optimal.feature.genes
  genes=Input@assays$RNA@var.features
  iG1=dubstepR.out[["corr.info"]]
  iG=as.data.frame(iG1[,-1])
  iG=as.data.frame(iG[order(iG[,1],decreasing = TRUE),])
  rownames(iG)=rownames(iG1)
  
  #newdata=data[,colnames(data) %in% genes ]
  return(iG)
  
}

data1=data1[,-caret::nearZeroVar(data1)]
Labels=metadata1$labels
VariableD=DUBStepRfun(dataGSE103334,Labels)
saveRDS(VariableD,"VariableGenes.rds")
