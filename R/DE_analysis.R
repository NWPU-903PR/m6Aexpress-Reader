##differential expression
DE_analysis <- function(gene_readscount,con_name1,con_name2,num_con1,num_con2,p_value=0.05){
  colnames(gene_readscount)[2:ncol(gene_readscount)] <- c(paste0(con_name1,1:num_con1), paste0(con_name2,1:num_con2))
  rownames(gene_readscount) <- NULL
  gene_readscount <- gene_readscount[which(rowSums(gene_readscount[,-1])>10),]
  #############check the DE infor
  gene_name <- as.character(gene_readscount$gene_name)
  gene_reads_count <- gene_readscount[,-1]
  rownames(gene_reads_count) <- gene_name
  coldata <- data.frame(condition = c(rep(con_name1,num_con1),rep(con_name2,num_con2)))
  dds <- DESeqDataSetFromMatrix(countData = gene_reads_count,
                                colData = coldata,
                                design = ~ condition)
  
  dds <-  DESeq(dds)
  dds <- estimateSizeFactors(dds)
  size_factor <- sizeFactors(dds)
  res <- results(dds)
  resLFC <- lfcShrink(dds, coef=2)
  
  DE_result <- as.data.frame(resLFC)
  DE_result <- data.frame(gene_name=gene_name,DE_result)
  rownames(DE_result) <- NULL
  sigDE_result <- DE_result[DE_result$pvalue<p_value,]
  sigDE_genecount <- gene_readscount[which(!is.na(match(gene_readscount$gene_name,sigDE_result$gene_name))),]
  sigDE_infor <- list(sigDE_genecount,sigDE_result,size_factor)
  names(sigDE_infor) <- c("sigDE_genecount","sigDE_geneinfor","size_factor")
  return(sigDE_infor)
}
