obtain_gene_readscount <- function(Input_data,annotation_file,isPairedEnd=F){
  gene_count <- featureCounts(Input_data,isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_name", 
                              annot.ext = annotation_file, isPairedEnd=F)
  gene_readscount <- as.data.frame(gene_count$counts)
  gene_names <- as.character(rownames(gene_readscount))
  gene_readscount <- data.frame(gene_name=gene_names,gene_readscount)
  rownames(gene_readscount) <- NULL
  return(gene_readscount)
}
