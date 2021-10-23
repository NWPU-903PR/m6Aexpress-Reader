bindsites_maplong_tr <- function(binding_sites,annotation_file,parclip=TRUE){
  ##obtain the longest transcript from transcript annotation file
  txdbfile <- GenomicFeatures::makeTxDbFromGFF(annotation_file)
  genes_txdb <- genes(txdbfile)
  exbytx_txdb <- exonsBy(txdbfile,by = "tx")
  isoform_ambiguity_method = "longest_tx"
  if(isoform_ambiguity_method == "longest_tx"){
    Longest_tx_table <- find_longest_transcript(exbytx_txdb,txdbfile)
    Kept_tx_indx <- Longest_tx_table$TXID[Longest_tx_table$longest]
    rm(Longest_tx_table)
  } else {
    Kept_tx_indx <- T
  }
  exbytx_txdb <- exbytx_txdb[Kept_tx_indx]
  exbytx_txdb <- exbytx_txdb[countOverlaps(exbytx_txdb,exbytx_txdb) == 1]
  if(parclip==TRUE){
    names(binding_sites) <- c("high_confi_TC","bind_cluster_infor")
    bind_cluster <- binding_sites$bind_cluster_infor
    highconfi_TC_GR <- binding_sites$high_confi_TC
    highconfi_TC <- as.data.frame(highconfi_TC_GR)
    target_sites_TC <-  GRanges(seqnames = as.character(highconfi_TC$seqnames),
                                IRanges(start = as.numeric(as.character(highconfi_TC$start)),
                                        width = 1),strand = as.character(highconfi_TC$strand))
    
    mcols(target_sites_TC) <- mcols(highconfi_TC_GR)
    bind_cluster_GR <- GRanges(seqnames = as.character(bind_cluster$seqnames),
                               IRanges(start = as.numeric(as.character(bind_cluster$start)),
                                       end = as.numeric(as.character(bind_cluster$end))),strand = as.character(bind_cluster$strand))
    ind_tr <- findOverlaps(target_sites_TC,exbytx_txdb,type = "within")
    select_label <- as.numeric(as.character(unique(ind_tr@from)))
    select_TC_GR <- target_sites_TC[select_label]
    ##name gene id
    geneoverlap_ID <- findOverlaps(bind_cluster_GR,genes_txdb,type="within")
    overlap_gene <- genes_txdb[unique(geneoverlap_ID@to)]
    overlap_genename <- as.character(overlap_gene$gene_id)
    add_genename_bindcluster <- data.frame()
    for (i in 1:length(overlap_genename)) {
      one_overlap <- findOverlaps(bind_cluster_GR,genes_txdb[genes_txdb$gene_id==overlap_genename[i]],type="within")
      one_bind_cluster <- cbind(overlap_genename[i], bind_cluster[unique(one_overlap@from),])
      colnames(one_bind_cluster)[1] <- "gene_name"
      add_genename_bindcluster <- rbind(add_genename_bindcluster,one_bind_cluster)
    }
    
    bindclusterGR <- GRanges(seqnames = as.character(add_genename_bindcluster$seqnames),
                             IRanges(start = as.numeric(as.character(add_genename_bindcluster$start)),
                                     end = as.numeric(as.character(add_genename_bindcluster$end))),
                             strand = as.character(add_genename_bindcluster$strand))
    mcols(bindclusterGR) <- add_genename_bindcluster[,c(1,8:ncol(add_genename_bindcluster))]
    bindclusterGR <- bindclusterGR[countOverlaps(bindclusterGR,bindclusterGR) == 1]
    new_bindcluster <- as.data.frame(bindclusterGR)
    selectcluster_ind <- findOverlaps(select_TC_GR,bindclusterGR,type = "within")
    lastTC_GR <- select_TC_GR[unique(selectcluster_ind@from)]
    select_bindcluster <- new_bindcluster[unique(selectcluster_ind@to),]
    select_clusterGR <- GRanges(seqnames = as.character(select_bindcluster$seqnames),
                                IRanges(start = as.numeric(as.character(select_bindcluster$start)),
                                        end = as.numeric(as.character(select_bindcluster$end))),
                                strand = as.character(select_bindcluster$strand),gene_name=as.character(select_bindcluster$gene_name))
    
    
    last_genename <- unique(as.character(select_bindcluster$gene_name))
    ##
    onegeneTC_GR <- function(gene_names,TC_GR,cluster_GR){
      last_overlap <- findOverlaps(TC_GR,cluster_GR[cluster_GR$gene_name==gene_names],type="within")
      onegene_TC <- TC_GR[unique(last_overlap@from)]
      mcols(onegene_TC) <- data.frame(mcols(onegene_TC),gene_name=gene_names)
      return(onegene_TC)
    }
    lastTC_GRList <- lapply(X=last_genename,FUN = onegeneTC_GR,TC_GR=lastTC_GR,cluster_GR=select_clusterGR)
    names(lastTC_GRList) <- last_genename
    last_TC_GRList <- GRangesList(lastTC_GRList)
    bind_TCsite_filter <- list(select_bindcluster,last_TC_GRList)
    names(bind_TCsite_filter) <- c("select_bind_cluster","select_TC_sites")
    return(bind_TCsite_filter)
  }
  if(parclip==FALSE){
    ##map to the longest transcript 
    ind_tr <- findOverlaps(binding_sites,exbytx_txdb,type = "within")
    select_label <- as.numeric(as.character(unique(ind_tr@from)))
    select_sites_GR <- binding_sites[select_label]
    ##name gene id
    geneoverlap_ID <- findOverlaps(select_sites_GR,genes_txdb,type="within")
    overlap_gene <- genes_txdb[unique(geneoverlap_ID@to)]
    overlap_genename <- as.character(overlap_gene$gene_id)
    ##add gene name to sites
    onegenesites_GR <- function(gene_names,site_GR,genes_txdb){
      last_overlap <- findOverlaps(site_GR,genes_txdb[genes_txdb$gene_id==gene_names],type="within")
      onegene_sites <- site_GR[unique(last_overlap@from)]
      mcols(onegene_sites) <- data.frame(mcols(onegene_sites),gene_name=gene_names)
      return(onegene_sites)
    }
    lastsites_GRList <- lapply(X=overlap_genename,FUN = onegenesites_GR,site_GR=select_sites_GR,genes_txdb=genes_txdb)
    names(lastsites_GRList) <- overlap_genename
    last_TC_GRList <- GRangesList(lastsites_GRList)
    return(last_TC_GRList)
  }
}
