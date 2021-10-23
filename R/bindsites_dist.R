##conversion TC sites
targetTC_infor <- function(target_TC,parclip=TRUE){
  if(parclip==TRUE){
    targetTC_infor <-  data.frame(seqnames=as.character(target_TC@seqnames),
                                  start=as.numeric(as.character(target_TC@ranges)),
                                  end=as.numeric(as.character(target_TC@ranges)),
                                  strand=as.character(target_TC@strand),
                                  coverage=as.numeric(as.character(target_TC$coverage)),
                                  count=as.numeric(as.character(target_TC$count)),
                                  rsf=as.numeric(as.character(target_TC$rsf)),
                                  gene_name= as.character(target_TC$gene_name))
  }
  if(parclip==FALSE){
    targetTC_infor <-  data.frame(seqnames=as.character(target_TC@seqnames),
                                start=as.numeric(as.character(target_TC@ranges)),
                                end=as.numeric(as.character(target_TC@ranges)),
                                strand=as.character(target_TC@strand),
                                score=as.numeric(target_TC$score),
                                gene_name= as.character(target_TC$gene_name))
  }
  return(targetTC_infor)
}

dist_fun <- function(overlap_bindsites_infor,bindgene_nobind_peakcenter,annotation_file,parclip=TRUE){
  txdbfile <- GenomicFeatures::makeTxDbFromGFF(annotation_file)
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
    overlap_bindsites <- overlap_bindsites_infor[[1]]
    bind_cluster_overlap <- overlap_bindsites_infor[[2]]
    bind_sites_data <- targetTC_infor(target_TC=overlap_bindsites,parclip=parclip)
    bind_cluster_GR <- GRanges(seqnames = as.character(bind_cluster_overlap$seqnames),
                               IRanges(start = as.numeric(as.character(bind_cluster_overlap$start)),
                                       end = as.numeric(as.character(bind_cluster_overlap$end))),
                               strand = as.character(bind_cluster_overlap$strand))
  }
  if(parclip==TRUE){
    bind_sites_data <- targetTC_infor(target_TC=overlap_bindsites_infor,parclip=parclip)
  }
  target_sites_se <- GRanges(seqnames = as.character(bindgene_nobind_peakcenter$seqnames),
                             IRanges(start = as.numeric(as.character(bindgene_nobind_peakcenter$start)),
                                     width = 1),strand = as.character(bindgene_nobind_peakcenter$strand))
  
  target_bindsites <-  GRanges(seqnames = as.character(bind_sites_data$seqnames),
                               IRanges(start = as.numeric(as.character(bind_sites_data$start)),
                                       width = 1),strand = as.character(bind_sites_data$strand))
  target_sites_map <- mapToTranscripts(target_sites_se, exbytx_txdb,ignore.strand=F)
  overlapbind_map <- mapToTranscripts(target_bindsites, exbytx_txdb,ignore.strand=F)
  
  overlap_hit_tx <- intersect(as.character(seqnames(target_sites_map)),
                              as.character(seqnames(overlapbind_map)))
  
  select_targetsites_map <- target_sites_map[which(!is.na(match(seqnames(target_sites_map),overlap_hit_tx)))]
  # select_targetsites_se <- target_sites_se[select_targetsites_map$xHits]
  
  
  compare_label <- list()
  abs_pos <- vector()
  norm_pos <- vector()
  mindist_TC_data <- data.frame()
  for (i in 1:length(select_targetsites_map)) {
    # site_label[[i]] <- which(!is.na(match(as.character(seqnames(target_sites_map)),overlap_hit[i])))
    # compare_label[[i]] <- which(!is.na(match(as.character(seqnames(compared_map)),overlap_hit[i])))
    compare_label[[i]] <- which(!is.na(match(as.character(seqnames(overlapbind_map)),as.character(seqnames(select_targetsites_map)[i]))))
    # select_target_sites_map <- target_sites_map[site_label]
    select_target_sites_map <- select_targetsites_map[i]
    select_compared_map <- overlapbind_map[compare_label[[i]]]
    ##obtain min distance to TC sites
    # site_pos <- start(target_sites_se[select_target_sites_map$xHits])
    site_pos <- start(select_target_sites_map)
    # compared_pos <- start(target_sites_TC[select_compared_map$xHits][which.max(overlap_TC[select_compared_map$xHits,]$rsf)])
    compared_pos <- start(select_compared_map[which.max(bind_sites_data[select_compared_map$xHits,]$score)])
    # tx_width <- as.numeric(sum(width(exbytx_txdb[which(!is.na(match(names(exbytx_txdb),overlap_hit[i])))])))
    # abs_pos[i] <- min(abs(site_pos-compared_pos))
    abs_pos[i] <- abs(site_pos-compared_pos)
    # mindist_TC <- target_sites_TC[select_compared_map[which.min(min(abs(site_pos-compared_pos)))]$xHits]
    if(parclip==TRUE){
      mindist_TC <- target_bindsites[select_compared_map$xHits][which.max(bind_sites_data[select_compared_map$xHits,]$rsf)]
      peak_num <- rep(paste0("peak",select_targetsites_map[i]$xHits),1)
      onemindist_TC <-data.frame(as.data.frame(mindist_TC),peak_num=peak_num)
      mindist_TC_data <- rbind(mindist_TC_data,onemindist_TC)  
      TC_bindcluster <- findOverlaps( mindist_TC, bindcluster_GR )
      oneselect_cluster <- bind_cluster_overlap[unique(TC_bindcluster@to),]
      peak_num2 <- rep(paste0("peak",select_targetsites_map[i]$xHits),nrow(oneselect_cluster))
      oneselect_bindcluster <- data.frame(oneselect_cluster,peak_num2)
      bind_cluster_select <-rbind(bind_cluster_select,oneselect_bindcluster)
    }
    if(parclip==FALSE){
      mindist_TC <- target_bindsites[select_compared_map$xHits][which.max(bind_sites_data[select_compared_map$xHits,]$score)]
      peak_num <- rep(paste0("peak",select_targetsites_map[i]$xHits),1)
      onemindist_TC <- data.frame(as.data.frame(mindist_TC),peak_num=peak_num)
      mindist_TC_data <- rbind(mindist_TC_data,onemindist_TC) 
    }
    ##obtain binding cluster
    # TC_bindcluster <- findOverlaps( target_sites_TC[select_compared_map$xHits], DFbindcluster_GR )
    
    # norm_pos[i] <- abs_pos[i]/tx_width
    norm_pos[i] <- abs_pos[i]
  }
  
  # site_num <- target_sites_map[site_label]$xHits
  if(parclip==TRUE){
    compared_num <- overlapbind_map[unlist(compare_label)]$xHits
    names(norm_pos) <- unique(as.character(bind_cluster_select$peak_num2))
    peaksites_compared_infor <- list(mindist_TC_data,bind_cluster_select,compared_num,norm_pos)
    names(peaksites_compared_infor) <- c("select_mindistTC_sites","select_bind_cluster",
                                         "select_compared_label",
                                         "norm_dist")
    return(peaksites_compared_infor)
  }
  
  if(parclip==FALSE){
    compared_num <- overlapbind_map[unlist(compare_label)]$xHits
    names(norm_pos) <- unique(as.character(mindist_TC_data$peak_num))
    peaksites_compared_infor <- list(mindist_TC_data,compared_num,norm_pos)
    names(peaksites_compared_infor) <- c("select_mindistTC_sites","select_compared_label",
                                         "norm_dist")
    return(peaksites_compared_infor)
  }
  
}

