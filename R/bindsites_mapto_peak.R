bindsites_mapto_peak <- function(peak_sites_infor,mapped_peak_GR,bind_sites,parclip=TRUE){
  consis_peak_infor <- peak_sites_infor
  map_peak_GR <- mapped_peak_GR[countOverlaps(mapped_peak_GR,mapped_peak_GR,type = "equal")==1]
  consis_GR <- GRanges(seqnames = as.character(consis_peak_infor$seqnames),
                             IRanges(start = as.numeric(as.character(consis_peak_infor$start)),
                                     end = as.numeric(as.character(consis_peak_infor$end))),
                                     strand = as.character(consis_peak_infor$strand),
                                     gene_name=consis_peak_infor$gene_name)
  ###
  new_consis_peak <- data.frame()
  rm_label <- vector()
  for (i in 1:length(map_peak_GR)) {
    start_tr <- findOverlaps(consis_GR,map_peak_GR[i],type = "start")
    end_tr <- findOverlaps(consis_GR,map_peak_GR[i],type = "end")
    consis_tr <- intersect(unique(start_tr@from),unique(end_tr@from))
    if(length(consis_tr)==0){
      rm_label[i] <- i
    }
    if(length(consis_tr)>0){
      rm_label[i] <- 0
    }
    one_peak <- consis_peak_infor[consis_tr,]
    # mcols(consis_peak_GR[i]) <- data.frame(mcols(consis_peak_GR[i]),gene_name=as.character(one_peak$gene_name))
    new_consis_peak <- rbind(new_consis_peak, one_peak)
  }
  consis_peak_GR <- consis_GR[which(rm_label==0)]
  names(consis_peak_GR) <- paste0("peak",1:length(consis_peak_GR))
  if(parclip==TRUE){
    highconfi_sub <- unlist(bind_sites$select_TC_sites)
    bind_cluster <- bind_sites$select_bind_cluster
    clip_genename <- as.character(bind_cluster$gene_name)
    chrs <- as.character(bind_cluster$seqnames)
    starts <- as.numeric(as.character(bind_cluster$start))
    ends <- as.numeric(as.character(bind_cluster$end))
    Strands <- as.character(bind_cluster$strand)
    cluster_GR <- GRanges(seqnames = chrs, ranges = IRanges(starts, ends),strand =Strands,names=clip_genename)
    ##overlap binding sites to peak
    ind_tr <- findOverlaps(highconfi_sub, consis_peak_GR, type = "within")
    highconfi_sub_overlap <- highconfi_sub[unique(ind_tr@from)]
    TC_cluster_tr <- findOverlaps(highconfi_sub_overlap,cluster_GR,type = "within")
    select_cluster <-bind_cluster[unique(TC_cluster_tr@to),] 
    peak_overlap_TC <- consis_peak_GR[unique(ind_tr@to)]
    ##obtain gene name given peak number order
    peak_label <- as.character(names(peak_overlap_TC))
    peak_number <- as.numeric(as.character(str_remove(peak_label,"peak")))
    overlap_TC_peaksites <- new_consis_peak[peak_number,]
    ##DF2bind gene
    peak_geneID <- unique(as.character(new_consis_peak$gene_name))
    # DFbind_geneID <- unique(as.character(peak_overlap_TC$gene_name))
    bind_geneID <- unique(as.character(overlap_TC_peaksites$gene_name))
    nonbind_geneID <- peak_geneID[which(is.na(match(peak_geneID,bind_geneID)))]
    nonbindgene_peak <- new_consis_peak[which(!is.na(match(new_consis_peak$gene_name,nonbind_geneID))),]
    bindgene_allpeak <- new_consis_peak[which(!is.na(match(new_consis_peak$gene_name,bind_geneID))),]
    # DF2bindgene_overlap_peak <- peak_overlap_TC
    bindgene_overlap_peak <- overlap_TC_peaksites
    bindgene_nonoverlap_peak <- bindgene_allpeak[which(is.na(match(bindgene_allpeak$start,bindgene_overlap_peak$start))),]
    
    consis_peaksites_infor <- list(bindgene_overlap_peak,bindgene_nonoverlap_peak,nonbindgene_peak)
    names(consis_peaksites_infor) <- c("bindgene_bind_peak","bindgene_nonbind_peak","nonbindgene_peak")
    ##binding sites overlap peak infor
    bindcluster_overlap_infor <- list(highconfi_sub_overlap,select_cluster)
    names(bindcluster_overlap_infor) <- c("highconfi_overlap_TC","overlap_cluster")
    ##combine peak and reader binding sites
    bindingsites_peak_overlap <- list(bindcluster_overlap_infor,consis_peaksites_infor)
    names(bindingsites_peak_overlap) <- c("binding_sites_overlap","peak_infor")
    return(bindingsites_peak_overlap)
  }
  if(parclip==FALSE){
    bind_sites_GRlist <- bind_sites
    bind_sites_GR <- unlist(bind_sites_GRlist)
    ##TC overlap in peak sites
    ind_tr <- findOverlaps(bind_sites_GR, consis_peak_GR, type = "within")
    bindsites_overlap <- bind_sites_GR[unique(ind_tr@from)]
    
    peak_overlap_site <- consis_peak_GR[unique(ind_tr@to)]
    ##obtain gene name given peak number order
    peak_label <- as.character(names(peak_overlap_site))
    peak_number <- as.numeric(as.character(str_remove(peak_label,"peak")))
    overlap_peaksites <- new_consis_peak[peak_number,]
    ##DF2bind gene
    peak_geneID <- unique(as.character(new_consis_peak$gene_name))
    # DFbind_geneID <- unique(as.character(peak_overlap_TC$gene_name))
    bind_geneID <- unique(as.character(overlap_peaksites$gene_name))
    nonbind_geneID <- peak_geneID[which(is.na(match(peak_geneID,bind_geneID)))]
    nonbindgene_peak <- new_consis_peak[which(!is.na(match(new_consis_peak$gene_name,nonbind_geneID))),]
    bindgene_allpeak <- new_consis_peak[which(!is.na(match(new_consis_peak$gene_name,bind_geneID))),]
    # DF2bindgene_overlap_peak <- peak_overlap_TC
    bindgene_overlap_peak <- overlap_peaksites
    bindgene_nonoverlap_peak <- bindgene_allpeak[which(is.na(match(bindgene_allpeak$start,
                                                                   bindgene_overlap_peak$start))),]
    
    consis_peaksites_infor <- list(bindgene_overlap_peak,bindgene_nonoverlap_peak,nonbindgene_peak)
    names(consis_peaksites_infor) <- c("bindgene_overlap_peak","bindgene_nonoverlap_peak","nonbindgene_peak")
    bindingsites_peak_overlap <- list(bindsites_overlap,consis_peaksites_infor)
    names(bindingsites_peak_overlap) <- c("binding_sites_overlap","peak_infor")
    return(bindingsites_peak_overlap)
  }
}
