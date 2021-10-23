add_peak_SNR <- function(bindgene_bindpeak,bindgene_nobindpeak,bindgene_nobind_peakdist_infor,overlapped_bindsites,parclip=TRUE){
  make_GR <- function(target_sites){
    target_GR <- GRanges(seqnames = as.character(target_sites$seqnames),
                         IRanges(start = as.numeric(as.character(target_sites$start)),
                                 end = as.numeric(as.character(target_sites$end))),
                         strand = as.character(target_sites$strand))
    return(target_GR)
  }
  bindgene_peaksites <- bindgene_bindpeak
  bindpeak_GR <- make_GR(target_sites=bindgene_peaksites)
  if(parclip=TRUE){
    overlap_TC <- overlapped_bindsites[[1]]
    bindclusters <- overlapped_bindsites[[2]]
    TC_label <- bindgene_nobind_peakdist_infor$select_compared_label
    peak_dist <- bindgene_nobind_peakdist_infor$norm_dist
    dist_TC <- overlap_TC[TC_label]
    ###
    bindclusters_GR <- make_GR(target_sites=bindclusters)
    distTC_bindcluster <- bindclusters[unique(findOverlaps(dist_TC,bindclusters_GR,type = "within")@to),]
    distTC_bindcluster_GR <-  make_GR(target_sites=distTC_bindcluster)

    indtr <- findOverlaps(dist_TC,bindpeak_GR,type = "within")
    select_bindpeak <- bindgene_peaksites[unique(indtr@to),]
    select_bindpeak_GR <- bindpeak_GR[unique(indtr@to)]
    ###
    peaksites_clusterbind <- data.frame()
    add_SNRinfor_peak <- data.frame()
    for (i in 1:nrow(select_bindpeak)){
      ind_tr <- findOverlaps(dist_TC,select_bindpeak_GR[i],type = "within")
      select_TC <- dist_TC[unique(ind_tr@from)]
      TC_cluster <- findOverlaps(select_TC,distTC_bindcluster_GR,type="within")
      slect_cluster <- distTC_bindcluster[unique(TC_cluster@to),]
      onepeak_SumLogOdds <- max(as.numeric(as.character(slect_cluster$SumLogOdds)))
      onepeak_RelLogOdds <- max(as.numeric(as.character(slect_cluster$RelLogOdds)))
      oneadd_SNRinforpeak <- data.frame(select_bindpeak[i,],SumLogOdds=onepeak_SumLogOdds,RelLogOdds=onepeak_RelLogOdds)
      add_SNRinfor_peak <- rbind(add_SNRinfor_peak,oneadd_SNRinforpeak)
    }
    ##bind gene all bind peak sites
    onlybindgene_sites <- bindgene_peaksites[which(is.na(match(bindgene_peaksites$gene_name,add_SNRinfor_peak$gene_name))),]
    onlybindgene_sites_GR <- make_GR(target_sites=onlybindgene_sites)
    onlybindgene_sites_addSNR <- data.frame()
    for (i in 1:nrow(onlybindgene_sites)) {
      onepeak_ind <- findOverlaps(overlap_TC,onlybindgene_sites_GR[i],type = "within")
      onepeak_TC <- overlap_TC[unique(onepeak_ind@from)]
      onepeak_cluster <- bindclusters[unique(findOverlaps(onepeak_TC,bindclusters_GR,type = "within")@to),]
      onepeak_SumLogOdds <- max(onepeak_cluster$SumLogOdds)
      onepeak_RelLogOdds <- max(onepeak_cluster$RelLogOdds)
      onepeak_addSNR <- data.frame(onlybindgene_sites[i,],SumLogOdds=onepeak_SumLogOdds,RelLogOdds=onepeak_RelLogOdds)
      onlybindgene_sites_addSNR <- rbind(onlybindgene_sites_addSNR, onepeak_addSNR)
    }
    
     bindgene_nobind_peak <- bindgene_nobindpeak
     bindgene_nobind_adddist <- data.frame()
     for (i in 1:length(peak_dist)) {
       onebindgene_nobind_adddist <- data.frame(bindgene_nobind_peak[as.numeric(as.character(str_remove(names(peak_dist),"peak")))[i],],dist=peak_dist[i])
       bindgene_nobind_adddist <- rbind(bindgene_nobind_adddist,onebindgene_nobind_adddist)
     }

    bindgene_SNRinfor_nobind_adddist <- list(onlybindgene_sites_addSNR,add_SNRinfor_peak,bindgene_nobind_adddist)
    names(bindgene_SNRinfor_nobind_adddist) <- c("onlybindgene_bindsites_addSNRinfor","bindgene_distsites_addSNRinfor","bindgene_nobind_adddist")
    return(bindgene_SNRinfor_nobind_adddist)
  }
  if(parclip=FALSE){
    overlap_bindsites <-overlapped_bindsites
    mindist_sites <- bindgene_nobind_peakdist_infor$select_mindistTC_sites
    peak_dist <- bindgene_nobind_peakdist_infor$norm_dist
    bindsites_dist_label <- bindgene_nobind_peakdist_infor$select_compared_label
    select_bindsites <- overlap_bindsites[(bindsites_dist_label)]
 
    indtr <- findOverlaps(alldist_bindsites_GR,bindpeak_GR,type = "within")
    select_bindpeak <- bindgene_peaksites[unique(indtr@to),]
    select_bindpeak_GR <- bindpeak_GR[unique(indtr@to)]
    add_SNR <- function(select_reader_bindpeak, select_bindpeak_GR,selectbind_site){
      add_SNRinfor_peak <- data.frame()
      for (i in 1:nrow(select_reader_bindpeak)) {
        onepeak_ind <- findOverlaps(selectbind_site,select_bindpeak_GR[i],type = "within")
        onepeak_bindsites <- selectbind_site[unique(onepeak_ind@from)]
        onepeak_bindSNR <- mean(onepeak_bindsites$score)
        onepeak_addSNR <- data.frame(select_reader_bindpeak[i,],SNR=onepeak_bindSNR)
        add_SNRinfor_peak <- rbind(add_SNRinfor_peak,onepeak_addSNR)
      }
    }

    add_SNRinfor_partpeak <- add_SNR(select_reader_bindpeak=select_bindpeak, select_bindpeak_GR=select_bindpeak_GR,selectbind_site=select_bindsites)
    ##only or all bind DM sites gene and add the SNR 
    onlybindgene_peak <- bindgene_peaksites[which(is.na(match(bindgene_peaksites$gene_name,add_SNRinfor_partpeak$gene_name))),]
    onlybindgene_sites_GR <- make_GR(target_sites = onlybindgene_peak)
    add_SNRinfor_allpeak <- add_SNR(select_reader_bindpeak=onlybindgene_peak, select_bindpeak_GR=onlybindgene_sites_GR,selectbind_site=overlap_bindsites)
    
    bindgene_nobind_peak <- bindgene_nobindpeak
    bindgene_nobind_adddist <- data.frame()
    for (i in 1:length(peak_dist)) {
      onebindgene_nobind_adddist <- data.frame(bindgene_nobind_peak[as.numeric(as.character(str_remove(names(peak_dist),"peak")))[i],],dist=peak_dist[i])
      bindgene_nobind_adddist <- rbind(bindgene_nobind_adddist,onebindgene_nobind_adddist)
    }
    bindgene_SNRinfor_nobind_adddist <- list(add_SNRinfor_allpeak,add_SNRinfor_partpeak,bindgene_nobind_adddist)
    names(bindgene_SNRinfor_nobind_adddist) <- c("onlybindgene_bindsites_addSNRinfor","bindgene_distsites_addSNRinfor","bindgene_nobind_adddist")
    return(bindgene_SNRinfor_nobind_adddist)
    
  }
}