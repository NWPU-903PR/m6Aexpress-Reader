###binding sites overlap filtered peak sites
mapped_filterpeak_bindsites <- function(overlap_bind_sites,bindsites_peak,parclip=TRUE){
  if(parclip==TRUE){
    TC_sites <- overlap_bind_sites$highconfi_overlap_TC
    bindcluster <- overlap_bind_sites$overlap_cluster
    bindcluster_GR <- GRanges(seqnames = as.character(bindcluster$seqnames),
                              IRanges(start = as.numeric(as.character(bindcluster$start)),
                                      end = as.numeric(as.character(bindcluster$end))),
                              strand=as.character(bindcluster$strand))
    
    bindsites_peak_GR <- GRanges(seqnames = as.character(bindsites_peak$seqnames),
                                 IRanges(start = as.numeric(as.character(bindsites_peak$start)),
                                         end = as.numeric(as.character(bindsites_peak$end))),
                                 strand=as.character(bindsites_peak$strand))
    
    ind_tr <- findOverlaps(TC_sites,bindsites_peak_GR,type = "within")
    TC_overlap_filter <- TC_sites[unique(ind_tr@from)]
    TC_cluster <- findOverlaps(TC_overlap_filter,bindcluster_GR,type="within")
    bindcluster_filter <- bindcluster[unique(TC_cluster@to),]
    ###
    bindsite_overlap_filterpeak_infor <- list(TC_overlap_filter,bindcluster_filter)
    names(bindsite_overlap_filterpeak_infor) <- c("bindsites_overlap_filterpeak","bindcluster_overlap_filterpeak")
    return(bindsite_overlap_filterpeak_infor)
  }
  if(parclip==FALSE){
    bindsites_peak_GR <- GRanges(seqnames = as.character(bindsites_peak$seqnames),
                                 IRanges(start = as.numeric(as.character(bindsites_peak$start)),
                                         end = as.numeric(as.character(bindsites_peak$end))),
                                 strand=as.character(bindsites_peak$strand))
    
    ind_tr <- findOverlaps(overlap_bind_sites,bindsites_peak_GR,type = "within")
    bindsite_overlap_filterpeak_infor <- overlap_bind_sites[unique(ind_tr@from)]
    return(bindsite_overlap_filterpeak_infor)
  }
}
