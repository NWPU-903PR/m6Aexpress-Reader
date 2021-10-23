.peak_site_filter <- function(m6Asite_readsinfor,filter_reads_num){
  methy_read <- m6Asite_readsinfor[,grep("IP",colnames(m6Asite_readsinfor))]
  CTL_methy_read <- methy_read[,1:(ncol(methy_read)/2)]
  Treated_methy_read <- methy_read[,-match(colnames(CTL_methy_read),colnames(methy_read))]
  colnames(CTL_methy_read) <- paste0("untreated_IP",1:ncol(CTL_methy_read))
  colnames(Treated_methy_read) <- paste0("treated_IP",1:ncol(Treated_methy_read))
  ####
  CTL_methy_selectlabel <- which(rowMeans(CTL_methy_read)>filter_reads_num)
  Treated_methy_selectlabel <- which(rowMeans(Treated_methy_read)>filter_reads_num)
  methy_selectlable <- intersect(CTL_methy_selectlabel,Treated_methy_selectlabel)
  ####
  unmethy_read <- m6Asite_readsinfor[,grep("Input", colnames(m6Asite_readsinfor))]
  CTL_unmethy_read <- as.data.frame(unmethy_read[,1:(ncol(unmethy_read)/2)])
  Treated_unmethy_read <- unmethy_read[,-(match(colnames(CTL_unmethy_read),colnames(unmethy_read)))]
  colnames(CTL_unmethy_read) <- paste0("untreated_Input",1:ncol(CTL_unmethy_read))
  colnames(Treated_unmethy_read) <- paste0("treated_Input", 1:ncol(Treated_unmethy_read))
  ####
  CTL_unmethy_selectlabel <- which(rowMeans(CTL_unmethy_read)>filter_reads_num)
  Treated_unmethy_selectlabel <- which(rowMeans(Treated_unmethy_read)>filter_reads_num)
  unmethy_selectlable <- intersect(CTL_unmethy_selectlabel,Treated_unmethy_selectlabel)
  ####
  select_label <- intersect(methy_selectlable,unmethy_selectlable)
  
  peak_reads_filter <- m6Asite_readsinfor[select_label,]
  return(peak_reads_filter)
}


filter_peak_readscount <- function(bind_nobindgene_peak,filter_reads_num=5){
  bindgene_bindpeak <- bind_nobindgene_peak$bindgene_overlap_peak
  bindgene_nonbindpeak <- bind_nobindgene_peak$bindgene_nonoverlap_peak
  nonbindgene_peak <- bind_nobindgene_peak$nonbindgene_peak
  bindgene_bindpeak_filter <- .peak_site_filter(m6Asite_readsinfor=bindgene_bindpeak,filter_reads_num=filter_reads_num)
  bindgene_nonbindpeak_filter <- .peak_site_filter(m6Asite_readsinfor=bindgene_nonbindpeak,filter_reads_num=filter_reads_num)
  nonbindgene_peak_filter <- .peak_site_filter(m6Asite_readsinfor=nonbindgene_peak,filter_reads_num=filter_reads_num)
  filtered_peaksites_infor <- list(bindgene_bindpeak_filter,bindgene_bindpeak_filter,bindgene_nonbindpeak_filter)
  names(filtered_peaksites_infor) <- c("bindgene_bindpeak_filter","bindgene_nonbindpeak_filter","nonbindgene_filter")
  return(filtered_peaksites_infor)
}