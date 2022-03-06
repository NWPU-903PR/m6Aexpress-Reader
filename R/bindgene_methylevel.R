#########quantify methylation lelvel 
.peak_methy_level <- function(IP_Input_read,size_factor){
  IP_site_read <- IP_Input_read[,grep("IP",colnames(IP_Input_read))]
  Input_site_read <- IP_Input_read[,(grep("Input",colnames(IP_Input_read)))]
  IP_Input <- cbind(IP_site_read, Input_site_read)
  IP_Input_norm <- as.data.frame(t(t(IP_Input)/size_factor)) 
  IP_norm_site <- IP_Input_norm[,1:(ncol(IP_Input_norm)/2)]
  Input_norm_site <- IP_Input_norm[,((ncol(IP_Input_norm)/2)+1):(ncol(IP_Input_norm))]
  methy_level <- log((IP_norm_site+0.01)/(Input_norm_site+0.01)) 
  methy_level_infor <- data.frame(IP_Input_read[,-c(grep("IP",colnames(IP_Input_read)),
                                                    (grep("Input",colnames(IP_Input_read))))], methy_level)
  for (i in 1:nrow(methy_level_infor)) {
    for(j in grep("IP",colnames(methy_level_infor))[1]:ncol(methy_level_infor)){
      if(methy_level_infor[i,j]<=0){
        methy_level_infor[i,j] <- NA
      }
    }
  }
  return(methy_level_infor)
}

.weight_methy_level <- function(methy_site_infor,library_size,add_dist=T){
  size_factor <- as.numeric(library_size/exp(mean(log(library_size))))
  norm_methy_level <- .peak_methy_level(methy_site_infor,size_factor)
  norm_methy_level <- na.omit(norm_methy_level)
  # distdecay <- exp(-as.numeric(as.character(norm_methy_level$norm_pos)))
  table_genename <- as.data.frame(table(as.character(norm_methy_level$gene_name)))
  colnames(table_genename) <- c("gene_name","freq")
  select_gene_name <- unique(as.character(norm_methy_level$gene_name))
  if(add_dist==T){
    colnames(norm_methy_level)[grep("dist",colnames(norm_methy_level))] <- "dist_TC"
    dist <- as.numeric(as.character(norm_methy_level$dist_TC))
    distdecay <- exp(-dist/round(quantile(dist,0.75)))
    peak_SNR <- exp(as.numeric(as.character(norm_methy_level$RelLogOdds)))
    # peak_SNR <- (as.numeric(as.character(norm_methy_level$RelLogOdds)))
    # peak_weight <- peak_SNR*distdecay
    peak_weight <- distdecay
    weighted_peak <- data.frame(as.character(norm_methy_level$gene_name),peak_weight)
    colnames(weighted_peak)[1] <- "gene_name"
  }else{
    peak_SNR <-  exp(as.numeric(as.character(norm_methy_level$RelLogOdds)))
    # peak_SNR <-  (as.numeric(as.character(norm_methy_level$RelLogOdds)))
    # peak_weight <- peak_SNR
    peak_weight <- 1
    weighted_peak <- data.frame(as.character(norm_methy_level$gene_name),peak_weight)
    colnames(weighted_peak)[1] <- "gene_name"
  }
  
  genelevel_weighted <- data.frame()
  for (i in 1:length(select_gene_name)) {
    add_weight_methy <- norm_methy_level[which(!is.na(match(norm_methy_level$gene_name, select_gene_name[i]))),grep("IP",colnames(norm_methy_level))[1]:ncol(norm_methy_level)]*
      weighted_peak[weighted_peak$gene_name==select_gene_name[i],-1]
    mean_site_methyweighted <- round(colMeans(add_weight_methy),2)
    # mean_site_methyweighted <- round(add_weight_methy[which.max(abs(rowMeans(add_weight_methy[,((ncol(add_weight_methy)/2)+1):(ncol(add_weight_methy))])-
    #                                                                   rowMeans(add_weight_methy[,1:(ncol(add_weight_methy)/2)]))),],2)
    genelevel_weighted <- rbind(genelevel_weighted, mean_site_methyweighted)
  }
  new_weightedlevel <- cbind(select_gene_name, genelevel_weighted)
  colnames(new_weightedlevel) <- c("gene_name", colnames(norm_methy_level)[grep("IP",colnames(norm_methy_level))[1]:ncol(norm_methy_level)])
  ##select methy-level
  select_weightedlevel <- new_weightedlevel[rowSums(new_weightedlevel[,-1])>0,]
  last_select_weightedlevel <- select_weightedlevel[rowSums(round(select_weightedlevel[,-1],2))>0,]
  last_select_gene <- na.omit(last_select_weightedlevel)
  last_gene_name <- as.character(last_select_gene$gene_name) 
  last_selectweightedlevel <- cbind(last_gene_name,round(last_select_gene[,-1],2))
  colnames(last_selectweightedlevel)[1] <- "gene_name"
  rownames(last_selectweightedlevel) <- NULL
  geneweightedlevel <- last_selectweightedlevel[,-1]
  meanmethy <- rowMeans(geneweightedlevel)
  select_label <- which(meanmethy>0.1)
  last_geneweightedlevel <- last_selectweightedlevel[select_label, ]
  return(last_geneweightedlevel)
}
###
bindgene_methylevel <- function(bindgene_peakSNR,library_size){
  onlybindgene_bindsites_SNR <- bindgene_peakSNR$onlybindgene_bindsites_addSNRinfor
  bindgene_distsites_SNR <- bindgene_peakSNR$bindgene_distsites_addSNRinfor
  bindgene_nobind_adddist <- bindgene_peakSNR$bindgene_nobind_adddist
  bindgene_bindpeak_nobindpeak_genename <- intersect(bindgene_distsites_SNR$gene_name,bindgene_nobind_adddist$gene_name)
  select_bindgene_distsites_SNR <- bindgene_distsites_SNR[which(!is.na(match(bindgene_distsites_SNR$gene_name,bindgene_bindpeak_nobindpeak_genename))),]
  select_bindgene_nobind_adddist <- bindgene_nobind_adddist[which(!is.na(match(bindgene_nobind_adddist$gene_name,bindgene_bindpeak_nobindpeak_genename))),]
  another_allbindgene_SNR <- bindgene_distsites_SNR[which(is.na(match(bindgene_distsites_SNR$gene_name,
                                                                      bindgene_bindpeak_nobindpeak_genename))),]
  bindgene_bindallpeak_SNR <- rbind(onlybindgene_bindsites_SNR,another_allbindgene_SNR)
  
  
  
  weightedpeak_genemethy_intensity_onlybindgene_sites <- .weight_methy_level(methy_site_infor=bindgene_bindallpeak_SNR,library_size=library_size,add_dist=F)
  weightedpeak_genemethy_intensity_bindgene_nobind_sites_dist <- .weight_methy_level(methy_site_infor=select_bindgene_nobind_adddist,library_size=library_size,add_dist=T)
  weightedpeak_genemethy_intensity_bindgene_bindother_sites <- .weight_methy_level(methy_site_infor=select_bindgene_distsites_SNR,library_size=library_size,add_dist=F)
  ######DFbind DM sites gene methylation intensity
  combine_genename <- intersect(weightedpeak_genemethy_intensity_bindgene_nobind_sites_dist$gene_name,weightedpeak_genemethy_intensity_bindgene_bindother_sites$gene_name)
  bindgenesites_genemethyintensity <- data.frame()
  for (i in 1:length(combine_genename)) {
    one_combined_gene <- weightedpeak_genemethy_intensity_bindgene_nobind_sites_dist[weightedpeak_genemethy_intensity_bindgene_nobind_sites_dist$gene_name==combine_genename[i],-1]+
      weightedpeak_genemethy_intensity_bindgene_bindother_sites[weightedpeak_genemethy_intensity_bindgene_bindother_sites$gene_name==combine_genename[i],-1]
    one_combined_gene <- data.frame(gene_name=combine_genename[i],one_combined_gene)
    bindgenesites_genemethyintensity <- rbind(bindgenesites_genemethyintensity,one_combined_gene)
  }
  bindgene_weightedmethyintensity <- list(weightedpeak_genemethy_intensity_onlybindgene_sites,bindgenesites_genemethyintensity)
  names(bindgene_weightedmethyintensity) <- c("onlybindgene_sites_methyintensity_weight","bindgenesites_methyintensity_weight")
  return(bindgene_weightedmethyintensity)
}
