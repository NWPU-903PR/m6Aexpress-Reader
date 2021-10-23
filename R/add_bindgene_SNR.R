###bind gene SNR
bindgene_SNR <- function(bindgene_peak_SNR_infor){
  bindgene_bindallsites_addSNRinfor <- bindgene_peak_SNR_infor$onlybindgene_bindsites_addSNRinfor
  bindgene_bindpartsites_addSNRinfor <- bindgene_peak_SNR_infor$bindgene_distsites_addSNRinfor
  bindgene_SNR <- function(bindgene_infor, gene_name){
    bindgenes_SNR <-data.frame()
    for (i in 1:length(gene_name)) {
      onebindgene_SNR <- log(sum(exp(bindgene_infor[bindgene_infor$gene_name==gene_name[i],]$SNR)))
      onebindgenes_SNR <- data.frame(gene_name=gene_name[i],SNR=onebindgene_SNR)
      bindgenes_SNR <- rbind(bindgenes_SNR,onebindgenes_SNR)
    }
    return(bindgenes_SNR)
  }
  bindgene_bindallsites_SNR <- bindgene_SNR(bindgene_infor=bindgene_bindallsites_addSNRinfor,
                                            gene_name = unique(as.character(bindgene_bindallsites_addSNRinfor$gene_name)))
  
  bindgene_partbind_SNR <- bindgene_SNR(bindgene_infor=bindgene_bindpartsites_addSNRinfor,
                                        gene_name = unique(as.character(bindgene_bindpartsites_addSNRinfor$gene_name)))
  bindgenes_SNR <- data.frame(rbind(bindgene_bindallsites_SNR,bindgene_partbind_SNR),
                              bind_case=c(rep("bind_allsites",nrow(bindgene_bindallsites_SNR)),
                                          rep("bind_partsites",nrow(bindgene_partbind_SNR))))
  bindgenes_SNR$SNR <- round(exp(bindgenes_SNR$SNR),2)
  
  return(bindgenes_SNR)
}






