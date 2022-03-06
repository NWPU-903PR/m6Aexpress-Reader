##obtain DM gene 
obtain_DM_geneinfor <- function(bindgene_methylevel,nobindgene_methylevel,p_value=0.05,num_cond1,num_cond2){
  all_bind_gene_methyintensity <- bindgene_methylevel$onlybindgene_sites_methyintensity_weight
  part_bindgene_methyintensity <- bindgene_methylevel$bindgenesites_methyintensity_weight
  bind_gene_methyintensity <- rbind(all_bind_gene_methyintensity,part_bindgene_methyintensity)
  nobindgene_methyintensity <- nobindgene_methylevel
  ###
  cases <- rep(c("bind_gene","no_bindgene"),c(nrow(bind_gene_methyintensity),nrow(nobindgene_methyintensity)))
  bind_nobind_gene_methyintensity <- rbind(bind_gene_methyintensity,nobindgene_methyintensity)
  bind_nobind_gene_methyintensity_infor <- data.frame(bind_nobind_gene_methyintensity,case=cases)
  methy_intensity <- bind_nobind_gene_methyintensity_infor[,-c(1,ncol(bind_nobind_gene_methyintensity_infor))]
  rownames(methy_intensity) <- as.character(bind_nobind_gene_methyintensity_infor$gene_name)
  design <- matrix(c(rep(1,(num_cond1+num_cond2)),rep(0,num_cond1),rep(1,num_cond2)),nrow = (num_cond1+num_cond2))
  fit <- lmFit(methy_intensity,design)
  fit <- eBayes(fit)
  diff_result <- topTable(fit,coef = 2,number = nrow(methy_intensity),genelist = as.character(rownames(methy_intensity)))
  selecth_gene <- as.character(diff_result[diff_result$P.Value<p_value,]$ID)
  # selecth_gene <- as.character(diff_result[diff_result$P.Value<0.05,]$ID)
  select_DM_gene_methy <- bind_nobind_gene_methyintensity_infor[which(!is.na(match(bind_nobind_gene_methyintensity_infor$gene_name,
                                                                             selecth_gene))),]
  select_DM_geneinfor <- diff_result[diff_result$P.Value<p_value,]
  rownames(select_DM_geneinfor) <- NULL
  colnames(select_DM_geneinfor)[grep("ID",colnames(select_DM_geneinfor))] <- "gene_name"
  sigDM_gene_methyinfor <- list(select_DM_gene_methy,select_DM_geneinfor)
  names(sigDM_gene_methyinfor) <- c("DM_gene_methylevel","DM_gene_methyinfor")
  return(sigDM_gene_methyinfor)
}
