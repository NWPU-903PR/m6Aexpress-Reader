##match DE-DM gene
obtain_candidate_gene_infor <- function(sig_DE_gene,sig_DM_gene,output_path){
  sigDE_genecount <- sig_DE_gene$sigDE_genecount
  sigDM_genemethy <- sig_DM_gene$DM_gene_methylevel
  sigDE_geneinfor <- sig_DE_gene$sigDE_geneinfor
  sigDM_geneinfor <- sig_DM_gene$DM_gene_methyinfor
  match_gene_name <- intersect(sigDE_genecount$gene_name, sigDM_genemethy$gene_name)
  candidate_gene_countmethy <- data.frame()
  candidate_gene_DEDMinfor <- data.frame()
  for (i in 1:length(match_gene_name)) {
    one_countmethy <- cbind(sigDE_genecount[which(!is.na(match(sigDE_genecount$gene_name,match_gene_name[i]))),],
                            sigDM_genemethy[which(!is.na(match(sigDM_genemethy$gene_name,match_gene_name[i]))),-c(1,ncol(sigDM_genemethy))])
    
    one_gene_DEDMinfor <- cbind(sigDE_geneinfor[which(!is.na(match(sigDE_geneinfor$gene_name,match_gene_name[i]))),c(1,3,(ncol(sigDE_geneinfor)-1):(ncol(sigDE_geneinfor)))],
                                sigDM_geneinfor[which(!is.na(match(sigDM_geneinfor$gene_name,match_gene_name[i]))),c(2,5:6)],sigDM_genemethy[which(!is.na(match(sigDM_genemethy$gene_name,match_gene_name[i]))),ncol(sigDM_genemethy)])
    
    colnames(one_gene_DEDMinfor)[2:ncol(one_gene_DEDMinfor)] <- c("DE_LFC","DE_pvalue","DE_padj","DM_LFC","DM_pvalue","DM_padj","case")
    candidate_gene_countmethy <- rbind(candidate_gene_countmethy,one_countmethy)
    candidate_gene_DEDMinfor <- rbind(candidate_gene_DEDMinfor,one_gene_DEDMinfor)
  }
  write.table(candidate_gene_countmethy, file = paste0(output_path,"/candidate_expr_methy.tab"), quote = FALSE, row.names = FALSE)
  writexl::write_xlsx(candidate_gene_DEDMinfor,path = paste0(output_path,"candidate_DEDM_infor.xlsx"))
  cand_gene_countmethy_path <- paste0(output_path,"/candidate_expr_methy.tab")
  cand_gene_DEDM_infor_path <- paste0(output_path,"candidate_DEDM_infor.xlsx")
  size_factor <- sig_DE_gene$size_factor
  output_file_infor <- list(cand_gene_countmethy_path,cand_gene_DEDM_infor_path,size_factor)
  names(output_file_infor) <- c("gene_countmety_path","DEDM_infor_path","size_factor")
  return(output_file_infor)
} 