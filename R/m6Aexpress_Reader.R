m6Aexpress_Reader_model <- function(candidate_gene_infor,bindgene_strength_infor,CUTOFF_TYPE,pvalue, FDR,out_dir=NA){
  py_code <- system.file("extdata", "R_runpython.py", package = "m6AexpressReader")
  source_python(py_code)
  fileNameCount=candidate_gene_infor$gene_countmety_path
  librarySizes=as.numeric(candidate_gene_infor$size_factor)
  out_result <- suppressMessages(try(fun_R_call(fileNameCount, librarySizes),silent=TRUE))
  j=0
  while((is.null(nrow(out_result)))&j<10){
    out_result <- suppressMessages(try(fun_R_call(fileNameCount, librarySizes),silent=TRUE))
    j <- j+1
  }
  if(j==10){
    print("There is no significant m6A regulated expression gene")
  }
  if(j<10){
    genecoutmethy <-read.table(fileNameCount,header = T)
    size_factor<-librarySizes
    exprmethyre <- out_result
    match_count_methy <- data.frame()
    for (i in 1:length(exprmethyre$Gene_ID)) {
      one_gene <- genecoutmethy[genecoutmethy$gene_name==as.character(exprmethyre$Gene_ID)[i],]
      match_count_methy <- rbind(match_count_methy,one_gene)
    }
    
    match_methy <- data.frame()
    for (i in 1:nrow(match_count_methy)) {
      one_methy <- exprmethyre[exprmethyre$Gene_ID==as.character(match_count_methy$gene_name)[i],]
      match_methy <- rbind(match_methy, one_methy)
    }
    
    ##match gene binding strength
    bindSNR <- vector()
    for (i in 1:nrow(match_count_methy)) {
      if(length(intersect(bindgene_strength_infor$gene_name,as.character(match_count_methy$gene_name)[i]))=0){
        bindSNR[i] <- 1
      }
      if(length(intersect(bindgene_strength_infor$gene_name,as.character(match_count_methy$gene_name)[i]))>0){
        onematch_SNR <- bindgene_strength_infor[bindgene_strength_infor$gene_name==as.character(match_count_methy$gene_name)[i],]
        bindSNR[i] <- (as.numeric(onematch_SNR$SNR))
      }
    }
    
    genecount <- match_count_methy[,2:((ncol(match_count_methy)+1)/2)]
    methyintensity <- match_count_methy[,(((ncol(match_count_methy)+1)/2)+1):ncol(match_count_methy)]
    
    alpha <- vector()
    beta <- vector()
    intersecpt <- vector()
    for (i in 1:nrow(match_count_methy)) {
      alpha[i] <- as.numeric(exprmethyre[exprmethyre$Gene_ID==as.character(match_count_methy$gene_name)[i],ncol(exprmethyre)])
      beta[i] <-  as.numeric(exprmethyre[exprmethyre$Gene_ID==as.character(match_count_methy$gene_name)[i],3])
      intersecpt[i] <-as.numeric(exprmethyre[exprmethyre$Gene_ID==as.character(match_count_methy$gene_name)[i],2])
    }
    
    beta_value <- cbind(intersecpt,beta)
    base_level <- rep(1, ncol(methyintensity))
    gene_name <- as.character(match_count_methy$gene_name)
    
    ##empirical prior estimate for beta
    ##use bet0 and bet1 prior with mu
    beta0_lower <- quantile(intersecpt,0.15)
    beta0_up <- quantile(intersecpt,0.85)
    beta1_lower <- quantile(beta,0.15)
    beta1_up <- quantile(beta,0.85)
    select_beta0 <- intersecpt[which((intersecpt>beta0_lower)&(intersecpt<beta0_up))]
    select_beta1 <- beta[which((beta>beta1_lower)&(beta<beta1_up))]
    nu=mean(select_beta1)
    iu=mean(select_beta0)
    var_beta<-var(select_beta1)
    var_IN <- var(select_beta0)
    lambda <- (1/var_beta)
    lambda_IN <- (1/var_IN)
    ##final estimate the beta value
    new_beta <- data.frame()
    for (i in 1:nrow(methyintensity)) {
      methy_gene <- t(methyintensity[i,])
      gene_matrix <- cbind( base_level,methy_gene)
      model_matrix <- matrix(cbind(gene_matrix[,1],gene_matrix[,2]),ncol = ncol(gene_matrix))
      mu <- as.numeric(size_factor*t(exp(gene_matrix%*% (beta_value[i,])))) 
      gene_count <- as.numeric(genecount[i,])
      disp <- alpha[i]
      
      beta_mat <- matrix(c(beta_value[i,1],beta_value[i,2]),ncol = ncol(beta_value))
      counts <- matrix(gene_count,nrow = length(gene_count))
      matrix_size_factor <-matrix(size_factor, ncol = length(size_factor))
      
      
      m = nrow(model_matrix)
      # ridge <- diag(lambda, nrow = ncol(model_matrix), ncol = ncol(model_matrix))
      beta_hat = t(beta_mat)
      ##get loglike hood function of beta
      beta_loglike <- function(beta) {
        mu_row <- as.numeric(t(matrix_size_factor)*exp(model_matrix%*%beta))
        count_row <- as.numeric(counts)
        nbinomLogLike <- (bindSNR[i])*sum(dnbinom(count_row,mu=mu_row,size=1/disp,log=TRUE))
        # nbinomLogLike <- (1)*sum(dnbinom(count_row,mu=mu_row,size=1/disp,log=TRUE))
        logPrior <- (1)*sum(dnorm(beta,mean=c(iu,nu),sd=c(sqrt(1/lambda_IN), sqrt(1/lambda)),log=TRUE))
        # logPrior <- (1)*sum(dnorm(beta[2],mean=c(0,0),sd=c(c1[[2]], c2[[2]]),log=TRUE))
        beta_LogLike <- -1*(nbinomLogLike+logPrior)
        
      }
      o <- optim(beta_hat, beta_loglike, method="L-BFGS-B",hessian = TRUE)
      one_new <- t(o$par)
      se.beta.hat <- sqrt(diag(solve(o$hessian)))
      ##wald test
      
      sigma <- solve(o$hessian)
      
      wdtest <- wald.test(b =as.numeric(one_new) , Sigma = sigma, Terms = 2)
      pvalue <- as.numeric(wdtest$result$chi2[3])
      one_data <- t(c(as.numeric(one_new),se.beta.hat, pvalue))
      new_beta <- rbind(new_beta, one_data)
    }
    ##get loglike hood function of beta
    gene_name <-as.character(match_methy$Gene_ID)
    adj_beta <- cbind(gene_name, new_beta)
    colnames(adj_beta) <- c("gene_name", "Beta0","Beta1","SE_Beta0","SE_Beat1", "pvalue")
    pvalues <- as.numeric(adj_beta$pvalue)
    padj <- p.adjust(pvalues, method = "BH")
    padj_beta <- cbind(adj_beta, padj)
    m6Aregexp_addDEDM_infor <- data.frame()
    DEDM_gene_infor <- readxl::read_xlsx(candidate_gene_infor$DEDM_infor_path)
    overlap_gene <- intersect(padj_beta$gene_name,DEDM_gene_infor$gene_name)
    for (i in 1:nrow(padj_beta)) {
      one_addDEDM <- cbind(padj_beta[which(!is.na(match(padj_beta$gene_name,overlap_gene[i]))),],
                           DEDM_gene_infor[which(!is.na(match(DEDM_gene_infor$gene_name,overlap_gene[i]))),c(grep("LFC",colnames(DEDM_gene_infor)),
                                                                                                             ncol(DEDM_gene_infor))])
      m6Aregexp_addDEDM_infor <- rbind(m6Aregexp_addDEDM_infor,one_addDEDM)
    }
    
    if(is.na(out_dir)){
      out_dir = dirname(Input_file[[1]])
    }
    if (CUTOFF_TYPE =="FDR") {
      select_m6Aregexp_gene <- m6Aregexp_addDEDM_infor[m6Aregexp_addDEDM_infor$padj<FDR,]
      write.table(select_m6Aregexp_gene,file=paste0(out_dir,sep="/","m6AexpressReader_result.xls"), sep="\t",row.names =FALSE,quote = FALSE)
      
    }
    if (CUTOFF_TYPE =="pvalue") {
      select_m6Aregexp_gene <- m6Aregexp_addDEDM_infor[m6Aregexp_addDEDM_infor$pvalue<pvalue,]
      write.table(select_m6Aregexp_gene,file=paste0(out_dir,sep="/", "m6AexpressReader_result.xls"), sep="\t",row.names =FALSE,quote = FALSE)
      
    }
    return(select_m6Aregexp_gene)
  }
  
}

