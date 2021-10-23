library(stringr)
library(wavClusteR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicFeatures)
parclip_reader_bindingsites <- function(par_bam,annotation_file){
  Bam <- readSortedBam(filename = par_bam)
  Bam <- readSortedBam(filename = DF2file )
  tagmd <- Bam$tag.MD
  subtagmd <- sub("\\^","_",tagmd)
  str_replacestrs <- str_replace_all(subtagmd,"_",NA_character_)
  select_label <- which(!is.na(str_replacestrs))
  countTable <- getAllSub( Bam[select_label], minCov = 10 )
  model <- fitMixtureModel(countTable, substitution = "TC")
  support <- getExpInterval(model, bayes = TRUE ) 
  highConfSub <- getHighConfSub( countTable, 
                                 support = support, 
                                 substitution = "TC" )
  coverage <- coverage(Bam[select_label])
  clusters <- getClusters( highConfSub = highConfSub,
                           coverage = coverage,
                           sortedBam = Bam[select_label],
                           cores = 10 )
  
  wavclusters <- filterClusters( clusters = clusters, 
                                 highConfSub = highConfSub,
                                 coverage = coverage, 
                                 model = model, 
                                 genome = Hsapiens, 
                                 refBase = "T", 
                                 minWidth = 12)
  
  txDB <- GenomicFeatures::makeTxDbFromGFF(annotation_file)
  genes_infor <- genes(txDB)
  wavcluster_dataframe <- data.frame(wavclusters)
  ind_tr <- findOverlaps(wavclusters, genes_infor, type = "within")
  select_genes <- genes_infor[unique(ind_tr@to)]
  select_geneID <- as.character(names(select_genes))
  genecluster <- data.frame()
  for (i in 1:length(select_geneID)) {
    one_genes <- select_genes[which(!is.na(match(names(select_genes),select_geneID[i])))]
    one_overlap_tr <- findOverlaps(wavclusters, one_genes, type = "within")
    one_genes_cluster <- wavclusters[unique(one_overlap_tr@from)]
    one_genescluster <- as.data.frame(one_genes_cluster)
    one_gene_ID <- rep(select_geneID[i],nrow(one_genescluster))
    one_genecluster <- data.frame(gene_ID=one_gene_ID,one_genescluster)
    genecluster <- rbind(genecluster,one_genecluster)
  }
  clip_result <- list(highConfSub,genecluster)
  names(clip_result) <- c("high_confi_TC","bind_cluster_infor")
  return(clip_result)
}


