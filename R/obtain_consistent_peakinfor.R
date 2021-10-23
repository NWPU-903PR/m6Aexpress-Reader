library(gdata)
library(rtracklayer)
library(AnnotationDbi)
library(org.Hs.eg.db)

obtain_consistent_peakinfor <- function(peak_file, peak_bed){
  read_peak <- import(peak_bed)
  read_peak <- as.data.frame(read_peak)
  peak_name <- as.character(read_peak$name)
  con_peak <- TRUE
  select_peak <-cbind( as.character(read_peak$seqnames), as.numeric(as.character(read_peak$start)), 
                       as.numeric(as.character(read_peak$end)),as.numeric(as.character(read_peak$width)),
                       as.character(read_peak$strand),peak_name)
  colnames(select_peak) <- c("seqnames", "start", "end", "width", "strand", "gene_name")
  select_peak <- as.data.frame(select_peak)
  peak <- peak_file[["PEAK"]]
  READS_COUNT <- peak_file[["READS_COUNT"]]
  ## get size factor
  reads_count <- READS_COUNT[,-((ncol(READS_COUNT)-1):ncol(READS_COUNT))]
  totalreads <- colSums(reads_count)
  ##judue the peak type: consistent peak or all peak
  if(con_peak==TRUE){
    peak_loci <- peak$loci2peak_consistent
    
  }else{
    peak_loci <- peak$loci2peak_merged
  }
  # peak_reads_count
  no_peak=length(peak_loci[,1])
  peak_reads_count = READS_COUNT[1:no_peak,]
  no_sample=length(peak_reads_count[1,])-2
  # cut the unnecessary information
  peak_reads_count = READS_COUNT[1:no_peak,1:no_sample]
  # count
  for (ipeak in 1:no_peak) {
    temp=peak_loci[ipeak,]
    temp2=colSums(READS_COUNT[temp[1]:temp[2],1:no_sample])
    peak_reads_count[ipeak,1:no_sample]=temp2
  }
  # remove the overlapping window effects
  peak_reads_count = round (peak_reads_count * 30 / 200);
  colnames(peak_reads_count) <- c(paste0("IP", (1:(no_sample/2))), paste0("Input", (1:(no_sample/2))))
  # get peak
  peak_report <- data.frame()
  for (i in 1:no_peak) {
    peak_row_id=peak_loci[i,]
    
    # batch id
    batch_id=unique(READS_COUNT$batch_id[peak_row_id])
    lg.p=min(peak$PW$log.p[peak_row_id[1]:peak_row_id[2]])/log(10)
    lg.fdr=min(peak$PW$log.fdr[peak_row_id[1]:peak_row_id[2]])/log(10)
    fold_enrchment=exp(max(peak$PW$log.fc[peak_row_id[1]:peak_row_id[2]]))
    
    # get sig digits
    lg.p=signif(lg.p, digits = 3)
    lg.fdr=signif(lg.fdr, digits = 3)
    fold_enrchment=signif(fold_enrchment, digits = 3)
    
    test_result <- data.frame(lg.p=lg.p,lg.fdr=lg.fdr, fold_enrchment=fold_enrchment)               
    peak_report=rbind(peak_report,test_result)
  }
  ## peak site read and log(fdr) 
  peak_site_reads <- cbind(select_peak, peak_reads_count, peak_report[,2:3])
  peak_site_reads <- peak_site_reads[which(!is.na(peak_site_reads$gene_name)), ]
  peak_site_infor <- list(peak_site_reads, totalreads)
  names(peak_site_infor) <- c("peak_sites_infor","library_size")
  return(peak_site_infor)
  
}
