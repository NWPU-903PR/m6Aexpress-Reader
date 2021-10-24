# Usage Example
## In Treated VS Control context
### *Peak calling for methylation sites in case-control context*
```r
library(exomePeak)
library(GenomicFeatures)
library(m6ALogisticModel)
f1 <- "./CTL_IP1.bam"
f2 <- "./CTL_IP2.bam"
f3 <- "./CTL_Input1.bam"
f4 <- "./CTL_Input2.bam"
f5 <- "./M3KD_IP1.bam"
f6 <- "./M3KD_IP2.bam"
f7 <- "./M3KD_Input1.bam"
f8 <- "./M3KD_Input2.bam"
IP_BAM <- c(f1,f2,f5,f6)
INPUT_BAM <- c(f3,f4,f7,f8)

GENE_ANNO_GTF = "./hg19_GTF/genes.gtf"
txdbfile <- GenomicFeatures::makeTxDbFromGFF(GENE_ANNO_GTF)
txdb <- txdbfile
result = exomepeak(TXDB=txdb, IP_BAM=IP_BAM, INPUT_BAM=INPUT_BAM, OUTPUT_DIR= "./exomePeak_calling/")
##obtain consistent peak sites information
load("./exomePeak_output/exomePeak.Rdata")
peak_file <- tmp_rs
consisten_peak <- "./exomePeak_output/con_peak.bed"
peak_site_infor <- obtain_consistent_peakinfor(peak_file, peak_bed=consisten_peak)
peak_site_filter <- peak_sites_filter(m6Asite_readsinfor=peak_site_infor,num_cond1=2,filter_reads_num=5)
##Mapping peak sites to the longest transcript 
#map_consist_peak_longTX <- map_peak_longTX(filepath=consisten_peak,annotation_file=GENE_ANNO_GTF)
map_consist_peak_longTX <- map_peak_longTX(filepath=consisten_peak,annotation_file=GENE_ANNO_GTF,peak_sites_infor=peak_site_filter)
```
### *Identify reader binding sites from CLIP-seq data*
#### *For PAR-CLIP-seq data*
```r
library(stringr)
library(wavClusteR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicFeatures)
library(m6ALogisticModel)
YTHDF2_binding <- "./YTHDF2_binding.bam"
obtain_reader_bindingsites <- parclip_reader_bindingsites(par_bam=YTHDF2_binding,annotation_file=GENE_ANNO_GTF)
##map to the longest transcript
bindsites_map_longestTX <- bindsites_maplong_tr(binding_sites=obtain_reader_bindingsites,annotation_file=GENE_ANNO_GTF,parclip=TRUE)
```
#### *For eCLIP-seq or ICLIP-seq data*
```sh
nohup pureclip -i ./IGF2BP1_eCLIP/IGF2BP1_rep1.bam -i ./IGF2BP1_eCLIP/IGF2BP1_rep2.bam -bai ./IGF2BP1_eCLIP/IGF2BP1_rep1.bam.bai -bai ./IGF2BP1_eCLIP/IGF2BP1_rep2.bam.bai  -g ./hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -o ./IGF2BP1_eCLIP/IGF2BP1_eCLIP.bed -nt 10 &

nohup pureclip -i ./IGF2BP3_eCLIP/IGF2BP3_rep1.bam -i ./IGF2BP3_eCLIP/IGF2BP3_rep2.bam -bai ./IGF2BP3_eCLIP/IGF2BP3_rep1.bam.bai -bai ./IGF2BP3_eCLIP/IGF2BP3_rep2.bam.bai  -g ./hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -o ./IGF2BP3_eCLIP/IGF2BP3_eCLIP.bed -nt 10 &
```
```r
##Merge IGF2BP1 and IGF2BP3 binding sites together
IGF2BP1_bindingsites <- "./IGF2BP1_eCLIP/IGF2BP1_eCLIP.bed"
IGF2BP3_bindingsites <- "./IGF2BP3_eCLIP/IGF2BP3_eCLIP.bed"
IGF2BPsbindingsites <- mergbinding_sites(one_bindingsites=IGF2BP1_bindingsites, two_bindingsites=IGF2BP3_bindingsites)
##map the IGF2BPs binding sites to the longest transcript
bindsites_map_longestTX <- bindsites_maplong_tr(binding_sites=IGF2BPsbindingsites,annotation_file=GENE_ANNO_GTF,parclip=FALSE)
```
### *Mapping reader binding sites to the consistent peak sites in the longest transcript*
```r
###mapping the binding sites from PAR-CLIP to the consistent peak sites
reader_peak_overlap <- bindsites_mapto_peak(peak_sites_infor=peak_site_infor,mapped_peak_GR=map_consist_peak_longTX,
                                              bind_sites=bindsites_map_longestTX,parclip=TRUE)
###mapping the binding sites from iCLIP/eCLIP to the consistent peak sites  
reader_peak_overlap <- bindsites_mapto_peak(peak_sites_infor=peak_site_infor,mapped_peak_GR=map_consist_peak_longTX,
                                              bind_sites=bindsites_map_longestTX,parclip=FALSE)
```
### *Obtain high condifident peak sites by remove lower reads count in peak sites, which are bind or no-bind by reader*
```r
reader_bindor_nobind_peak <- reader_peak_overlap$consis_peak_infor
peaksites_filter <- bindornobind_gene_peakfilter(bind_nobindgene_peak=reader_bindor_nobind_peak,filter_reads_num=5)
##binding sites overlapp to filtered peak sites with reader binding sites
bindgene_bindpeak <- peaksites_filter$bindgene_bindpeak_filter
overlapped_bindsites <- reader_peak_overlap$binding_sites_overlap
##reader binding sites from parclip-seq data
bindsites_overlap_filterpeak <- mapped_filterpeak_bindsites(overlap_bind_sites=overlapped_bindsites,bindsites_peak=bindgene_bindpeak,parclip=TRUE)
##reader binding sites from eCLIP/iCLIP data
bindsites_overlap_filterpeak <- mapped_filterpeak_bindsites(overlap_bind_sites=overlapped_bindsites,bindsites_peak=bindgene_bindpeak,parclip=FALSE)
```
### *Obtain the distance between peak and reader binding sites or stop codon*
```r
##Get peak center
bindgene_nonbind_peaksite <- peaksites_filter$bindgene_nonbindpeak_filter
bindgene_nobindsite_peakcenter <- findpeakcenter(targetpeaks=bindgene_nonbind_peaksite,annotation_file=GENE_ANNO_GTF,maplongtx_peak=bindsites_map_longestTX)
nobindgene_peaksite <- peaksites_filter$nonbindgene_filter
nobindgene_peakcenter <- findpeakcenter(targetpeaks=nobindgene_peaksite,annotation_file=GENE_ANNO_GTF,maplongtx_peak=bindsites_map_longestTX)
##Obtain the min distance information to binding sites (single base)
###For PAR-CLIP-seq data
bindgene_nobind_peakdist <- dist_fun(overlap_bindsites_infor=bindsites_overlap_filterpeak,bindgene_nobind_peakcenter=bindgene_nobindsite_peakcenter,
                                      annotation_file=GENE_ANNO_GTF,parclip=TRUE)
###For eCLIP/iCLIP data
bindgene_nobind_peakdist <- dist_fun(overlap_bindsites_infor=bindsites_overlap_filterpeak,bindgene_nobind_peakcenter=bindgene_nobindsite_peakcenter,
                                      annotation_file=GENE_ANNO_GTF,parclip=FALSE)

##Obtain the distance to stop codon                                
nobind_gene_dist_stopcodon <- dist_stopcodon(target_peakcenter=nobindgene_peakcenter,annotation_file=GENE_ANNO_GTF)

```
### *Add binding signal stregnth for the bound peak sites*
```r
##Add reader binding singal strength for bound peak sites and add the distance information to reader binding gene, whose peak without bingding
bindgene_bindpeak <- peaksites_filter$bindgene_bindpeak_filter
bindgene_nobindpeak <- peaksites_filter$bindgene_nonbindpeak_filter
add_binding_strength_dist <- add_peak_SNR(bindgene_bindpeak=bindgene_bindpeak,
                                           bindgene_nobindpeak=bindgene_nobindpeak,
                                           bindgene_nobind_peakdist_infor=bindgene_nobind_peakdist,
                                           overlapped_bindsites=bindsites_overlap_filterpeak,parclip=TRUE)
```
### *Obtain the methylation level of reader binding gene or no binding gene*
```r
## No bind gene methylation level
nobindgene_methylevel <- nobindgene_gene_methy_level_distdecay(methy_site_infor=nobindgene_peaksite,library_size=peak_site_infor[[2]],
                                                               peak_dist_stopcodon=nobind_gene_dist_stopcodon)
##bind gene methylation level
bindgene_methylevel <- bindgene_methylevel (bindgene_peakSNR=add_binding_strength_dist,library_size=peak_site_infor[[2]])
```
### *Obtain gene reads count and differential expression gene*
```r
library(Rsubread)
library(DESeq2)
##Obtian gene reads count
f1 <- "./CTL_Input1.bam"
f2 <- "./CTL_Input2.bam"
f3 <- "./M3KD_Input1.bam"
f4 <- "./M3KD_Input2.bam"
Input_data <- c(f1,f2,f3,f4)
gene_readscount <- obtain_gene_readscount(Input_data=Input_data,annotation_file=GENE_ANNO_GTF,isPairedEnd=F)
##Obtain DE gene information
DE_infor <- DE_analysis(gene_readscount=gene_readscount,con_name1="Control",con_name2="METTL3KD",num_con1=2,num_con2=2,p_value=0.05)
```
### *Obtain differential methylation gene*
```r
library(limma)
DM_geneinfor <- obtain_DM_geneinfor(bindgene_methylevel=bindgene_methylevel,nobindgene_methylevel=nobindgene_methylevel,p_value=0.05,num_cond1=2,num_cond2=2)
```
### *Obtain candidate gene*
```r
candidate_gene_infor <- obtain_candidate_gene_infor(sig_DE_gene=DE_infor,sig_DM_gene=DM_geneinfor,output_path="./m6AexpressReader/")
```
### *Predicated m6A regulated expression gene by m6AexpressReader model*
```r
## Obtain the reader binding strength for reader binding gene
gene_bind_strength <- bindgene_SNR(bindgene_peak_SNR_infor=add_binding_strength_dist)
## Predicate the m6A-reg-exp gene by m6Aexpress-Reader model
m6Aregexpgene_m6AexpressReader <- m6Aexpress_Reader_model()
```
