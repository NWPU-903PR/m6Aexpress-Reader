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
load("./exomePeak_calling/exomePeak_output/exomePeak.Rdata")
peak_file <- tmp_rs
consisten_peak <- "./exomePeak_calling/exomePeak_output/con_peak.bed"
peak_site_infor <- peak_infor(peak_file, peak_bed=consisten_peak,con_peak=TRUE)
##Mapping peak sites to the longest transcript 
map_consist_peak_longTX <- map_peak_longTX(filepath=consisten_peak,annotation_file=GENE_ANNO_GTF)
```
### *Identify reader binding sites from CLIP-seq data*
#### *For PAR-CLIP-seq data*
```r
library(stringr)
library(wavClusteR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicFeatures)
library(m6ALogisticModel)
reader_binding <- "./reader_binding.bam"
obtain_reader_bindingsites <- reader_bindingsites(par_bam=YTHDF2_binding,annotation_file=GENE_ANNO_GTF)
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
reader_peak_overlap <- bindsites_overlap_peak(peak_sites_infor=peak_site_infor,mapped_peak_GR=map_consist_peak_longTX,
                                              bind_sites=bindsites_map_longestTX,parclip=TRUE)
###mapping the binding sites from iCLIP/eCLIP to the consistent peak sites  
reader_peak_overlap <- bindsites_overlap_peak(peak_sites_infor=peak_site_infor,mapped_peak_GR=map_consist_peak_longTX,
                                              bind_sites=bindsites_map_longestTX,parclip=FALSE)
```
### *Obtain high condifident peak sites by remove lower reads count in peak sites, which are bind or no-bind by reader*
```r
reader_bindor_nobind_peak <- reader_peak_overlap$consis_peak_infor
peaksites_filter <- bindornobind_gene_peakfilter(bind_nobindgene_peak=reader_bindor_nobind_peak,filter_reads_num=5)
```
### *Obtain the distance between peak and reader binding sites or stop codon
```r
##Get peak center
bindgene_nonbind_peaksite <- peaksites_filter$bindgene_nonbindpeak_filter
bindgene_nobindsite_peakcenter <- findpeakcenter(targetpeaks=bindgene_nonbind_peaksite,annotation_file=GENE_ANNO_GTF,maplongtx_peak=bindsites_map_longestTX)
nobindgene_peaksite <- peaksites_filter$nonbindgene_filter
nobindgene_peakcenter <- findpeakcenter(targetpeaks=nobindgene_peaksite,annotation_file=GENE_ANNO_GTF,maplongtx_peak=bindsites_map_longestTX)

```
