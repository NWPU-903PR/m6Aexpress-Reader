# Usage Example
## Step by Step Analysis
### For two (case-control) conditions
#### *Peak calling for methylation sites in case-control context*
```r
library(exomePeak)
library(GenomicFeatures)
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

```
