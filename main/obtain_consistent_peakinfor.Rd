\name{obtain_consistent_peakinfor}
\alias{obtain_consistent_peakinfor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain consistent peak sites information from exomePeak
}
\description{
Obatin each peak sites information including peak sites range, reads count of each peak.
}
\usage{
obtain_consistent_peakinfor(peak_file, peak_bed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{peak_file}{
a Rdata format file including the peak sites information obtained from \code{exomePeak} R package.
}
 \item{peak_bed}{
a bed format file including the consistent peak sites range across two conditions obtained from \code{exomePeak} R package.
}

}
\details{
For peak calling purpose, it requires the IP and input control samples:
An IP sample is the aligned BAM file from the immunoprecipitated sample using RNA modification antibodies such as anti-m6A; 
The input control sample is the aligned BAM file from the total RNAseq shotgun sequencing.
}
\value{
By default, \code{Get_peak_sites} will output results both 

1. A folder named "exomePeak_output" by default included peaks and consist peaks BED/XLS files on diskunder the specified directory (default: current working directory) and a GRangesList object under the R environment, including information generated in the peak calling process can be accessed with the "mcols" command.
2. A list object under current R environment, one includes the consistent peak sites information (e.g. "gene name","seqname","start","end","strand", reads count in each peak of each input samples), another included the libray sizes factor to normalize methylation level of each peak site in further analysis.
}
\references{
Meng, Jia, Xiaodong Cui, Manjeet K. Rao, Yidong Chen, and Yufei Huang. "Exome-based analysis for RNA epigenome sequencing data." Bioinformatics 29, no. 12 (2013): 1565-1567.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{exomepeak}}
}
\examples{
\dontrun{
## Peak calling for case-control context
f1 <- system.file("extdata", "IP1.bam", package="m6Aexpress")
f2 <- system.file("extdata", "IP2.bam", package="m6Aexpress")
f3 <- system.file("extdata", "IP3.bam", package="m6Aexpress")
f4 <- system.file("extdata", "IP4.bam", package="m6Aexpress")
f5 <- system.file("extdata", "Input1.bam", package="m6Aexpress")
f6 <- system.file("extdata", "Input2.bam", package="m6Aexpress")
f7 <- system.file("extdata", "Input3.bam", package="m6Aexpress")
f8 <- system.file("extdata", "Input4.bam", package="m6Aexpress")
IP_BAM <- c(f1,f2,f3,f4)
INPUT_BAM <- c(f5,f6,f7,f8)

# Input the annotation file
gtf <- system.file("extdata", "hg19toy.gtf", package="m6Aexpress")
# Get peak sites information
Get_peak_infor <- Get_peakinfor(IP_BAM, INPUT_BAM,GENE_ANNO_GTF=gtf)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
