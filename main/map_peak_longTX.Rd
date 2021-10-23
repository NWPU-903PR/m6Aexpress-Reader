\name{map_peak_longTX {m6AexpressReader}}
\alias{map_peak_longTX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Map peak sites to the longest transcript of each gene
}
\description{
  Map peak sites to the longest transcirpt of each gene based on transcriptome coordinates co to avoid isoform amibiguity.
}
\usage{
  map_peak_longTX(filepath, annotation_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filepath}{
    a character string giving the peak file path (in .bed format) obtained by \code{exomePeak} output.
  }
  \item{annotation_file}{
    A character string giving name of a user-provided annotation file in GTF format.
  }
  
}

\value{
  Results will return peak sites mapped to the longest transcript of genes in CompressedGRangesList format.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
  
  \seealso{
    \code{\link{find_longest_transcript}}
  }
\examples{
  \dontrun{
    GENE_ANNO_GTF = "./hg19_GTF/genes.gtf"
    consisten_peak <- "./exomePeak_output/con_peak.bed"
    map_consist_peak_longTX <- map_peak_longTX(filepath=consisten_peak,annotation_file=GENE_ANNO_GT)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
