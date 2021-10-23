\name{parclip_reader_bindingsites}
\alias{parclip_reader_bindingsites}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain the reader binding sites from PAR-CLIP-seq data
}
\description{
Obtain reader binding sites in single base and binding cluster from PAR-CLIP-seq data.
}
\usage{
parclip_reader_bindingsites(par_bam, annotation_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par_bam}{
A character string giving file path of a user-provided PAR-CLIP-seq BAM file.
}
  \item{annotation_file}{
A character string giving name of a user-provided annotation file in GTF format.
}
}

\value{
Results will return a list including the reader binding sites in signle base and reader binding cluster with reader binding strength
  \item{high_confi_TC }{The high condifident T to C site, which denotes the reader binding site in single base}
  \item{bind_cluster_infor }{reader binding cluster information, which include the reader binding strength}
}
\references{
Sievers C, Schlumpf T, Sawarkar R, Comoglio F, Paro R (2012). "Mixture models and wavelet transforms reveal high confidence RNA-protein interaction sites in MOV10 PAR-CLIP data." 
\emph{Nucleic Acids Research}, 40(20), e160 doi: \url{10.1093/nar/gks697}.
\url{http://www.ncbi.nlm.nih.gov/pubmed/?term=sievers+comoglio}.
}
\author{
Teng Zhang <tengzhang156@126.com>
}

\examples{
  \dontrun{
  YTHDF2_binding <- "./YTHDF2_binding.bam"
  GENE_ANNO_GTF = "./hg19_GTF/genes.gtf"
  obtain_reader_bindingsites <- parclip_reader_bindingsites(par_bam=YTHDF2_binding,annotation_file=GENE_ANNO_GTF)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
