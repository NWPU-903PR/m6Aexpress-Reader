\name{bindsites_maplong_tr}
\alias{bindsites_maplong_tr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Map the reader binding site to longest transcirpt
}
\description{
Map reader binding site to the longest transcirpt of each gene based on transcriptome coordinates co to avoid isoform amibiguity.
}
\usage{
  bindsites_maplong_tr(binding_sites,annotation_file,parclip=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{binding_sites}{
    If the reader binding sites are obtained by \code{wavClusteR} from PAR-CLIP-seq data, the list file \code{binding_sites} includes the reader binding sites in GRange format and reader binding cluster with reader binding streng in dataframe format. If the reader binding sites are obtained by \code{PureCLIP} from eCLIP/iCLIP data, \code{binding_sites} includes the reader binding sites with binding strength in GRanges format.
  }
  \item{annotation_file}{
    A character string giving name of a user-provided annotation file in GTF format.
  }
  \item{parclip}{
    logical indicating whether the reader binding sites from PAR-CLIP-seq data (\code{parclip=TRUE}) or from eCLIP/iCLIP (\code{parclip=FALSE}) 
  }
}

\value{
  If \code{parclip=TRUE}, the result will return a list, which includes the reader binding sites mapped to the longest transcript in GRangesList format and its corresponding reader binding cluster in dataframe format. If \code{parclip=FALSE}, the result will return a GRangesList, which includes the reader binding sites with reader binding strength.
}

\author{
  Teng Zhang <tengzhang156@126.com>
}

\examples{
  \dontrun{
bindsites_map_longestTX <- bindsites_maplong_tr(binding_sites=obtain_reader_bindingsites,annotation_file=GENE_ANNO_GTF,parclip=TRUE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
