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
Results will return a list including the peak sites information in data.frame format and library sizes (total reads count in each sample) in vector format
  \item{peak_sites_infor}{A dataframe including the consistent peak sites across two conditions}
  \item{library_sizes}{A integer vector, which specifices the total reads count of each sample}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{exomepeak}}
}
\examples{
\dontrun{
load("./exomePeak_output/exomePeak.Rdata")
peak_file <- tmp_rs
consisten_peak <- "./exomePeak_output/con_peak.bed"
peak_site_infor <- obtain_consistent_peakinfor(peak_file, peak_bed=consisten_peak)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
