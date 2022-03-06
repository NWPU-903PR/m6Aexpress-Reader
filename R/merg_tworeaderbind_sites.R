mergbinding_sites <- function(one_bindingsites, two_bindingsites){
  IGF2BP1_clipsite <- read.table(one_bindingsites,header=F)
  IGF2BP1_clipsite <- read.table(f1,header=F)
  colnames(IGF2BP1_clipsite) <- c("seqname","start","end","state","score","strand","all_score")
  IGF2BP3_clipsite <- read.table(two_bindingsites,header = F)
  colnames(IGF2BP3_clipsite) <- c("seqname","start","end","state","score","strand","all_score")
  IGF2BP3_clipsite <- IGF2BP3_clipsite[which(is.na(match(IGF2BP3_clipsite$seqname,"chrM"))),]
  IGF2BP1_GR <- GRanges(seqnames = as.character(IGF2BP1_clipsite$seqname),
                        IRanges(start = as.numeric(as.character(IGF2BP1_clipsite$start)),
                                width = 1),strand = as.character(IGF2BP1_clipsite$strand),
                        score=as.numeric(as.character(IGF2BP1_clipsite$score)))
  IGF2BP3_GR <- GRanges(seqnames = as.character(IGF2BP3_clipsite$seqname),
                        IRanges(start = as.numeric(as.character(IGF2BP3_clipsite$start)),
                                width = 1),strand = as.character(IGF2BP3_clipsite$strand),
                        score=as.numeric(as.character(IGF2BP3_clipsite$score)))
  
  IGF2BP_ind <- findOverlaps(IGF2BP1_GR,IGF2BP3_GR,type = "equal")
  IGF2BP1_GR_overlap <- IGF2BP1_GR[unique(IGF2BP_ind@from)]
  IGF2BP3_GR_overlap <- IGF2BP3_GR[unique(IGF2BP_ind@to)]
  IGF2BP1_score <- exp(as.numeric(as.character(IGF2BP1_GR_overlap$score)))
  IGF2BP3_score <- exp(as.numeric(as.character(IGF2BP3_GR_overlap$score)))
  scores <- data.frame(IGF2BP1_score=IGF2BP1_score,IGF2BP3_score=IGF2BP3_score)
  Score <- rowSums(scores)
  score <- log(Score)
  IGF2BP13_common <-GRanges(seqnames = as.character(IGF2BP1_GR_overlap@seqnames),
                            IRanges(start = as.numeric(as.character(start(IGF2BP1_GR_overlap@ranges))),
                                    width = 1),strand = as.character(IGF2BP1_GR_overlap@strand),
                            score=score)
  
  IGF2BP1_GR_unique <- IGF2BP1_GR[-unique(IGF2BP_ind@from)]
  # IGF2BP1_GR_unique$score <- exp(as.numeric(as.character(IGF2BP1_GR_unique$score)))
  IGF2BP3_GR_unique <- IGF2BP3_GR[-unique(IGF2BP_ind@to)]
  # IGF2BP3_GR_unique$score <- exp(as.numeric(as.character(IGF2BP3_GR_unique$score)))
  IGF2BP13_clipsite <- rbind(as.data.frame(IGF2BP13_common),
                             as.data.frame(IGF2BP1_GR_unique),
                             as.data.frame(IGF2BP3_GR_unique))
  
  target_sites_GR <-  GRanges(seqnames = as.character(IGF2BP13_clipsite$seqname),
                              IRanges(start = as.numeric(as.character(IGF2BP13_clipsite$start)),
                                      width = 1),strand = as.character(IGF2BP13_clipsite$strand),
                              score=as.numeric(as.character(IGF2BP13_clipsite$score)))
  return(target_sites_GR)
}





