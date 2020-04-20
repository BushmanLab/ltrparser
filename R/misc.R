library(stringr)
library(merf)

extract_initial_soft_clip <- function(cigar, revStrand=FALSE) {
  if(revStrand) {
    result <- gregexpr("[0-9]+S$",cigar)[[1]]
    returnVal <- as.integer(substr(cigar, result[1], result[1] + attr(result, "match.length") - 2))
  } else {
    result <- gregexpr("^[0-9]+S",cigar)[[1]]
    returnVal <- as.integer(substr(cigar, 1, attr(result, "match.length") - 1))
  }
  return(returnVal)
}

identify_candidates <- function(bowtieResultDf, primerLength, minLTRlength, maxLTRlength) {
  bowtieResultDf %>%
    group_by(seq, cigar, revStrand) %>%
    summarize(count=n()) %>%
    rowwise() %>%
    mutate(initialSoftClip=extract_initial_soft_clip(cigar, revStrand)) %>%
    filter(initialSoftClip >= minLTRlength & initialSoftClip <= maxLTRlength) %>%
    mutate(primer=substr(seq, 1, primerLength),
           LTR=substr(seq, primerLength + 1, initialSoftClip)) %>%
    as.data.frame() %>%
    group_by(primer, LTR) %>%
    summarize(count=sum(count)) %>%
    ungroup() %>%
    as.data.frame()
}

adjust_candidates <- function(candidates, terminalSeq) {
  if(is.null(candidates)) return(NULL)
  if(nrow(candidates)==0) return(candidates)
  adjustedCandidatesDf <- do.call(rbind, lapply(candidates$LTR, function(LTR) {
    if(str_sub(LTR, -2)==terminalSeq) {
      data.frame(LTR=LTR, adjustedLTR=LTR)
    } else {
      if(str_detect(LTR, terminalSeq)) {
        truncatedLTRlength <- max(stringr::str_locate_all(LTR, terminalSeq)[[1]][,2])
        truncatedLTR <- substr(LTR,1,truncatedLTRlength)
      } else {
        truncatedLTR <- NULL
      }
      if(nchar(terminalSeq)==2 & str_sub(LTR, -1)==substr(terminalSeq, 1, 1)) {
        extendedLTR <- paste0(LTR, substr(terminalSeq, 2, 2))
      } else {
        extendedLTR <- paste0(LTR, terminalSeq)
      }
      data.frame(LTR=LTR, adjustedLTR=c(truncatedLTR, extendedLTR))
    }
  }))
  return(candidates %>% left_join(adjustedCandidatesDf, by="LTR") %>% select(-LTR) %>% dplyr::rename(LTR=adjustedLTR))
}

filter_candidates <- function(candidates) {
  if(is.null(candidates)) return(NULL)
  if(nrow(candidates)==0) return(candidates)
  return(candidates %>%
           filter(count > 1 & count > sum(count)/20) %>%
           arrange(-count) %>%
           mutate(qname=as.character(as.integer(factor(LTR)))))
}

get_unique_ltrs <- function(candidatesDf) {
  df <- candidatesDf %>%
    distinct(LTR, qname)
  ltrs <- df$LTR
  names(ltrs) <- df$qname
  return(ltrs)
}

update_candidates <- function(candidatesDf, bowtieResult) {
  tryCatch(
    {
      bowtieResultDf <- parse_bowtie_results(bowtieResult)
      outputTable <- candidatesDf %>%
        left_join(bowtieResultDf %>% select(qname, MD), by="qname") %>%
        mutate(MD=gsub("MD:Z:","",MD),
               LANL=!is.na(MD)) %>%
        ungroup()
    },
    error=function(cond) {
      candidatesDf %>%
        mutate(MD=NA,
               LANL=FALSE) %>%
        ungroup()
    })
}

parse_bowtie_results <- function(bowtieResult) {
  bowtieDfHeader <- c("qname","flag","rname","pos","mapq","cigar","rnext","pnext","tlen","seq","qual",
                      "AS","XN","XM","XO","XG","NM","MD","YT")
  bowtieResultDf <- do.call(rbind, strsplit(bowtieResult, "\t")) %>%
    as.data.frame()
  colnames(bowtieResultDf) <- bowtieDfHeader
  bowtieResultDf$revStrand <- sapply(bowtieResultDf$flag, function(flag) {
    intToBits(as.numeric(flag))[5]==1
  })
  bowtieResultDf <- bowtieResultDf %>%
    mutate(seq=if_else(revStrand, revComp(seq), seq))
}

filter_seqs <- function(seqsDf, candidates, primerLength) {
  do.call(rbind,
          lapply(paste0(candidates$primer, ":", candidates$LTR),
                 function(primerLTRpair) {
                   junk <- strsplit(primerLTRpair, split = ":")[[1]]
                   primer <- junk[1]
                   LTR <- junk[2]
                   seqsDf %>%
                     filter(substr(seq, 1, primerLength)==primer) %>%
                     mutate(primer=primer,
                            LTR=LTR,
                            LTRdist=stringdist(substr(seq, primerLength + 1, primerLength + nchar(LTR)), LTR)) %>%
                     filter(LTRdist <= 2)
                 }))
}




get_closest_ltr_name <- function(ltr, db) {
  ltr_distances <- stringdistmatrix(ltr, db, useNames = TRUE)
  ltr_closest_name <- names(db)[which.min(ltr_distances)]
  ltr_closest_seq <- db[ltr_closest_name]
  ltr_closest_dist <- ltr_distances[1,ltr_closest_seq]
  
  rc_distances <- stringdistmatrix(revComp(ltr), db, useNames = TRUE)
  rc_closest_name <- names(db)[which.min(rc_distances)]
  rc_closest_seq <- db[rc_closest_name]
  rc_closest_dist <- rc_distances[1,rc_closest_seq]

  if(rc_closest_dist < ltr_closest_dist) {
    return(rc_closest_name)
  } else {
    return(ltr_closest_name)
  }
}

