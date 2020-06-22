library(stringr)
library(merf)

read_config_file <- function(configFile) {
  config <- read.table(configFile, sep = ":",
                       col.names = c("Parameters","Value"))
  
  params <- config$Value
  names(params) <- config$Parameters
  params <- as.list(params)
  
  params[["Primer_Length"]] <- as.numeric(params[["Primer_Length"]])
  params[["Min_LTR_Length"]] <- as.numeric(params[["Min_LTR_Length"]])
  params[["Max_LTR_Length"]] <- as.numeric(params[["Max_LTR_Length"]])
  params[["LTR_Max_Dist"]] <- as.numeric(params[["LTR_Max_Dist"]])
  params[["Min_Proportion"]] <- as.numeric(params[["Min_Proportion"]])
  params[["Map_To_LANL"]] <- as.logical(params[["Map_To_LANL"]])
  
  return(params)
}

get_linker_match_string <- function(seq, expectedLinkerRegexVector) {
  returnString <- unlist(sapply(expectedLinkerRegexVector, function(expectedLinkerRegex) {
    regExResult <- gregexpr(expectedLinkerRegex, seq)[[1]]
    if(regExResult[1] > 0) return(paste0(regExResult[1], ":", regExResult[1] + attr(regExResult, "match.length")[1] - 1))
    return("")
  }))
  returnString <- returnString[returnString != ""]
  if(length(returnString)==0) return("0:0")
  return(returnString)
}

extract_initial_soft_clip <- function(cigar, revStrand=FALSE) {
  if(revStrand) {
    result <- gregexpr("[0-9]+S$",cigar)[[1]]
    returnVal <- as.integer(substr(cigar, result[1], result[1] + attr(result, "match.length") - 2))
  } else {
    result <- gregexpr("^[0-9]+S",cigar)[[1]]
    returnVal <- as.integer(substr(cigar, 1, attr(result, "match.length") - 1))
  }
  if(is.na(returnVal) | is.null(returnVal)) returnVal <- 0
  return(returnVal)
}

extract_match_length <- function(cigar) {
  result <- gregexpr("[0-9]+M",cigar)[[1]]
  if(length(result)==0) return(0)
  returnVal <- sum(sapply(seq(1, length(result)), function(i) {
    matchStart <- result[i]
    matchLength <- attr(result, "match.length")[i]
    as.integer(substr(cigar, matchStart, matchStart + matchLength - 2))
  }))
  return(returnVal)
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
    rowwise() %>%
    mutate(seq=if_else(revStrand, revComp(seq), seq),
           initialSoftClip=extract_initial_soft_clip(cigar, revStrand),
           matchLength=extract_match_length(cigar)) %>%
    ungroup() %>%
    as.data.frame()
  return(bowtieResultDf)
}

summarize_seqs <- function(seqsDf) {
  if(!"human" %in% colnames(seqsDf)) {
    seqsDf <- mutate(seqsDf, human=FALSE)
  }
  if(!"hasLinker" %in% colnames(seqsDf)) {
    seqsDf <- mutate(seqsDf, hasLinker=FALSE)
  }
  summaryStats <- seqsDf %>%
    summarize(Total=sum(count),
              Distinct=n(),
              Max=max(count),
              TotalLinker=sum(count[hasLinker]),
              DistinctLinker=sum(hasLinker),
              MaxLinker=ifelse(TotalLinker > 0, max(count[hasLinker]), 0),
              TotalHuman=sum(count[human]),
              DistinctHuman=sum(human),
              MaxHuman=ifelse(DistinctHuman > 0, max(count[human]), 0))
  print(paste0("Processed ", summaryStats$Total, " reads."))
  return(summaryStats)
}

identify_candidates <- function(bowtieResultDf, primerLength, minLTRlength, maxLTRlength) {
  bowtieResultDf %>%
    group_by(seq, cigar, initialSoftClip) %>%
    summarize(count=n()) %>%
    filter(initialSoftClip >= primerLength + minLTRlength &
             initialSoftClip <= primerLength + maxLTRlength) %>%
    mutate(primer=substr(seq, 1, primerLength),
           LTR=substr(seq, primerLength + 1, initialSoftClip)) %>%
    as.data.frame() %>%
    group_by(primer, LTR) %>%
    summarize(count=sum(count)) %>%
    ungroup() %>%
    as.data.frame()
}


adjust_candidates <- function(candidates, minLTRlength, maxLTRlength, terminalSeq) {
  if(is.null(candidates)) return(NULL)
  if(nrow(candidates)==0) return(candidates)
  adjustedCandidatesDf <- do.call(rbind, lapply(unique(candidates$LTR), function(LTR) {
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
      data.frame(LTR=LTR, adjustedLTR=c(truncatedLTR, extendedLTR)) %>%
        filter(nchar(adjustedLTR) >= minLTRlength & nchar(adjustedLTR) <= maxLTRlength)
    }
  }))
  return(candidates %>%
           left_join(adjustedCandidatesDf, by="LTR") %>%
           filter(!is.na(adjustedLTR)) %>%
           select(-LTR) %>%
           dplyr::rename(LTR=adjustedLTR) %>%
           distinct() %>%
           group_by(primer, LTR) %>%
           summarize(count=sum(count)) %>%
           ungroup() %>%
           as.data.frame())
}

filter_candidates <- function(candidates, minCount = 10, minProportion = 0.01) {
  if(is.null(candidates)) return(NULL)
  if(nrow(candidates)==0) return(candidates)
  return(candidates %>%
           filter(count > minCount | count > sum(count) * minProportion) %>%
           arrange(-count) %>%
           distinct() %>%
           head(10))
}

preliminary_candidates <- function(hostBowtieResultDf, primerLength, minLTRlength, maxLTRlength, terminalSeq) {
  # primerLength <- config[["Primer_Length"]]
  # minLTRlength <- config[["Min_LTR_Length"]]
  # maxLTRlength <- config[["Max_LTR_Length"]]
  # terminalSeq <- config[["Terminal_Seq"]]
  candidates <- identify_candidates(hostBowtieResultDf, primerLength, minLTRlength, maxLTRlength)
  candidates <- adjust_candidates(candidates, minLTRlength, maxLTRlength, terminalSeq)
  candidates <- filter_candidates(candidates)
  return(candidates)
}

recover_seqs <- function(seqsDf, candidates, expectedLTR, LTRmaxDist = 2, terminalSeq) {
  process_pair <- function(pair) {
    primer <- str_split(pair, ":")[[1]][1]
    LTR <- str_split(pair, ":")[[1]][2]
    primerLength <- nchar(primer)
    LTRlength <- nchar(LTR)
    seqsDf %>%
      filter(substr(seq, 1, primerLength)==primer) %>%
      mutate(LTRdist=stringdist(substr(seq, primerLength + 1, primerLength + LTRlength), LTR)) %>%
      filter(LTRdist <= LTRmaxDist & substr(seq,
                                            primerLength + LTRlength - nchar(terminalSeq) + 1,
                                            primerLength + LTRlength) == terminalSeq) %>%
      mutate(primer=primer,
             LTR=LTR,
             expectedLTR=LTR==expectedLTR) %>%
      as.data.frame()
  }
  recoveredSeqsDf <- do.call(rbind, lapply(unite(candidates, pair, c("primer","LTR"), sep = ":") %>% pull(pair), process_pair))
}

segment_seqs <- function(recoveredSeqsDf, hostBowtieResultDf, expectedLinkerRegexVector) {
  result <- recoveredSeqsDf %>%
    left_join(hostBowtieResultDf %>%
                group_by(seq, matchLength, initialSoftClip) %>%
                summarize(count=n()) %>%
                group_by(seq) %>%
                arrange(-count, -matchLength) %>%
                mutate(row=seq(1:n())) %>%
                filter(row==1) %>%
                ungroup() %>%
                select(-count, -row), by="seq") %>%
    rowwise() %>%
    mutate(LTRend=nchar(primer)+nchar(LTR),
           LTRsegment=substr(seq, nchar(primer)+1, LTRend),
           matchLength=if_else(is.na(matchLength), as.integer(0), matchLength),
           originalSoftClip=initialSoftClip,
           initialSoftClip=if_else(is.na(initialSoftClip), 0, as.numeric(initialSoftClip)),
           matchEnd=initialSoftClip + matchLength,
           initialSoftClip=max(initialSoftClip, nchar(primer) + nchar(LTR)),
           matchStart=if_else(matchEnd > initialSoftClip, initialSoftClip + 1,0),
           matchLength=if_else(matchStart > 0, matchEnd - matchStart + 1, 0),
           human=matchLength > 0,
           remainder_u5=if_else(human & initialSoftClip > LTRend,
                                substr(seq, LTRend + 1, initialSoftClip),
                                ""),
           matchSeq=substr(seq, matchStart, matchEnd),
           remainder_u3_start=max(matchEnd+1, LTRend+1),
           remainder_u3=substr(seq, remainder_u3_start, nchar(seq)),
           linkerMatchString=get_linker_match_string(remainder_u3, expectedLinkerRegexVector)) %>%
    separate(linkerMatchString, c("linkerMatchStart", "linkerMatchEnd"), sep = ":") %>%
    mutate(linkerMatchStart=as.numeric(linkerMatchStart),
           linkerMatchEnd=as.numeric(linkerMatchEnd),
           hasLinker=linkerMatchStart > 0,
           linkerSeq=if_else(hasLinker, remainder_u3, ""),
           remainder_u3=if_else(hasLinker, "", remainder_u3)) %>%
    ungroup() %>%
    as.data.frame()
  return(result)
}

summarize_candidates <- function(segmentedSeqsDf, seqsDf, hostBowtieResultDf, expectedLinkerRegexVector,
                                 minCount = 100, minProportion = 0.05) {
  missingSeqsDf <- seqsDf %>%
    filter(!seq %in% segmentedSeqsDf$seq) %>%
    rowwise() %>%
    mutate(human=seq %in% hostBowtieResultDf$seq,
           linkerMatchString=get_linker_match_string(seq, expectedLinkerRegexVector)) %>%
    separate(linkerMatchString, c("linkerMatchStart", "linkerMatchEnd"), sep = ":") %>%
    mutate(linkerMatchStart=as.numeric(linkerMatchStart),
           hasLinker=linkerMatchStart > 0) %>%
    select(seq, count, human, hasLinker) %>%
    ungroup() %>%
    as.data.frame()
  presentSeqsDf <- seqsDf %>%
    inner_join(segmentedSeqsDf %>%
                 group_by(seq) %>%
                 summarize(human=any(human), hasLinker=any(hasLinker)) %>%
                 ungroup() %>%
                 as.data.frame(), by="seq")
  readSummary <- summarize_seqs(rbind(missingSeqsDf, presentSeqsDf)) %>%
    mutate(primer="", LTR="", expectedLTR="", id=0)
    
  candidateSummary <- segmentedSeqsDf %>%
    group_by(primer, LTR, expectedLTR) %>%
    summarize(Total=sum(count),
              Distinct=n(),
              Max=max(count),
              TotalLinker=sum(count[hasLinker]),
              DistinctLinker=sum(hasLinker),
              MaxLinker=ifelse(TotalLinker > 0, max(count[hasLinker]), 0),
              TotalHuman=sum(count[human]),
              DistinctHuman=sum(human),
              MaxHuman=ifelse(DistinctHuman > 0, max(count[human]), 0)) %>%
    ungroup()
  if(nrow(candidateSummary %>% filter(Total >= minCount | Total >= readSummary$Total * minProportion)) == 0) {
    return(rbind(candidateSummary %>%
                   arrange(-Total) %>%
                   head() %>%
                   mutate(id=seq(1:n())) %>%
                   as.data.frame(),
                 readSummary))
    }
  return(candidateSummary %>%
           filter(Total >= minCount | Total >= readSummary$Total * minProportion) %>%
           arrange(-Total) %>%
           mutate(id=seq(1:n())) %>%
           as.data.frame() %>%
           rbind(readSummary))
}

name_candidates <- function(candidates, minProportion = 0.05) {
  if(is.null(candidates)) return(NULL)
  if(nrow(candidates)==0) return(candidates)
  return(candidates %>%
           mutate(qname=as.character(as.integer(factor(LTR)))) %>%
           distinct())
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

create_display_seq <-
  function(primer, LTRsegment, remainder_u5, matchSeq, linkerSeq, remainder_u3,
           maxPrimerLength=nchar(primer), maxLTRlength=nchar(LTRsegment),
           maxRemainderU5Length=nchar(remainder_u5), maxMatchLength=nchar(matchSeq),
           maxLinkerLength=nchar(linkerSeq), maxRemainderU3Length=nchar(remainder_u3),
           buffer=5) {
    return(paste0(c(primer, rep("N", maxPrimerLength - nchar(primer) + buffer),
                    substr(LTRsegment, 1, nchar(LTRsegment) - 2), rep("N", maxLTRlength - nchar(LTRsegment) + buffer),
                    substr(LTRsegment, nchar(LTRsegment)-1, nchar(LTRsegment)), rep("N", buffer),
                    remainder_u5, rep("N", maxRemainderU5Length - nchar(remainder_u5) + buffer),
                    matchSeq, rep("N", maxMatchLength - nchar(matchSeq) + buffer),
                    linkerSeq, rep("N", maxLinkerLength - nchar(linkerSeq) + buffer),
                    remainder_u3,  rep("N", maxRemainderU3Length - nchar(remainder_u3))),
                  collapse = ""))
  }
