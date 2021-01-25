options(stringsAsFactors = FALSE,
        scipen = 99,
        width = 999,
        knitr.table.format = "latex",
        tinytex.clean = TRUE)

source("/home/kevin/dev/ltrparser/R/misc.R")

## Not run:
#library(dnaplotr)
#libsToLoad <- c("dplyr", "tidyr", "argparse", "ShortRead", "dnar", "stringdist", "rmarkdown")
libsToLoad <- c("dplyr", "tidyr", "ShortRead", "dnar", "stringdist", "rmarkdown", "ggplot2", "msa")
nullObj <- lapply(libsToLoad, get_package, silent = TRUE)


# Main code body

parse_ltrs <- function(inputFastq, outputCsv, configFile, expectedLTR = "", expectedLinker = "", description = "") {
  
  inputFastq <- try_filepath(inputFastq, "input")
  outputCsv <- try_dirpath(outputCsv, "output")
  configFile <- try_filepath(configFile, "config")

  expectedLinkerRegexVector <- gsub("N", ".", revComp(expectedLinker))
  print(expectedLinkerRegexVector)
  
  outputPdf <- gsub(pattern = "csv", replacement = "pdf", x = outputCsv)
  
  config <- read_config_file(configFile)

  expectedLTR <- paste0(expectedLTR, config[["Terminal_Seq"]], collapse = "")
  
  fastq <- ShortRead::readFastq(inputFastq)
  
  #fastq.trimmed <- trimTailw(fastq, 2, '<', 5)
  #seqsDf <- data.frame(seq=as.character(sread(fastq.trimmed))) %>%
  seqsDf <- data.frame(seq=as.character(sread(fastq))) %>%
    group_by(seq) %>%
    summarize(count=n()) %>%
    as.data.frame()
  
  if(nrow(seqsDf)==0) {
    print(paste0("No sequencing reads found for file: ", inputFastq))
    return(0)
  }
  
  hostBowtieResult <-
    system(paste0("bowtie2 --quiet -p 16 ", config[["Host_Bowtie_Options"]],
                  " -x ", config[["Host_Genome_Bowtie_Index"]], " -U ",
                  inputFastq, " | samtools view -S -F4 - "), intern = TRUE)
  
  if(length(hostBowtieResult)==0) {
    print(paste0("No reads mapped to human genome for file: ", inputFastq))
    seqsDf <- seqsDf %>%
      rowwise() %>%
      mutate(linkerMatchString=get_linker_match_string(seq, expectedLinkerRegexVector)) %>%
      separate(linkerMatchString, c("linkerMatchStart", "linkerMatchEnd"), sep = ":") %>%
      mutate(linkerMatchStart=as.numeric(linkerMatchStart),
             linkerMatchEnd=as.numeric(linkerMatchEnd),
             hasLinker=linkerMatchStart > 0)
    summaryStats <- summarize_seqs(seqsDf)
    render_output(config[["RMD_Script"]],
                  outputPdf,
                  params = list(inputFastq=inputFastq,
                                configFile=configFile,
                                seqsDf=seqsDf,
                                summaryStats=summaryStats))
    return(0)
  }
  
  hostBowtieResultDf <- parse_bowtie_results(hostBowtieResult)
  candidates <- preliminary_candidates(seqsDf,
                                       hostBowtieResultDf,
                                       config[["Primer_Length"]],
                                       config[["Min_LTR_Length"]],
                                       config[["Max_LTR_Length"]],
                                       config[["Terminal_Seq"]])
  
  mappingLengthDistribution <- hostBowtieResultDf %>%
    group_by(matchLength) %>%
    summarize(count=n()) %>%
    ungroup() %>%
    as.data.frame()
  
  if(nrow(candidates)==0) {
    print(paste0("No common primer-LTR pairs found for file: ", inputFastq))
    summaryStats <- summarize_seqs(seqsDf)
    render_output(config[["RMD_Script"]],
                  outputPdf,
                  params = list(inputFastq=inputFastq,
                                configFile=configFile,
                                seqsDf=seqsDf,
                                summaryStats=summaryStats,
                                mappingLengthDistribution=mappingLengthDistribution))
    return(0)
  }
  
  recoveredSeqsDf <- recover_seqs(seqsDf, candidates, expectedLTR, config[["LTR_Max_Dist"]], config[["Terminal_Seq"]])
  
  if(nrow(recoveredSeqsDf)==0) {
    print(paste0("Common primer-LTR pairs somehow not actually common: ", inputFastq))
    summaryStats <- summarize_seqs(seqsDf)
    render_output(config[["RMD_Script"]],
                  outputPdf,
                  params = list(inputFastq=inputFastq,
                                configFile=configFile,
                                seqsDf=seqsDf,
                                summaryStats=summaryStats,
                                mappingLengthDistribution=mappingLengthDistribution))
    return(0)
  }
  
  segmentedSeqsDf <- segment_seqs(recoveredSeqsDf, hostBowtieResultDf, expectedLinkerRegexVector)
  
  candidateSummary <- summarize_candidates(segmentedSeqsDf, seqsDf, hostBowtieResultDf, expectedLinkerRegexVector)
  
  if(!config[["Map_To_LANL"]]) {
    # This block runs for gene therapy samples or other samples where we don't want to map against LANL
    write.csv(candidates, outputCsv, quote = FALSE, row.names = FALSE)
    render_output(config[["RMD_Script"]],
                  outputPdf,
                  params = list(inputFastq=inputFastq,
                                configFile=configFile,
                                seqsDf=seqsDf,
                                mappingLengthDistribution=mappingLengthDistribution,
                                candidates=candidates,
                                candidateSummary=candidateSummary,
                                segmentedSeqsDf=segmentedSeqsDf))
    return(0)
  }
  
  namedCandidates <- name_candidates(candidates)
  
  ltrsToAnalyze <- get_unique_ltrs(namedCandidates)
  
  tempLTRfile <- tempfile()
  writeFasta(object = ShortRead(DNAStringSet(ltrsToAnalyze), id = BStringSet(names(ltrsToAnalyze))),
             file = tempLTRfile, compress = FALSE)
  
  viralBowtieResult <-
    system(paste0("bowtie2 --quiet -p 16 -f ", config[["Viral_Bowtie_Options"]],
                  " -x ", config[["Viral_Genome_Bowtie_Index"]], " -U ",
                  tempLTRfile, " | samtools view -S -F4 - "), intern = TRUE)
  
  candidates <- update_candidates(namedCandidates, viralBowtieResult)
  
  candidateSummary <- candidateSummary %>%
    left_join(select(candidates, primer, LTR, LANL), by=c("primer","LTR"))
  
  write.csv(candidates, outputCsv, quote = FALSE, row.names = FALSE)
  
  print("End of script, calling script to generate PDF.")
  render_output(config[["RMD_Script"]],
                outputPdf,
                params = list(inputFastq=inputFastq,
                              configFile=configFile,
                              expectedLTR=expectedLTR,
                              expectedLinker=expectedLinker,
                              seqsDf=seqsDf,
                              mappingLengthDistribution=mappingLengthDistribution,
                              candidates=candidates,
                              candidateSummary=candidateSummary,
                              segmentedSeqsDf=segmentedSeqsDf,
                              description=description))
  
  
  msaSeqsDf <- candidateSummary %>% filter(nchar(LTR) > 0) %>% arrange(-id) %>%
    mutate(id=if_else(expectedLTR=="TRUE", paste0(id, " (HXB2)"), as.character(id)))
  msaSeqs <- msaSeqsDf$LTR
  names(msaSeqs) <- msaSeqsDf$id
  if(!expectedLTR %in% msaSeqs) {
    names(expectedLTR) <- "HXB2"
    msaSeqs <- c(msaSeqs, expectedLTR)
  }
  MSA <- msa::msa(DNAStringSet(msaSeqs), order = "input")
  
  msatex = gsub("R2.pdf", "R2.msa.tex", outputPdf)
  msapdf = gsub("R2.pdf", "R2.msa.pdf", outputPdf)
  tmppdf = gsub("R2.pdf", "R2.tmp.pdf", outputPdf)
  msaPrettyPrint(MSA, output = "tex", file = msatex, verbose = FALSE, askForOverwrite = FALSE,
                 paperWidth = 8.5, paperHeight = 11)
  wd <- getwd()
  setwd(dirname(outputPdf))
  tools::texi2pdf(msatex, clean = TRUE)
  setwd(wd)
  
  pdftools::pdf_combine(input = c(outputPdf, msapdf), output = tmppdf)
  system(paste0("cp ", tmppdf, " ", outputPdf))
  system(paste0("rm ", tmppdf))
  system(paste0("rm ", msatex))
  system(paste0("rm ", msapdf))

  return(0)
}

