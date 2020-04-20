options(stringsAsFactors = FALSE,
        scipen = 99,
        width = 999,
        knitr.table.format = "latex",
        tinytex.clean = TRUE)

source("/home/kevin/dev/ltrparser/R/misc.R")

## Not run:
#library(dnaplotr)
libsToLoad <- c("dplyr", "tidyr", "argparse", "ShortRead", "dnar", "stringdist") # etc.
nullObj <- lapply(libsToLoad, get_package, silent = TRUE)



# Gather command line arguments
parser <- ArgumentParser(
  usage = "Rscript ltrParser.R -c [config file]"
)

parser$add_argument(
  "-i", "--input", nargs = 1, type = "character"
)

parser$add_argument(
  "-o", "--output", nargs = 1, type = "character"
)

parser$add_argument(
  "-c", "--configFile", nargs = 1, type = "character"
)

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

# Main code body

inputFastq <- try_filepath(args$input, "input")
outputCsv <- try_dirpath(args$output, "output")
configFile <- try_filepath(args$configFile, "config")

outputPdf <- gsub(pattern = ".csv", replacement = ".pdf", x = outputCsv)

config <- read.table(configFile, sep = ":",
                     col.names = c("Parameters","Value"))
params <- config$Value
names(params) <- config$Parameters

ltrDatabaseFile <- params["LTR_Database"]
primerLength <- as.numeric(params["Primer_Length"])
minLTRlength <- as.numeric(params["Min_LTR_Length"])
maxLTRlength <- as.numeric(params["Max_LTR_Length"])
hostGenomeBowtieIndex <- params["Host_Genome_Bowtie_Index"]
hostBowtieOptions <- params["Host_Bowtie_Options"]
viralGenomeBowtieIndex <- params["Viral_Genome_Bowtie_Index"]
viralBowtieOptions <- params["Viral_Bowtie_Options"]
terminalSeq <- params["Terminal_Seq"]
rmdScript <- params["RMD_Script"]

ltrFasta <- readFasta(ltrDatabaseFile)
ltrDatabase <- as.character(sread(ltrFasta))
names(ltrDatabase) <- as.character(ShortRead::id(ltrFasta))

seqsDf <- data.frame(seq=as.character(ShortRead::sread(ShortRead::readFastq(inputFastq)))) %>%
  group_by(seq) %>%
  summarize(count=n()) %>%
  ungroup() %>%
  as.data.frame()

if(nrow(seqsDf)==0) {
  print(paste0("No sequencing reads found for file: ", inputFastq))
  q()
}


summaryStats <- data.frame("Total number of reads in sample"=sum(seqsDf$count),
                           "Number of distinct reads in sample"=nrow(seqsDf),
                           "Number of copies of most common read"=max(seqsDf$count))

hostBowtieResult <-
  system(paste0("bowtie2 --quiet -p 16 ", hostBowtieOptions, " -x ", hostGenomeBowtieIndex, " -U ",
                inputFastq, " | samtools view -S -F4 - "), intern = TRUE)

if(length(hostBowtieResult)==0) {
  print(paste0("No reads mapped to human genome for file: ", inputFastq))
  q()
}

hostBowtieResultDf <- parse_bowtie_results(hostBowtieResult)

seqsDf <- seqsDf %>%
  mutate(seq %in% hostBowtieResultDf$seq)
candidates <- identify_candidates(hostBowtieResultDf, primerLength, minLTRlength, maxLTRlength)
adjustedCandidates <- adjust_candidates(candidates, terminalSeq)
commonCandidates <- filter_candidates(adjustedCandidates)

if(nrow(commonCandidates)==0) {
  print(paste0("No common primer-LTR pairs found for file: ", inputFastq))
  render_output(rmdScript, outputPdf)
  q()
}

ltrsToAnalyze <- get_unique_ltrs(commonCandidates)

tempLTRfile <- tempfile()
writeFasta(object = ShortRead(DNAStringSet(ltrsToAnalyze), id = BStringSet(names(ltrsToAnalyze))),
           file = tempLTRfile, compress = FALSE)

viralBowtieResult <-
  system(paste0("bowtie2 --quiet -p 16 -f ", viralBowtieOptions," -x ", viralGenomeBowtieIndex, " -U ",
                tempLTRfile, " | samtools view -S -F4 - "), intern = TRUE)

outputTable <- update_candidates(commonCandidates, viralBowtieResult)
  
write.csv(outputTable, outputCsv, quote = FALSE, row.names = FALSE)

filteredSeqsDf <- filter_seqs(seqsDf, commonCandidates, primerLength)

render_output(rmdScript, outputPdf)

q()

