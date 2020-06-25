source("ltr-parser.R")
configFile <- "../config/hiv.config.tsv"
testMetadataFile <- "../testdata/testdata.tsv"
testMetadata <- read.table(file = testMetadataFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
for(testId in testMetadata$TestID) {
  print(Sys.time())
  metadataRow <- testMetadata %>% filter(TestID==testId)
  inputFastq <- metadataRow$fastq
  outputCsv <- gsub("fastq.gz", "csv", inputFastq)
  expectedLTR <- metadataRow$LTR
  expectedLinker <- metadataRow$Linker
  description = paste0(metadataRow$Project, " project, sample: ", metadataRow$Sample, " (", metadataRow$Region, ")")
  parse_ltrs(inputFastq, outputCsv, configFile, expectedLTR, expectedLinker, description)
  print(Sys.time())
}


