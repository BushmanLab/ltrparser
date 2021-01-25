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


source("/home/kevin/dev/ltrparser/R/ltr-parser.R")
configFile <- "/home/kevin/dev/ltrparser/config/genetherapy.config.tsv"
testMetadataFile <- "/home/kevin/dev/ltrparser/testdata/testdata2.tsv"
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

# Jones

source("/home/kevin/dev/ltrparser/R/ltr-parser.R")
configFile <- "/home/kevin/dev/ltrparser/config/hiv.config.tsv"
testMetadataFile <- "/home/kevin/dev/ltrparser/testdata/jones_20200928.tsv"
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

# Wistar

source("/home/kevin/dev/ltrparser/R/ltr-parser.R")
configFile <- "/home/kevin/dev/ltrparser/config/hiv.config.tsv"
testMetadataFile <- "/home/kevin/dev/ltrparser/testdata/wistar_20201004.tsv"
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

source("/home/kevin/dev/ltrparser/R/ltr-parser.R")
configFile <- "/home/kevin/dev/ltrparser/config/hiv.config.tsv"
testMetadataFile <- "/home/kevin/dev/ltrparser/testdata/wistar_20201117.tsv"
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

source("/home/kevin/dev/ltrparser/R/ltr-parser.R")
configFile <- "/home/kevin/dev/ltrparser/config/hiv.config.tsv"
testMetadataFile <- "/home/kevin/dev/ltrparser/testdata/persaud_20201120.tsv"
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

source("/home/kevin/dev/ltrparser/R/ltr-parser.R")
configFile <- "/home/kevin/dev/ltrparser/config/hiv.config.GTSP3795.tsv"
testMetadataFile <- "/home/kevin/dev/ltrparser/testdata/persaud_20201123.tsv"
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

source("/home/kevin/dev/ltrparser/R/ltr-parser.R")
configFile <- "/home/kevin/dev/ltrparser/config/hiv.config.GTSP3793.tsv"
testMetadataFile <- "/home/kevin/dev/ltrparser/testdata/persaud_20201124.tsv"
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

source("/home/kevin/dev/ltrparser/R/ltr-parser.R")
configFile <- "/home/kevin/dev/ltrparser/config/hiv.config.GTSP3788.tsv"
testMetadataFile <- "/home/kevin/dev/ltrparser/testdata/persaud_20201125.tsv"
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


source("/home/kevin/dev/ltrparser/R/ltr-parser.R")
configFile <- "/home/kevin/dev/ltrparser/config/hiv.config.tsv"
testMetadataFile <- "/home/kevin/dev/ltrparser/testdata/blankson_20201211.tsv"
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

source("/home/kevin/dev/ltrparser/R/ltr-parser.R")
configFile <- "/home/kevin/dev/ltrparser/config/hiv.config.tsv"
testMetadataFile <- "/home/kevin/dev/ltrparser/testdata/scott_20201211.tsv"
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

source("/home/kevin/dev/ltrparser/R/ltr-parser.R")
configFile <- "/home/kevin/dev/ltrparser/config/hiv.config.tsv"
testMetadataFile <- "/home/kevin/dev/ltrparser/testdata/blankson_20210115.tsv"
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

source("/home/kevin/dev/ltrparser/R/ltr-parser.R")
configFile <- "/home/kevin/dev/ltrparser/config/hiv.config.tsv"
testMetadataFile <- "/home/kevin/dev/ltrparser/testdata/blankson_ES24.tsv"
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
