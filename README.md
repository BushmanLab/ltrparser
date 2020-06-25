# ltrparser
Functions to identify and parse LTR sequences from HIV integration site and gene therapy samples.

# Setup
Check the `config/hiv.config.tsv` file to make sure Bowtie indexes are in the correct location.
You will need to set "Host" index for all samples, but only set "Viral" index if you want to map against an HIV database (e.g. for HIV Integration site data).

# Running
From the project directory, run:
`Rscript R/test.R`
