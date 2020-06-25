# ltrparser
Functions to identify and parse LTR sequences from HIV integration site and gene therapy samples.

# Setup
Check the `config/hiv.config.tsv` file to make sure Bowtie indexes are in the correct location.
You will need to set "Host" index for all samples, but only set "Viral" index if you want to map against an HIV database (e.g. for HIV Integration site data).

# Running Test
From the project directory, run:

`cd R/`

`Rscript test.R`

This will run the script on samples described in the file `testdata/testdata.tsv` using the parameters in the `config/hiv.config.tsv` config file.

CSV and PDF files created will be placed in `testdata/demo/csv` and `testdata/demo/pdf`, respectively (they will also already be there upon download but the script will recreate them).
