source("/home/kevin/dev/ltrparser/R/misc.R")
library(testthat)

test_that("soft_clip", {
  expect_equal(extract_initial_soft_clip("17S"),17)
  expect_equal(extract_initial_soft_clip("17S44M", FALSE),17)
  expect_equal(extract_initial_soft_clip("17S44M22S", TRUE),22)
  
  expect_true(is.na(extract_initial_soft_clip("100M", TRUE)))
  expect_true(is.na(extract_initial_soft_clip("100M")))
  expect_true(is.na(extract_initial_soft_clip(NA)))
})

test_that("identify_candidates", {
  
})

test_that("adjust_candidates", {
  
})

test_that("filter_candidates", {
  
})

test_that("unique_ltrs", {
  
})

test_that("update_candidates", {
  
})

test_that("parse_bowtie", {
  
})

test_that("filter_seqs", {
  
})