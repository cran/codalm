## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(codalm)
library(ggtern)

## ----load_data----------------------------------------------------------------
data("WhiteCells", package = 'ggtern')
image <- subset(WhiteCells, Experiment == "ImageAnalysis")
image_mat <- as.matrix(image[,c("G", "L", "M")])
microscopic <- subset(WhiteCells, Experiment == "MicroscopicInspection")
microscopic_mat <- as.matrix(microscopic[,c("G", "L", "M")])

image_mat  <- image_mat  / rowSums(image_mat)
microscopic_mat <- microscopic_mat / rowSums(microscopic_mat)

## ----estimate_B---------------------------------------------------------------
B_est <- codalm(y = microscopic_mat, x = image_mat)
B_est

## ----bootstrap----------------------------------------------------------------
B_ci <- codalm_ci(y = microscopic_mat, x = image_mat, nboot = 50,
                   conf = .95)
B_ci$ci_L
B_ci$ci_U

## ----bootstrap_parallel-------------------------------------------------------
ncores <- 2
Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE = "quiet")
B_ci_parallel <- codalm_ci(y = microscopic_mat, x = image_mat, nboot = 50,
                   conf = .95, parallel = TRUE, ncpus = ncores, strategy = 'multisession')
identical(B_ci$ci_L, B_ci_parallel$ci_L)
identical(B_ci$ci_U, B_ci_parallel$ci_U)

## ----independence_test--------------------------------------------------------
set.seed(123)
x <- gtools::rdirichlet(100, rep(1, 3))
y <- gtools::rdirichlet(100, rep(1, 3))
indep_test_pval <- codalm_indep_test(y = y, x = x,
                                     nperms = 100, init.seed = 123)
indep_test_pval

## ----independence_test_parallel-----------------------------------------------
indep_test_pval_parallel <- codalm_indep_test(y = y, x = x,
                                              nperms = 100, init.seed = 123,
                                              parallel = TRUE, ncpus = ncores,
                                              strategy = 'multisession')
indep_test_pval_parallel

