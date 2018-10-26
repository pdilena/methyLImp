## methyLImp - linear imputation model for continuous variables
## 
## This R script contains the functions to introduce
## missing values in the data and calulate the
## performance statistics.
## 
## Author: P. Di Lena, pietro.dilena@unibo.it

# ROOT MEAN SQUARE ERROR
rmse <- function(o,p) {sqrt(mean((o-p)^2))}
# MEAN ABSOLUTE ERROR
mae  <- function(o,p) {mean(abs(o-p))}
# PEARSON CORRELATION COEFFICIENT
pcc  <- function(o,p) {cor(o,p)}
# MEAN ABSOLUTE PERCENTAGE ERROR
mape <- function(o,p) {100*mean(abs((o-p)/o))}

# Calculates the performances statistics
# with respect to the above four metrics
gen_stat <- function(original, imputed) {
  # Remove too small values for MAPE metric
  ori_mape <- original[abs(original) > 10^-4]
  imp_mape <-  imputed[abs(original) > 10^-4]

  a <- round(rmse(original,imputed),digits=3)
  b <- round(mae(original,imputed),digits=3)
  c <- round(pcc(original,imputed),digits=3)
  d <- round(mape(ori_mape,imp_mape),digits=3)
  stat <- c(a,b,c,d)
  names(stat) <- c("RMSE","MAE","PCC","MAPE")
  return(stat)
}

# Introduces a random amount of missing values
# in a single column of the data. 
gen_randNA <- function(dat, col = 1, frac = 0.1) {
  n <- nrow(dat)
  NA_tot <- logical(n)
  NA_tot[sample(n,round(n*frac))] <- TRUE

  dat[NA_tot,col] <- NA
  return(dat)
}
