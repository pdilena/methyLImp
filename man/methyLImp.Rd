\name{methyLImp}
\alias{methyLImp}
\title{
Linear regression model for beta-values imputation
}
\description{
Imputes missing DNA methylation data
}
\usage{
methyLImp(dat, min = NULL, max = NULL, max.sv = NULL, col.list = NULL) 
}
\arguments{
  \item{dat}{
Matrix of betas or M-values
}
  \item{min}{
Minimum value of the distribution
}  
  \item{max}{
Maximum value of the distribution
}
  \item{max.sv}{
Max number of singular values to be used in the psuedoinvers matrix.
The parameter is used to reduce the matrix rank.
}
  \item{col.list}{
Convert col.list, if any, from names to numbers
}
}
%\details{}
\value{
Returns imputed values for the missing values.
}
\references{
Di Lena P et al: Missing value estimation methods for DNA methylation data. submitted to Bioinformatics
}
\author{
Pietro Di Lena
}
% \note{}

\examples{
## Load the methyLImp dataset, containing no missing value
data(gse64495)
summary(gse64495)

## Artificially introduce 10% missing values in the first sample
## with the gen_randNA function
set.seed(50)
samp <- 1
frac <- 0.1
gse64495.mis <- gen_randNA(gse64495,samp,frac)
summary(gse64495.mis)

## Impute the missing values with the methyLImp routine.
## Note that variables need to be on the columns and
## samples on the rows.
gse64495.imp <- methyLImp(t(gse64495.mis),min=0,max=0)
gse64495.imp <- t(gse64495.imp)

## Compare imputed and original values
miss <- is.na(gse64495.mis[,samp])
orig <- gse64495[miss,samp]
pred <- gse64495.imp[miss,samp]
gen_stat(orig,pred)
}