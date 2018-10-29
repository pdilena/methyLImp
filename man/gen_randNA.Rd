\name{gen_randNA}
\alias{gen_randNA}
\title{
Generation of missing values at random locations in a data array
}
\description{
Introduces a random amount of missing values in a single column of the data array. 
}
\usage{
gen_randNA(dat, col, frac) 
}
\arguments{
  \item{dat}{
Matrix of betas or M-values. Variables are on the rows and samples on the columns.
}
  \item{col}{
Column for missing data insertion.
}  
  \item{frac}{
Fraction of missing data to be inserted.
}
}
%\details{}
\value{
Returns input data array with the missing values inserted.
}
\references{

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
}
