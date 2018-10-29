\name{gen_stat}
\alias{gen_stat}
\title{
Calculation of imputation performance statistics
}
\description{
Calculates the performances statistics of an imputed methylation array with respect to four metrics: Root Mean Square Error, Mean Absolute Error, Pearson Correlation Coefficient and Mean Absolute Percentage Error. 
}
\usage{
gen_stat(original, imputed) 
}
\arguments{
  \item{original}{
Vector of original beta or M-values.
}
  \item{imputed}{
Vector of imputed beta or M-values.
}  
}
%\details{}
\value{
Returns a detailed summary of the performance statistics of the imputed array.
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
## Impute the missing values with the methyLImp routine methyLImp
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
