\name{methyLImp}
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
  \item{data}{
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
methyLImp(data,min=0,max=1)
}