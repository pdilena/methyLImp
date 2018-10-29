\name{inv.plogit}
\alias{inv.plogit}
\title{
Moore-Penrose generalized inverse of a matrix
}
\description{
Computes the Moore-Penrose generalized inverse of a matrix. Allows rank
reduction of the generalized inverse.
 
This function is directly taken from MASS package (code on GPLv3 license)
and modified in order to include the rank reduction option. The added code
for rank reduction is commented in the implementation.
}
\usage{
pinvr(X, max.sv = min(dim(X)), tol = sqrt(.Machine$double.eps))
}
\arguments{
  \item{X}{
Real or complex valued matrix to be transformed  
}
  \item{max.sv}{
Maximum for single value decomposition
}  
  \item{tol}{
Tolerance
}  
}
%\details{}
\value{
Returns the inverse logit transformed values.
}
\references{
Di Lena P et al: Missing value estimation methods for DNA methylation data. submitted to Bioinformatics
}
\author{
Pietro Di Lena
}
% \note{}
