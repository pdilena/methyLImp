\name{pinvr}
\alias{pinvr}
\title{
Moore-Penrose generalized inverse of a matrix
}
\description{
Computes the Moore-Penrose generalized inverse of a matrix. Allows rank
reduction of the generalized inverse.
 
This function is directly taken from MASS package
and modified in order to include the rank reduction option. 
}
\usage{
pinvr(X, max.sv = min(dim(X)), tol = sqrt(.Machine$double.eps))
}
\arguments{
  \item{X}{
Real or complex valued matrix to be inverted.
}
  \item{max.sv}{
Maximum number of singlular values to be used in the single value decomposition.
}  
  \item{tol}{
A relative tolerance to detect zero singular values.
}  
}
%\details{}
\value{
Returns the Moore-Penrose pseuo-inverse of the matrix X. 
}
\author{
Pietro Di Lena
}
% \note{}
