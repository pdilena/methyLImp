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
Matrix of betas or M-values
}
  \item{col}{
Column for missing data insertion
}  
  \item{frac}{
Fraction of missing data to be inserted
}
}
%\details{}
\value{
Returns input data array with the missing values inserted.
}
\references{
Di Lena P et al: Missing value estimation methods for DNA methylation data. submitted to Bioinformatics
}
\author{
Pietro Di Lena
}
% \note{}

\examples{
}