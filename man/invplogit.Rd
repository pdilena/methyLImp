\name{inv.plogit}
\alias{inv.plogit}
\title{
Inverse pseudo-logit function 
}
\description{
Inverse of the pseudo-logit function bounded to [-X,X], where
X depends on the double machine precision.
}
\usage{
inv.plogit(x, min=0, max=1)
}
\arguments{
  \item{x}{
Data to be transformed  
}
  \item{min}{
minimum value for bound
}  
  \item{max}{
maximum value for bond
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
