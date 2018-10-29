\name{plogit}
\alias{plogit}
\title{
Pseudo-logit function 
}
\description{
Pseudo-logit function bounded to [-X,X], where
X depends on the double machine precision.
}
\usage{
plogit(x, min=0, max=1)
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
Returns the logit transformed values.
}
\references{
Di Lena P et al: Missing value estimation methods for DNA methylation data. submitted to Bioinformatics
}
\author{
Pietro Di Lena
}
% \note{}
