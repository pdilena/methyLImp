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
Numeric value in the domain [min,max].
}
  \item{min}{
Minimum range value.
}  
  \item{max}{
Maximum range value.
}  
}
%\details{}
\value{
Returns the logit transformed values.
}
\author{
Pietro Di Lena
}
% \note{}
\examples{
# Display system's double machine precision
print(.Machine$double.neg.eps)
# Without logit function
print(log(0))
# With logit function
print(log(plogit(0))
}