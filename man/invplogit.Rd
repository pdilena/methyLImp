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
Numeric value in the domain [-inf,+inf]  
}
  \item{min}{
Minimum range value 
}  
  \item{max}{
Maximum range value
}  
}
%\details{}
\value{
Returns the inverse logit transformed values.
}
\author{
Pietro Di Lena
}
% \note{}
\examples{
# Minimum value of the pseudo logit function
x <- plogit(0)

# Both return 0
inv.plogit(x-1)
inv.plogit(x)
}