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
p <-0
# Logistic function at 0
log(p/(1-p))

# Logistic function at approx 10^-16
p <- .Machine$double.neg.eps
log(p/(1-p))

# Pseudo logistic function at 0 is equivalent
# to logistic function at 10^-16
plogit(0)
}