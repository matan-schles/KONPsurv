\name{KONPSURV-package}
\alias{KONPSURV-package}
\alias{KONPSURV}
\docType{package}
\title{
 KONP Tests for Testing the Equality of K Distributions for Right-Censored Data
}
\description{
 An implementation of the K-sample omnibus non-proportional hazrds (KONP) tests, which are non-parametric permutation tests, for testing the equality of K distributions for right-censored data.
}
\details{
  The package contains one function:

  
  \code{\link{konp_test}}: the permutation test for equality of K distributions with right-censored data.

}
\author{
Author and Maintainer: Matan Schlesinger \email{matan.schles@gmail.com}

Author: Malka Gorfine \email{gorfinem@post.tau.ac.il}

}
\references{
This will contain a a reference to the article when it will be published
}
\keyword{ package }

\examples{
  \donttest{
## Generate some data to preform the test
set.seed(1)
n <- 50
time <- rexp(n)
delta <- sample(c(0,1),n,TRUE)
trt <- c(rep(0,25),rep(1,25))

konp_test(time,delta,trt,n.perm=10^3)  
  }
}
