\name{KONPsurv-package}
\alias{KONPsurv-package}
\alias{KONPsurv}
\docType{package}
\title{
 KONP Tests for Testing the Equality of K Distributions for Right-Censored Data
}
\description{
An implementation of the K-sample omnibus non-proportional hazrds (KONP) tests.

The KONP tests are powerful non-parametric tests for comparing K (>=2) hazard functions based on right-censored data. These tests are consistent against any differences between the hazard functions of the groups.
The KONP tests are often more powerful than other existing tests, especially under non-proportional hazard functions.
}
\details{
  The package contains one function:

  
  \code{\link{konp_test}}: non-parametric tests for equality of K distributions using right-censored data.

}
\author{
Author and Maintainer: Matan Schlesinger \email{matan.schles@gmail.com}

Author: Malka Gorfine \email{gorfinem@tauex.tau.ac.il}

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
status <- sample(c(0,1),n,TRUE)
group <- c(rep(0,25),rep(1,25))

konp_test(time,status,group,n_perm=10^3)  
  }
}
