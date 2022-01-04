\name{KONPsurv-package}
\alias{KONPsurv-package}
\alias{KONPsurv}
\docType{package}
\title{
 KONP Tests for Testing the Equality of K Distributions for Right-Censored Data
}
\description{
An implementation of the \eqn{K}-sample omnibus non-proportional hazrds (KONP) tests.

The KONP tests are powerful non-parametric tests for comparing \eqn{K} (>=2) hazard functions based on right-censored data. These tests are consistent against any differences between the hazard functions of the groups.
The KONP tests are often more powerful than other existing tests, especially under non-proportional hazard functions.
}
\details{
  The package contains one function:

  
  \code{\link{konp_test}}: non-parametric tests for equality of \eqn{K} distributions using right-censored data.

}
\author{
Author and Maintainer: Matan Schlesinger \email{matan.schles@gmail.com}

Author: Malka Gorfine \email{gorfinem@tauex.tau.ac.il}

}
\references{
Gorfine, M., Schlesinger, M., & Hsu, L. (2020). \eqn{K}-sample omnibus non-proportional hazards tests based on right-censored data. Statistical Methods in Medical Research, 29(10), 2830–2850.
\doi{10.1177/0962280220907355}


}
\keyword{ package }

\examples{
  \donttest{
# gastric cancer data  
data(gastric)

konp_test(gastric$time, gastric$status, gastric$group, n_perm=10^3)
  }
}
