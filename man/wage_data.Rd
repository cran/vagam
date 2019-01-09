\name{wage_data}
\alias{wage_data}
\docType{data}
\title{
Union membership data set
}
\description{
1985 North American population survey containing information on union membership and various worker's attributes. 
}
\usage{data("wage_data")}
\format{
  A data frame with 534 observations on the following 11 variables.
  \describe{
    \item{\code{education}}{a numeric vector.}
    \item{\code{south}}{a factor with levels \code{no} \code{yes}.}
    \item{\code{gender}}{a factor with levels \code{female} \code{male}.}
    \item{\code{experience}}{a numeric vector.}
    \item{\code{union}}{a factor with levels \code{member} \code{not_member}.}
    \item{\code{wage}}{a numeric vector.}
    \item{\code{age}}{a numeric vector.}
    \item{\code{race}}{a factor with levels \code{Hispanic} \code{Other} \code{White}.}
    \item{\code{occupation}}{a factor with levels \code{Clerical} \code{Management} \code{Other} \code{Professional} \code{Sales} \code{Service}.}
    \item{\code{section}}{a factor with levels \code{Construction} \code{Manufacturing} \code{Other}.}
    \item{\code{marital}}{a factor with levels \code{Married} \code{Unmarried}.}
  }
}
\details{
The data consist of \eqn{n=534} observations, with the response being a Bernoulli variable of whether they were a member of union (1 = yes; 0 = no), and six covariates: gender (1 = female, 0 = male), race (1 = white; 0 = other), an indicator variable for whether the worker lives in the south (1 = yes; 0 = no), age in years, hourly wage, and number of years in education. 

One of the aims of the survey is to uncover associations between workers' characteristics and their probability of union membership. The dataset is used in Ruppert et al., (2003) and Hui et al. (2018), among others, to illustrate the application of Semiparametric regression, as it is believed that union membership may vary non-linearly with the three continuous variables (age, wage, education). 
}
\source{
http://mldata.org/repository/data/viewslug/statlib-20050214-cps_85_wages/
}
\references{
\itemize{
\item Berndt, E. (1991). The Practice of Econometrics: Classic and Contemporary. Addison-Wesley Publishing Company, Reading, Massachusetts.

\item Hui, F. K. C., You, C., Shang, H. L., and Mueller, S. (2018). Semiparametric regression using variational approximations, \emph{Journal of the American Statistical Association}, \bold{forthcoming}.

\item Ruppert, D., Wand, M. P., and Carroll, R. (2003). Semiparametric Regression. Cambridge University Press, New York.
}}
\examples{
data(wage_data)

## Please see examples in the help file for the vagam function.
}
\keyword{datasets}
