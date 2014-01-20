\name{hotellingT}
\alias{hotellingT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Performs hotellingT test in screening step
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hotellingT(set, perm, d, l)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{set}{
%%     ~~Describe \code{set} here~~
}
  \item{perm}{
%%     ~~Describe \code{perm} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{l}{
%%     ~~Describe \code{l} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (set, perm, d, l) 
{
    hotelling <- function(d1, d2) {
        k <- ncol(d1)
        n1 <- nrow(d1)
        n2 <- nrow(d2)
        if (k > n1 | k > n2) {
            mysolve <- ginv
        }
        else {
            mysolve <- solve
        }
        xbar1 <- apply(d1, 2, mean)
        xbar2 <- apply(d2, 2, mean)
        dbar <- xbar2 - xbar1
        v <- ((n1 - 1) * var(d1) + (n2 - 1) * var(d2))/(n1 + 
            n2 - 2)
        t2 <- n1 * n2 * dbar \%*\% mysolve(v) \%*\% dbar/(n1 + n2)
        f <- (n1 + n2 - k - 1) * t2/((n1 + n2 - 2) * k)
        return(1 - pf(f, k, n1 + n2 - k - 1))
    }
    if (is.factor(l)) {
        ls <- levels(l)
    }
    else {
        ls <- unique(l)
    }
    if (length(ls) > 2) 
        stop("Too many groups!")
    d1 <- d[set, l == ls[1]]
    d2 <- d[set, l == ls[2]]
    hotelling(t(d1), t(d2))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
