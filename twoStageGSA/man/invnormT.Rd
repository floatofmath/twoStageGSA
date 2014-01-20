\name{invnormT}
\alias{invnormT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  provides inverse normal combination test for screening step
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
invnormT(set, perm, d, l)
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
    p <- apply(d[set, ], 1, function(y) t.test(y ~ l)$p.value)
    q <- sum(qnorm(1 - p))/sqrt(length(p))
    return(1 - pnorm(q))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
