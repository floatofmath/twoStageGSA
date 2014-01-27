\name{screening}
\alias{screening}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Performs the screening step
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
screening(data, labels, genes = rownames(data), genesets, B = 0, q = 0.05, settest = "globalTest", min = 6, max = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{labels}{
%%     ~~Describe \code{labels} here~~
}
\item{genes}{
  vector of gene identifiers specifying the gene each row of data
  corresponds to and that match to the identifyers used to define gene
  sets. If not given rownames are assumed to be gene identifyers.
  }
  \item{genesets}{
%%     ~~Describe \code{genesets} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{q}{
%%     ~~Describe \code{q} here~~
}
  \item{settest}{
%%     ~~Describe \code{settest} here~~
}
  \item{min}{
%%     ~~Describe \code{min} here~~
}
  \item{max}{
%%     ~~Describe \code{max} here~~
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
function (data, labels, genesets, B = 0, q = 0.05, settest = "globalTest", 
    min = 6, max = 100) 
{
    toSmall <- which(sapply(genesets, function(set) sum(rownames(data) \%in\% 
        set) <= min))
    toLarge <- which(sapply(genesets, function(set) sum(rownames(data) \%in\% 
        set) >= max))
    genesets <- genesets[-c(toSmall, toLarge)]
    inSet <- unlist(sapply(genesets, function(set) which(rownames(data) \%in\% 
        set)))
    impW <- tabulate(inSet)/length(inSet)
    inSet <- sort(unique(inSet))
    data <- data[inSet, ]
    S <- length(genesets)
    n <- max(sapply(genesets, length))
    G <- nrow(data)
    m <- ncol(data)
    pmrpp <- sapply(genesets, match.fun(settest), perm = B, d = data, 
        l = labels)
    names(pmrpp) <- names(genesets)
    q.vals <- p.adjust(pmrpp, method = "fdr")
    sigSets <- names(genesets)[which(q.vals <= q)]
    return(list(sigSets = sigSets, qVal = q.vals))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
