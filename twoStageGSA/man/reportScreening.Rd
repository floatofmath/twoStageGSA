\name{reportScreening}
\alias{reportScreening}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Produces human readable summary report
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
reportScreening(scr,pw.annotation)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scr}{
%%     ~~Describe \code{scr} here~~
}
\item{pw.annotation}{
  Annotation environment that maps pathway identifiers to some
  annotation
}
\item{pw.kv}{
  Character value of the form c('KEY','VALUE') specifying which key- and
  valuetypes should be used in the corresponding annotation database
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
function (scr) 
{
    if (length(scr$sigSets) < 1) 
        return(NULL)
    if (length(scr$sigSets) == 1) {
        tout <- c(scr$sigSets, round(scr$qVal[scr$sigSets], 3), 
            unlist(mget(sub("hsa", "", scr$sigSets), KEGGPATHID2NAME, 
                ifnotfound = NA)))
        names(tout) <- c("Gene Set ID", "Gene Set Description", 
            "FDR adj. q-Value")
        return(tout)
    }
    tout <- cbind(scr$sigSets, unlist(mget(sub("hsa", "", scr$sigSets), 
        KEGGPATHID2NAME, ifnotfound = NA)), round(scr$qVal[scr$sigSets], 
        3))
    tout <- tout[order(tout[, 2]), ]
    colnames(tout) <- c("Gene Set ID", "Gene Set Description", 
        "FDR adj. q-Value")
    o <- order(tout[, 3])
    rownames(tout) <- 1:nrow(tout)
    tout[o, ]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
