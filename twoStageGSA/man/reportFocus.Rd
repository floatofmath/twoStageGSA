\name{reportFocus}
\alias{reportFocus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Produces a nice report
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
reportFocus(scr, pval, fc, q, linebreak = 75,pw.annotation,gene.annotation,pw.kv,gene.kv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scr}{
%%     ~~Describe \code{scr} here~~
}
  \item{pval}{
%%     ~~Describe \code{pval} here~~
}
  \item{fc}{
%%     ~~Describe \code{fc} here~~
}
  \item{q}{
%%     ~~Describe \code{q} here~~
}
\item{linebreak}{
%%     ~~Describe \code{linebreak} here~~
}
\item{pw.annotation}{
  Annotation environment that maps pathway identifiers to some
  annotation
}
\item{gene.annotation}{
Annotation environment that maps gene identifiers to some
  annotation
}
\item{pw.kv}{
  Character value of the form c('KEY','VALUE') specifying which key- and
  valuetypes should be used in the corresponding annotation database
  }
\item{gene.kv}{
  Character value of the form c('KEY','VALUE') specifying which key- and
  valuetypes should be used in the corresponding annotation database
w}
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
function (scr, pval, fc, q, linebreak = 75) 
{
    parser <- function(s) {
        if (!all(s >= q)) {
            idx <- which(s < q)
            ans <- paste(mget(paste(names(s)[idx], "at", sep = "_"), 
                hugene10stv1hsentrezgSYMBOL), " (", round(s[idx], 
                3), ", ", c("down", "NA", "up")[2 + sign(fc[paste(names(s)[idx], 
                "at", sep = "_"), ])], ")", sep = "", collapse = "; ")
            if (nchar(ans) > linebreak) {
                splt <- c(seq(1, nchar(ans), linebreak), nchar(ans))
                ans <- paste(strwrap(ans, linebreak), collapse = "\n")
            }
        }
        else {
            return(NULL)
        }
    }
    if (length(scr$sigSets) < 1) 
        return(NULL)
    if (length(scr$sigSets) == 1) {
        return(parser(unlist(pval)))
    }
    out <- sapply(pval, parser)
    names(out) <- scr$sigSets
    out <- out[!sapply(out, is.null)]
    desc <- unlist(mget(sub("hsa", "", names(out)), KEGGPATHID2NAME))
    outtab <- cbind(names(out), desc, unlist(out))
    rownames(outtab) <- 1:length(out)
    colnames(outtab) <- c("Gene Set ID", "Gene Set Description", 
        paste("Considerable genes with local p-value <", q))
    return(outtab)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
