\name{focus}
\alias{focus}
\title{Perform focus step of two-stage GSA
}
\description{
  Performs focus step of two-stage GSA. In this step single genes are
  tested at significance levels resulting from the screening
  step, which only tests gene sets. 
}
\usage{
focus(data, labels, genesets, sigSets, B = 100, test = "wilcoxon", side = "abs", adjust = "WY", reuse.p = F)
}
\arguments{
  \item{data}{
    Matrix containing expression values. Genes are in rows, samples in columns.
}
  \item{labels}{
    Vector containing gene identifiers that correspond to those in genesets
}
  \item{genesets}{
    List of genesets, given as vectors of gene identifiers
}
  \item{sigSets}{
    Output from screening step
}
  \item{B}{
    Number of resampling iterations
}
  \item{test}{
    name of test function, if reuse.p is FALSE a test function that
    given if reuse.p is TRUE a function that given takes a formula
    argument (x~labels) and returns an object with $p.value slot. Vector
    of gene identifiers returns the corresponding unadjusted 
    p-values 
}
  \item{side}{
    should one or two-sided tests be performed
}
  \item{adjust}{
    multiplicity adjustment to use
}
  \item{reuse.p}{
    Sometimes it is desired to re-use the p-values of a previous
    differential gene expression analysis. If TRUE test should be a
    vector of unadjusted p-values
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
function (data, labels, genesets, sigSets, B = 100, test = "wilcoxon", 
    side = "abs", adjust = "WY", reuse.p = F) 
{
    focusWY <- function(idx) {
        idx <- which(rownames(data) \%in\% genesets[[idx]])
        out <- mt.minP(data[idx, ], labels, test = test, B = B, 
            side = side)
        ans <- out$adjp
        names(out) <- out$index
    }
    focusBonfH <- function(idx) {
        testfun <- match.fun(test)
        idx <- rownames(data)[which(rownames(data) \%in\% genesets[[idx]])]
        if (reuse.p) {
            return(testfun(idx))
        }
        out <- apply(data[idx, ], 1, testfun)
        ans <- p.adjust(out, method = "holm")
        return(ans)
    }
    focusTest <- function(idx) {
        switch(adjust, WY = focusWY(idx), holm = focusBonfH(idx))
    }
    pWY <- lapply(sigSets, focusTest)
    return(pWY)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
