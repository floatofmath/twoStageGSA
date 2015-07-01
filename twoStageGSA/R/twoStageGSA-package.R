

#' performs global test in screening step %% ~~function to do ... ~~
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param set %% ~~Describe \code{set} here~~
#' @param perm %% ~~Describe \code{perm} here~~
#' @param d %% ~~Describe \code{d} here~~
#' @param l %% ~~Describe \code{l} here~~
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ...
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#' 
#' ## The function is currently defined as
#' function (set, perm, d, l) 
#' {
#'     out <- gt(l, t(d[set, ]), model = "logistic", permutations = perm)
#'     p.value(out)
#'   }
#' 
NULL





#' Performs hotellingT test in screening step
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param set %% ~~Describe \code{set} here~~
#' @param perm %% ~~Describe \code{perm} here~~
#' @param d %% ~~Describe \code{d} here~~
#' @param l %% ~~Describe \code{l} here~~
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ...
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#' 
#' ## The function is currently defined as
#' function (set, perm, d, l) 
#' {
#'     hotelling <- function(d1, d2) {
#'         k <- ncol(d1)
#'         n1 <- nrow(d1)
#'         n2 <- nrow(d2)
#'         if (k > n1 | k > n2) {
#'             mysolve <- ginv
#'         }
#'         else {
#'             mysolve <- solve
#'         }
#'         xbar1 <- apply(d1, 2, mean)
#'         xbar2 <- apply(d2, 2, mean)
#'         dbar <- xbar2 - xbar1
#'         v <- ((n1 - 1) * var(d1) + (n2 - 1) * var(d2))/(n1 + 
#'             n2 - 2)
#'         t2 <- n1 * n2 * dbar %*% mysolve(v) %*% dbar/(n1 + n2)
#'         f <- (n1 + n2 - k - 1) * t2/((n1 + n2 - 2) * k)
#'         return(1 - pf(f, k, n1 + n2 - k - 1))
#'     }
#'     if (is.factor(l)) {
#'         ls <- levels(l)
#'     }
#'     else {
#'         ls <- unique(l)
#'     }
#'     if (length(ls) > 2) 
#'         stop("Too many groups!")
#'     d1 <- d[set, l == ls[1]]
#'     d2 <- d[set, l == ls[2]]
#'     hotelling(t(d1), t(d2))
#'   }
#' 
NULL





#' provides inverse normal combination test for screening step
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param set %% ~~Describe \code{set} here~~
#' @param perm %% ~~Describe \code{perm} here~~
#' @param d %% ~~Describe \code{d} here~~
#' @param l %% ~~Describe \code{l} here~~
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ...
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#' 
#' ## The function is currently defined as
#' function (set, perm, d, l) 
#' {
#'     p <- apply(d[set, ], 1, function(y) t.test(y ~ l)$p.value)
#'     q <- sum(qnorm(1 - p))/sqrt(length(p))
#'     return(1 - pnorm(q))
#'   }
#' 
NULL





#' Provides nettletons test for screening step
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param set %% ~~Describe \code{set} here~~
#' @param perm %% ~~Describe \code{perm} here~~
#' @param d %% ~~Describe \code{d} here~~
#' @param l %% ~~Describe \code{l} here~~
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ...
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#' 
#' ## The function is currently defined as
#' function (set, perm, d, l) 
#' {
#'     mrpp(l, e.com(t(d[set, ])), nperm = perm)
#'   }
#' 
NULL





#' What the package does (short line) ~~ package title ~~
#' 
#' More about what it does (maybe more than one line) ~~ A concise (1-5 lines)
#' description of the package ~~
#' 
#' \tabular{ll}{ Package: \tab twoStageGSA\cr Type: \tab Package\cr Version:
#' \tab 1.0\cr Date: \tab 2014-01-16\cr License: \tab What license is it
#' under?\cr } ~~ An overview of how to use the package, including the most
#' important ~~ ~~ functions ~~
#' 
#' @name twoStageGSA-package
#' @aliases twoStageGSA-package twoStageGSA
#' @docType package
#' @author Who wrote it
#' 
#' Maintainer: Who to complain to <yourfault@@somewhere.net> ~~ The author
#' and/or maintainer of the package ~~
#' @seealso ~~ Optional links to other man pages, e.g. ~~ ~~
#' \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
#' @references ~~ Literature or other references for background information ~~
#' @keywords package
#' @examples
#' 
#' ~~ simple examples of the most important functions ~~
#' 
NULL



