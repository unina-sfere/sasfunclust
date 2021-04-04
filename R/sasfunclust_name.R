


#' @title Sparse and smooth functional data clustering
#' @details
#'
#'\tabular{ll}{
#'Package: \tab sasfunclust\cr
#'Type: \tab Package\cr
#'Version: \tab `r packageVersion("sasfunclust")` \cr
#'Date: \tab  `r Sys.Date()` \cr
#'License: \tab `r packageDescription("sasfunclust", fields="License")`\cr
#'}
#'
#'
#'
#'
#'
#' @author Fabio Centofanti, Antonio Lepore, Biagio Palumbo
#' @references
#' Centofanti, F., Lepore, A., & Palumbo, B. (2021).
#' Sparse and Smooth Functional Data Clustering.
#' \emph{arXiv preprint arXiv:2103.15224}.
#'
#' @seealso \code{\link{sasfclust}},  \code{\link{sasfclust_cv}}
#' @examples
#' \donttest{
#' library(sasfunclust)
#' train<-simulate_data("Scenario I",n_i=20,var_e = 1,var_b = 0.5^2)
#' lambda_s_seq=10^seq(-4,-3)
#' lambda_l_seq=10^seq(-1,0)
#' G_seq=2
#' mod_cv<-sasfclust_cv(X=train$X,grid=train$grid,G_seq=G_seq,
#' lambda_l_seq = lambda_l_seq,lambda_s_seq =lambda_s_seq,maxit = 5,K_fold = 2,q=10)
#' plot(mod_cv)
#'
#' mod<-sasfclust(X=train$X,grid=train$grid,G=mod_cv$G_opt,
#'  lambda_l = mod_cv$lambda_l_opt,lambda_s =mod_cv$lambda_s_opt,maxit = 5,q=10)
#'
#' print(mod$clus$classes)
#' plot(mod)
#'}
"_PACKAGE"
#> [1] "_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
