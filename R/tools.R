#' Tridiagonal solver
#' 
#' Converted into a R code from the original code of Gordon Bonan: Bonan, G.
#' (2019). Climate Change and Terrestrial Ecosystem Modeling. Cambridge:
#' Cambridge University Press. doi:10.1017/9781107339217
#' 
#' @description 
#' % Solve for U given the set of equations R * U = D, where U is a vector
#' of length N, D is a vector of length N, and R is an N x N tridiagonal
#' matrix defined by the vectors A, B, C each of length N. A(1) and
#' C(N) are undefined and are not referenced.
#'
#'     |B(1) C(1) ...  ...  ...                     |
#'     |A(2) B(2) C(2) ...  ...                     |
#' R = |     A(3) B(3) C(3) ...                     |
#'     |                    ... A(N-1) B(N-1) C(N-1)|
#'     |                    ... ...    A(N)   B(N)  |
#'
#' The system of equations is written as:
#'
#'    A_i * U_i-1 + B_i * U_i + C_i * U_i+1 = D_i
#'
#' for i = 1 to N. The solution is found by rewriting the
#' equations so that:
#'
#'    U_i = F_i - E_i * U_i+1
#' 
#' @param a See description.
#' @param b See description.
#' @param c See description.
#' @param d See description.
#' @param n See description.
#' 
#' @return Solution U
#' @export
f.tridiagonal.solver=function(a, b, c, d, n){
  e=rep(NA,n-1)
  e[1] = c[1] / b[1]
  
  for (i in 2:(n-1)){
    e[i] = c[i] / (b[i] - a[i] * e[i-1])
  } 

  f=rep(NA,n)
  f[1] = d[1] / b[1]
  
  for (i in 2:(n)){
    f[i] = (d[i] - a[i] * f[i-1]) / (b[i] - a[i] * e[i-1])
  }
  
  u=rep(NA,n)
  u[n] = f[n];
  
  for (i in seq(n-1,1,-1)){
    u[i] = f[i] - e[i] * u[i+1];
  }
  return(u)
}


listk <- function(...) {
  # get variable names from input expressions
  cols <- as.list(substitute(list(...)))[-1]
  vars <- names(cols)
  Id_noname <- if (is.null(vars)) seq_along(cols) else which(vars == "")

  if (length(Id_noname) > 0) {
    vars[Id_noname] <- sapply(cols[Id_noname], deparse)
  }
  # ifelse(is.null(vars), Id_noname <- seq_along(cols), Id_noname <- which(vars == ""))
  x <- setNames(list(...), vars)
  return(x)
}

clamp <- function(x, lims = c(0, 1), fill.na = FALSE) {
  if (fill.na) {
    x[x < lims[1]] <- NA_real_
    x[x > lims[2]] <- NA_real_
  } else {
    x[x < lims[1]] <- lims[1]
    x[x > lims[2]] <- lims[2]
  }
  x
}
