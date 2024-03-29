% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gtr_solve.R
\name{gtr_solve}
\alias{gtr_solve}
\title{GTR Solver Function}
\usage{
gtr_solve(
  TYX,
  tYX,
  Kscore = 0,
  gam = 0,
  egam = 0,
  lam = 0,
  lam.vec = NULL,
  maxit = 200,
  algor = "ECOS",
  reltol = 1e-04,
  feastol = 1e-04,
  abstol = 1e-04,
  quiet = F,
  zeros = NULL,
  doprimal = FALSE,
  btarg = 0,
  silent = F,
  nXs = NULL,
  nYS = NULL,
  weights = 1,
  cval = 0.1,
  pen = NULL,
  beta2 = F,
  Xs = NULL,
  sYgrid = NULL,
  bounded = F,
  Cbound = 1e+06,
  threshold = 1e-05
)
}
\arguments{
\item{threshold}{}
}
\description{
This function solves the optimisation problem. Used to be Xpert8d_ss.
}
