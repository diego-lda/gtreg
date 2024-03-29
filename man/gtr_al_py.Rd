% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grt_al_py.R
\name{gtr_al_py}
\alias{gtr_al_py}
\title{GTR Adaptive Lasso Function.}
\usage{
gtr_al_py(
  y,
  x,
  X.type = NULL,
  gam = 0,
  lam.vec = rep(1, 6),
  xdf = 0,
  ydf = 4,
  yorder = 4,
  nyg.max = 101,
  Ysing = FALSE,
  y_knots = NULL,
  ng.qgm = 101,
  addxint = T,
  info = NULL,
  yorth = FALSE,
  xorth = FALSE,
  reltol = 0.001,
  feastol = 0.001,
  abstol = 0.001,
  tol.res = 0.1,
  cval = 0.1,
  lb = 1,
  threshold = 1e-05,
  e0mode = F,
  constrained = F,
  AL = F,
  doprimal = F,
  fac = 1,
  maxit = 90,
  model = NULL,
  bounded = F,
  beta2 = F,
  Cbound = Inf,
  ugrid = seq(0.1, 0.9, by = 0.1),
  gam.grid = log_space(-3, 3, 20),
  algor = "ECOS",
  coord.bare = NULL,
  coord.spline = NULL,
  coord.tensor = NULL,
  easy = T,
  method = "BIC",
  delta.ok = F,
  parallel = F,
  eta.check = T
)
}
\arguments{
\item{parallel}{}
}
\description{
This is a wrap-up function for adaptive lasso GTR/spline-spline. It allows for a Lasso first-step.
}
