#' ---
#' title: "Singularity test of life-table spline implementation"
#' author: "Jonas Schöley"
#' date: "August 17, 2017"
#' ---

#' We fit a spline to the log-transformed survival. In order to avoid taking the
#' log of 0 and to avoid fitting a Hyman filtered spline to log(y) = 0 (the
#' spline performs bad if y = 0), we use the value transformation $g(y) =
#' \log(y+2)$, with the inverse $g^{-1}(y) = exp(y)-2$. We transform the domain
#' as well using $q(x) = \log(x+1)$ with inverse $q^{-1}(x)=1/(x-1)$ in order to
#' capture the rapid decline in hazard of death right right after birth observed
#' in infant mortality.

FitLTSpline <- function (df, x, lx) {
  require(rlang)
  x = enquo(x); lx = enquo(lx)

  logloglxFun <- splinefun(x = log(df[,quo_name(x)]+1),
                           y = log(df[,quo_name(lx)]+2),
                           method = 'hyman')

  fnct <- function (x, type = 'lx') {
    lxFun <- function (x) { exp(logloglxFun(log(x+1)))-2 }
    dxFun <- function (x) { -lxFun(x) * (1/(x+1)) * logloglxFun(log(x+1), deriv = 1) }
    hxFun <- function (x) {
      hx = dxFun(x) / lxFun(x)
      # set hazard to 0 if lx = 0
      ifelse(is.nan(hx), 0, hx)
    }
    ddxhxFun <- function (x) { numDeriv::grad(hxFun, x) }
    switch(type, lx = lxFun(x), dx = dxFun(x), hx = hxFun(x), ddxhx = ddxhxFun(x))
  }

  return(fnct)
}

domain = seq(0, 10, 0.1)

##'# Test 1: No deaths

df = data.frame(x = 0:10, lx = 1)

plot(df)
lines(x = domain, FitLTSpline(df, x, lx)(domain, type = 'lx'))
plot(x = domain, FitLTSpline(df, x, lx)(domain, type = 'dx'), type = 'l')
plot(x = domain, FitLTSpline(df, x, lx)(domain, type = 'hx'), type = 'l')
plot(x = domain, FitLTSpline(df, x, lx)(domain, type = 'ddxhx'), type = 'l')

##'# Test 2: All deaths at second observed time point

df = data.frame(x = 0:10, lx = c(1, rep(0, 10)))

plot(df)
lines(x = domain, FitLTSpline(df, x, lx)(domain, type = 'lx'))
plot(x = domain, FitLTSpline(df, x, lx)(domain, type = 'dx'), type = 'l')
plot(x = domain, FitLTSpline(df, x, lx)(domain, type = 'hx'), type = 'l')
plot(x = domain, FitLTSpline(df, x, lx)(domain, type = 'ddxhx'), type = 'l')

##'# Test 3: All deaths at last observed time point

df = data.frame(x = 0:10, lx = c(rep(1, 10), 0))

plot(df)
lines(x = domain, FitLTSpline(df, x, lx)(domain, type = 'lx'))
plot(x = domain, FitLTSpline(df, x, lx)(domain, type = 'dx'), type = 'l')
plot(x = domain, FitLTSpline(df, x, lx)(domain, type = 'hx'), type = 'l')
plot(x = domain, FitLTSpline(df, x, lx)(domain, type = 'ddxhx'), type = 'l')


##'# Test 4: All deaths at single point mid-domain

df = data.frame(x = 0:10, lx = c(1,1,1,1,1,0,0,0,0,0,0))

plot(df)
lines(x = domain, FitLTSpline(df, x, lx)(domain, type = 'lx'))
plot(x = domain, FitLTSpline(df, x, lx)(domain, type = 'dx'), type = 'l')
plot(x = domain, FitLTSpline(df, x, lx)(domain, type = 'hx'), type = 'l')
plot(x = domain, FitLTSpline(df, x, lx)(domain, type = 'ddxhx'), type = 'l')