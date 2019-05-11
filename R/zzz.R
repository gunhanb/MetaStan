.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
}

if(getRversion() >= "3.1.0") utils::suppressForeignCheck("dose")
if(getRversion() >= "2.15.1")  utils::globalVariables("dose")
