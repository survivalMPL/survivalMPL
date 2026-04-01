# zzz.R — package hooks
#
# .onLoad() is called after all R files in R/ are sourced.  Using it here
# guarantees that every basis-*.R file has already defined its spec list by
# the time register_basis() is called — no Collate: ordering needed.

.onLoad <- function(libname, pkgname) {
  # Phase 3: basis registrations will be added here, e.g.
  #   register_basis(.msplines_spec)
  #   register_basis(.uniform_spec)
  #   register_basis(.gaussian_spec)
  #   register_basis(.epanechnikov_spec)
}
