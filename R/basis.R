#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# basis.R — Registry infrastructure for basis types
#
# Defines:
#   .basis_registry   private environment holding all registered basis entries
#   register_basis()  add a new basis type to the registry
#   get_basis()       look up by canonical name or alias
#   list_bases()      enumerate all registered bases
#
# Thin dispatchers (drop-in replacements for old if/else helpers):
#   resolve_basis_name()   replaces basis.name_mpl()
#   compute_knots()        replaces knots_mpl()
#   compute_basis_matrix() replaces basis_mpl()
#   compute_penalty()      replaces penalty_mpl()
#   compute_penalty_order() replaces penalty.order_mpl()
#   basis_label()          replaces if/else label chains in plot/summary
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Private registry environment — not exported
.basis_registry <- new.env(parent = emptyenv())


# ── Shared helpers ────────────────────────────────────────────────────────────

#' Compute the Alpha knot sequence from event times.
#' Called by every basis-specific knots_fn.
#' @keywords internal
.compute_alpha_knots <- function(control, events) {
  if (control$n.knots[2] == 0L) {
    quantile(events, seq(0, 1, length.out = control$n.knots[1] + 2L))
  } else {
    Alpha1 <- quantile(events,
                       seq(0, control$range.quant[2],
                           length.out = control$n.knots[1] + 1L))
    Alpha2 <- seq(quantile(events, control$range.quant[2]),
                  range(events)[2L],
                  length = control$n.knots[2] + 2L)
    c(Alpha1, Alpha2[-1L])
  }
}


#' Register a New Basis Type
#'
#' @param spec A list with the following named elements:
#'   \describe{
#'     \item{name}{Canonical name string (e.g. \code{"msplines"}).}
#'     \item{aliases}{Character vector of accepted short names.}
#'     \item{knots_fn}{Function \code{(control, events)} returning a named list
#'       with at least \code{m}, \code{Alpha}, \code{Delta}.}
#'     \item{matrix_fn}{Function \code{(x, knots, order, which)} returning an
#'       \eqn{n \times m} matrix (density or cumulative basis).}
#'     \item{penalty_fn}{Function \code{(control, knots)} returning an
#'       \eqn{m \times m} roughness penalty matrix.}
#'     \item{label}{Display string used in plots and summaries.}
#'     \item{default_n_knots}{Length-2 numeric vector (default for n.knots).}
#'     \item{penalty_order_fn}{Function \code{(p, order)} returning the penalty
#'       order as an integer.}
#'   }
#' @return The registered entry (invisibly).
#' @keywords internal
register_basis <- function(spec) {
  stopifnot(
    is.list(spec),
    is.function(spec$knots_fn),
    is.function(spec$matrix_fn),
    is.function(spec$penalty_fn),
    is.character(spec$name),
    nchar(spec$name) > 0
  )

  # Collision check
  for (key in c(spec$name, spec$aliases)) {
    existing <- .basis_registry[[key]]
    if (!is.null(existing) && existing$name != spec$name) {
      stop(sprintf(
        "Alias '%s' is already registered by basis '%s'.",
        key, existing$name
      ), call. = FALSE)
    }
  }

  .basis_registry[[spec$name]] <- spec
  for (a in spec$aliases) {
    .basis_registry[[a]] <- spec
  }

  invisible(spec)
}


#' Look Up a Registered Basis by Name or Alias
#'
#' @param name Character string — canonical name or any registered alias.
#' @return The registered basis entry (a list).
#' @keywords internal
get_basis <- function(name) {
  key   <- tolower(trimws(name))
  entry <- .basis_registry[[key]]
  if (is.null(entry)) {
    avail <- list_bases()
    stop(sprintf(
      "Unknown basis '%s'. Available bases: %s",
      name, paste(avail$name, collapse = ", ")
    ), call. = FALSE)
  }
  entry
}


#' List All Registered Bases
#'
#' @return A \code{data.frame} with columns \code{name} and \code{label}.
#' @export
list_bases <- function() {
  all_keys <- ls(.basis_registry)
  if (length(all_keys) == 0L) {
    return(data.frame(name = character(0), label = character(0),
                      stringsAsFactors = FALSE))
  }
  entries  <- lapply(all_keys, function(k) .basis_registry[[k]])
  names_v  <- vapply(entries, `[[`, character(1), "name")
  idx      <- !duplicated(names_v)
  data.frame(
    name  = names_v[idx],
    label = vapply(entries[idx], `[[`, character(1), "label"),
    stringsAsFactors = FALSE
  )
}


# ── Thin dispatchers ──────────────────────────────────────────────────────────

#' Resolve a basis alias to its canonical name.
#' Replaces \code{basis.name_mpl()}.
#' @keywords internal
resolve_basis_name <- function(name) {
  get_basis(name)$name
}

#' Compute the knot sequence for a given basis.
#' Replaces \code{knots_mpl()}.
#' @keywords internal
compute_knots <- function(control, events) {
  get_basis(control$basis)$knots_fn(control, events)
}

#' Compute the density or cumulative basis matrix.
#' Replaces \code{basis_mpl()}.
#' @keywords internal
compute_basis_matrix <- function(x, knots, basis_name, order, which = 1) {
  get_basis(basis_name)$matrix_fn(x, knots, order, which)
}

#' Compute the roughness penalty matrix.
#' Replaces \code{penalty_mpl()}.
#' @keywords internal
compute_penalty <- function(control, knots) {
  get_basis(control$basis)$penalty_fn(control, knots)
}

#' Compute the penalty order for a basis.
#' Replaces \code{penalty.order_mpl()}.
#' @keywords internal
compute_penalty_order <- function(basis_name, p, order) {
  get_basis(basis_name)$penalty_order_fn(p, order)
}

#' Return the display label for a basis.
#' Used by plot and summary methods.
#' @keywords internal
basis_label <- function(basis_name) {
  get_basis(basis_name)$label
}
