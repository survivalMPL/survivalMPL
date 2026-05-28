# Register a New Basis Type

Register a New Basis Type

## Usage

``` r
register_basis(spec)
```

## Arguments

- spec:

  A list with the following named elements:

  name

  :   Canonical name string (e.g. `"msplines"`).

  aliases

  :   Character vector of accepted short names.

  knots_fn

  :   Function `(control, events)` returning a named list with at least
      `m`, `Alpha`, `Delta`.

  matrix_fn

  :   Function `(x, knots, order, which)` returning an \\n \times m\\
      matrix (density or cumulative basis).

  penalty_fn

  :   Function `(control, knots)` returning an \\m \times m\\ roughness
      penalty matrix.

  label

  :   Display string used in plots and summaries.

  default_n_knots

  :   Length-2 numeric vector (default for n.knots).

  penalty_order_fn

  :   Function `(p, order)` returning the penalty order as an integer.

## Value

The registered entry (invisibly).
