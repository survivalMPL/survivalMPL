# MPL Proportional Hazards Regression Object

Objects returned by \[coxph_mpl()\] represent proportional hazards
models fitted via maximum penalised likelihood. They provide methods for
`print`, `summary`, `plot`, `residuals`, and `predict`.

## Format

A list with the following elements:

- coef:

  List with regression estimates (`Beta`) of length \\p\\ and baseline
  hazard estimates (`Theta`) of length \\m\\.

- se:

  List with standard-error matrices for `Beta` and `Theta` across
  inference methods.

- covar:

  List of covariance matrices for the available inference methods.

- ploglik:

  Length-2 vector with penalised log-likelihood components.

- iter:

  Length-3 vector with iteration counts for smoothing and parameter
  updates.

- knots:

  List with basis parameters: `m`, `Alpha`, `Delta`, and (for Gaussian
  bases) `Sigma`.

- control:

  Control settings as returned by \[coxph_mpl.control()\].

- dim:

  List with `n`, number of events, ties, number of covariates (`p`), and
  number of baseline parameters (`m`).

- call:

  The matched call.

- data:

  List with outcome times, censoring indicators, and design matrix `X`.

## Details

All components listed below must be present in a valid `"coxph_mpl"`
object.

## See also

\[coxph_mpl()\], \[summary.coxph_mpl()\], \[coef.coxph_mpl()\],
\[plot.coxph_mpl()\], \[residuals.coxph_mpl()\], \[predict.coxph_mpl()\]
