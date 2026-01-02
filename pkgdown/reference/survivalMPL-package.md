# Penalised Maximum Likelihood for Survival Analysis Models

Simultaneously estimate regression coefficients and provide a smooth
non-parametric baseline hazard for proportional hazards Cox models using
maximum penalised likelihood (MPL). Right, left, and interval censoring
are supported.

## Details

Optimisation combines Newton updates for regression parameters with a
multiplicative scheme for baseline hazards while enforcing
non-negativity (Ma, Couturier, Heritier and Marschner, 2021). Inference
is available for regression coefficients, baseline hazard, cumulative
hazard, and survival functions, as well as their predictions.

This software is provided without warranties or guarantees of any kind.

## References

Ma, J., Couturier, D.-L., Heritier, S., and Marschner, I.C. (2021).
Penalized likelihood estimation of the proportional hazards model for
survival data with interval censoring. *International Journal of
Biostatistics*, doi:10.1515/ijb-2020-0104.

Ma, J., Heritier, S., and Lo, S. (2014). On the Maximum Penalised
Likelihood Approach for Proportional Hazard Models with Right Censored
Survival Data. *Computational Statistics and Data Analysis* **74**,
142-156.

Ma, J. (2010). Positively constrained multiplicative iterative algorithm
for maximum penalised likelihood tomographic reconstruction. *IEEE
Transactions On Signal Processing* **57**, 181-192.

## See also

Useful links:

- <https://CRAN.R-project.org/package=survivalMPL>

## Author

**Maintainer**: Dominique-Laurent Couturier
<dominique.couturier@mrc-bsu.cam.ac.uk>
([ORCID](https://orcid.org/0000-0001-5774-5036))

Authors:

- Jun Ma

- Stephane Heritier

- Maurizio Manuguerra

- Serigne Lo
