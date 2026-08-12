# Changelog

## survivalMPL 0.2-4.9003 (development)

- **`hiroshima` dataset rebuilt, breaking change.** It is now one row
  per subject (86,611 subjects), reconstructed from the RERF LSS14
  grouped person-time table, rather than one pseudo-subject per
  grouped-table cell (53,782 rows). The pseudo-subject version was not
  usable for analysis: it produced a *negative* colon-dose coefficient
  for all-cause mortality, an artefact of high-dose cells having more
  person-years and hence a longer mean attained age. Changes for
  existing code:
  - `entry` is now the attained age at the 1950-10-01 start of
    follow-up, the age at which a subject actually enters the risk set,
    rather than a cell mean age at exposure. Age at exposure is
    available as the new `agex` column (exactly `entry - 5.15`, so the
    two are collinear).
  - `dose` is now in **Gy**, not mGy.
  - `status` is a per-subject death indicator (50,620 deaths), not a
    “any death in this cell” indicator.
  - reproduces the source table’s subjects, deaths and colon-cancer
    deaths exactly, and its person-years to within 0.33%.

## survivalMPL 0.2-4 (2025-03-23)

CRAN release: 2025-03-24

- changed maintainer email address
- added contributor

## survivalMPL 0.2-3 (2022-11-21)

CRAN release: 2022-11-21

- sorted minor typos in documentation

## survivalMPL 0.2-2 (2022-11-21)

CRAN release: 2022-11-20

- extended the example section of `coxph_mpl`
- amended code to allow models without covariates
- removed a duplicated item in `coxph_mpl.object`

## survivalMPL 0.2-1 (2021-11-14)

CRAN release: 2021-11-15

- change of maintainer (D-L Couturier)
- minor fixes (warning messages from `try` and `model.matrix`)

## survivalMPL 0.2 (2017-12-11)

CRAN release: 2017-12-11

- interval censoring implementation

## survivalMPL 0.1.2 (2017-10-13)

CRAN release: 2017-10-13

- change of maintainer (Maurizio Manuguerra)
