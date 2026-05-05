# Development Workflow

## Branch conventions

- `main` — stable, CRAN-ready
- `test-dev` — current development branch

## Commit discipline

- Run `devtools::test()` after every change — snapshots must match before committing
- No logic changes in the same commit as a file split or rename
- No new abstractions beyond what the task requires
- No comments unless the WHY is non-obvious (naming should be self-explanatory)
- No backwards-compatibility shims for removed code — delete cleanly

## Typical session

```r
devtools::load_all()    # after any R/ change
devtools::test()        # confirm green
devtools::check()       # before pushing or tagging
```

## File locations

| What | Where |
|------|-------|
| Main fitting engine | `R/coxph.r` |
| Basis registry | `R/basis.R` |
| Package dataset | `data/bcos2.rda` + `R/data-bcos2.R` |
| Test snapshots | `tests/testthat/_snaps/` |
| Book PDF | `old_sources/Penalized_likelihood_based_survival_data_analysis_using_R_05_04_24.pdf` |
| Scratch / investigation scripts | `non-infinite_issue.R`, `book_sim.R` (root, untracked) |
