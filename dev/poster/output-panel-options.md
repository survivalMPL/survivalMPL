# Replacing the "AND WHAT IT PRINTS" panel

Nothing here has been applied to the poster.  All numbers are real, taken from
the cached bcos2 msplines fit (`dev/poster/cache/fits.rds`).

## What is wrong with the current panel

Fourteen lines, of which:

| lines | content | does it sell the package? |
|---|---|---|
| 1 | title banner | no |
| 3 | log-likelihood, lambda, convergence | actively unhelpful, see below |
| 4 | dataset name, n, events, censored | no, this is the data description |
| 3 | the coefficient table | yes, but any Cox implementation prints this |
| 1 | elision comment | neutral |

Two of those lines invite doubt rather than confidence:

* `Convergence : Yes (10 + 212463 iter.)` reads as "this thing grinds".
* `Estimated smoothing value : 2997128305` is a huge unexplained number.

And `Number of events : 0 ( 0%)` only reads as a strength once the annotation
underneath explains it; on its own it looks like a broken fit.

---

## Option 1 (recommended) - what you can ask the fitted object for

Same console idiom as the panel above it, but every line is an output the
partial likelihood cannot hand you from these data without inventing event
times first.

```
> summary(fit)$Beta
                 Estimate Std. Error z-value Pr(>|z|)
treatmentRadChem  0.89355    0.29414 3.03787  0.00238

> predict(fit, type = "survival", i = 47, time = c(12, 24, 36, 48))
  time survival     se    low   high
    12   0.7630 0.0510 0.6611 0.8649
    24   0.4499 0.0574 0.3351 0.5648
    36   0.2273 0.0461 0.1351 0.3194
    48   0.0976 0.0353 0.0269 0.1683

> predict(fit, type = "risk", i = 47, time = 24)$risk
[1] 0.05034
```

Note underneath:

> Absolute survival at any time you name, with a standard error, and an
> instantaneous hazard rate to go with it - from 94 women, not one of whom has
> an observed event time.

Nine lines rather than fourteen, and it answers the question a trialist
actually arrives with ("what is the chance this patient is still event-free at
two years?") rather than describing the dataset.

`i = 47` is the first RadChem row of `bcos2`; if the index looks cryptic on the
poster, `i` can be dropped and the call left at mean covariates.

---

## Option 2 - side by side against the usual workaround

More directly promotional, but read the caveat before choosing it.

```
                               coxph() on midpoints   coxph_mpl() on intervals
treatment, log HR                0.906  (SE 0.285)      0.894  (SE 0.294)
S(24) for radio + chemo          0.443  (0.327-0.600)   0.450  (0.335-0.565)
h0(24)                           not estimable          0.0206 (SE 0.0031)
event times invented             56 of 94               none
```

**Caveat, and it matters:** the treatment effect is essentially identical
(0.906 against 0.894).  A sceptical reader will conclude the workaround was
fine, which is the opposite of the intended message.  What genuinely differs is
the third and fourth rows - a hazard *rate* the partial likelihood cannot
produce, and 56 fabricated event times it needs in order to run at all.  If you
take this option, consider dropping the first row entirely and letting the
panel be about h0(t) and the invented times.

---

## Option 3 - keep a console block, but curate it

If you want to keep the recognisable `summary()` output, cut it to the four
lines that carry information and drop the bookkeeping:

```
Cox Proportional Hazards Model Fit Using MPL

Number of obs.   : 94        Number of events : 0 ( 0%)
Estimated smoothing value : 2997128305   (chosen by marginal likelihood)

                 Estimate Std. Error z-value Pr(>|z|)
treatmentRadChem  0.89355    0.29414 3.03787  0.00238
```

This is the smallest change from what is on the poster now, but it is still
mostly a data description, so it inherits the original complaint.

---

## Recommendation

Option 1.  It is the only one of the three where every line is something the
package gives you and the alternatives do not, and it reads as a capability
rather than as a dataset summary.  Option 2 works as well but only if the
coefficient row comes out.

---

# FINAL: combined call / output / interpretation

Agreed shape, three tiers in one panel, replacing both the current
"A COMPLETE CALL" and "AND WHAT IT PRINTS" boxes.

An earlier draft quoted a hazard ratio of 2.44 (95% CI 1.37 to 4.35) in the
interpretation without that number appearing anywhere in the output above it -
the reader had to exponentiate 0.89355 in their head.  The reason was a gap in
the package: `coxph_mpl` had no `vcov` method, so `confint(fit)` failed with
`no applicable method for 'vcov'`.  `R/vcov.R` now supplies both, so the
hazard ratio and its interval can be printed by a call the audience already
knows.

```
A WORKED EXAMPLE                    bcos2 · breast cosmesis after radiotherapy

1  THE CALL
   library(survivalMPL)
   fit <- coxph_mpl(Surv(left, right, type = "interval2") ~ treatment,
                    data = bcos2, basis = "msplines")

2  WHAT IT RETURNS
   > summary(fit)$Beta
                    Estimate Std. Error z-value Pr(>|z|)
   treatmentRadChem  0.89355    0.29414 3.03787  0.00238

   > round(exp(cbind(HR = coef(fit), confint(fit))), 2)
                      HR 2.5 % 97.5 %
   treatmentRadChem 2.44  1.37   4.35

   > predict(fit, type = "survival", i = 47, time = c(24, 36))
     time survival     se    low   high
       24   0.4499 0.0574 0.3351 0.5648
       36   0.2273 0.0461 0.1351 0.3194

3  WHAT IT MEANS
   Adding chemotherapy raises the hazard of cosmetic deterioration by a
   factor of 2.44 (95% CI 1.37 to 4.35).  For a woman given both, the fitted
   probability of remaining free of deterioration is 45% (34 to 56) at two
   years and 23% (14 to 32) at three.  Every one of the 94 women is interval-,
   left- or right-censored: not one observed event time, and none invented.
```

Every number in step 3 now appears in step 2: 2.44 / 1.37 / 4.35 from the
`confint` line, 45% and 23% from the `predict` rows.

## Knock-on changes elsewhere on the poster

* Stat band: "6 S3 methods for fitted objects" becomes **8** (coef, confint,
  plot, predict, print, residuals, summary, vcov).
* SMOOTHING & METHODS chip row: adding `vcov()` and `confint()` takes the row
  past the right edge of its panel, so it needs either a second row of chips or
  slightly narrower ones.  Not done - flagging it as a layout decision.
