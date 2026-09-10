# ICTMC 2026 poster

Three build passes live here.

* **Pass 1** (`05-build-poster.R`) turned the original hand-laid-out deck in
  `source/poster-source.pptx` into `ICTMC-2026-survivalMPL-poster.pptx`.
* **Pass 2** (`06-update-v3.R`) works on `ICTMC-2026-survivalMPL-poster_v3.pptx`
  — the deck after the layout was reworked by hand — and writes
  `..._v4a.pptx` and `..._v4b.pptx`, which differ only in the FIG 2 image.

* **Pass 3** (`07-update-v5.R`) works on `ICTMC-2026-survivalMPL-poster_v5.pptx`
and writes `..._v6.pptx`: it merges the call and output panels into one
three-tier "A WORKED EXAMPLE" panel, takes the S3 method count from 6 to 8, and
rebuilds the chip row with `vcov()` and `confint()` at 19.5pt so all eight fit
on one line.

No pass ever edits its own input, so all three can be re-run at will.

## Rebuilding

```sh
Rscript dev/poster/01-fits.R                  # ~5 min, caches to cache/fits.rds
Rscript dev/poster/01b-fits-lambda.R          # ~6 min, adds the fixed-lambda fits
POSTER_FONT_DIR=/path/to/ibm-plex-ttf \
  Rscript dev/poster/03-fig2-variants.R
Rscript dev/poster/06-update-v3.R           # -> v4a, v4b
Rscript dev/poster/07-update-v5.R           # -> v6
```

`POSTER_FONT_DIR` should hold `IBMPlexSans-{Regular,SemiBold,Italic}.ttf` and
`IBMPlexMono-{Regular,SemiBold}.ttf` (from <https://github.com/IBM/plex>).
Without it the figures fall back to the platform sans/mono and no longer match
the poster's typography.

## Files

| file | what it does |
|---|---|
| `_theme.R` | palette, fonts, axis style, `baseline()` for h0/S0 at x = 0 |
| `_pptx.R` | small DrawingML writers and the slide-surgery helpers |
| `01-fits.R` | the four melanoma fits (one per basis) and the bcos2 fit |
| `01b-fits-lambda.R` | the same model at fixed smoothing values, for FIG 2b |
| `02-fig1-what-is-observed.R` | FIG 1 (unused in v3 onwards; kept for reuse) |
| `03-fig2-variants.R` | FIG 2 in both versions |
| `04-fig3-basis-shapes.R` | the psi_u(t) strip, cropped per basis card in v3 |
| `05-build-poster.R` | pass 1 |
| `06-update-v3.R` | pass 2 |

Figure sizes in the scripts are the sizes the pictures occupy on the A0 sheet,
so a point of type in a figure is a point of type on the poster.  FIG 2 is
16.04 x 4.40 in in both versions, which is the size the picture already has on
the v3 deck — swapping versions is therefore a byte swap of
`ppt/media/image1.png`, with no geometry change.

## A note on lambda

The smoothing value is not on a dataset-independent scale, and imposing one
from a cold start is not the same as reaching it through the outer loop.

* `melanoma` selects lambda ~ 1e-5.  Fixing lambda at anything from 1e-5 upward
  runs to the iteration cap; by lambda = 10 the penalised log-likelihood has
  collapsed to about -1e12.  Melanoma therefore has no room to illustrate the
  penalty at all.
* `bcos2` selects lambda ~ 3e9.  Fixed fits are well behaved up to about 1e4
  (df 13.0 -> 6.4), degrade quietly to 1e7 (ploglik -173 against -144 for the
  warm-started automatic fit at a *larger* lambda), and collapse by 1e8.

FIG 2b therefore shows lambda = 0 and lambda = 1e4 from `01b-fits-lambda.R`
plus the automatic fit itself for the heavily-penalised curve, so every line in
the panel is one the package would really return.  This gap between cold-start
and warm-started fits looks like a robustness issue in the multiplicative
update rather than a plotting problem.

## If a rebuild seems to do nothing

PowerPoint holds an exclusive lock on an open deck.  `06-update-v3.R` writes to
a scratch file and copies it into place, and falls back to `dev/poster/out/`
with a warning if the target cannot be written - it will not silently leave the
previous build on disk.

## Still to do by hand

* Two QR codes are placeholders: the "Poster QR" box in the header and the
  blank square in the CRAN panel.  Drop images in; the boxes are already sized.
* `IBM Plex Sans`, `IBM Plex Mono` and `Archivo` must be installed on whatever
  machine exports the final PDF, or PowerPoint substitutes and the line breaks
  move.
