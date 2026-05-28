# pkgdown Configuration

## Key facts

- Config: `_pkgdown.yml`
- Output directory: `pkgdown/` (`destination: pkgdown`)
- Homepage source: `index.md` at package root — **not** `README.md`
  (using `index.md` avoids a badge-stripping artifact that injects blank lines into code blocks)
- `README.md` is for GitHub only — keep badges and `r` language tags there

## The root .md rule

pkgdown builds **every `.md` file it finds at the package root** as a page.
Only these three are intentional:

| File | Purpose |
|------|---------|
| `index.md` | pkgdown homepage |
| `README.md` | GitHub repository page |
| `NEWS.md` | Changelog / news page |

Never leave other `.md` files at the root — they will be built and published.
Development docs go in `dev/`, historical docs in `old_sources/`.

## Navbar

Tutorials dropdown contains: Getting Started, Interval Censoring, Basis Functions, coxph vs coxph_mpl.
Defined in `_pkgdown.yml` under `navbar.components.tutorials`.

## Build commands

```r
pkgdown::build_home()      # homepage only
pkgdown::build_articles()  # vignettes only
pkgdown::build_site()      # full rebuild
```

## Blank-lines issue — RESOLVED

Each `build_site()` run was accumulating blank lines in code blocks in `pkgdown/index.html`.

**Root causes (both fixed):**
1. `pkgdown/index.html` was committed at an older pkgdown HTML format; each new build produced
   a different structure that accumulated empty `<span></span>` tags.
2. `README.md` had a double blank line (`\r\n\r\n\r\n`) inside the Quick Start code block,
   between the closing `)` and `summary(fit_lung)`. Fixed by direct binary patch (the Edit tool
   could not match CRLF sequences line-by-line).

**Standard workflow going forward:**
```r
pkgdown::clean_site()   # flush stale output — run at least once after any pkgdown version change
pkgdown::build_site()
```
