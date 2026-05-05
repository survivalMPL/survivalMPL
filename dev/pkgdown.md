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

Tutorials dropdown contains: Getting Started, Interval Censoring.
Defined in `_pkgdown.yml` under `navbar.components.tutorials`.

## Build commands

```r
pkgdown::build_home()      # homepage only
pkgdown::build_articles()  # vignettes only
pkgdown::build_site()      # full rebuild
```
