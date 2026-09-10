## Fifth pass over the poster: three edits to the likelihood-term column of the
## OBSERVATION SCHEMES SUPPORTED panel.
##
##   in  : ICTMC-2026-survivalMPL-poster_v7.pptx   (never edited)
##   out : ICTMC-2026-survivalMPL-poster_v8.pptx
##
## 1. Exact event: f(t) becomes h(t) S(t).  Both are correct - R/coxph.r:716
##    scores exact events as the derivative of log h(t) + log S(t) = log f(t) -
##    but f was the only symbol in the column not written in S, so a reader had
##    to supply f = -S' themselves.  h(t) S(t) is self-consistent with the other
##    four rows and puts the hazard on the poster in the one row that needs it,
##    which is the quantity MPL estimates and the partial likelihood profiles
##    away.
##
## 2 and 3. Left truncation is not a fifth contribution, it is a divisor applied
##    to whichever of the four rows above describes the subject: the sample is
##    conditioned on T > entry, so every term is divided by S(entry).  In the
##    code this is H(t) - H(entry) applied to the cumulative term only
##    (R/coxph.r:236), which is the same operation since
##    exp{-mu[H(t)-H(e)]} = S(t)/S(e).  Sitting in a column of five standalone
##    terms it read as a fifth one, so it now carries a leading ellipsis and is
##    set in the muted grey used for annotations rather than the blue used for
##    the four real terms.

source(file.path("dev", "poster", "_pptx.R"))

SRC <- "ICTMC-2026-survivalMPL-poster_v7.pptx"
OUT <- "ICTMC-2026-survivalMPL-poster_v8.pptx"

BLUE <- "025894"
MUTE <- "606369"
ELL  <- "\u2026"   # horizontal ellipsis
DIV  <- "\u00f7"   # division sign, as written on the deck

stopifnot(file.exists(SRC))

work <- file.path(tempdir(), "poster-v8")
unlink(work, recursive = TRUE)
dir.create(work, recursive = TRUE)
utils::unzip(SRC, exdir = work)

slide_path <- file.path(work, "ppt", "slides", "slide1.xml")
xml <- paste(readLines(slide_path, warn = FALSE, encoding = "UTF-8"),
             collapse = "\n")
s <- slide_split(xml)
sh <- s$shapes
message("v7 has ", length(sh), " shapes")

expect <- function(i, txt) {
  got <- paste(shape_text(sh[[i]]), collapse = "")
  if (!identical(got, txt)) {
    stop(sprintf("v7 shape %d is '%s', expected '%s'", i, got, txt))
  }
  invisible(TRUE)
}

## 1. -------------------------------------------------- exact event term ----
expect(22, "f(t)")
sh[[22]] <- set_run_text(sh[[22]], 1, "h(t) S(t)")

## 2 and 3. ------------------------------------------ left truncation term ----
expect(50, paste0(DIV, " S(entry)"))
sh[[50]] <- set_run_text(sh[[50]], 1, paste0(ELL, " ", DIV, " S(entry)"))
## Recolour just this run: the four real contributions keep the blue.
n_before <- length(gregexpr(BLUE, sh[[50]], fixed = TRUE)[[1]])
stopifnot(n_before == 1L)
sh[[50]] <- sub(BLUE, MUTE, sh[[50]], fixed = TRUE)

s$shapes <- sh
con <- file(slide_path, open = "wb")
writeLines(enc2utf8(slide_join(s)), con, useBytes = TRUE)
close(con)

## PowerPoint holds an exclusive lock on an open deck: write to a scratch file
## and copy it into place, and say so plainly if that is not possible.
tmp_zip <- file.path(tempdir(), "v8.zip")
if (file.exists(tmp_zip)) unlink(tmp_zip)
old <- setwd(work)
files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
rc <- utils::zip(tmp_zip, files, flags = "-q -X")
setwd(old)
if (rc != 0 || !file.exists(tmp_zip)) stop("zip failed")

if (!suppressWarnings(file.copy(tmp_zip, OUT, overwrite = TRUE))) {
  alt_dir <- file.path("dev", "poster", "out")
  dir.create(alt_dir, recursive = TRUE, showWarnings = FALSE)
  OUT <- file.path(alt_dir, basename(OUT))
  if (!file.copy(tmp_zip, OUT, overwrite = TRUE)) stop("could not write ", OUT)
  warning("target locked; wrote ", OUT, " instead", call. = FALSE)
}
message("wrote ", OUT, " (", round(file.size(OUT) / 1024), " KB)")
