## Fourth pass over the poster.
##
##   in  : ICTMC-2026-survivalMPL-poster_v6.pptx   (never edited)
##   out : ICTMC-2026-survivalMPL-poster_v7.pptx
##
## One change, in five parts: FIG 2 becomes the taller, more legible version
## from 03c-fig2-readable.R.
##
##   1. ppt/media/image1.png is replaced with fig2-readable.png.
##   2. The picture box grows from 4.27 to 5.10 in.  It has to match the file's
##      aspect exactly or PowerPoint stretches the image, which would undo the
##      whole point of the exercise.
##   3. The height comes from two places.  The SMOOTHING & METHODS panel has
##      0.35 in of dead space below its chip row, so its box shrinks; and the
##      IN PRACTICE label with the closing paragraph move down 0.22 in, which is
##      the slack between the paragraph and the footer rule.
##   4. The caption moves up with the panel above it.
##   5. The two paper-coloured rectangles that masked the old figure's in-plot
##      titles are deleted: the new figure has no subtitles, so they would now
##      sit over the data.

source(file.path("dev", "poster", "_pptx.R"))

SRC <- "ICTMC-2026-survivalMPL-poster_v6.pptx"
OUT <- "ICTMC-2026-survivalMPL-poster_v7.pptx"
FIG <- file.path("dev", "poster", "figures", "fig2-readable.png")

stopifnot(file.exists(SRC), file.exists(FIG))

## Geometry, in inches.  Vertical budget below the basis cards:
##   SMOOTHING panel   31.51 -> 34.29   (was 34.52)
##   caption           34.38 -> 34.74
##   FIG 2             34.80 -> 39.90
##   IN PRACTICE       39.96
##   closing paragraph 40.56 -> 43.24   (footer rule at 43.35)
SMOOTH_H  <- 2.78
CAP_Y     <- 34.38
FIG_X     <- 15.82; FIG_Y <- 34.80; FIG_W <- 15.56; FIG_H <- 5.10
SHIFT     <- 0.22

work <- file.path(tempdir(), "poster-v7")
unlink(work, recursive = TRUE)
dir.create(work, recursive = TRUE)
utils::unzip(SRC, exdir = work)

slide_path <- file.path(work, "ppt", "slides", "slide1.xml")
xml <- paste(readLines(slide_path, warn = FALSE, encoding = "UTF-8"),
             collapse = "\n")
s <- slide_split(xml)
sh <- s$shapes
message("v6 has ", length(sh), " shapes")

expect <- function(i, txt) {
  got <- paste(shape_text(sh[[i]]), collapse = "")
  if (!grepl(txt, got, fixed = TRUE)) {
    stop(sprintf("v6 shape %d is '%s', expected '%s'", i, got, txt))
  }
  invisible(TRUE)
}
expect(102, "Four bases, automatic smoothing")
expect(130, "SMOOTHING")
expect(132, "By combining a rigorous penalised likelihood")
expect(152, "IN PRACTICE")

## The picture and the two masks carry no text, so check them by name.
shape_name <- function(i) sub('^.*name="([^"]*)".*$', "\\1", sh[[i]])
stopifnot(shape_name(151) == "Fig 2",
          shape_name(153) == "Rectangle 28",
          shape_name(154) == "Rectangle 29")

sh[[129]] <- set_xfrm(sh[[129]], h = SMOOTH_H)          # SMOOTHING panel box
sh[[102]] <- set_xfrm(sh[[102]], y = CAP_Y)             # FIG 2 caption
sh[[151]] <- set_xfrm(sh[[151]], x = FIG_X, y = FIG_Y,  # the picture itself
                      w = FIG_W, h = FIG_H)
sh[[152]] <- set_xfrm(sh[[152]], y = 39.74 + SHIFT)     # IN PRACTICE label
sh[[132]] <- set_xfrm(sh[[132]], y = 40.34 + SHIFT)     # closing paragraph
sh[[153]] <- ""                                         # masks over old titles
sh[[154]] <- ""

s$shapes <- sh
con <- file(slide_path, open = "wb")
writeLines(enc2utf8(slide_join(s)), con, useBytes = TRUE)
close(con)

## Resolve the media file from the picture's own relationship rather than
## assuming a filename.  PowerPoint renumbers media on save: in v3 the FIG 2
## picture was image1.png, but by v5 it had become image2.png and image1.png
## was the FIG 3 basis strip.  Hard-coding the name replaced the wrong file and
## silently broke the four basis-card crops while leaving FIG 2 untouched.
rels_path <- file.path(work, "ppt", "slides", "_rels", "slide1.xml.rels")
rels <- paste(readLines(rels_path, warn = FALSE), collapse = "")
rid <- regmatches(sh[[151]], regexpr('r:embed="[^"]+"', sh[[151]]))
rid <- sub('r:embed="([^"]+)"', "\\1", rid)
stopifnot(length(rid) == 1L, nzchar(rid))
pat <- sprintf('Id="%s"[^>]*Target="\\.\\./([^"]+)"', rid)
target <- sub(pat, "\\1", regmatches(rels, regexpr(pat, rels)))
stopifnot(length(target) == 1L, grepl("^media/", target))
message("FIG 2 picture uses ", rid, " -> ", target)
stopifnot(file.copy(FIG, file.path(work, "ppt", target), overwrite = TRUE))

## PowerPoint holds an exclusive lock on an open deck: write to a scratch file
## and copy it into place, and say so plainly if that is not possible.
tmp_zip <- file.path(tempdir(), "v7.zip")
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
