## Third pass over the poster, applied to the hand-edited v5 deck.
##
##   in  : ICTMC-2026-survivalMPL-poster_v5.pptx   (never edited)
##   out : ICTMC-2026-survivalMPL-poster_v6.pptx
##
## Changes:
##   1. the "A COMPLETE CALL" and "AND WHAT IT PRINTS" panels merge into one
##      three-tier "A WORKED EXAMPLE" panel: the call, its output, and what the
##      output means, so a reader can enter at whichever depth suits them;
##   2. the stat band's S3 method count goes 6 -> 8, and the chip row and the
##      Post-processing stage gain vcov() and confint(), which R/vcov.R now
##      provides;
##   3. the chip text drops from 22.5pt to 19.5pt, because eight chips at the
##      old size overrun the right edge of their panel.
##
## Every number in tier 3 is printed by a call in tier 2 - that is the point of
## the panel, so do not quote a figure here that the output above it does not
## show.  All output is verbatim from the cached bcos2 msplines fit.

source(file.path("dev", "poster", "_pptx.R"))

SRC <- "ICTMC-2026-survivalMPL-poster_v5.pptx"
OUT <- "ICTMC-2026-survivalMPL-poster_v6.pptx"

INK  <- "121417"; INK2 <- "303338"; MUTE <- "606369"
BLUE <- "025894"; NAVY <- "003A68"; CHIP <- "D9EAFC"
RULE <- "D3D8E0"; PAPER <- "FDFCF9"; PANEL <- "F0F6FC"
MDOT <- "\u00b7"

stopifnot(file.exists(SRC))

## ---------------------------------------------------------------- geometry --
## Left column: the two panels being merged run 28.91 to 38.36, and
## DOCUMENTATION starts at 39.27.
EX_X <- 1.1458; EX_Y <- 28.91; EX_W <- 13.9486; EX_H <- 9.60
PAD  <- 0.48
BODY_W <- 13.20

## Chip row, inside the SMOOTHING & METHODS panel (x 15.80 to 31.35).
CHIP_X <- 16.1900; CHIP_Y <- 33.5600; CHIP_H <- 0.6146
CHIP_SZ <- 1950
CHIP_GAP <- 0.13
## 19.5pt IBM Plex Mono is about 0.1625 in per character; 0.42 in of padding
## keeps the same optical inset the hand-made chips had.
chip_w <- function(txt) nchar(txt) * 0.1625 + 0.42

## ------------------------------------------------------------------ unpack --
work <- file.path(tempdir(), "poster-v6")
unlink(work, recursive = TRUE)
dir.create(work, recursive = TRUE)
utils::unzip(SRC, exdir = work)

slide_path <- file.path(work, "ppt", "slides", "slide1.xml")
xml <- paste(readLines(slide_path, warn = FALSE, encoding = "UTF-8"),
             collapse = "\n")
s <- slide_split(xml)
sh <- s$shapes
message("v5 has ", length(sh), " shapes")

expect <- function(i, txt) {
  got <- paste(shape_text(sh[[i]]), collapse = "")
  if (!grepl(txt, got, fixed = TRUE)) {
    stop(sprintf("v5 shape %d is '%s', expected '%s'", i, got, txt))
  }
  invisible(TRUE)
}

## 1. ---------------------------------------------- the worked-example panel --
expect(129, "");                 # code panel background
expect(131, "A COMPLETE CALL")
expect(169, "AND WHAT IT PRINTS")
for (i in c(129:132, 167:171)) sh[[i]] <- ""

step <- function(n, label) {
  para(run(sprintf("%d  ", n), sz = 1950, bold = TRUE, col = BLUE,
           face = MONO_FACE),
       run(label, sz = 1950, col = BLUE, face = MONO_FACE, spc = 273))
}
code <- function(txt, col = INK) {
  para(run(txt, sz = 1800, col = col, face = MONO_FACE), lnSpc = 1.14)
}

## Tier 2 is verbatim console output.  exp(cbind(HR = coef(fit), confint(fit)))
## is what puts 2.44 / 1.37 / 4.35 on the poster; before R/vcov.R existed,
## confint() errored and those numbers had to be worked out by the reader.
out_lines <- c(
  "> summary(fit)$Beta",
  "                 Estimate Std. Error z-value Pr(>|z|)",
  "treatmentRadChem  0.89355    0.29414 3.03787  0.00238",
  "",
  "> round(exp(cbind(HR = coef(fit), confint(fit))), 2)",
  "                   HR 2.5 % 97.5 %",
  "treatmentRadChem 2.44  1.37   4.35",
  "",
  "> predict(fit, type = \"survival\", i = 47, time = c(24, 36))",
  "  time survival     se    low   high",
  "    24   0.4499 0.0574 0.3351 0.5648",
  "    36   0.2273 0.0461 0.1351 0.3194")

extra <- character()
addx <- function(...) extra <<- c(extra, ...)

addx(sp_rect(EX_X, EX_Y, EX_W, EX_H, fill = PANEL, name = "Example panel"))
addx(sp_rect(EX_X, EX_Y, 0.125, EX_H, fill = BLUE, name = "Example rule"))
addx(sp_text(EX_X + PAD, EX_Y + 0.24, 8.0, 0.36,
             para(run("A WORKED EXAMPLE", sz = 1950, col = BLUE,
                      face = MONO_FACE, spc = 273)),
             name = "Example label"))
addx(sp_text(EX_X + EX_W - 0.42 - 8.6, EX_Y + 0.26, 8.6, 0.34,
             para(run(paste0("bcos2 ", MDOT,
                             " breast cosmesis after radiotherapy"),
                      sz = 1650, col = MUTE, face = MONO_FACE),
                  algn = "r"),
             name = "Example dataset"))

addx(sp_text(EX_X + PAD, EX_Y + 0.80, 8.0, 0.36, step(1, "THE CALL"),
             name = "Step 1 label"))
addx(sp_text(EX_X + PAD, EX_Y + 1.24, BODY_W, 1.10, c(
  code("library(survivalMPL)"),
  code("fit <- coxph_mpl(Surv(left, right, type = \"interval2\") ~ treatment,"),
  code("                 data = bcos2, basis = \"msplines\")")),
  name = "Step 1 body"))

addx(sp_text(EX_X + PAD, EX_Y + 2.52, 8.0, 0.36, step(2, "WHAT IT RETURNS"),
             name = "Step 2 label"))
addx(sp_text(EX_X + PAD, EX_Y + 2.96, BODY_W, 4.20,
             lapply(out_lines, code),
             name = "Step 2 body"))

addx(sp_text(EX_X + PAD, EX_Y + 7.28, 8.0, 0.36, step(3, "WHAT IT MEANS"),
             name = "Step 3 label"))
addx(sp_text(EX_X + PAD, EX_Y + 7.72, BODY_W, 2.10,
             para(run(paste0(
               "Adding chemotherapy raises the hazard of cosmetic ",
               "deterioration by a factor of 2.44 (95% CI 1.37 to 4.35). For ",
               "a woman given both, the fitted probability of remaining free ",
               "of deterioration is 45% (34 to 56) at two years and 23% ",
               "(14 to 32) at three. Not one of the 94 women has an observed ",
               "event time, and none was invented."),
               sz = 2400, col = INK2), lnSpc = 1.15),
             name = "Step 3 body"))

## 2. ---------------------------------------------- the S3 method count and row --
## coef, confint, plot, predict, print, residuals, summary, vcov.
expect(63, "6"); sh[[63]] <- set_run_text(sh[[63]], 1, "8")
expect(99, "print ")
sh[[99]] <- set_run_text(sh[[99]], 1, paste(
  c("print", "summary", "predict", "residuals", "plot", "vcov", "confint"),
  collapse = paste0(" ", MDOT, " ")))

## The chip row is rebuilt rather than appended to, so all eight chips end up
## the same size at the smaller type.
expect(137, "coef()")
expect(147, "print()")
for (i in 136:147) sh[[i]] <- ""

chips <- c("coef()", "summary()", "predict()", "residuals()",
           "plot()", "print()", "vcov()", "confint()")
cx <- CHIP_X
for (ch in chips) {
  w <- chip_w(ch)
  addx(sp_rect(cx, CHIP_Y, w, CHIP_H, fill = CHIP, name = paste("Chip", ch)))
  addx(sp_text(cx, CHIP_Y + 0.09, w, 0.45,
               para(run(ch, sz = CHIP_SZ, col = NAVY, face = MONO_FACE),
                    algn = "ctr"),
               name = paste("Chip label", ch)))
  cx <- cx + w + CHIP_GAP
}
message(sprintf("chip row runs %.2f to %.2f in (panel edge 31.35)",
                CHIP_X, cx - CHIP_GAP))
if (cx - CHIP_GAP > 31.35) warning("chip row overruns its panel", call. = FALSE)

## ------------------------------------------------------------------ rebuild --
s$shapes <- sh
xml <- slide_join(s)
xml <- sub("</p:spTree>", paste0(paste0(extra, collapse = ""), "</p:spTree>"),
           xml, fixed = TRUE)

con <- file(slide_path, open = "wb")
writeLines(enc2utf8(xml), con, useBytes = TRUE)
close(con)

## PowerPoint keeps an exclusive lock on an open deck, so write to a scratch
## file and copy it into place, and say so plainly if that is not possible.
tmp_zip <- file.path(tempdir(), "v6.zip")
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
