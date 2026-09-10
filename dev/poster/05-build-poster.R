## Build the ICTMC 2026 poster from the hand-laid-out source deck.
##
##   in  : dev/poster/source/poster-source.pptx   (never edited)
##   out : ICTMC-2026-survivalMPL-poster.pptx     (regenerated every run)
##
## What this does, in order:
##   1. strips the yellow "fix this wording" highlights,
##   2. corrects the counts and claims that the code does not support,
##   3. rebuilds the observation-schemes diagram with a legend and a fifth row,
##   4. drops in the three figures, a worked call, and the repository URL.
##
## Everything is expressed in inches from the top-left of the A0 sheet
## (33.11 x 46.81 in), matching how the original shapes were placed.

source(file.path("dev", "poster", "_pptx.R"))

SRC <- file.path("dev", "poster", "source", "poster-source.pptx")
OUT <- "ICTMC-2026-survivalMPL-poster.pptx"
FIG <- file.path("dev", "poster", "figures")

stopifnot(file.exists(SRC))
for (f in c("fig1-what-is-observed.png", "fig2-four-bases.png",
            "fig3-basis-shapes.png")) {
  stopifnot(file.exists(file.path(FIG, f)))
}

## ---------------------------------------------------------------- palette --
INK    <- "121417"; INK2 <- "303338"; MUTE <- "606369"
BLUE   <- "025894"; NAVY <- "003A68"; BLUE_LT <- "B0C7DD"
ORANGE <- "EB883B"; ORANGE_D <- "903300"
RULE   <- "D3D8E0"; RULE2 <- "CCD1D9"; PAPER <- "FDFCF9"; PANEL <- "F0F6FC"

## Characters that must not depend on the encoding of this file.
EMD <- "\u2014"; MDOT <- "\u00b7"; MINUS <- "\u2212"; DIV <- "\u00f7"
LAM <- "\u03bb"; THT <- "\u03b8"; GE <- "\u2265"; TIMES <- "\u00d7"
SUB1 <- "\u2081"; SUB2 <- "\u2082"

## ------------------------------------------------------------------ unpack --
work <- file.path(tempdir(), "poster-build")
unlink(work, recursive = TRUE)
dir.create(work, recursive = TRUE)
utils::unzip(SRC, exdir = work)

slide_path <- file.path(work, "ppt", "slides", "slide1.xml")
xml <- paste(readLines(slide_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

## 1. --------------------------------------------------- drop the highlights --
n_hl <- length(gregexpr("<a:highlight>", xml, fixed = TRUE)[[1]])
## Stale autofit line-squeeze values left over from the original wording; they
## make PowerPoint re-draw the last line of a box on export.
xml <- gsub(' lnSpcReduction="[0-9]+"', "", xml)
xml <- gsub("<a:highlight>.*?</a:highlight>", "", xml, perl = TRUE)
message("removed ", n_hl, " highlight runs")

s <- slide_split(xml)
sh <- s$shapes
message("slide has ", length(sh), " shapes")

## Guard against editing the wrong shape if the source deck is ever re-laid out.
expect <- function(i, txt) {
  got <- paste(shape_text(sh[[i]]), collapse = "")
  if (!grepl(txt, got, fixed = TRUE)) {
    stop(sprintf("shape %d is '%s', expected to contain '%s'", i, got, txt))
  }
  invisible(TRUE)
}

## 2. ------------------------------------------------- corrections to claims --

## THE PROBLEM: the likelihood also takes exact events, and entry= adds delayed
## entry, so say so.
expect(10, "Most survival software")
sh[[10]] <- set_run_text(sh[[10]], 3, paste0(
  " fits a Cox model to exact, right-, left- and interval-censored data in ",
  "one likelihood, with optional delayed entry, and returns a smooth ",
  "estimate of the baseline hazard rather than a step function."))

## Stat band.  Four observation schemes plus left truncation is five; there are
## six S3 generics for coxph_mpl objects and eight articles on the site.
expect(40, "4"); sh[[40]] <- set_run_text(sh[[40]], 1, "5")
expect(41, "censoring schemes")
sh[[41]] <- set_run_text(sh[[41]], 1, "observation schemes, one likelihood")
expect(43, "5"); sh[[43]] <- set_run_text(sh[[43]], 1, "6")
expect(44, "S3 method families")
sh[[44]] <- set_run_text(sh[[44]], 1, "S3 methods for fitted objects")
expect(45, "7"); sh[[45]] <- set_run_text(sh[[45]], 1, "8")
## The bullet list should name all eight articles.
expect(121, "Basis comparison")
sh[[121]] <- set_xfrm(sh[[121]], w = 6.0)
sh[[121]] <- set_run_text(sh[[121]], 1, "Bases and control parameters")
expect(46, "worked tutorials")
sh[[46]] <- set_run_text(sh[[46]], 1, "worked tutorials on the docs site")

## Architecture stages: name what each one actually takes or returns.
expect(55, "response encoding")
sh[[55]] <- set_run_text(sh[[55]], 1,
  paste0("formula ", MDOT, ' Surv(type = "interval2") ', MDOT, " entry="))
## Smoothing is chosen by a marginal likelihood (Laplace approximation), not by
## REML - see coxph_mpl.control() and the control-parameters article.
expect(63, "REML")
sh[[63]] <- set_run_text(sh[[63]], 1,
  paste0("penalised likelihood ", MDOT, " non-negativity ", MDOT, " marginal-likelihood ", LAM))
sh[[63]] <- set_run_text(sh[[63]], 2, "")
expect(79, "summary")
sh[[79]] <- set_run_text(sh[[79]], 1,
  paste0("print ", MDOT, " summary ", MDOT, " predict ", MDOT, " residuals ", MDOT, " plot"))

## EXTENSION POINT: the registry holds bases, not censoring schemes, and the
## memory win is a specific one - the Hessian blocks no longer build n x n
## diagonal matrices, which is what makes the 86,611-subject fit feasible.
expect(83, "basis registry")
sh[[83]] <- set_run_text(sh[[83]], 1, paste0(
  "A basis is one file: three functions, registered at load time. A fifth ",
  "family is a new file, not an edit to "))
sh[[83]] <- set_run_text(sh[[83]], 2, "coxph_mpl()")
sh[[83]] <- set_run_text(sh[[83]], 3, paste0(
  ". The same separation let the Hessian assembly drop its n ", TIMES, " n ",
  "matrices ", EMD, " the reason the 86,611-subject ", "hiroshima",
  " fit runs at all ", EMD, " with estimates unchanged."))

## BASELINE HAZARD BASES strapline.
expect(85, "pick how it is represented")
sh[[85]] <- set_run_text(sh[[85]], 1,
  "The baseline hazard is a non-negative mixture of basis ")
sh[[85]] <- set_run_text(sh[[85]], 2, "functions.")

## Basis cards.  "uniform" is the package default; "msplines" is the
## recommendation, not the default.
expect(88, "Piecewise-constant")
sh[[88]] <- set_run_text(sh[[88]], 1, paste0(
  "Piecewise-constant steps. The package default, and the only basis that ",
  "supports delayed entry."))
expect(91, "Non-negative spline")
sh[[91]] <- set_run_text(sh[[91]], 1, paste0(
  "Cubic M-splines, non-negative by construction. Closed-form roughness ",
  "penalty; the usual recommendation."))
expect(94, "Kernel basis")
sh[[94]] <- set_run_text(sh[[94]], 1, paste0(
  "Truncated Gaussian kernels. The smoothest estimates; first- or ",
  "second-order penalty."))
expect(97, "Compact-support kernel")
sh[[97]] <- set_run_text(sh[[97]], 1, paste0(
  "Compact parabolic kernels. Local support like M-splines; second-order ",
  "penalty only."))

## SMOOTHING & METHODS body.
expect(101, "restricted maximum likelihood")
sh[[101]] <- set_run_text(sh[[101]], 1, paste0(
  LAM, " is selected automatically by marginal likelihood, or held fixed at a ",
  "value you supply. Non-negativity of ", THT, " is a constraint inside the ",
  "optimiser, not a correction applied afterwards."))
sh[[101]] <- set_run_text(sh[[101]], 2, "")

## Datasets: there are four, and only two of them are real.
expect(126, "Three real datasets")
sh[[126]] <- set_run_text(sh[[126]], 1, "Four datasets ship with the package: ")
sh[[126]] <- set_run_text(sh[[126]], 2, "bcos2 ")
sh[[126]] <- set_run_text(sh[[126]], 3, "(94 women, breast cosmesis), ")
sh[[126]] <- set_run_text(sh[[126]], 4, "hiv ")
sh[[126]] <- set_run_text(sh[[126]], 5, "and ")
sh[[126]] <- set_run_text(sh[[126]], 6, "melanoma ")
sh[[126]] <- set_run_text(sh[[126]], 7, "(300 simulated subjects each), and ")
sh[[126]] <- insert_run(sh[[126]], after = 7, template = 6, text = "hiroshima ")
sh[[126]] <- insert_run(sh[[126]], after = 8, template = 7, text = paste0(
  "(86,611 subjects, delayed entry). Every example runs as written."))

## WHY IT MATTERS.
expect(133, "Rigorous penalised likelihood")
sh[[133]] <- set_xfrm(sh[[133]], h = 1.60)
sh[[133]] <- gsub('<a:lnSpc><a:spcPct val="109091"/></a:lnSpc>', "", sh[[133]], fixed = TRUE)
sh[[133]] <- set_run_text(sh[[133]], 1, paste0(
  "Scheduled visits make interval censoring the norm. survivalMPL fits the ",
  "intervals directly ", EMD, " no imputation, nothing discarded ", EMD,
  " and returns a smooth baseline hazard with standard errors."))

## 3. --------------------------------- rebuild the observation-schemes panel --
## The old panel had four rows headed "CENSORING SCHEMES SUPPORTED", which both
## missed exact events and filed left truncation under censoring.  It also gave
## the same blue bar two different meanings from row to row.  One grammar now:
##   blue   = under observation and event-free
##   orange = the event is known to lie in this window
##   tick   = an assessment visit
##   dot    = an exactly recorded event time

X0 <- 17.9845; L <- 13.4738; X1 <- X0 + L
BOX_X <- 17.4741; BOX_Y <- 8.7126; BOX_W <- 14.4946; BOX_H <- 7.60

TITLE_Y <- 9.1813
LEG_Y   <- 9.70
ROW_Y   <- 10.32
PITCH   <- 1.155

AX_DY   <- 0.84      # axis line, relative to the row's label top
BAR_H   <- 0.052
BAND_DY <- 0.795; BAND_H <- 0.14
TICK_DY <- 0.66; TICK_H <- 0.34; TICK_W <- 0.055
ANN_DY  <- 0.40; ANN_H  <- 0.32
DOT     <- 0.185

at <- function(f) X0 + f * L

panel <- character()
add <- function(...) panel <<- c(panel, ...)

add(sp_text(X0, TITLE_Y, 14.82, 0.396,
            para(run("OBSERVATION SCHEMES SUPPORTED", sz = 1950, col = MUTE,
                     face = MONO_FACE, spc = 273)),
            name = "Panel title"))

## legend
leg <- list(
  list(lab = "under observation", kind = "bar",  col = BLUE),
  list(lab = "event lies in here", kind = "band", col = ORANGE),
  list(lab = "assessment visit",   kind = "tick", col = ORANGE_D),
  list(lab = "exact event time",   kind = "dot",  col = ORANGE_D)
)
for (i in seq_along(leg)) {
  lx <- X0 + (i - 1) * 3.37
  e <- leg[[i]]
  if (e$kind == "bar")  add(sp_rect(lx, LEG_Y + 0.175, 0.46, BAR_H, fill = e$col))
  if (e$kind == "band") add(sp_rect(lx, LEG_Y + 0.13, 0.46, 0.13, fill = e$col))
  if (e$kind == "tick") add(sp_rect(lx + 0.2, LEG_Y + 0.06, TICK_W, 0.28, fill = e$col))
  if (e$kind == "dot")  add(sp_rect(lx + 0.14, LEG_Y + 0.11, DOT, DOT,
                                    fill = e$col, geom = "ellipse"))
  add(sp_text(lx + 0.62, LEG_Y, 2.70, 0.40,
              para(run(e$lab, sz = 1700, col = MUTE, face = MONO_FACE)),
              name = "Legend label"))
}

## rows: label, likelihood term, and the marks on the timeline
rows <- list(
  list(lab = "Exact event", term = "f(t)",
       blue = c(0, 0.52), dot = 0.52,
       ann = list(list(f = 0.52, txt = "event time recorded"))),
  list(lab = "Right censored", term = "S(t)",
       blue = c(0, 0.52), band = c(0.52, 1), tick = 0.52,
       ann = list(list(f = 0.52, txt = "last contact"))),
  list(lab = "Left censored", term = paste0("1 ", MINUS, " S(t)"),
       band = c(0, 0.28), tick = 0.28,
       ann = list(list(f = 0.28, txt = "first assessment"))),
  list(lab = "Interval censored",
       term = paste0("S(t", SUB1, ") ", MINUS, " S(t", SUB2, ")"),
       blue = c(0, 0.36), band = c(0.36, 0.70), tick = c(0.36, 0.70),
       ann = list(list(f = 0.36, txt = "visit k"),
                  list(f = 0.70, txt = "visit k+1"))),
  list(lab = "Left truncated", term = paste0(DIV, " S(entry)"),
       blue = c(0.24, 1), tick = 0.24,
       ann = list(list(f = 0.24, txt = "delayed entry")))
)

for (k in seq_along(rows)) {
  r  <- rows[[k]]
  y0 <- ROW_Y + (k - 1) * PITCH

  add(sp_text(X0, y0, 9.0, 0.46,
              para(run(r$lab, sz = 2550, bold = TRUE, col = INK)),
              name = "Row label"))
  add(sp_text(X1 - 5.6, y0 + 0.03, 5.6, 0.42,
              para(run(r$term, sz = 2050, col = BLUE, face = MONO_FACE),
                   algn = "r"),
              name = "Likelihood term"))

  ## timescale, then the marks on top of it
  add(sp_rect(X0, y0 + AX_DY, L, BAR_H, fill = RULE2, name = "Timescale"))
  if (!is.null(r$grey))
    add(sp_rect(at(r$grey[1]), y0 + AX_DY, diff(r$grey) * L, BAR_H, fill = RULE))
  if (!is.null(r$blue))
    add(sp_rect(at(r$blue[1]), y0 + AX_DY, diff(r$blue) * L, BAR_H, fill = BLUE))
  if (!is.null(r$band))
    add(sp_rect(at(r$band[1]), y0 + BAND_DY, diff(r$band) * L, BAND_H, fill = ORANGE))
  for (f in r$tick)
    add(sp_rect(at(f) - TICK_W / 2, y0 + TICK_DY, TICK_W, TICK_H, fill = ORANGE_D))
  if (!is.null(r$dot))
    add(sp_rect(at(r$dot) - DOT / 2, y0 + AX_DY + BAR_H / 2 - DOT / 2, DOT, DOT,
                fill = ORANGE_D, geom = "ellipse"))
  for (a in r$ann)
    add(sp_text(at(a$f) - 1.35, y0 + ANN_DY, 2.70, ANN_H,
                para(run(a$txt, sz = 1650, col = MUTE, face = MONO_FACE),
                     algn = "ctr"),
                name = "Row annotation"))
}

## A footnote: truncation is orthogonal to the four censoring types, and the
## entry= implementation is basis-limited.
add(sp_text(X0, ROW_Y + 5 * PITCH - 0.20, L, 0.34,
            para(run(paste0("Left truncation combines with any of the four; ",
                            "entry= currently requires the uniform basis."),
                     sz = 1650, col = MUTE, face = MONO_FACE)),
            name = "Panel footnote"))

## The old panel is shapes 12..34; shape 11 is the box behind it.
expect(11, "")
expect(12, "CENSORING SCHEMES SUPPORTED")
expect(34, "delayed entry")
sh[[11]] <- set_xfrm(sh[[11]], y = BOX_Y, h = BOX_H)
sh[[12]] <- paste0(panel, collapse = "")
for (i in 13:34) sh[[i]] <- ""

## 4. ------------------------------------------------------ figures and code --
## FIG 2's caption already exists (it was written for a figure that was never
## placed); move it to sit above the figure it describes.
expect(98, "FIG 1")
sh[[98]] <- set_xfrm(sh[[98]], x = 15.9277, y = 25.58, w = 16.04, h = 0.36)
sh[[98]] <- set_run_text(sh[[98]], 1, paste0(
  "FIG 2  Same melanoma model, four bases. Smoothing chosen automatically for each."))

extra <- character()
addx <- function(...) extra <<- c(extra, ...)

cap <- function(x, y, w, txt) {
  sp_text(x, y, w, 0.36,
          para(run(txt, sz = 1950, col = MUTE, face = MONO_FACE, spc = 60)),
          name = "Figure caption")
}

## FIG 1, under THE PROBLEM
addx(cap(1.1458, 13.36, 15.40, paste0(
  "FIG 1  bcos2, the breast-cosmesis trial: what is recorded, and what comes back.")))
addx(sp_pic(1.1458, 13.78, 15.40, 2.55, "rIdFig1", name = "Fig 1"))

## FIG 2, under the basis cards (caption is shape 98, moved above)
addx(sp_pic(15.9277, 26.02, 16.04, 4.40, "rIdFig2", name = "Fig 2"))

## FIG 3, the strip above WHY IT MATTERS
addx(cap(1.1458, 38.96, 24.80, paste0(
  "FIG 3  The four basis families: psi_u(t), drawn with 8 knots for legibility.")))
addx(sp_pic(1.1458, 39.38, 24.80, 2.35, "rIdFig3", name = "Fig 3"))

## A worked call, in the gap under EXTENSION POINT.
CODE_X <- 1.1458; CODE_Y <- 31.58; CODE_W <- 13.9486; CODE_H <- 2.46
addx(sp_rect(CODE_X, CODE_Y, CODE_W, CODE_H, fill = PANEL, name = "Code panel"))
addx(sp_rect(CODE_X, CODE_Y, 0.125, CODE_H, fill = BLUE, name = "Code rule"))
addx(sp_text(CODE_X + 0.48, CODE_Y + 0.20, 12.0, 0.36,
             para(run("A COMPLETE CALL", sz = 1950, col = BLUE,
                      face = MONO_FACE, spc = 273)),
             name = "Code label"))

code_line <- function(txt, col = INK) para(run(txt, sz = 1800, col = col,
                                               face = MONO_FACE), lnSpc = 1.14)
addx(sp_text(CODE_X + 0.48, CODE_Y + 0.64, 13.2, 1.7, c(
  code_line("library(survivalMPL)"),
  code_line('fit <- coxph_mpl(Surv(left, right, type = "interval2") ~ treatment,'),
  code_line('                 data = bcos2, basis = "msplines")'),
  code_line("summary(fit)   # coefficients, lambda, degrees of freedom", col = MUTE),
  code_line("plot(fit)      # bases, hazard, cumulative hazard, survival", col = MUTE)),
  name = "Code body"))

## One more method chip, so the row matches the six generics.
addx(sp_rect(26.7923, 32.6460, 1.8333, 0.6146, fill = "D9EAFC",
             name = "Chip print"))
addx(sp_text(27.0006, 32.7502, 1.5000, 0.4479,
             para(run("print()", sz = 2100, col = BLUE, face = MONO_FACE),
                  algn = "ctr"),
             name = "Chip print label"))

## The poster had no address on it anywhere.
addx(sp_text(12.20, 4.62, 16.0, 0.52,
             para(run("github.com/survivalMPL/survivalMPL", sz = 2400,
                      col = BLUE_LT, face = MONO_FACE)),
             name = "Repo URL"))

## ------------------------------------------------------------------ rebuild --
s$shapes <- sh
xml <- slide_join(s)
xml <- sub("</p:spTree>", paste0(paste0(extra, collapse = ""), "</p:spTree>"),
           xml, fixed = TRUE)

con <- file(slide_path, open = "wb")
writeLines(enc2utf8(xml), con, useBytes = TRUE)
close(con)

## media + relationships + content types
dir.create(file.path(work, "ppt", "media"), showWarnings = FALSE)
figs <- c(rIdFig1 = "fig1-what-is-observed.png",
          rIdFig2 = "fig2-four-bases.png",
          rIdFig3 = "fig3-basis-shapes.png")
for (i in seq_along(figs)) {
  file.copy(file.path(FIG, figs[i]),
            file.path(work, "ppt", "media", sprintf("image%d.png", i)),
            overwrite = TRUE)
}

rels_path <- file.path(work, "ppt", "slides", "_rels", "slide1.xml.rels")
rels <- paste(readLines(rels_path, warn = FALSE), collapse = "")
new_rels <- paste0(sprintf(paste0('<Relationship Id="%s" Type="http://schemas.openxmlformats.org/',
                                  'officeDocument/2006/relationships/image" Target="../media/image%d.png"/>'),
                           names(figs), seq_along(figs)), collapse = "")
rels <- sub("</Relationships>", paste0(new_rels, "</Relationships>"), rels, fixed = TRUE)
writeLines(rels, rels_path)

ct_path <- file.path(work, "[Content_Types].xml")
ct <- paste(readLines(ct_path, warn = FALSE), collapse = "")
if (!grepl('Extension="png"', ct, fixed = TRUE)) {
  ct <- sub("<Default Extension=\"rels\"",
            "<Default Extension=\"png\" ContentType=\"image/png\"/><Default Extension=\"rels\"",
            ct, fixed = TRUE)
}
writeLines(ct, ct_path)

## ------------------------------------------------------------------- rezip --
if (file.exists(OUT)) unlink(OUT)
old <- setwd(work)
files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
utils::zip(file.path(old, OUT), files, flags = "-q -X")
setwd(old)

message("wrote ", OUT, " (", round(file.size(OUT) / 1024), " KB)")
