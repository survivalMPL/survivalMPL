## Second pass over the poster, applied to the hand-edited v3 deck.
##
##   in  : ICTMC-2026-survivalMPL-poster_v3.pptx   (never edited)
##   out : ICTMC-2026-survivalMPL-poster_v4a.pptx  (FIG 2 = two panels)
##         ICTMC-2026-survivalMPL-poster_v4b.pptx  (FIG 2 = three, incl. lambda)
##
## Adds, on top of v3:
##   1. two more paragraphs under THE PROBLEM, filling the gap above the stats,
##   2. an output panel under A COMPLETE CALL showing what summary() prints,
##   3. a replacement FIG 2 image (the picture keeps v3's position and size,
##      only ppt/media/image1.png changes).
##
## v3's own layout is left alone otherwise.  Shape indices below refer to v3.

source(file.path("dev", "poster", "_pptx.R"))

SRC <- "ICTMC-2026-survivalMPL-poster_v3.pptx"
FIG <- file.path("dev", "poster", "figures")
## Each variant supplies its FIG 2 image, plus replacement caption text where
## the picture no longer matches the caption v3 already carries.
VARIANTS <- list(
  v4a = list(png = "fig2a-two-panels.png", caption = NULL),
  v4b = list(png = "fig2b-with-lambda.png",
             caption = paste0(
               "Four bases on melanoma, smoothing automatic; ",
               "right-hand panel varies lambda on bcos2.")))

stopifnot(file.exists(SRC))

INK  <- "121417"; INK2 <- "303338"; MUTE <- "606369"
BLUE <- "025894"; RULE <- "D3D8E0"; PAPER <- "FDFCF9"
EMD  <- "\u2014"; BETA <- "\u03b2"; SUB0 <- "\u2080"; NDASH <- "\u2013"

## ---------------------------------------------------------------- geometry --
## THE PROBLEM: body text ends at 13.28, the stats band starts at 16.51.
PROB_X <- 1.1458; PROB_Y <- 13.42; PROB_W <- 14.484

## Left column below the code panel (which ends at 31.37) is empty down to
## DOCUMENTATION at 39.27.
OUT_X <- 1.1458; OUT_Y <- 31.55; OUT_W <- 13.947; OUT_H <- 7.35
PAD   <- 0.48

body_para <- function(txt, first = FALSE) {
  para(run(txt, sz = 2700, col = INK2), lnSpc = 1.103,
       space_before = if (first) NULL else 14)
}

console <- function(txt, col = INK) {
  para(run(txt, sz = 1800, col = col, face = MONO_FACE), lnSpc = 1.14)
}

## ------------------------------------------------------------------- build --
build <- function(variant, fig2_png, out_file, caption = NULL) {
  .next_id_reset <- NULL  # ids are unique per file; a fresh run is fine

  work <- file.path(tempdir(), paste0("poster-", variant))
  unlink(work, recursive = TRUE)
  dir.create(work, recursive = TRUE)
  utils::unzip(SRC, exdir = work)

  slide_path <- file.path(work, "ppt", "slides", "slide1.xml")
  xml <- paste(readLines(slide_path, warn = FALSE, encoding = "UTF-8"),
               collapse = "\n")
  s <- slide_split(xml)
  sh <- s$shapes

  expect <- function(i, txt) {
    got <- paste(shape_text(sh[[i]]), collapse = "")
    if (!grepl(txt, got, fixed = TRUE)) {
      stop(sprintf("v3 shape %d is '%s', expected '%s'", i, got, txt))
    }
    invisible(TRUE)
  }
  expect(10, "Most survival software")   # THE PROBLEM body
  expect(133, "A COMPLETE CALL")         # code panel label
  expect(103, "Same model, four bases")  # FIG 2 caption

  if (!is.null(caption)) sh[[103]] <- set_run_text(sh[[103]], 1, caption)

  extra <- character()
  addx <- function(...) extra <<- c(extra, ...)

  ## 1. ------------------------------------------------- THE PROBLEM, cont. --
  addx(sp_text(PROB_X, PROB_Y, PROB_W, 3.00, c(
    body_para(paste0(
      "Two workarounds are common. Imputing a single time ", EMD,
      " the midpoint, say ", EMD,
      " treats a guess as data: standard errors come out too small, and the ",
      "bias grows with the gap between visits."),
      first = TRUE),
    body_para(paste0(
      "MPL instead keeps every interval as recorded and fits the baseline ",
      "hazard jointly with ", BETA, ", so both come back with standard errors."))),
    name = "Problem continued"))

  ## 2. ---------------------------------------------- what the call prints --
  ## Verbatim from summary(fit) on the cached bcos2 msplines fit; the baseline
  ## parameter block is elided and the elision is marked as a comment.
  out_lines <- c(
    "Cox Proportional Hazards Model Fit Using MPL",
    "",
    "Penalized log-likelihood  :  -144.0955",
    "Estimated smoothing value :  2997128305",
    "Convergence               :  Yes (10 + 212463 iter.)",
    "",
    "Data             : bcos2",
    "Number of obs.   : 94",
    "Number of events :  0 (  0%)",
    "Number of cens.  : 94 (100%)",
    "",
    "Regression parameters : Surv(left, right, type = \"interval2\") ~ treatment",
    "                 Estimate Std. Error z-value Pr(>|z|)",
    "treatmentRadChem  0.89355    0.29414  3.0379 0.002383 **")

  addx(sp_rect(OUT_X, OUT_Y, OUT_W, OUT_H, fill = PAPER, line = RULE,
               name = "Output panel"))
  addx(sp_rect(OUT_X, OUT_Y, 0.125, OUT_H, fill = BLUE, name = "Output rule"))
  addx(sp_text(OUT_X + PAD, OUT_Y + 0.22, 12.0, 0.36,
               para(run("AND WHAT IT PRINTS", sz = 1950, col = BLUE,
                        face = MONO_FACE, spc = 273)),
               name = "Output label"))
  addx(sp_text(OUT_X + PAD, OUT_Y + 0.70, 13.2, 5.35,
               c(lapply(out_lines, console),
                 list(console(""),
                      console(paste0("# ", EMD,
                                     " then the 13 M-spline baseline parameters"),
                              col = MUTE))),
               name = "Output body"))
  addx(sp_rect(OUT_X + PAD, OUT_Y + 6.30, 12.6, 0.021, fill = RULE,
               name = "Output divider"))
  addx(sp_text(OUT_X + PAD, OUT_Y + 6.46, 12.9, 0.90,
               para(run(paste0(
                 "All 94 women are censored ", EMD,
                 " no exact event time anywhere ", EMD,
                 " yet the fit returns HR 2.44 (95% CI 1.37", NDASH,
                 "4.35) for added chemotherapy, plus a baseline hazard."),
                 sz = 2100, col = INK2), lnSpc = 1.15),
               name = "Output note"))

  ## ------------------------------------------------------------- assemble --
  s$shapes <- sh
  xml <- slide_join(s)
  xml <- sub("</p:spTree>", paste0(paste0(extra, collapse = ""), "</p:spTree>"),
             xml, fixed = TRUE)

  con <- file(slide_path, open = "wb")
  writeLines(enc2utf8(xml), con, useBytes = TRUE)
  close(con)

  ## 3. ------------------------------------------------------- swap FIG 2 ----
  ## rId3 -> ppt/media/image1.png is the FIG 2 picture, uncropped, and both
  ## replacements are drawn at the same 16.04 x 4.40 in, so only bytes change.
  stopifnot(file.exists(file.path(FIG, fig2_png)))
  file.copy(file.path(FIG, fig2_png),
            file.path(work, "ppt", "media", "image1.png"), overwrite = TRUE)

  ## PowerPoint keeps an exclusive lock on an open deck, so a rebuild while the
  ## previous output is on screen silently leaves the old file in place.  Write
  ## to a scratch name first, then move it into position, and say plainly when
  ## that is not possible rather than reporting the stale file's size as a win.
  tmp_zip <- file.path(tempdir(), paste0(variant, ".zip"))
  if (file.exists(tmp_zip)) unlink(tmp_zip)
  old <- setwd(work)
  files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
  rc <- utils::zip(tmp_zip, files, flags = "-q -X")
  setwd(old)
  if (rc != 0 || !file.exists(tmp_zip)) stop("zip failed for ", variant)

  ok <- suppressWarnings(file.copy(tmp_zip, out_file, overwrite = TRUE))
  if (!ok) {
    alt_dir <- file.path("dev", "poster", "out")
    dir.create(alt_dir, recursive = TRUE, showWarnings = FALSE)
    alt <- file.path(alt_dir, basename(out_file))
    if (!file.copy(tmp_zip, alt, overwrite = TRUE)) stop("could not write ", alt)
    warning(sprintf("%s is locked (open in PowerPoint?); wrote %s instead",
                    out_file, alt), call. = FALSE)
    out_file <- alt
  }
  message("wrote ", out_file, " (", round(file.size(out_file) / 1024), " KB)")
}

for (v in names(VARIANTS)) {
  build(v, VARIANTS[[v]]$png,
        sprintf("ICTMC-2026-survivalMPL-poster_%s.pptx", v),
        caption = VARIANTS[[v]]$caption)
}
