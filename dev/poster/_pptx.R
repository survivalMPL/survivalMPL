## Minimal DrawingML writers, enough to build the shapes this poster needs.
## Everything here emits the same flavour of markup PowerPoint itself wrote for
## the hand-made shapes already on the slide, so edited and generated shapes are
## indistinguishable once the file is reopened.

EMU_PER_IN <- 914400
emu <- function(inch) sprintf("%.0f", inch * EMU_PER_IN)

## Shape ids only have to be unique within the slide; the hand-made shapes stop
## well below 900.
.next_id <- local({
  i <- 900
  function() {
    i <<- i + 1
    i
  }
})

xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

## Typeface block, repeated for latin / east-asian / complex-script as in the
## original file - PowerPoint drops the font otherwise on some platforms.
.font <- function(face) {
  if (is.null(face)) return("")
  sprintf(paste0('<a:latin typeface="%s" pitchFamily="34" charset="0"/>',
                 '<a:ea typeface="%s" pitchFamily="34" charset="-122"/>',
                 '<a:cs typeface="%s" pitchFamily="34" charset="-120"/>'),
          face, face, face)
}

SANS_FACE <- "IBM Plex Sans"
MONO_FACE <- "IBM Plex Mono"
DISP_FACE <- "Archivo"

## One run of text.  `spc` is letter-spacing in 1/100 pt, as the small-caps
## section labels on the poster use it.
run <- function(text, sz = 2550, bold = FALSE, col = "121417",
                face = SANS_FACE, spc = NULL, italic = FALSE) {
  sprintf('<a:r><a:rPr lang="en-GB" sz="%d"%s%s kern="0"%s dirty="0"><a:solidFill><a:srgbClr val="%s"/></a:solidFill>%s</a:rPr><a:t>%s</a:t></a:r>',
          sz,
          if (bold) ' b="1"' else "",
          if (italic) ' i="1"' else "",
          if (is.null(spc)) "" else sprintf(' spc="%d"', spc),
          col, .font(face), xml_escape(text))
}

## One paragraph.  `lnSpc` is a line-height multiplier.
para <- function(..., algn = "l", lnSpc = NULL, space_before = NULL) {
  runs <- c(...)
  sprintf('<a:p><a:pPr marL="0" indent="0" algn="%s">%s%s<a:buNone/></a:pPr>%s</a:p>',
          algn,
          if (is.null(lnSpc)) "" else
            sprintf('<a:lnSpc><a:spcPct val="%d"/></a:lnSpc>', round(lnSpc * 100000)),
          if (is.null(space_before)) "" else
            sprintf('<a:spcBef><a:spcPts val="%d"/></a:spcBef>', round(space_before * 100)),
          paste0(runs, collapse = ""))
}

.xfrm <- function(x, y, w, h) {
  sprintf('<a:xfrm><a:off x="%s" y="%s"/><a:ext cx="%s" cy="%s"/></a:xfrm>',
          emu(x), emu(y), emu(w), emu(h))
}

## A text box.  Insets match the ones the poster already uses.
sp_text <- function(x, y, w, h, paras, name = "Text", anchor = "t",
                    inset = 0.0278) {
  ins <- emu(inset)
  sprintf(paste0('<p:sp><p:nvSpPr><p:cNvPr id="%d" name="%s"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr>',
                 '<p:spPr>%s<a:prstGeom prst="rect"><a:avLst/></a:prstGeom><a:noFill/><a:ln/></p:spPr>',
                 '<p:txBody><a:bodyPr wrap="square" lIns="%s" tIns="%s" rIns="%s" bIns="%s" rtlCol="0" anchor="%s"><a:normAutofit/></a:bodyPr><a:lstStyle/>%s</p:txBody></p:sp>'),
          .next_id(), name, .xfrm(x, y, w, h), ins, ins, ins, ins, anchor,
          paste0(paras, collapse = ""))
}

## A filled rectangle, optionally outlined.  Also used for the rules and the
## bars in the observation-schemes diagram.
sp_rect <- function(x, y, w, h, fill = NULL, line = NULL, lw = 0.0104,
                    name = "Shape", geom = "rect", radius = NULL) {
  g <- if (is.null(radius)) {
    sprintf('<a:prstGeom prst="%s"><a:avLst/></a:prstGeom>', geom)
  } else {
    sprintf('<a:prstGeom prst="roundRect"><a:avLst><a:gd name="adj" fmla="val %d"/></a:avLst></a:prstGeom>',
            round(radius * 100000))
  }
  sprintf(paste0('<p:sp><p:nvSpPr><p:cNvPr id="%d" name="%s"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr>',
                 '<p:spPr>%s%s%s%s</p:spPr>',
                 '<p:txBody><a:bodyPr/><a:lstStyle/><a:p><a:endParaRPr lang="en-GB"/></a:p></p:txBody></p:sp>'),
          .next_id(), name, .xfrm(x, y, w, h), g,
          if (is.null(fill)) "<a:noFill/>" else
            sprintf('<a:solidFill><a:srgbClr val="%s"/></a:solidFill>', fill),
          if (is.null(line)) "<a:ln/>" else
            sprintf('<a:ln w="%s"><a:solidFill><a:srgbClr val="%s"/></a:solidFill></a:ln>',
                    emu(lw), line))
}

## A picture.  `rid` is the relationship id added to slide1.xml.rels.
sp_pic <- function(x, y, w, h, rid, name = "Picture") {
  sprintf(paste0('<p:pic><p:nvPicPr><p:cNvPr id="%d" name="%s"/><p:cNvPicPr>',
                 '<a:picLocks noChangeAspect="1"/></p:cNvPicPr><p:nvPr/></p:nvPicPr>',
                 '<p:blipFill><a:blip r:embed="%s"/><a:stretch><a:fillRect/></a:stretch></p:blipFill>',
                 '<p:spPr>%s<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr></p:pic>'),
          .next_id(), name, rid, .xfrm(x, y, w, h))
}

## ---------------------------------------------------------------------------
## Slide surgery.  The slide is a flat list of shapes inside <p:spTree>, so it
## can be split into head / shapes / tail, edited by index, and glued back.

slide_split <- function(xml) {
  pat <- "<p:(sp|pic|graphicFrame|cxnSp)>.*?</p:\\1>"
  m <- gregexpr(pat, xml, perl = TRUE)[[1]]
  stopifnot(m[1] > 0)
  starts <- as.integer(m)
  ends <- starts + attr(m, "match.length") - 1L
  shapes <- substring(xml, starts, ends)
  gaps <- character(length(shapes))
  for (i in seq_along(shapes)) {
    gaps[i] <- if (i == 1) "" else substring(xml, ends[i - 1] + 1, starts[i] - 1)
  }
  out <- list(head = substring(xml, 1, starts[1] - 1),
              shapes = shapes, gaps = gaps,
              tail = substring(xml, ends[length(ends)] + 1, nchar(xml)))
  stopifnot(identical(slide_join(out), xml))
  out
}

slide_join <- function(s) {
  paste0(s$head, paste0(s$gaps, s$shapes, collapse = ""), s$tail)
}

## Replace the text of the n-th <a:t> inside one shape.
set_run_text <- function(shape, n, text) {
  m <- gregexpr("<a:t>.*?</a:t>", shape, perl = TRUE)[[1]]
  stopifnot(length(m) >= n, m[1] > 0)
  st <- as.integer(m)[n]
  en <- st + attr(m, "match.length")[n] - 1L
  paste0(substring(shape, 1, st - 1),
         "<a:t>", xml_escape(text), "</a:t>",
         substring(shape, en + 1))
}

shape_text <- function(shape) {
  m <- regmatches(shape, gregexpr("<a:t>.*?</a:t>", shape, perl = TRUE))[[1]]
  gsub("</?a:t>", "", m)
}

## Move / resize a shape in place.
set_xfrm <- function(shape, x = NULL, y = NULL, w = NULL, h = NULL) {
  off <- regmatches(shape, regexpr('<a:off x="-?\\d+" y="-?\\d+"/>', shape, perl = TRUE))
  ext <- regmatches(shape, regexpr('<a:ext cx="\\d+" cy="\\d+"/>', shape, perl = TRUE))
  cur <- as.numeric(regmatches(off, gregexpr("-?\\d+", off))[[1]])
  cue <- as.numeric(regmatches(ext, gregexpr("\\d+", ext))[[1]])
  nx <- if (is.null(x)) cur[1] else as.numeric(emu(x))
  ny <- if (is.null(y)) cur[2] else as.numeric(emu(y))
  nw <- if (is.null(w)) cue[1] else as.numeric(emu(w))
  nh <- if (is.null(h)) cue[2] else as.numeric(emu(h))
  shape <- sub('<a:off x="-?\\d+" y="-?\\d+"/>',
               sprintf('<a:off x="%.0f" y="%.0f"/>', nx, ny), shape, perl = TRUE)
  sub('<a:ext cx="\\d+" cy="\\d+"/>',
      sprintf('<a:ext cx="%.0f" cy="%.0f"/>', nw, nh), shape, perl = TRUE)
}

## Insert a new run into a shape's paragraph, copying the formatting of an
## existing run.  Used where a rewritten sentence needs one more code-styled
## dataset name than the original had.
insert_run <- function(shape, after, template, text) {
  m <- gregexpr("<a:r>.*?</a:r>", shape, perl = TRUE)[[1]]
  stopifnot(m[1] > 0, length(m) >= max(after, template))
  st <- as.integer(m); en <- st + attr(m, "match.length") - 1L
  tpl <- substring(shape, st[template], en[template])
  tpl <- sub("<a:t>.*?</a:t>", paste0("<a:t>", xml_escape(text), "</a:t>"),
             tpl, perl = TRUE)
  paste0(substring(shape, 1, en[after]), tpl, substring(shape, en[after] + 1))
}
