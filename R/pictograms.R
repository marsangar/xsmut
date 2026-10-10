# Species pictograms (silhouettes) for the plot panels.
#
# Silhouettes come from PhyloPic (https://www.phylopic.org) through the optional
# rphylopic package; images that need no attribution and have no ShareAlike
# clause are preferred. Credits for whatever is used are returned so they can be
# cited. Users can also supply their own image (PNG path, raster array or grob).
# Every failure (package missing, offline, no image) degrades silently to "no
# pictogram" so plotting never breaks.

.pictogram_cache <- new.env(parent = emptyenv())

.phylopic_lookup <- function(name) {
  key <- tolower(name)
  if (exists(key, envir = .pictogram_cache, inherits = FALSE)) return(get(key, envir = .pictogram_cache))
  res <- NULL
  # only look up proper binomials ("Macaca mulatta"), never free-text labels
  if (grepl("^[A-Z][a-z]+ [a-z]+", name) && requireNamespace("rphylopic", quietly = TRUE)) {
    res <- tryCatch({
      uuid <- tryCatch(rphylopic::get_uuid(name = name, n = 1, filter = c("by", "sa")),
                       error = function(e) NULL)
      if (is.null(uuid) || !length(uuid) || is.na(uuid[1]))
        uuid <- rphylopic::get_uuid(name = name, n = 1)
      img <- rphylopic::get_phylopic(uuid = uuid[1], format = "raster", height = 256)
      credit <- tryCatch(rphylopic::get_attribution(uuid = uuid[1], text = TRUE),
                         error = function(e) NA_character_)
      list(img = img, uuid = uuid[1], credit = paste(credit, collapse = " "))
    }, error = function(e) NULL)
  }
  assign(key, res, envir = .pictogram_cache)
  res
}

# x: NULL (look up `name` on PhyloPic), FALSE (none), a PNG path, a raster/array, or a grob
.as_pictogram_grob <- function(x = NULL, name = NULL, colour = "grey25", alpha = 0.9) {
  if (isFALSE(x)) return(NULL)
  credit <- NA_character_
  if (is.null(x)) {
    if (is.null(name)) return(NULL)
    hit <- .phylopic_lookup(name)
    if (is.null(hit)) return(NULL)
    x <- hit$img; credit <- hit$credit
    x <- tryCatch(rphylopic::recolor_phylopic(x, alpha = alpha, fill = colour), error = function(e) x)
  }
  if (is.character(x) && length(x) == 1 && file.exists(x)) {
    if (!requireNamespace("png", quietly = TRUE)) return(NULL)
    x <- png::readPNG(x)
  }
  g <- if (inherits(x, "grob")) x else
    tryCatch(grid::rasterGrob(x, x = grid::unit(1, "npc"), hjust = 1, interpolate = TRUE),
             error = function(e) NULL)
  if (is.null(g)) return(NULL)
  attr(g, "credit") <- credit
  g
}

# Place a pictogram in the outer corner of a panel (beside its title strip).
.add_pictogram <- function(p, grob, flip = FALSE, size = 0.14) {
  if (is.null(grob)) return(p)
  if (flip) {
    p + patchwork::inset_element(grob, left = 1 - size * 2.5, bottom = 0, right = 1, top = size,
                                 align_to = "full", clip = FALSE, on_top = TRUE)
  } else {
    p + patchwork::inset_element(grob, left = 1 - size * 2.5, bottom = 1 - size, right = 1, top = 1,
                                 align_to = "full", clip = FALSE, on_top = TRUE)
  }
}
