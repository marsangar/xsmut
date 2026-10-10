#' @import ggplot2
NULL

# Piecewise-linear transform that keeps exons at native scale and squeezes each
# intron to a fixed width (fraction of total exon length). Returns a function.
.make_transform <- function(exons, compress_introns = TRUE, intron_frac = 0.04) {
  ex <- as.data.table(exons)[, .(start = min(start), end = max(end)), by = exon_rank][order(start)]
  if (!compress_introns || nrow(ex) < 2) return(list(f = identity, exons = ex))
  intron_w <- intron_frac * sum(ex$end - ex$start + 1)
  # cumulative plot offsets
  ex[, plot_start := 0]
  for (i in seq_len(nrow(ex))) {
    ex$plot_start[i] <- if (i == 1) ex$start[1] else ex$plot_start[i - 1] + (ex$end[i - 1] - ex$start[i - 1] + 1) + intron_w
  }
  ex[, plot_end := plot_start + (end - start)]
  f <- function(x) {
    out <- rep(NA_real_, length(x))
    for (i in seq_len(nrow(ex))) {
      in_ex <- !is.na(x) & x >= ex$start[i] & x <= ex$end[i]
      out[in_ex] <- ex$plot_start[i] + (x[in_ex] - ex$start[i])
      if (i < nrow(ex)) {
        in_int <- !is.na(x) & x > ex$end[i] & x < ex$start[i + 1]
        frac <- (x[in_int] - ex$end[i]) / (ex$start[i + 1] - ex$end[i])
        out[in_int] <- ex$plot_end[i] + frac * intron_w
      }
    }
    # outside gene: clamp linearly beyond the ends
    lo <- !is.na(x) & x < ex$start[1]; out[lo] <- ex$plot_start[1] - (ex$start[1] - x[lo])
    hi <- !is.na(x) & x > ex$end[nrow(ex)]; out[hi] <- ex$plot_end[nrow(ex)] + (x[hi] - ex$end[nrow(ex)])
    out
  }
  list(f = f, exons = ex)
}

.consequence_class <- function(x) {
  x <- tolower(as.character(x))
  data.table::fcase(
    grepl("missense", x), "Missense",
    grepl("nonsense|stop_gained|stop gained", x), "Nonsense",
    grepl("frameshift", x), "Frameshift",
    grepl("inframe|in frame|deletion|insertion", x), "In-frame indel",
    grepl("splice", x), "Splice",
    grepl("synonymous|silent", x), "Synonymous",
    is.na(x) | x == "", "Unknown",
    default = "Other"
  )
}

.pal <- c(Missense = "#2E86AB", Nonsense = "#C73E1D", Frameshift = "#F18F01",
          "In-frame indel" = "#7B2D8E", Splice = "#3B8B5A", Synonymous = "#9E9E9E",
          Other = "#5E5E5E", Unknown = "#BDBDBD")

# named presets for the hotspot halo; any R colour is also accepted
.hl_presets <- c(gold = "#E0A800", red = "#E8000B")

.resolve_colour <- function(x) {
  if (is.null(x) || !nzchar(x)) return(.hl_presets[["gold"]])
  if (tolower(x) %in% names(.hl_presets)) .hl_presets[[tolower(x)]] else x
}

# darker shade of a colour, for text drawn in the highlight colour
.darken <- function(col, f = 0.6) {
  rgb <- grDevices::col2rgb(col)[, 1] * f
  grDevices::rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
}

# "macaca_mulatta" -> "Macaca mulatta" (genus capitalised, epithet(s) lower case)
.species_display_name <- function(x) {
  s <- gsub("_", " ", tolower(trimws(as.character(x))))
  paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}

# integer-only axis breaks, without the 0 baseline (it touches the gene model)
.int_breaks <- function(lims) {
  lims <- lims[is.finite(lims)]
  if (!length(lims)) return(numeric())
  b <- pretty(lims, n = 4)
  b <- b[abs(b - round(b)) < 1e-8 & b >= min(lims) - 1e-8 & b <= max(lims) + 1e-8]
  if (length(b) < 2) b <- seq(ceiling(min(lims)), floor(max(lims)))
  b[b != 0]
}

.fmt_n <- function(x) format(x, big.mark = ",", trim = TRUE)

# Small one-paragraph COSMIC summary for a gene from cosmic_gene_listed()
.gene_info_text <- function(info, max_tumour_types = 8, width = 150) {
  if (is.null(info) || !nrow(info)) return(NULL)
  info <- as.data.table(info)[1]
  get <- function(col) if (col %in% names(info)) as.character(info[[col]]) else NA_character_
  ok <- function(v) !is.na(v) && nzchar(v)
  parts <- character()
  role <- get("role")
  census_used <- isTRUE(as.logical(get("census_supplied")))
  in_cgc <- isTRUE(as.logical(get("in_census")))
  # CGC tier is deliberately not shown (saves header space); flag only absence
  if (census_used && !in_cgc) parts <- c(parts, "Not in COSMIC Cancer Gene Census")
  if (ok(role)) parts <- c(parts, paste0("Role: ", role))
  tt <- get("tumour_types_somatic")
  if (ok(tt)) {
    tts <- trimws(strsplit(tt, ",")[[1]])
    tts <- tts[nzchar(tts)]
    if (length(tts) > max_tumour_types) tts <- c(tts[seq_len(max_tumour_types)], "...")
    parts <- c(parts, paste0("Somatic tumour types: ", paste(tts, collapse = ", ")))
  }
  if (!length(parts)) return(NULL)
  paste(strwrap(paste(parts, collapse = "  \u00b7  "), width = width), collapse = "\n")
}

.lollipop <- function(dt, tf, title, subtitle = NULL, flip = FALSE, min_label_n = 3,
                      hl_x = NULL, hl_col = "#E0A800", base_size = 11,
                      title_face = "bold") {
  # dt: pos, aa_change, consequence, n, [hl (logical)]
  dt <- as.data.table(dt)
  small <- base_size * 0.75
  # panel titles: top panel -> above; flipped (bottom) panel -> at its base
  title_layer <- if (flip) {
    list(labs(title = NULL, subtitle = NULL, x = title, caption = subtitle),
         theme(axis.title.x = element_text(face = title_face, size = base_size, hjust = 0,
                                           margin = margin(t = 3)),
               plot.caption = element_text(size = small, colour = "grey35", hjust = 0,
                                           margin = margin(t = 1)),
               plot.caption.position = "panel",
               plot.margin = margin(0, 5.5, 5.5, 5.5)))
  } else {
    list(labs(title = title, subtitle = subtitle, x = NULL),
         theme(plot.title = element_text(face = title_face, size = base_size, margin = margin(b = 1)),
               plot.subtitle = element_text(size = small, colour = "grey35", margin = margin(b = 3)),
               plot.margin = margin(5.5, 5.5, 0, 5.5)))
  }
  if (!nrow(dt)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "no mappable mutations",
                               colour = "grey50", size = small / .pt) +
             theme_void(base_size = base_size) + title_layer)
  }
  if (!"hl" %in% names(dt)) dt[, hl := FALSE]
  dt[is.na(hl), hl := FALSE]
  dt[, class := .consequence_class(consequence)]
  dt[, x := tf(pos)]
  dt[, ymax := if (flip) -n else n]
  lab    <- dt[n >= min_label_n & !hl]
  lab_hl <- dt[hl == TRUE]

  p <- ggplot(dt)
  if (length(hl_x)) {
    p <- p + geom_vline(xintercept = tf(hl_x), colour = hl_col, alpha = 0.5,
                        linewidth = 0.35, linetype = "dashed")
  }
  p <- p +
    geom_segment(aes(x = x, xend = x, y = 0, yend = ymax), colour = "grey55", linewidth = 0.3)
  if (nrow(lab_hl)) {
    # halo drawn underneath the point: thick border of the same mapped size
    p <- p + geom_point(data = lab_hl, aes(x = x, y = ymax, size = n), shape = 21,
                        fill = NA, colour = hl_col, stroke = 2.4)
  }
  p <- p +
    geom_point(aes(x = x, y = ymax, colour = class, size = n), alpha = 0.9) +
    scale_colour_manual(values = .pal, drop = TRUE, name = "Consequence") +
    scale_size_area(max_size = 5, guide = "none") +
    scale_y_continuous(breaks = .int_breaks, labels = function(v) .fmt_n(abs(v)),
                       expand = if (flip) expansion(mult = c(0.22, 0)) else expansion(mult = c(0, 0.22))) +
    labs(y = "Samples") +
    theme_classic(base_size = base_size) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
          axis.title.y = element_text(size = small, colour = "grey25"),
          axis.text.y = element_text(size = small, colour = "grey25"),
          axis.line.y = element_line(linewidth = 0.3, colour = "grey40"),
          axis.ticks.y = element_line(linewidth = 0.3, colour = "grey40"),
          legend.position = "right",
          legend.title = element_text(size = small, face = "bold"),
          legend.text = element_text(size = small)) +
    title_layer
  vj <- if (flip) 1.7 else -0.9
  lab_size <- base_size * 0.62 / .pt
  if (nrow(lab)) {
    p <- p + geom_text(data = lab, aes(x = x, y = ymax, label = aa_change),
                       vjust = vj, size = lab_size, colour = "grey15", check_overlap = TRUE)
  }
  if (nrow(lab_hl)) {
    p <- p + geom_text(data = lab_hl, aes(x = x, y = ymax, label = aa_change),
                       vjust = vj, size = lab_size * 1.05, fontface = "bold", colour = .darken(hl_col, 0.8))
  }
  p
}

.gene_track <- function(model, tf, show_exon_numbers = TRUE, intron_linewidth = 1.1,
                        hl_x = NULL, hl_col = "#E0A800", base_size = 11) {
  ex <- as.data.table(model$exons)
  ex[, `:=`(px = tf(start), pxend = tf(end))]
  ex[, h := ifelse(feature == "CDS", 0.42, 0.22)]
  merged <- ex[, .(px = min(px), pxend = max(pxend), start = min(start), end = max(end)), by = exon_rank][order(px)]
  strand_lab <- if (model$gene$strand == 1) "forward (+) strand" else "reverse (-) strand"
  small <- base_size * 0.62
  locus <- sprintf("%s  \u00b7  chr%s:%s-%s  \u00b7  %s  \u00b7  %s  \u00b7  %s",
                   model$gene$symbol, model$gene$chrom, .fmt_n(model$gene$start),
                   .fmt_n(model$gene$end), model$transcript$id, model$gene$assembly, strand_lab)
  p <- ggplot()
  if (length(hl_x)) {
    p <- p + geom_vline(xintercept = tf(hl_x), colour = hl_col, alpha = 0.5,
                        linewidth = 0.35, linetype = "dashed")
  }
  p <- p +
    annotate("segment", x = min(merged$px), xend = max(merged$pxend), y = 0, yend = 0,
             colour = "grey35", linewidth = intron_linewidth, lineend = "butt") +
    geom_rect(data = ex, aes(xmin = px, xmax = pxend, ymin = -h, ymax = h, fill = feature), colour = NA) +
    scale_fill_manual(values = c(CDS = "#1F3A5F", UTR = "#9DB4D0"), name = NULL) +
    annotate("text", x = -Inf, y = 0.78, label = locus, hjust = -0.01, vjust = 1,
             size = small / .pt, colour = "grey30") +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = base_size) +
    theme(legend.position = "right", legend.text = element_text(size = base_size * 0.75),
          plot.margin = margin(0, 5.5, 0, 5.5))
  if (show_exon_numbers) {
    p <- p + geom_text(data = merged, aes(x = (px + pxend) / 2, y = -0.66, label = exon_rank),
                       size = small / .pt, colour = "grey40")
  }
  p
}

#' Plot user (species) mutations and human COSMIC mutations on a human gene model
#'
#' Publication-oriented layout: mutation panels sit flush against the gene
#' model, the lower panel carries its title at its base, sample axes use
#' integer breaks, and summary counts are shown in small text under each title.
#'
#' @param model `gene_model` from [get_gene_model()] (human)
#' @param cosmic data.table from [cosmic_gene_mutations()] (needs `start`,
#'   `aa_change`, `consequence`)
#' @param species_mut data.table from [map_mutations_to_human()] (needs
#'   `human_pos`; optional `aa_change`, `consequence`, `sample_id`)
#' @param species_label label for the species panel, e.g. "Macaca mulatta (Mmul_10)".
#'   A trailing "(assembly)" is moved to the small-print line.
#' @param species_panel where the user mutations go: "top" (default; COSMIC is
#'   then drawn below the gene model) or "bottom" (COSMIC on top)
#' @param gene_info optional one-row table from [cosmic_gene_listed()]; its
#'   role (oncogene / TSG / fusion) and somatic tumour types are shown in the
#'   figure header
#' @param compress_introns squeeze introns to a fixed width so exons dominate
#' @param coding_only drop COSMIC rows without a protein change
#' @param min_label_n label COSMIC mutations recurring in at least this many samples
#' @param highlight_hotspots highlight user mutations at/near the top COSMIC
#'   hotspot residues (see [match_cosmic_hotspots()])
#' @param hotspot_top_n,hotspot_window,hotspot_min_samples passed to
#'   [match_cosmic_hotspots()]
#' @param hotspot_matches optional precomputed result of [match_cosmic_hotspots()]
#' @param highlight_colour halo colour for matches: "gold" (default), "red"
#'   (more visible on slides) or any R colour
#' @param intron_linewidth line width of the intron backbone in the gene track
#' @param base_size base font size (pt); small print is 75% of it
#' @param italic_species draw the species name in italics
#' @param max_tumour_types truncate the COSMIC somatic tumour-type list
#' @param panel_heights relative heights c(upper, gene model, lower); default
#'   gives COSMIC more room than the species panel
#' @param legend_position "right", "bottom" or "none"
#' @param show_pictograms draw a species silhouette in the corner of each
#'   mutation panel (PhyloPic via the optional rphylopic package; needs internet
#'   on first use per session, silently skipped otherwise)
#' @param pictograms optional overrides, `list(species = ..., human = ...)`;
#'   each a PNG path, raster array, grob, or `FALSE` to hide that one
#' @param pictogram_species scientific name used to look up the species
#'   silhouette; defaults to the name in `species_label`
#' @param pictogram_size height of the pictogram as a fraction of its panel
#' @return a patchwork object (ggplot); use `print()` or `ggsave()`. If
#'   highlighting is on, the match table is in `attr(p, "hotspot_matches")`.
#' @export
plot_gene_mutations <- function(model, cosmic, species_mut, species_label = "Species",
                                species_panel = c("top", "bottom"), gene_info = NULL,
                                compress_introns = TRUE, coding_only = TRUE, min_label_n = 3,
                                highlight_hotspots = FALSE, hotspot_top_n = 10,
                                hotspot_window = 3, hotspot_min_samples = 2,
                                hotspot_matches = NULL, highlight_colour = "gold",
                                intron_linewidth = 1.1, base_size = 11, italic_species = TRUE,
                                max_tumour_types = 8, panel_heights = NULL,
                                legend_position = "right", show_pictograms = TRUE,
                                pictograms = list(), pictogram_species = NULL,
                                pictogram_size = 0.14) {
  species_panel <- match.arg(species_panel)
  hl_col <- .resolve_colour(highlight_colour)
  small <- base_size * 0.75
  tr <- .make_transform(model$exons, compress_introns)
  tf <- tr$f

  cs <- as.data.table(cosmic)[!is.na(start)]
  if (coding_only) cs <- cs[!is.na(aa_pos)]
  cs_agg <- cs[, .(n = .N), by = .(pos = start, aa_change, consequence)]

  sp_all <- as.data.table(species_mut)
  sp <- sp_all[!is.na(human_pos)]
  if (!"aa_change" %in% names(sp)) sp[, aa_change := paste0("p.", aa_pos)]
  if (!"consequence" %in% names(sp)) sp[, consequence := NA_character_]
  sp_agg <- sp[, .(n = .N), by = .(pos = human_pos, aa_change, consequence)]

  # --- hotspot matching ------------------------------------------------------
  hl_x <- NULL; hs_caption <- NULL
  if (isTRUE(highlight_hotspots)) {
    if (is.null(hotspot_matches)) {
      hotspot_matches <- match_cosmic_hotspots(species_mut, cosmic, top_n = hotspot_top_n,
                                               window = hotspot_window,
                                               min_samples = hotspot_min_samples)
    }
    m <- as.data.table(hotspot_matches)
    if (nrow(m)) {
      sp_agg[, hl := paste(pos, aa_change) %in% paste(m$human_pos, m$aa_change)]
      hs <- attr(hotspot_matches, "hotspots")
      hit_hs <- if (!is.null(hs)) hs[aa_pos %in% m$hotspot_aa_pos] else NULL
      cs_agg[, hl := !is.null(hit_hs) & paste(pos, aa_change) %in% paste(hit_hs$start, hit_hs$top_change)]
      hl_x <- unique(m$human_pos)
    }
    hs_caption <- sprintf(
      "Highlighted: %d mutation(s) within \u00b1%d aa of the top %d COSMIC hotspot residues (min. %d samples each).",
      nrow(m), hotspot_window, hotspot_top_n, hotspot_min_samples)
  }

  # --- titles and small-print summaries --------------------------------------
  n_samples <- function(d) if ("sample_id" %in% names(d)) data.table::uniqueN(d$sample_id[!is.na(d$sample_id)]) else NA_integer_
  sp_name <- species_label; sp_asm <- NULL
  if (grepl("^.+ \\(.+\\)$", species_label)) {
    sp_name <- sub("^(.+) \\((.+)\\)$", "\\1", species_label)
    sp_asm  <- sub("^(.+) \\((.+)\\)$", "\\2", species_label)
  }
  n_unmapped <- sum(is.na(sp_all$human_pos))
  ns_sp <- n_samples(sp_all)
  sp_sub <- paste(c(
    sp_asm,
    sprintf("%s mutations mapped to human%s", .fmt_n(nrow(sp)),
            if (n_unmapped) sprintf(" (%s unmapped)", .fmt_n(n_unmapped)) else ""),
    if (!is.na(ns_sp)) sprintf("%s samples", .fmt_n(ns_sp))
  ), collapse = "  \u00b7  ")
  ns_cs <- n_samples(cs)
  cs_sub <- paste(c(
    sprintf("%s distinct mutations", .fmt_n(nrow(cs_agg))),
    sprintf("%s mutation calls", .fmt_n(nrow(cs))),
    if (!is.na(ns_cs)) sprintf("%s samples", .fmt_n(ns_cs))
  ), collapse = "  \u00b7  ")
  cs_title <- sprintf("Human COSMIC (%s)", model$gene$assembly)

  # --- panels ----------------------------------------------------------------
  sp_top <- species_panel == "top"
  p_sp <- .lollipop(sp_agg, tf, sp_name, sp_sub, flip = !sp_top, min_label_n = 1,
                    hl_x = hl_x, hl_col = hl_col, base_size = base_size,
                    title_face = if (italic_species) "bold.italic" else "bold")
  p_cs <- .lollipop(cs_agg, tf, cs_title, cs_sub, flip = sp_top, min_label_n = min_label_n,
                    hl_x = hl_x, hl_col = hl_col, base_size = base_size)
  mid <- .gene_track(model, tf, intron_linewidth = intron_linewidth,
                     hl_x = hl_x, hl_col = hl_col, base_size = base_size)

  xr <- range(c(tf(model$exons$start), tf(model$exons$end), tf(cs_agg$pos), tf(sp_agg$pos)), na.rm = TRUE)
  pad <- diff(xr) * 0.02
  lock <- function(p) p + coord_cartesian(xlim = xr + c(-pad, pad), clip = "off")
  p_sp <- lock(p_sp); p_cs <- lock(p_cs)

  # --- pictograms (species silhouette on its panel, human on COSMIC) ----------
  credits <- character()
  if (isTRUE(show_pictograms)) {
    g_sp <- .as_pictogram_grob(pictograms$species,
                               name = if (is.null(pictogram_species)) sp_name else pictogram_species)
    g_hs <- .as_pictogram_grob(pictograms$human, name = "Homo sapiens")
    for (g in list(g_sp, g_hs)) {
      cr <- attr(g, "credit")
      if (!is.null(cr) && !is.na(cr) && nzchar(cr)) credits <- c(credits, cr)
    }
    p_sp <- .add_pictogram(p_sp, g_sp, flip = !sp_top, size = pictogram_size)
    p_cs <- .add_pictogram(p_cs, g_hs, flip = sp_top, size = pictogram_size)
  }
  mid <- mid + coord_cartesian(xlim = xr + c(-pad, pad), ylim = c(-0.8, 0.8), clip = "off")

  panels <- if (sp_top) list(p_sp, mid, p_cs) else list(p_cs, mid, p_sp)
  if (is.null(panel_heights)) panel_heights <- if (sp_top) c(3, 0.8, 4) else c(4, 0.8, 3)

  # --- header: gene, COSMIC role / tier / tumour types ------------------------
  gname <- if (!is.null(gene_info) && "name" %in% names(gene_info) && !is.na(gene_info$name[1]) &&
               nzchar(gene_info$name[1])) sprintf("%s: %s", model$gene$symbol, gene_info$name[1]) else model$gene$symbol
  info_txt <- .gene_info_text(gene_info, max_tumour_types = max_tumour_types)

  out <- patchwork::wrap_plots(panels, ncol = 1, heights = panel_heights) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = gname, subtitle = info_txt, caption = hs_caption,
      theme = theme(plot.title = element_text(face = "bold", size = base_size * 1.25),
                    plot.subtitle = element_text(size = small, colour = "grey30", lineheight = 1.1),
                    plot.caption = element_text(size = small, colour = .darken(hl_col, 0.75), hjust = 0))) &
    theme(legend.position = legend_position)
  attr(out, "hotspot_matches") <- hotspot_matches
  attr(out, "pictogram_credits") <- unique(credits)
  out
}
