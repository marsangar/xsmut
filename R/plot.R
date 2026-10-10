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

# "macaca_mulatta" -> "Macaca mulatta" (genus capitalised, epithet(s) lower case)
.species_display_name <- function(x) {
  s <- gsub("_", " ", tolower(trimws(as.character(x))))
  paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}

.lollipop <- function(dt, tf, title, flip = FALSE, min_label_n = 3,
                      hl_x = NULL, hl_col = "#E0A800") {
  # dt: pos, aa_change, consequence, n, [hl (logical)]
  dt <- as.data.table(dt)
  if (!nrow(dt)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "no mappable mutations", colour = "grey50") +
             theme_void() + labs(title = title))
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
    p <- p + geom_vline(xintercept = tf(hl_x), colour = hl_col, alpha = 0.45,
                        linewidth = 0.4, linetype = "dashed")
  }
  p <- p +
    geom_segment(aes(x = x, xend = x, y = 0, yend = ymax), colour = "grey60", linewidth = 0.3)
  if (nrow(lab_hl)) {
    # halo drawn underneath the point: thick border of the same mapped size
    p <- p + geom_point(data = lab_hl, aes(x = x, y = ymax, size = n), shape = 21,
                        fill = NA, colour = hl_col, stroke = 2.4)
  }
  p <- p +
    geom_point(aes(x = x, y = ymax, colour = class, size = n), alpha = 0.85) +
    scale_colour_manual(values = .pal, drop = TRUE, name = NULL) +
    scale_size_area(max_size = 6, guide = "none") +
    labs(title = title, y = "samples", x = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "right", plot.title = element_text(face = "bold"))
  vj <- if (flip) 1.6 else -0.8
  if (nrow(lab)) {
    p <- p + geom_text(data = lab, aes(x = x, y = ymax, label = aa_change),
                       vjust = vj, size = 2.7, check_overlap = TRUE)
  }
  if (nrow(lab_hl)) {
    p <- p + geom_text(data = lab_hl, aes(x = x, y = ymax, label = aa_change),
                       vjust = vj, size = 2.9, fontface = "bold", colour = "#7A5600")
  }
  if (flip) p <- p + scale_y_continuous(labels = function(v) abs(v))
  p
}

.gene_track <- function(model, tf, show_exon_numbers = TRUE, intron_linewidth = 1.1,
                        hl_x = NULL, hl_col = "#E0A800") {
  ex <- as.data.table(model$exons)
  ex[, `:=`(px = tf(start), pxend = tf(end))]
  ex[, h := ifelse(feature == "CDS", 0.5, 0.25)]
  merged <- ex[, .(px = min(px), pxend = max(pxend), start = min(start), end = max(end)), by = exon_rank][order(px)]
  strand_lab <- if (model$gene$strand == 1) "5' → 3' (+)" else "3' ← 5' (−)"
  p <- ggplot()
  if (length(hl_x)) {
    p <- p + geom_vline(xintercept = tf(hl_x), colour = hl_col, alpha = 0.45,
                        linewidth = 0.4, linetype = "dashed")
  }
  p <- p +
    annotate("segment", x = min(merged$px), xend = max(merged$pxend), y = 0, yend = 0,
             colour = "grey40", linewidth = intron_linewidth, lineend = "butt") +
    geom_rect(data = ex, aes(xmin = px, xmax = pxend, ymin = -h, ymax = h, fill = feature), colour = NA) +
    scale_fill_manual(values = c(CDS = "#1F3A5F", UTR = "#9DB4D0"), name = NULL) +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(x = sprintf("%s  chr%s:%s-%s (%s, %s)  %s", model$gene$symbol, model$gene$chrom,
                     format(model$gene$start, big.mark = ","), format(model$gene$end, big.mark = ","),
                     model$transcript$id, model$gene$assembly, strand_lab), y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
          legend.position = "right")
  if (show_exon_numbers) {
    p <- p + geom_text(data = merged, aes(x = (px + pxend) / 2, y = -0.75, label = exon_rank), size = 2.6, colour = "grey30")
  }
  p
}

#' Plot user (species) mutations and human COSMIC mutations on a human gene model
#'
#' @param model `gene_model` from [get_gene_model()] (human)
#' @param cosmic data.table from [cosmic_gene_mutations()] (needs `start`,
#'   `aa_change`, `consequence`)
#' @param species_mut data.table from [map_mutations_to_human()] (needs
#'   `human_pos`; optional `aa_change`, `consequence`)
#' @param species_label label for the species panel, e.g. "Macaca mulatta (Mmul_10)"
#' @param species_panel where the user mutations go: "top" (default; COSMIC is
#'   then drawn below the gene model) or "bottom" (COSMIC on top)
#' @param compress_introns squeeze introns to a fixed width so exons dominate
#' @param coding_only drop COSMIC rows without a protein change
#' @param min_label_n label COSMIC mutations recurring in at least this many samples
#' @param highlight_hotspots highlight user mutations at/near the top COSMIC
#'   hotspot residues (see [match_cosmic_hotspots()])
#' @param hotspot_top_n,hotspot_window,hotspot_min_samples passed to
#'   [match_cosmic_hotspots()]
#' @param hotspot_matches optional precomputed result of [match_cosmic_hotspots()]
#' @param intron_linewidth line width of the intron backbone in the gene track
#' @param highlight_colour colour of halos / guide lines for matches
#' @return a patchwork object (ggplot); use `print()` or `ggsave()`. If
#'   highlighting is on, the match table is in `attr(p, "hotspot_matches")`.
#' @export
plot_gene_mutations <- function(model, cosmic, species_mut, species_label = "Species",
                                species_panel = c("top", "bottom"),
                                compress_introns = TRUE, coding_only = TRUE, min_label_n = 3,
                                highlight_hotspots = FALSE, hotspot_top_n = 10,
                                hotspot_window = 3, hotspot_min_samples = 2,
                                hotspot_matches = NULL, intron_linewidth = 1.1,
                                highlight_colour = "#E0A800") {
  species_panel <- match.arg(species_panel)
  tr <- .make_transform(model$exons, compress_introns)
  tf <- tr$f

  cs <- as.data.table(cosmic)[!is.na(start)]
  if (coding_only) cs <- cs[!is.na(aa_pos)]
  cs_agg <- cs[, .(n = .N), by = .(pos = start, aa_change, consequence)]

  sp <- as.data.table(species_mut)[!is.na(human_pos)]
  if (!"aa_change" %in% names(sp)) sp[, aa_change := paste0("p.", aa_pos)]
  if (!"consequence" %in% names(sp)) sp[, consequence := NA_character_]
  sp_agg <- sp[, .(n = .N), by = .(pos = human_pos, aa_change, consequence)]

  # --- hotspot matching ------------------------------------------------------
  hl_x <- NULL; caption <- NULL
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
    caption <- sprintf(
      "Highlighted: %d user mutation(s) within ±%d aa of the top %d COSMIC hotspot residues (≥%d samples).",
      nrow(m), hotspot_window, hotspot_top_n, hotspot_min_samples)
  }

  # --- panels ----------------------------------------------------------------
  n_unmapped <- sum(is.na(as.data.table(species_mut)$human_pos))
  sp_title <- sprintf("%s mutations (n = %d mapped%s)", species_label, nrow(sp),
                      if (n_unmapped) sprintf(", %d unmapped", n_unmapped) else "")
  cs_title <- sprintf("Human COSMIC: %s (%d mutations, %d samples)",
                      model$gene$symbol, nrow(cs_agg), nrow(cs))
  sp_top <- species_panel == "top"

  p_sp <- .lollipop(sp_agg, tf, sp_title, flip = !sp_top, min_label_n = 1,
                    hl_x = hl_x, hl_col = highlight_colour)
  p_cs <- .lollipop(cs_agg, tf, cs_title, flip = sp_top, min_label_n = min_label_n,
                    hl_x = hl_x, hl_col = highlight_colour)
  mid <- .gene_track(model, tf, intron_linewidth = intron_linewidth,
                     hl_x = hl_x, hl_col = highlight_colour)

  xr <- range(c(tf(model$exons$start), tf(model$exons$end), tf(cs_agg$pos), tf(sp_agg$pos)), na.rm = TRUE)
  pad <- diff(xr) * 0.02
  lock <- function(p) p + coord_cartesian(xlim = xr + c(-pad, pad), clip = "off")
  p_sp <- lock(p_sp); p_cs <- lock(p_cs)
  mid <- mid + coord_cartesian(xlim = xr + c(-pad, pad), ylim = c(-1, 1))

  panels  <- if (sp_top) list(p_sp, mid, p_cs) else list(p_cs, mid, p_sp)
  heights <- if (sp_top) c(3, 1, 4) else c(4, 1, 3)

  out <- patchwork::wrap_plots(panels, ncol = 1, heights = heights) +
    patchwork::plot_layout(guides = "collect")
  if (!is.null(caption)) {
    out <- out + patchwork::plot_annotation(caption = caption,
                                            theme = theme(plot.caption = element_text(colour = "#7A5600")))
  }
  attr(out, "hotspot_matches") <- hotspot_matches
  out
}
