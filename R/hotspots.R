# Match user (species) mutations against the most recurrent human COSMIC residues.
#
# "Most important" = the top-N COSMIC residues of the gene ranked by the number
# of samples carrying a coding, non-synonymous change at that residue (the
# classic hotspot definition, e.g. TP53 R175, R248, R273; KRAS G12, G13, Q61).
# A user mutation matches when its human-mapped residue lies within `window` aa.

# "p.R175H" -> "H", "p.Q136*" -> "*", "p.P72fs*3" / "p.72fs" -> "fs", "p.E285_K286del" -> "del"
.aa_suffix <- function(x) {
  x <- as.character(x)
  s <- sub("^p\\.\\(?([A-Z][a-z]{0,2}|\\*)?\\d+(_[A-Z][a-z]{0,2}\\d+)?", "", x)
  s <- sub("\\)$", "", s)
  s[grepl("fs", s)] <- "fs"
  s[grepl("^del", s)] <- "del"
  s[grepl("^(ins|dup)", s)] <- "ins"
  s[is.na(x)] <- NA_character_
  s
}

#' Top recurrent COSMIC residues for a gene
#'
#' @param cosmic data.table from [cosmic_gene_mutations()]
#' @param top_n number of residues to keep
#' @param min_samples minimum samples for a residue to qualify
#' @return data.table: rank, aa_pos, n_samples, top_change, top_change_n,
#'   start (human genomic position of the top change), changes (all changes)
#' @export
cosmic_hotspots <- function(cosmic, top_n = 10, min_samples = 2) {
  cs <- as.data.table(cosmic)
  cs <- cs[!is.na(aa_pos) & !grepl("synonymous|silent|unknown|intronic", tolower(consequence))]
  if (!nrow(cs)) {
    return(data.table(rank = integer(), aa_pos = integer(), n_samples = integer(),
                      top_change = character(), top_change_n = integer(),
                      start = integer(), changes = character()))
  }
  by_change <- cs[, .(n = .N, start = start[1]), by = .(aa_pos, aa_change)][order(aa_pos, -n)]
  hs <- by_change[, .(n_samples = sum(n), top_change = aa_change[1], top_change_n = n[1],
                      start = start[1], changes = paste(unique(aa_change), collapse = ";")),
                  by = aa_pos]
  hs <- hs[n_samples >= min_samples][order(-n_samples, aa_pos)]
  hs <- hs[seq_len(min(top_n, nrow(hs)))]
  hs[, rank := seq_len(.N)]
  data.table::setcolorder(hs, c("rank", "aa_pos", "n_samples", "top_change", "top_change_n", "start", "changes"))
  hs[]
}

#' Find user mutations at or near the top COSMIC hotspots
#'
#' @param species_mut data.table from [map_mutations_to_human()]
#' @param cosmic data.table from [cosmic_gene_mutations()]
#' @param top_n,min_samples see [cosmic_hotspots()]
#' @param window maximum distance in human residues for a "close" match (0 = same residue only)
#' @return data.table, one row per matched user mutation (nearest hotspot kept):
#'   user columns, hotspot columns and `distance_aa`, `match_type`
#'   ("identical change", "same residue", "within N aa"). Attribute `hotspots`.
#' @export
match_cosmic_hotspots <- function(species_mut, cosmic, top_n = 10, window = 3, min_samples = 2) {
  hs <- cosmic_hotspots(cosmic, top_n = top_n, min_samples = min_samples)
  sp <- as.data.table(species_mut)[!is.na(human_aa_pos)]
  if (!"aa_change" %in% names(sp)) sp[, aa_change := paste0("p.", aa_pos)]
  if (!"sample_id" %in% names(sp)) sp[, sample_id := NA_character_]
  if (!"consequence" %in% names(sp)) sp[, consequence := NA_character_]

  empty <- data.table(sample_id = character(), aa_change = character(), consequence = character(),
                      human_aa_pos = integer(), human_pos = integer(), hotspot_rank = integer(),
                      hotspot_aa_pos = integer(), hotspot_change = character(),
                      hotspot_n_samples = integer(), distance_aa = integer(), match_type = character())
  if (!nrow(hs) || !nrow(sp)) { attr(empty, "hotspots") <- hs; return(empty) }

  sp[, .row := .I]
  pairs <- data.table::CJ(.row = sp$.row, .hs = seq_len(nrow(hs)))
  pairs[, distance_aa := abs(sp$human_aa_pos[.row] - hs$aa_pos[.hs])]
  pairs <- pairs[distance_aa <= window][order(.row, distance_aa, .hs)][, .SD[1], by = .row]
  if (!nrow(pairs)) { attr(empty, "hotspots") <- hs; return(empty) }

  out <- sp[pairs$.row, .(sample_id, aa_change, consequence, human_aa_pos, human_pos)]
  h <- hs[pairs$.hs]
  out[, `:=`(hotspot_rank = h$rank, hotspot_aa_pos = h$aa_pos, hotspot_change = h$top_change,
             hotspot_n_samples = h$n_samples, distance_aa = as.integer(pairs$distance_aa))]

  # identical change: same residue and same resulting aa / indel class as any COSMIC change there
  user_sfx <- .aa_suffix(out$aa_change)
  same_change <- mapply(function(s, ch) !is.na(s) && nzchar(s) && s %in% .aa_suffix(strsplit(ch, ";")[[1]]),
                        user_sfx, h$changes)
  out[, match_type := data.table::fcase(
    distance_aa == 0L & same_change, "identical change",
    distance_aa == 0L, "same residue",
    default = paste0("within ", distance_aa, " aa")
  )]
  attr(out, "hotspots") <- hs
  out[order(hotspot_rank, distance_aa)]
}
