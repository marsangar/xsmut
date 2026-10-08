# Map mutations annotated in species X onto human reference coordinates.
#
# Strategy (species-agnostic, works for any Ensembl species):
#   species aa position  --(ortholog protein alignment)-->  human aa position
#                        --(Ensembl /map/translation)-->    human genomic position
#
# An optional genomic route (chrom/pos + UCSC chain file) is provided via
# rtracklayer::liftOver for non-coding variants.

# Build a lookup: position in seq B (ungapped) -> position in seq A (ungapped),
# from two gapped, equal-length alignment strings.
.aln_map <- function(gapped_a, gapped_b) {
  a <- strsplit(gapped_a, "")[[1]]
  b <- strsplit(gapped_b, "")[[1]]
  stopifnot(length(a) == length(b))
  ia <- cumsum(a != "-"); ib <- cumsum(b != "-")
  keep <- b != "-"
  map <- rep(NA_integer_, sum(keep))
  map[ib[keep]] <- ifelse(a[keep] != "-", ia[keep], NA_integer_)
  map
}

.align_proteins <- function(seq_human, seq_target) {
  if (requireNamespace("pwalign", quietly = TRUE)) {
    aln <- pwalign::pairwiseAlignment(Biostrings::AAString(seq_target), Biostrings::AAString(seq_human),
                                      substitutionMatrix = "BLOSUM62", gapOpening = 10,
                                      gapExtension = 0.5, type = "global")
    list(target = as.character(pwalign::alignedPattern(aln)),
         human  = as.character(pwalign::alignedSubject(aln)))
  } else {
    aln <- Biostrings::pairwiseAlignment(Biostrings::AAString(seq_target), Biostrings::AAString(seq_human),
                                         substitutionMatrix = "BLOSUM62", gapOpening = 10,
                                         gapExtension = 0.5, type = "global")
    list(target = as.character(Biostrings::alignedPattern(aln)),
         human  = as.character(Biostrings::alignedSubject(aln)))
  }
}

# residue -> genomic coordinate of the first codon base, for a human protein
.residue_to_genome <- function(protein_id, protein_length, host) {
  blocks <- ens_map_translation(protein_id, 1, protein_length, host = host)
  if (!nrow(blocks)) stop("Ensembl returned no genomic mapping for ", protein_id)
  strand <- blocks$strand[1]
  blocks <- if (strand == 1) blocks[order(start)] else blocks[order(-start)]
  # expand CDS bases in transcript order
  bases <- unlist(lapply(seq_len(nrow(blocks)), function(i)
    if (strand == 1) seq(blocks$start[i], blocks$end[i]) else seq(blocks$end[i], blocks$start[i])))
  function(aa) {
    idx <- 3L * aa - 2L
    out <- rep(NA_integer_, length(aa))
    ok <- !is.na(idx) & idx >= 1 & idx <= length(bases)
    out[ok] <- bases[idx[ok]]
    out
  }
}

#' Map species mutations for one gene onto human coordinates
#'
#' @param mutations data.frame with at least `gene` and one of `aa_pos`
#'   (integer protein position in the species protein) or `aa_change` (HGVS,
#'   e.g. "p.R175H"). Optional: `sample_id`, `consequence`, `chrom`, `pos`,
#'   `ref`, `alt`.
#' @param gene human gene symbol
#' @param species Ensembl production name of the source species (see
#'   [ens_resolve_species()])
#' @param species_protein_id optional Ensembl protein ID the mutations were
#'   annotated against; if it differs from the ortholog Ensembl picks, a fresh
#'   pairwise alignment is computed.
#' @param human_assembly "GRCh38" or "GRCh37"
#' @param chain optional path to a UCSC chain file (species -> human) enabling
#'   genomic liftOver for rows lacking a protein position (requires rtracklayer)
#' @return data.table: input columns plus `human_aa_pos`, `human_chrom`,
#'   `human_pos`, `map_method`, `map_note`. Attributes: `ortholog`, `alignment`.
#' @export
map_mutations_to_human <- function(mutations, gene, species, species_protein_id = NULL,
                                   human_assembly = c("GRCh38", "GRCh37"), chain = NULL) {
  host <- ens_host(match.arg(human_assembly))
  mt <- as.data.table(mutations)
  if (is_dndscv(mt)) mt <- as_xsmut(mt)
  if (!"gene" %in% names(mt)) stop("`mutations` needs a `gene` column")
  g <- toupper(gene)
  mt <- mt[toupper(mt$gene) == g]
  if (!"aa_pos" %in% names(mt)) {
    mt[, aa_pos := if ("aa_change" %in% names(mt)) parse_aa_pos(aa_change) else NA_integer_]
  }
  mt[, aa_pos := as.integer(aa_pos)]

  orth <- get_ortholog(gene, species, host = host)
  human_pid <- orth$human$protein_id
  
  # one isoform for all rows (argument) / per-row isoforms (column, e.g. dNdScv pid) / Ensembl's pick
  if (!is.null(species_protein_id)) {
    mt[, species_protein_id := species_protein_id]
  } else if (!"species_protein_id" %in% names(mt)) {
    mt[, species_protein_id := orth$target$protein_id]
  }
  mt[is.na(species_protein_id) | species_protein_id == "", species_protein_id := orth$target$protein_id]
  
  align_one <- function(pid) {
    if (pid == orth$target$protein_id && !is.null(orth$human$seq) && !is.null(orth$target$seq)) {
      aln <- list(human = orth$human$seq, target = orth$target$seq); note <- "Ensembl ortholog alignment"
    } else {
      seq_h <- ens_protein_seq(human_pid, host = host)
      seq_t <- tryCatch(ens_protein_seq(pid, host = host), error = function(e) NULL)
      if (is.null(seq_t)) return(list(aln = NULL, t2h = integer(), note = paste("unknown protein", pid)))
      aln <- .align_proteins(seq_h, seq_t)
      note <- if (pid == orth$target$protein_id) "pairwise alignment (pwalign)" else sprintf("pairwise alignment to %s", pid)
    }
    list(aln = aln, t2h = .aln_map(aln$human, aln$target), note = note)
  }
  pids <- unique(mt$species_protein_id)
  alns <- setNames(lapply(pids, align_one), pids)
  
  human_len <- if (!is.null(orth$human$seq)) sum(strsplit(orth$human$seq, "")[[1]] != "-")
  else nchar(ens_protein_seq(human_pid, host = host))
  r2g    <- .residue_to_genome(human_pid, human_len, host = host)
  hchrom <- ens_map_translation(human_pid, 1, 1, host = host)$chrom[1]
  
  mt[, `:=`(human_aa_pos = NA_integer_, map_note = NA_character_)]
  for (pid in pids) {
    a <- alns[[pid]]
    i <- which(mt$species_protein_id == pid)
    if (is.null(a$aln)) { mt[i, map_note := a$note]; next }
    ok <- i[!is.na(mt$aa_pos[i]) & mt$aa_pos[i] >= 1 & mt$aa_pos[i] <= length(a$t2h)]
    mt[ok, human_aa_pos := a$t2h[aa_pos]]
    mt[i, map_note := data.table::fcase(
      !is.na(human_aa_pos), a$note,
      is.na(aa_pos), "no protein position",
      aa_pos < 1 | aa_pos > length(a$t2h), "aa_pos outside protein",
      default = "residue aligned to gap in human")]
  }
  mt[, human_pos := r2g(human_aa_pos)]
  mt[, human_chrom := ifelse(is.na(human_pos), NA_character_, hchrom)]
  mt[, map_method := ifelse(!is.na(human_pos), "protein", NA_character_)]
  # --- optional genomic liftOver for unmapped rows ---------------------------
  if (!is.null(chain) && all(c("chrom", "pos") %in% names(mt)) && any(is.na(mt$human_pos))) {
    if (!requireNamespace("rtracklayer", quietly = TRUE))
      stop("rtracklayer is required for chain-file liftOver")
    todo <- which(is.na(mt$human_pos) & !is.na(mt$pos))
    if (length(todo)) {
      gr <- GenomicRanges::GRanges(paste0("chr", sub("^chr", "", mt$chrom[todo])),
                                   IRanges::IRanges(mt$pos[todo], mt$pos[todo]))
      ch <- rtracklayer::import.chain(chain)
      lo <- rtracklayer::liftOver(gr, ch)
      n <- lengths(lo)
      hit <- which(n == 1)
      if (length(hit)) {
        lifted <- unlist(lo[hit])
        mt[todo[hit], `:=`(human_chrom = sub("^chr", "", as.character(GenomicRanges::seqnames(lifted))),
                           human_pos = as.integer(GenomicRanges::start(lifted)),
                           map_method = "liftOver", map_note = "chain-file liftOver")]
      }
    }
  }

  attr(mt, "ortholog") <- orth
  attr(mt, "alignment") <- lapply(alns, `[[`, "aln")   # one per species protein id  mt[]
}
