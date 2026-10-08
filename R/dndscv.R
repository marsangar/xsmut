# Adapter for dNdScv `annotmuts` output (dndscv()$annotmuts), e.g.
#
#   sampleID chr         pos       ref mut gene strand ref_cod mut_cod ref3_cod mut3_cod aachange ntchange codonsub impact    pid
#   MQD0001d NC_041754.1 95090953  A   G   RIT1 1      A       G       TAT      TGT      Y75C     A224G    TAT>TGT  Missense  ENSMMUP00000024831
#   MQD0001d NC_041755.1 108560386 CA  C   BAP1 -1     .       .       .        .        .        459-459-delfrshift .  no-SNV ENSMMUP00000022338
#
# Every row carries the Ensembl protein the call was annotated against (`pid`), so
# the mapper aligns exactly that isoform to the human ortholog instead of guessing.

.dndscv_cols <- c("sampleID", "chr", "pos", "ref", "mut", "gene", "aachange", "ntchange", "impact", "pid")

#' Is this a dNdScv `annotmuts` table?
#' @param x data.frame
#' @export
is_dndscv <- function(x) {
  all(c("sampleID", "gene", "aachange", "impact", "pid") %in% names(x))
}

.dndscv_impact <- function(impact, ntchange) {
  impact <- as.character(impact); nt <- as.character(ntchange)
  data.table::fcase(
    impact == "Missense",            "missense_variant",
    impact == "Nonsense",            "stop_gained",
    impact == "Synonymous",          "synonymous_variant",
    impact == "Essential_Splice",    "splice_site_variant",
    impact == "Stop_loss",           "stop_lost",
    impact == "no-SNV" & grepl("frshift", nt), "frameshift_variant",
    impact == "no-SNV" & grepl("del", nt),     "inframe_deletion",
    impact == "no-SNV" & grepl("ins", nt),     "inframe_insertion",
    impact == "no-SNV",              "indel",
    default = tolower(impact)
  )
}

# "A224G" -> c.224A>G ; "459-459-delfrshift" -> c.459del (frameshift) ; "." -> NA
.dndscv_cds <- function(ntchange) {
  nt <- as.character(ntchange)
  snv <- grepl("^[ACGT]\\d+[ACGT]$", nt)
  out <- rep(NA_character_, length(nt))
  out[snv] <- sub("^([ACGT])(\\d+)([ACGT])$", "c.\\2\\1>\\3", nt[snv])
  ind <- grepl("^\\d+-\\d+-(del|ins)", nt)
  if (any(ind)) {
    st <- sub("^(\\d+)-(\\d+)-.*$", "\\1", nt[ind]); en <- sub("^(\\d+)-(\\d+)-.*$", "\\2", nt[ind])
    op <- ifelse(grepl("del", nt[ind]), "del", "ins")
    out[ind] <- ifelse(st == en, paste0("c.", st, op), paste0("c.", st, "_", en, op))
  }
  out
}

# first CDS base touched by the change -> protein residue (ceiling(cds/3))
.dndscv_cds_pos <- function(ntchange) {
  nt <- as.character(ntchange)
  p <- suppressWarnings(as.integer(sub("^[ACGT]?(\\d+).*$", "\\1", nt)))
  p[!grepl("^[ACGT]?\\d+", nt)] <- NA_integer_
  p
}

#' Convert a dNdScv `annotmuts` table to the xsmut mutation schema
#'
#' Produces the columns the mapper and plots expect while keeping every original
#' column (prefixed `dndscv_` where it would clash). Protein positions come from
#' `aachange` for SNVs and from the CDS coordinate in `ntchange` for indels, so
#' frameshifts are placed at the first affected codon. Chromosome names (NCBI
#' accessions in the example) are passed through untouched; the mapper normalises
#' them with [normalize_chrom()].
#'
#' @param x data.frame from `dndscv()$annotmuts`
#' @return data.table with sample_id, chrom, pos, ref, alt, gene, strand,
#'   aa_change, aa_pos, cds_change, consequence, species_protein_id, impact
#' @export
as_xsmut <- function(x) {
  d <- as.data.table(x)
  miss <- setdiff(.dndscv_cols, names(d))
  if (length(miss)) stop("Not a dNdScv annotmuts table; missing: ", paste(miss, collapse = ", "))
  aach <- as.character(d$aachange)
  has_aa <- !is.na(aach) & !aach %in% c(".", "-", "")
  cds_pos <- .dndscv_cds_pos(d$ntchange)
  out <- data.table(
    sample_id          = as.character(d$sampleID),
    chrom              = as.character(d$chr),
    pos                = as.integer(d$pos),
    ref                = as.character(d$ref),
    alt                = as.character(d$mut),
    gene               = as.character(d$gene),
    strand             = if ("strand" %in% names(d)) as.integer(d$strand) else NA_integer_,
    aa_change          = ifelse(has_aa, paste0("p.", aach), NA_character_),
    aa_pos             = NA_integer_,
    cds_change         = .dndscv_cds(d$ntchange),
    consequence        = .dndscv_impact(d$impact, d$ntchange),
    species_protein_id = as.character(d$pid),
    impact             = as.character(d$impact)
  )
  out[, aa_pos := parse_aa_pos(aa_change)]
  # indels: fall back to the codon containing the first changed CDS base
  out[, .cds_pos := cds_pos]
  out[is.na(aa_pos) & !is.na(.cds_pos) & impact != "Essential_Splice",
      aa_pos := as.integer(ceiling(.cds_pos / 3))]
  out[, .cds_pos := NULL]
  # label indels in the plot by their start codon
  out[is.na(aa_change) & !is.na(aa_pos),
      aa_change := paste0("p.", aa_pos, data.table::fcase(
        consequence == "frameshift_variant", "fs",
        consequence == "inframe_deletion", "del",
        consequence == "inframe_insertion", "ins",
        default = "indel"))]
  # keep extra dndscv columns for reference
  extra <- setdiff(names(d), c(.dndscv_cols, "strand"))
  for (col in extra) data.table::set(out, j = paste0("dndscv_", col), value = d[[col]])
  out[]
}

#' Read a dNdScv annotmuts file (tsv/csv/rds) and convert it
#' @param path file path
#' @export
read_dndscv <- function(path) {
  x <- if (grepl("\\.rds$", path, ignore.case = TRUE)) readRDS(path) else fread(path)
  as_xsmut(x)
}
