#' Run the full cross-species comparison for one gene
#'
#' @param gene human gene symbol
#' @param mutations data.frame of species mutations (see [map_mutations_to_human()])
#' @param species species name in any form accepted by [ens_resolve_species()]
#' @param assembly expected species assembly (e.g. "Mmul_10"); checked against Ensembl
#' @param cosmic COSMIC mutation table from [cosmic_load()] or a file path
#' @param census_path optional Cancer Gene Census file
#' @param human_assembly "GRCh38" or "GRCh37" (must match your COSMIC download)
#' @param species_panel "top" (default) or "bottom": where user mutations are drawn
#' @param highlight_hotspots highlight user mutations at/near top COSMIC hotspots
#' @param hotspot_top_n number of top COSMIC residues treated as hotspots
#' @param hotspot_window max distance (human aa) for a close match; 0 = same residue
#' @param hotspot_min_samples minimum COSMIC samples for a residue to count
#' @param ... passed to [plot_gene_mutations()]
#' @return list: gene_check, species, species_name, assembly, model,
#'   cosmic_mutations, species_mutations, hotspot_matches (NULL unless
#'   highlighting), plot
#' @export
xsmut_pipeline <- function(gene, mutations, species, assembly = NULL, cosmic,
                           census_path = NULL, human_assembly = c("GRCh38", "GRCh37"),
                           species_panel = c("top", "bottom"),
                           highlight_hotspots = FALSE, hotspot_top_n = 10,
                           hotspot_window = 3, hotspot_min_samples = 2, ...) {
  human_assembly <- match.arg(human_assembly)
  species_panel <- match.arg(species_panel)
  host <- ens_host(human_assembly)
  sp <- ens_resolve_species(species, host = host)
  asm <- ens_assembly(sp, assembly, host = host)

  if (is.character(cosmic)) cosmic <- cosmic_load(cosmic, genes = gene)
  cos_mut <- cosmic_gene_mutations(gene, cosmic)
  check <- cosmic_gene_listed(gene, census_path = census_path, mutations = cosmic)
  if (!isTRUE(check$in_census) && check$n_cosmic_mutations == 0)
    warning(gene, " has no COSMIC entries in the supplied files")

  model <- get_gene_model(gene, species = "homo_sapiens", human_assembly = human_assembly)
  if (is_dndscv(mutations)) mutations <- as_xsmut(mutations)
  sp_mut <- map_mutations_to_human(mutations, gene, sp, human_assembly = human_assembly)

  matches <- if (isTRUE(highlight_hotspots)) {
    match_cosmic_hotspots(sp_mut, cos_mut, top_n = hotspot_top_n, window = hotspot_window,
                          min_samples = hotspot_min_samples)
  } else NULL

  sp_name <- .species_display_name(sp)          # "Macaca mulatta"
  label <- sprintf("%s (%s)", sp_name, if (!is.null(assembly)) assembly else asm$assembly_name)
  p <- plot_gene_mutations(model, cos_mut, sp_mut, species_label = label,
                           species_panel = species_panel,
                           highlight_hotspots = highlight_hotspots,
                           hotspot_top_n = hotspot_top_n, hotspot_window = hotspot_window,
                           hotspot_min_samples = hotspot_min_samples,
                           hotspot_matches = matches, ...)

  list(gene_check = check, species = sp, species_name = sp_name, assembly = asm, model = model,
       cosmic_mutations = cos_mut, species_mutations = sp_mut, hotspot_matches = matches,
       plot = p)
}

#' Launch the Shiny app
#' @param ... passed to [shiny::runApp()]
#' @export
run_app <- function(...) {
  shiny::runApp(system.file("app", package = "xsmut"), ...)
}
