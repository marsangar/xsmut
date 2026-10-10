test_that("species names keep binomial capitalisation", {
  expect_equal(xsmut:::.species_display_name("macaca_mulatta"), "Macaca mulatta")
  expect_equal(xsmut:::.species_display_name("canis_lupus_familiaris"), "Canis lupus familiaris")
})

test_that("aa suffix parsing", {
  expect_equal(xsmut:::.aa_suffix(c("p.R175H", "p.Q136*", "p.P72fs*3", "p.72fs", "p.E285_K286del")),
               c("H", "*", "fs", "fs", "del"))
})

test_that("hotspot matching: identical, same residue, nearby, none", {
  cos <- data.table::data.table(
    aa_change = c(rep("p.R175H", 10), rep("p.R175G", 2), rep("p.R248Q", 6), "p.V10A"),
    aa_pos = c(rep(175L, 12), rep(248L, 6), 10L),
    consequence = "missense_variant",
    start = c(rep(7675088L, 12), rep(7674220L, 6), 7676000L))
  sp <- data.table::data.table(
    sample_id = c("a", "b", "c", "d"),
    aa_change = c("p.R175H", "p.R175C", "p.T250M", "p.A50T"),
    consequence = "missense_variant",
    human_aa_pos = c(175L, 175L, 250L, 50L),
    human_pos = c(7675088L, 7675088L, 7674214L, 7676500L))
  m <- match_cosmic_hotspots(sp, cos, top_n = 2, window = 3)
  expect_equal(nrow(m), 3)
  expect_equal(m[sample_id == "a", match_type], "identical change")
  expect_equal(m[sample_id == "b", match_type], "same residue")
  expect_equal(m[sample_id == "c", match_type], "within 2 aa")
  expect_false("d" %in% m$sample_id)
})
