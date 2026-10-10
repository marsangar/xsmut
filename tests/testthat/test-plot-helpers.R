test_that("sample axes use integer breaks and drop the 0 baseline", {
  expect_equal(xsmut:::.int_breaks(c(0, 1.22)), 1)
  expect_equal(xsmut:::.int_breaks(c(-1.22, 0)), -1)
  b <- xsmut:::.int_breaks(c(0, 7.3))
  expect_true(all(b == round(b)) && !0 %in% b)
})

test_that("halo colour presets", {
  expect_equal(xsmut:::.resolve_colour("red"), "#E8000B")
  expect_equal(xsmut:::.resolve_colour("gold"), "#E0A800")
  expect_equal(xsmut:::.resolve_colour("#123456"), "#123456")
})

test_that("COSMIC gene header text", {
  info <- data.table::data.table(gene = "TP53", in_census = TRUE, tier = "1",
                                 role = "oncogene, TSG, fusion", census_supplied = TRUE,
                                 tumour_types_somatic = "breast, colorectal, lung, sarcoma")
  txt <- xsmut:::.gene_info_text(info, max_tumour_types = 2)
  expect_false(grepl("tier", txt))
  expect_match(txt, "Role: oncogene, TSG, fusion")
  expect_match(txt, "breast, colorectal, \\.\\.\\.")
  expect_match(xsmut:::.gene_info_text(data.table::data.table(in_census = FALSE, census_supplied = TRUE)),
               "Not in COSMIC")
  expect_null(xsmut:::.gene_info_text(data.table::data.table(in_census = TRUE, census_supplied = FALSE)))
})

test_that("species HGVS changes are rewritten in human numbering", {
  # species residue k -> human residue k + 1 (one extra residue at the human N-terminus)
  t2h <- 1:10 + 1L
  hseq <- "MPYRSTEKAQW"          # human residues 1..11
  out <- xsmut:::.to_human_hgvs(c("p.Y2C", "p.2fs", "p.R3_S4del", "p.Q9*", NA, "p.Arg3His"),
                                aa_pos = c(2L, 2L, 3L, 9L, 5L, 3L), t2h = t2h, human_seq = hseq)
  expect_equal(out, c("p.Y3C", "p.Y3fs", "p.R4_S5del", "p.Q10*", "p.T6", NA))
  # residue aligned to a human gap -> NA
  expect_true(is.na(xsmut:::.to_human_hgvs("p.A2V", 2L, c(1L, NA, 2L), "MA")))
})
