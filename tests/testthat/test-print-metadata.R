test_that("print_metadata lists every key with aligned values", {
  obj <- set_metadata(make_isolated_seed_fixture()$obj, "author", "Jane Doe")
  out <- cli::cli_fmt(print_metadata(obj))
  expect_true(any(grepl("^author +: Jane Doe$", out)))
  expect_true(any(grepl("^network_cutoff +: 0$", out)))
  key_lines <- grep(" : ", out, value = TRUE)
  expect_length(key_lines, length(names(obj@metadata)))
  expect_length(unique(regexpr(" : ", key_lines, fixed = TRUE)), 1L)
})
