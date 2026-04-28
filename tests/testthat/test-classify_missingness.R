test_that("classify_missingness rejects non-SummarizedExperiment input", {
  expect_error(classify_missingness(matrix(1:6, nrow = 2)))
  expect_error(classify_missingness(data.frame(a = 1:3)))
  expect_error(classify_missingness(NULL))
})

test_that("classify_missingness errors when no missing values are present", {
  counts <- matrix(
    seq_len(12),
    nrow = 4,
    dimnames = list(
      paste0("p", seq_len(4)),
      paste0("s", seq_len(3))
    )
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(intensity = counts),
    colData = S4Vectors::DataFrame(condition = c("A", "A", "B"))
  )
  expect_error(
    classify_missingness(se),
    regexp = "missing"
  )
})
