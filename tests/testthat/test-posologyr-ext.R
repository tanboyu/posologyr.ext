test_that("poso_ev, %then_poso%, and poso_add_obs work cleanly", {
  ev1 <- poso_ev(amt = 1000, dur = 1, cmt = 1, covariates = list(WT = 70))
  ev2 <- poso_ev(amt = 500, dur = 1, ii = 12, addl = 2, cmt = 1)
  plan <- (ev1 %then_poso% ev2) |>
    poso_add_obs(time = c(1.5, 11.5), dv = c(24.0, 8.0), is_tdm = TRUE)
  out <- as_poso_data(plan)
  expect_true(all(c("ID", "TIME", "DV", "AMT", "EVID", "WT") %in% names(out)))
  expect_true(all(out$WT == 70))
})

test_that("as_poso_map_tidy handles unequal row counts safely", {
  df_ip <- data.frame(time = c(0, 1, 2, 3), Cc = c(0, 20, 15, 11))
  df_pr <- data.frame(time = c(0, 1.5, 3), Cc = c(0, 16, 10))
  m_tidy <- as_poso_map_tidy(df_ip, df_pr, var_name = "Cc")
  expect_s3_class(m_tidy, "poso_map_tidy")
  expect_equal(nrow(m_tidy$indiv_df), 4L)
  expect_true("pred" %in% names(m_tidy$indiv_df))
  g <- plot(m_tidy, target_range = c(10, 20))
  expect_s3_class(g, "ggplot")
})

test_that("as_poso_sim and poso_plot_bundle work", {
  sim_df <- data.frame(
    id = rep(1:5, each = 4),
    time = rep(c(0, 1, 6, 12), times = 5),
    Cc = rnorm(20, mean = 15, sd = 2)
  )
  s <- as_poso_sim(sim_df, vars = "Cc")
  bundle <- poso_plot_bundle(prior_sim = s, include = "prior", target_range = c(10, 20))
  expect_s3_class(bundle, "poso_plot_bundle")
})
