# ==============================================================================
# 完整示例：万古霉素（Vancomycin）MIPD 全流程（整合 posologyr + posologyr.ext）
# ==============================================================================

library(rxode2)
library(posologyr)
library(posologyr.ext)
library(ggplot2)

mod_vancomycin_2cmt <- function() {
  ini({
    THETA_CL <- 4.50
    THETA_V1 <- 38.0
    THETA_Q  <- 3.20
    THETA_V2 <- 45.0

    ETA_CL ~ 0.09
    ETA_V1 ~ 0.0625
    ETA_Q  ~ 0.16
    ETA_V2 ~ 0.16

    add_sd  <- 1.20
    prop_sd <- 0.15
  })
  model({
    TVCL <- THETA_CL * (ClCr / 100)^0.75
    TVV1 <- THETA_V1 * (WT / 70)
    TVQ  <- THETA_Q  * (WT / 70)^0.75
    TVV2 <- THETA_V2 * (WT / 70)

    CL <- TVCL * exp(ETA_CL)
    V1 <- TVV1 * exp(ETA_V1)
    Q  <- TVQ  * exp(ETA_Q)
    V2 <- TVV2 * exp(ETA_V2)

    Cc <- centr / V1
    d/dt(centr)  <- -(CL / V1) * centr - (Q / V1) * centr + (Q / V2) * periph
    d/dt(periph) <-  (Q / V1) * centr - (Q / V2) * periph

    Cc ~ add(add_sd) + prop(prop_sd) + combined1()
  })
}

# 1. 事件与观察录入
patient_ev <- poso_ev(
  time = 0, amt = 1500, dur = 1.5, cmt = 1,
  covariates = list(WT = 65, ClCr = 55)
) %then_poso%
  poso_ev(time = 0, amt = 1000, dur = 1.0, ii = 12, addl = 2, cmt = 1) |>
  poso_add_obs(time = c(2.5, 35.5), dv = c(31.4, 18.6), cmt = 1, is_tdm = TRUE)

tdm_data <- as_poso_data(patient_ev, keep_meta = FALSE)

# 2. 多线程先验模拟（默认自动启用多线程，可通过 n_cores 指定）
prior_grid <- patient_ev |>
  poso_add_obs(time = seq(0, 48, by = 0.25), cmt = 1, is_tdm = FALSE)

sim_prior <- poso_sim_tidy(
  model   = mod_vancomycin_2cmt,
  events  = prior_grid,
  n_sim   = 300,
  vars    = "Cc",
  seed    = 2026,
  n_cores = 4
)

# 3. MAP 个体估计与一键生成 IPRED + PRED 轨迹对象
fit_map <- poso_estim_map(dat = tdm_data, prior_model = mod_vancomycin_2cmt)
map_obj <- poso_map_tidy(
  fit_map  = fit_map,
  model    = mod_vancomycin_2cmt,
  events   = prior_grid,
  var_name = "Cc"
)

# 4. 剂量优化（注意 poso_dose_conc 不接收 addl 参数）
opt_dose <- poso_dose_conc(
  dat         = tdm_data,
  prior_model = mod_vancomycin_2cmt,
  tdm         = TRUE,
  target_conc = 12.0,
  duration    = 1.0,
  time_dose   = 48,
  time_c      = 60,
  inter       = 12
)
recommended_dose <- round(opt_dose$dose / 50) * 50

# 5. 前瞻性方案对比（安全对齐不同时间网格行数，避免长度不一致报错）
full_plan <- (
  poso_ev(time = 0, amt = 1500, dur = 1.5, cmt = 1, covariates = list(WT = 65, ClCr = 55)) %then_poso%
    poso_ev(time = 0, amt = 1000, dur = 1.0, ii = 12, addl = 2, cmt = 1) %then_poso%
    poso_ev(time = 0, amt = recommended_dose, dur = 1.0, ii = 12, addl = 3, cmt = 1)
) |>
  poso_add_obs(time = seq(0, 96, by = 0.25), cmt = 1, is_tdm = FALSE) |>
  poso_add_obs(time = c(2.5, 35.5), dv = c(31.4, 18.6), cmt = 1, is_tdm = TRUE)

unadj_plan <- (
  poso_ev(time = 0, amt = 1500, dur = 1.5, cmt = 1, covariates = list(WT = 65, ClCr = 55)) %then_poso%
    poso_ev(time = 0, amt = 1000, dur = 1.0, ii = 12, addl = 6, cmt = 1)
) |>
  poso_add_obs(time = seq(0, 96, by = 0.25), cmt = 1, is_tdm = FALSE)

rx_mod <- rxode2::rxode2(mod_vancomycin_2cmt)
ipred_params <- fit_map$indiv_param[1, , drop = FALSE]
sol_adj   <- as.data.frame(rxode2::rxSolve(rx_mod, params = ipred_params, events = as_poso_data(full_plan, keep_meta = FALSE)))
sol_unadj <- as.data.frame(rxode2::rxSolve(rx_mod, params = ipred_params, events = as_poso_data(unadj_plan, keep_meta = FALSE)))

forecast_obj <- as_poso_map_tidy(
  ipred_df = sol_adj,
  pred_df  = sol_unadj,
  var_name = "Cc",
  obs_df   = data.frame(time = c(2.5, 35.5), dv = c(31.4, 18.6))
)

# 6. 默认输出组合图；也可通过 which = "map" 或 include = c("prior", "forecast") 自定义
bundle <- poso_plot_bundle(
  prior_sim     = sim_prior,
  map_tidy      = map_obj,
  forecast_tidy = forecast_obj,
  target_range  = c(10, 15)
)

plot(bundle)                       # 默认输出全部组合图
plot(bundle, which = c("map", "forecast")) # 按需自定义输出指定图
