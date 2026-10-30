#Loading packages
#---------------------------------------------------------
library(ggplot2)
library(gridExtra)
library(lubridate)
library(dplyr)
library(purrr)
library(tidyr)
library(sf)
library(terra)
library(tidyterra)
library(brms)

#Selecting numbers not in the exponential format
#---------------------------------------------------------
options(scipen = 999)

#Plotting crossings
#---------------------------------------------------------
load(file = "~/Downloads/ReproducibleDataset_20261029.Rdata")
head(d)

#---------------------------------------------------------
#Modeling step lengths
#---------------------------------------------------------
{
  #Histogram
  #---------------------------------------------------------
  ggplot() +
    geom_histogram(data = d, aes(x = sl), fill = "steelblue", alpha = 0.8) +
    theme_bw()
  ggplot() +
    geom_histogram(data = d, aes(x = log(sl)), fill = "steelblue", alpha = 0.8) +
    theme_bw()
  
  nrow(d |> dplyr::filter(sl <= 300))/nrow(d) #5% of steps within 300m
  
  #Choosing a distribution
  #---------------------------------------------------------
  mod.gamma <- brm(sl ~ 1,
                   data = d, 
                   family = Gamma(link = "log"),
                   chains = 4, cores = 4,
                   iter = 5000, warmup = 2000, thin = 10,
                   control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  mod.gev <- brm(sl ~ 1,
                 data = d, 
                 family = gen_extreme_value(link = "identity", link_sigma = "log", link_xi = "log1p"),
                 chains = 4, cores = 4,
                 iter = 5000, warmup = 2000, thin = 10,
                 control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  mod.skn <- brm(sl ~ 1,
                 data = d, 
                 family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"),
                 chains = 4, cores = 4,
                 iter = 5000, warmup = 2000, thin = 10,
                 control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  mod.shln <- brm(sl ~ 1,
                  data = d, 
                  family = shifted_lognormal(link = "identity", link_sigma = "log", link_ndt = "log"),
                  chains = 4, cores = 4,
                  iter = 5000, warmup = 2000, thin = 10,
                  control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  loo(mod.gamma) #-6886.6 ± 14.2
  loo(mod.gev) #-6772.2 ± 16.8
  loo(mod.skn) #-6775.4 ± 17.3
  loo(mod.shln) #-7177.9 ± 24.5
  
  mod.null <- brm(sl.std ~ 1,
                  data = d, 
                  family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"),
                  chains = 4, cores = 4,
                  iter = 5000, warmup = 2000, thin = 10,
                  control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  #Posterior predictive checks
  #---------------------------------------------------------
  pp_check(mod.gamma)
  pp_check(mod.gev)
  pp_check(mod.skn)
  pp_check(mod.shln)
  
  ggsave(
    plot = grid.arrange(
      pp_check(mod.skn),
      pp_check(mod.skn, type = "ecdf_overlay"), ncol = 2),
    filename = "~/Downloads/PPchecK_SL_20260809.pdf",
    width = 25, height = 12, units = "cm", dpi = 300)
  
  mod.null <- mod.skn
  
  #Model selection
  #---------------------------------------------------------
  mod.a.lm <- brm(sl.std ~ dist.front.abs.std + cross + dist.front.abs.std * cross,
                  data = d, 
                  family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"),
                  chains = 4, cores = 4,
                  iter = 5000, warmup = 2000, thin = 10,
                  control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  plot(conditional_effects(mod.a.lm))
  loo(mod.null) #-1012.7 ± 17.3
  loo(mod.a.lm) #-961.0 ± 18.9
  
  ggplot() +
    geom_point(aes(x = d$dist.front.abs, y = residuals(mod.a.lm)[,1])) +
    theme_bw()
  
  mod.a.lm.distr.int <- brm(bf(
    sl.std ~ dist.front.abs.std + cross + dist.front.abs.std * cross,
    sigma ~ dist.front.abs.std + cross + dist.front.abs.std * cross),
    data = d, 
    family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"),
    chains = 4, cores = 4,
    iter = 5000, warmup = 2000, thin = 10,
    control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  loo(mod.null) #-1012.7 ± 17.3
  loo(mod.a.lm) #-961.0 ± 18.9
  loo(mod.a.lm.distr.int) # -938.6 ± 18.6
  
  ggplot() +
    geom_point(aes(x = d$dist.front.abs, y = residuals(mod.a.lm.distr.int)[,1])) +
    theme_bw()
  plot(conditional_effects(mod.a.lm.distr.int), points = TRUE)
  
  mod.a.k3 <- brm(bf(
    sl.std ~ gp(dist.front.abs.std, by = cross, c = 5/4, k = 3),
    sigma ~ dist.front.abs.std + cross + dist.front.abs.std * cross),
    data = d, 
    family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"),
    chains = 4, cores = 4,
    iter = 5000, warmup = 2000, thin = 10,
    control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  loo(mod.null) #-1012.7 ± 17.3
  loo(mod.a.lm) #-961.0 ± 18.9
  loo(mod.a.lm.distr.int) # -938.6 ± 18.6
  loo(mod.a.k3) # -911.3 ± 19.5
  
  ggplot() +
    geom_point(aes(x = d$dist.front.abs, y = residuals(mod.a.k3)[,1])) +
    theme_bw()
  plot(conditional_effects(mod.a.k3), points = TRUE)
  
  mod.a.k5 <- brm(bf(
    sl.std ~ gp(dist.front.abs.std, by = cross, c = 5/4, k = 5),
    sigma ~ dist.front.abs.std + cross + dist.front.abs.std * cross),
    data = d, 
    family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"),
    chains = 4, cores = 4,
    iter = 5000, warmup = 2000, thin = 10,
    control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  loo(mod.null) #-1012.7 ± 17.3
  loo(mod.a.lm) #-961.0 ± 18.9
  loo(mod.a.lm.distr.int) # -938.6 ± 18.6
  loo(mod.a.k3) # -911.3 ± 19.5
  loo(mod.a.k5) # -904.1 ± 19.5
  ggplot() +
    geom_point(aes(x = d$dist.front.abs, y = residuals(mod.a.k5)[,1])) +
    theme_bw()
  plot(conditional_effects(mod.a.k5), points = TRUE)
  
  mod.a.k7 <- brm(bf(
    sl.std ~ gp(dist.front.abs.std, by = cross, c = 5/4, k = 7),
    sigma ~ dist.front.abs.std + cross + dist.front.abs.std * cross),
    data = d, 
    family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"),
    chains = 4, cores = 4,
    iter = 5000, warmup = 2000, thin = 10,
    control = list(adapt_delta = 0.99, max_treedepth = 11))
  loo(mod.null) #-1012.7 ± 17.3
  loo(mod.a.lm) #-961.0 ± 18.9
  loo(mod.a.lm.distr.int) # -938.6 ± 18.6
  loo(mod.a.k3) # -911.3 ± 19.5
  loo(mod.a.k5) # -904.1 ± 19.5
  loo(mod.a.k7) # -902.6 19.2
  ggplot() +
    geom_point(aes(x = d$dist.front.abs, y = residuals(mod.a.k7)[,1])) +
    theme_bw()
  plot(conditional_effects(mod.a.k7), points = TRUE)
  
  mod.a <- mod.a.k5
  
  mod.b <- brm(bf(
    sl.std ~ gp(dist.front.abs.std, by = cross, c = 5/4, k = 5),
    sigma ~ dist.front.abs.std + cross),
    data = d, 
    family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"),
    chains = 4, cores = 4,
    iter = 5000, warmup = 2000, thin = 10,
    control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  #bestmod.sl <- mod.b
  #save(bestmod.sl, file = "~/Downloads/BestMod_SL_20260810.Rdata")
  
  loo(mod.null) #-1012.7 ± 17.3
  loo(mod.a) # -904.1 ± 19.2
  loo(mod.b) # -902.4 ± 19.5
  ggplot() +
    geom_point(aes(x = d$dist.front.abs, y = residuals(mod.b)[,1])) +
    theme_bw()
  plot(conditional_effects(mod.a), points = TRUE)
  
  mod.c <- brm(bf(
    sl.std ~ gp(dist.front.abs.std, by = cross, c = 5/4, k = 5),
    sigma ~ dist.front.abs.std),
    data = d, 
    family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"),
    chains = 4, cores = 4,
    iter = 5000, warmup = 2000, thin = 10,
    control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  loo(mod.null) #-1012.7 ± 17.3
  loo(mod.a) # -904.1 ± 19.2
  loo(mod.b) # -902.4 ± 19.5
  loo(mod.c) # -917.5 ± 20.2
  ggplot() +
    geom_point(aes(x = d$dist.front.abs, y = residuals(mod.c)[,1])) +
    theme_bw()
  plot(conditional_effects(mod.c), points = TRUE)
  
  mod.d <- brm(bf(
    sl.std ~ gp(dist.front.abs.std, c = 5/4, k = 5),
    sigma ~ dist.front.abs.std + cross),
    data = d, 
    family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"),
    chains = 4, cores = 4,
    iter = 5000, warmup = 2000, thin = 10,
    control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  loo(mod.null) #-1012.7 ± 17.3
  loo(mod.a) # -904.1 ± 19.2
  loo(mod.b) # -902.4 ± 19.5
  loo(mod.c) # -917.5 ± 20.2
  loo(mod.d) # -930.7 ± 19.1
  
  #Plotting ELPDs
  #---------------------------------------------------------
  d.loo <- setNames(as.data.frame(cbind(c("mod.null", "mod.a", "mod.b", "mod.c", "mod.d"),
                                        c(-1012.7, -904.1, -902.4, -917.5, -930.7),
                                        c(17.3, 19.2, 19.5, 20.2, 19.1))),
                    c("model", "elpd", "se"))
  
  d.loo <- d.loo |> mutate(elpd = as.numeric(as.character(elpd)),
                           se = as.numeric(as.character(se)),
                           model = factor(model, levels = rev(c("mod.null", "mod.a", "mod.b", "mod.c", "mod.d"))),
                           lower = (elpd-se), upper = (elpd + se))
  ggsave(
    plot = 
      ggplot() +
      geom_errorbar(data = d.loo, aes(x = model, ymin = lower, ymax = upper, width = 0.05)) +
      geom_point(data = d.loo, aes(x = model, y = elpd), shape = 21, size = 2, fill = "white") +
      geom_point(data = d.loo, aes(x = model, y = elpd), shape = 21, size = 2, fill = "steelblue", alpha = 0.6) +
      scale_x_discrete(labels = rev(c(
        "~ 1",
        "mean ~ gp(distance, by = crossing, c = 5/4, k = 7),\nsigma ~ distance + crossing + distance * crossing)",
        "mean ~ gp(distance, by = crossing, c = 5/4, k = 5),\nsigma ~ distance + crossing\n(Best candidate model)",
        "mean ~ gp(distance, by = crossing, c = 5/4, k = 5),\nsigma ~ distance",
        "mean ~ gp(distance, c = 5/4, k = 5),\nsigma ~ distance + crossing"
      ))) +
      coord_flip() +
      labs(y = "ELPD") +
      theme_bw() +
      theme(axis.title.y = element_blank()),
    filename = "~/Downloads/ModelSelection_SL_20260810.pdf",
    width = 25, height = 12, units = "cm", dpi = 300)
  
  #Inspecting the best candidate model
  #---------------------------------------------------------
  bestmod <- mod.b
  bayes_R2(mod.b) #23.1% variation in the response
  acf(residuals(mod.b)[,1]) #Significant temporal autocorrelation
  
  ggsave(
    plot = grid.arrange(
      ggplot() +
        geom_point(aes(x = predict(mod.b)[,1], y = residuals(mod.b)[,1]), color = "steelblue", alpha = 0.8) +
        labs(x = "Fitted values", y = "Model residuals") +
        theme_bw(),
      ggplot() +
        geom_point(aes(x = d$dist.front.abs/1000, y = residuals(mod.b)[,1]), color = "steelblue", alpha = 0.8) +
        labs(x = "Distance from the front line (km)", y = "Model residuals") +
        theme_bw(), ncol = 2),
    filename = "~/Downloads/Residuals_SL_20260810.pdf",
    width = 25, height = 12, units = "cm", dpi = 300)
  
  #Plotting marginal effects
  #---------------------------------------------------------
  df.condeff <- conditional_effects(bestmod, effects = "dist.front.abs.std",
                                    conditions = data.frame("cross"= c("1st", "2nd", "3rd")))$dist.front.abs.std
  head(df.condeff)
  df.condeff <- df.condeff |> dplyr::select(effect1__, cond__, estimate__, se__, lower__, upper__) |>
    dplyr::rename("dist.front.abs.std" = "effect1__", "cross" = "cond__", "est" = "estimate__",
                  "se" = "se__", "lower" = "lower__", "upper" = "upper__")
  
  ggsave(
    plot = grid.arrange(
      ggplot() +
        geom_ribbon(data = df.condeff |> dplyr::filter(cross == 1), aes(x = dist.front.abs.std, ymin = lower, ymax = upper),
                    linewidth = 0.1, fill = "firebrick", color = "black", alpha = 0.2) +
        geom_point(data = d |> dplyr::filter(cross == "1st"), aes(x = dist.front.abs.std, y = sl.std),
                   size = 0.4, color = "steelblue", alpha = 0.8) +
        geom_line(data = df.condeff |> dplyr::filter(cross == 1), aes(x = dist.front.abs.std, y = est), linetype = "dashed",
                  color = "firebrick", linewidth = 0.7) +
        labs(x = "Distance from the front line (Z-score)", y = "Step length (Z-score)") +
        theme_bw(),
      
      ggplot() +
        geom_ribbon(data = df.condeff |> dplyr::filter(cross == 2), aes(x = dist.front.abs.std, ymin = lower, ymax = upper),
                    linewidth = 0.1, fill = "firebrick", color = "black", alpha = 0.2) +
        geom_point(data = d |> dplyr::filter(cross == "2nd"), aes(x = dist.front.abs.std, y = sl.std),
                   size = 0.4, color = "steelblue", alpha = 0.8) +
        geom_line(data = df.condeff |> dplyr::filter(cross == 2), aes(x = dist.front.abs.std, y = est), linetype = "dashed",
                  color = "firebrick", linewidth = 0.7) +
        labs(x = "Distance from the front line (Z-score)", y = "Step length (Z-score)") +
        scale_x_continuous(limits = c(-1.5, 1.5)) +
        theme_bw(),
      
      ggplot() +
        geom_ribbon(data = df.condeff |> dplyr::filter(cross == 3), aes(x = dist.front.abs.std, ymin = lower, ymax = upper),
                    linewidth = 0.1, fill = "firebrick", color = "black", alpha = 0.2) +
        geom_point(data = d |> dplyr::filter(cross == "3rd"), aes(x = dist.front.abs.std, y = sl.std),
                   size = 0.4, color = "steelblue", alpha = 0.8) +
        geom_line(data = df.condeff |> dplyr::filter(cross == 3), aes(x = dist.front.abs.std, y = est), linetype = "dashed",
                  color = "firebrick", linewidth = 0.7) +
        labs(x = "Distance from the front line (Z-score)", y = "Step length (Z-score)") +
        scale_x_continuous(limits = c(-1.5, 1.5)) +
        theme_bw(), ncol = 3),
    filename = "~/Downloads/CondEff_SL_20260810.pdf",
    width = 30, height = 10, units = "cm", dpi = 300)
}

#---------------------------------------------------------
#Modeling flying height
#---------------------------------------------------------
{
  dh <- d[which(!is.na(d$h)),]
  
  #Histogram
  #---------------------------------------------------------
  ggplot() +
    geom_histogram(data = dh, aes(x = h), fill = "steelblue", alpha = 0.8) +
    theme_bw()
  ggplot() +
    geom_histogram(data = dh, aes(x = log(h)), fill = "steelblue", alpha = 0.8) +
    theme_bw()
  
  #Choosing a distribution
  #---------------------------------------------------------
  mod.gev <- brm(h.std ~ 1,
                 data = dh, 
                 family = gen_extreme_value(link = "identity", link_sigma = "log", link_xi = "log1p"),
                 chains = 4, cores = 4,
                 iter = 5000, warmup = 2000, thin = 10,
                 control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  mod.skn <- brm(h.std ~ 1,
                 data = dh, 
                 family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"),
                 chains = 4, cores = 4,
                 iter = 5000, warmup = 2000, thin = 10,
                 control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  loo(mod.gev) #-869.9 ± 18.8
  loo(mod.skn) #-866.8 ± 18.1
  mod.null <- mod.skn
  
  #Posterior predictive checks
  #---------------------------------------------------------
  pp_check(mod.skn)
  ggsave(
    plot = grid.arrange(
      pp_check(mod.skn),
      pp_check(mod.skn, type = "ecdf_overlay"), ncol = 2),
    filename = "~/Downloads/PPchecK_Height_20260809.pdf",
    width = 25, height = 12, units = "cm", dpi = 300)
  
  #Model selection
  #---------------------------------------------------------
  mod.a.lm <- brm(h.std ~ dist.front.abs.std + cross + dist.front.abs.std * cross,
                  data = dh, 
                  family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"),
                  chains = 4, cores = 4,
                  iter = 5000, warmup = 2000, thin = 10,
                  control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  plot(conditional_effects(mod.a.lm))
  loo(mod.null) #-866.8 ± 18.1
  loo(mod.a.lm) #-829.1 ± 39.3
  
  ggplot() +
    geom_point(aes(x = dh$dist.front.abs, y = residuals(mod.a.lm)[,1])) +
    theme_bw()
  
  plot(conditional_effects(mod.a.lm), points = TRUE)
  
  mod.a <- mod.a.lm
  
  mod.b <- brm(h.std ~ dist.front.abs.std + cross,
               data = dh, 
               family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"),
               chains = 4, cores = 4,
               iter = 5000, warmup = 2000, thin = 10,
               control = list(adapt_delta = 0.99, max_treedepth = 10))
  loo(mod.null) #-866.8 ± 18.1
  loo(mod.a.lm) #-829.1 ± 39.3
  loo(mod.b) #-833.6 ± 19.4
  
  mod.c <- brm(h.std ~ dist.front.abs.std,
               data = dh, 
               family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"),
               chains = 4, cores = 4,
               iter = 5000, warmup = 2000, thin = 10,
               control = list(adapt_delta = 0.99, max_treedepth = 10))
  loo(mod.null) #-866.8 ± 18.1
  loo(mod.a.lm) #-829.1 ± 39.3
  loo(mod.b) #-833.6 ± 19.4
  loo(mod.c) #-865.0 ± 18.2
  
  mod.d <- brm(h.std ~ cross,
               data = dh, 
               family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"),
               chains = 4, cores = 4,
               iter = 5000, warmup = 2000, thin = 10,
               control = list(adapt_delta = 0.99, max_treedepth = 10))
  loo(mod.null) #-866.8 ± 18.1
  loo(mod.a.lm) #-829.1 ± 39.3
  loo(mod.b) #-833.6 ± 19.4
  loo(mod.c) #-865.0 ± 18.2
  loo(mod.d) #-833.2 ± 19.4
  bayes_R2(mod.d) #10% of variability in flying height
  
  #Best candidate model
  #---------------------------------------------------------
  bestmod <- mod.d
  grid.arrange(
    ggplot() +
      geom_point(aes(x = predict(mod.b)[,1], y = residuals(mod.b)[,1]), color = "steelblue", alpha = 0.8) +
      labs(x = "Fitted values", y = "Model residuals") +
      theme_bw(),
    ggplot() +
      geom_boxplot(aes(x = dh$cross, y = residuals(mod.b)[,1])) +
      labs(x = "Crossings of Ukraine", y = "Model residuals") +
      theme_bw(), ncol = 2)
  
  #Plotting ELPD and conditional effects
  #---------------------------------------------------------
  d.loo <- setNames(as.data.frame(cbind(c("mod.null", "mod.a", "mod.b", "mod.c", "mod.d"),
                                        c(-866.8, -829.1, -833.6, -865.0, -833.2 ),
                                        c(18.1, 39.3, 19.4, 18.2, 19.4))),
                    c("model", "elpd", "se"))
  
  d.loo <- d.loo |> mutate(elpd = as.numeric(as.character(elpd)),
                           se = as.numeric(as.character(se)),
                           model = factor(model, levels = rev(c("mod.null", "mod.a", "mod.b", "mod.c", "mod.d"))),
                           lower = (elpd-se), upper = (elpd + se))
  
  bestmod <- mod.d
  df.condeff <- conditional_effects(bestmod, effects = "cross")$cross
  
  df.condeff <- df.condeff |> dplyr::select(cross, estimate__, se__, lower__, upper__) |>
    dplyr::rename("est" = "estimate__",
                  "se" = "se__", "lower" = "lower__", "upper" = "upper__")
  
  ggsave(
    plot = 
      grid.arrange(
        ggplot() +
          geom_errorbar(data = d.loo, aes(x = model, ymin = lower, ymax = upper, width = 0.05)) +
          geom_point(data = d.loo, aes(x = model, y = elpd), shape = 21, size = 2, fill = "white") +
          geom_point(data = d.loo, aes(x = model, y = elpd), shape = 21, size = 2, fill = "steelblue", alpha = 0.6) +
          scale_x_discrete(labels = rev(c(
            "~ 1",
            "mean ~ distance + crossing +\ndistance * crossing",
            "mean ~ distance + crossing",
            "mean ~ distance",
            "mean ~ crossing"
          ))) +
          coord_flip() +
          labs(y = "ELPD") +
          theme_bw() +
          theme(axis.title.y = element_blank()),
        ggplot() +
          geom_errorbar(data = df.condeff, aes(x = cross, ymin = lower, ymax = upper), width = 0.1) +
          geom_point(data = df.condeff, aes(x = cross, y = est), shape = 21, size = 2, fill = "white") +
          geom_point(data = df.condeff, aes(x = cross, y = est), shape = 21, size = 2, fill = "steelblue", alpha = 0.6) +
          labs(x = "Crossing", y = "Estimated flying height (Z-scores)") +
          theme_bw()
        , ncol = 2),
    filename = "~/Downloads/ModelSelection_Height_20260810.pdf",
    width = 35, height = 12, units = "cm", dpi = 300)
}

#---------------------------------------------------------
#Modeling turning angles
#---------------------------------------------------------
{
  #Histogram
  #---------------------------------------------------------
  ggplot() +
    geom_histogram(data = d, aes(x = ta), fill = "steelblue", alpha = 0.8) +
    scale_x_continuous(breaks = seq(from = -180, to = 180, by = 40)) +
    theme_bw()
  
  #Testing a null model
  #---------------------------------------------------------
  mod.null.vm <- brm(ta.rad ~ 1,
                     data = d, 
                     family = von_mises(),
                     chains = 4, cores = 4,
                     iter = 5000, warmup = 2000, thin = 10,
                     control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  mod.null.gaus <- brm(ta| trunc(lb = -180, ub = 180) ~ 1,
                       data = d, 
                       family = "Gaussian",
                       chains = 4, cores = 4,
                       iter = 5000, warmup = 2000, thin = 10,
                       control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  loo(mod.null.vm) #-842.9 ± 29.7
  loo(mod.null.gaus) #-3807.1 ± 30.3
  
  #Posterior predictive checks
  #---------------------------------------------------------
  ggsave(
    plot = grid.arrange(
      pp_check(mod.null.gaus),
      pp_check(mod.null.gaus, type = "ecdf_overlay"), ncol = 2),
    filename = "~/Downloads/PPchecK_TA_20260809.pdf",
    width = 25, height = 12, units = "cm", dpi = 300)
  
  grid.arrange(
    pp_check(mod.null.gaus),
    pp_check(mod.null.gaus, type = "ecdf_overlay"), ncol = 2) #pp-checks for a Gaussian model are not the best
  
  mod.null <- mod.null.vm
  
  #Model selection
  #---------------------------------------------------------
  mod.a.lm <- brm(ta.rad ~ dist.front.abs.std + cross + dist.front.abs.std * cross,
                  data = d, 
                  family = von_mises(link = "tan_half", link_kappa = "log"), #kappa is the inverse of variance
                  chains = 4, cores = 4,
                  iter = 5000, warmup = 2000, thin = 10,
                  control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  plot(d$dist.front.abs.std, residuals(mod.a.lm)[,1])
  plot(conditional_effects(mod.a.lm))
  loo(mod.null.vm) #-842.9 ± 29.7
  loo(mod.a.lm) #-846.5 ± 29.4 No effect on the mean
  
  mod.a.lm.distr.int <- brm(bf(
    ta.rad ~ dist.front.abs.std + cross + dist.front.abs.std * cross,
    kappa ~ dist.front.abs.std + cross + dist.front.abs.std * cross),
    data = d, 
    family = von_mises(link = "tan_half", link_kappa = "log"), #kappa is the inverse of variance
    chains = 4, cores = 4,
    iter = 5000, warmup = 2000, thin = 10,
    control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  loo(mod.null.vm) #-842.9 ± 29.7
  loo(mod.a.lm) #-846.3 ± 29.4
  loo(mod.a.lm.distr.int) #-841.8 ± 30.7
  
  mod.b <- brm(bf(
    ta.rad ~ dist.front.abs.std + cross ,
    kappa ~ dist.front.abs.std + cross + dist.front.abs.std * cross),
    data = d, 
    family = von_mises(link = "tan_half", link_kappa = "log"), #kappa is the inverse of variance
    chains = 4, cores = 4,
    iter = 5000, warmup = 2000, thin = 10,
    control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  loo(mod.null.vm) #-842.9 ± 29.7
  loo(mod.a.lm) #-846.3 ± 29.4
  loo(mod.a.lm.distr.int) #-841.8 ± 30.7
  loo(mod.b) #-840.7 ± 30.7
  mod.b$formula
  plot(conditional_effects(mod.b), points = TRUE)
  
  mod.c <- brm(bf(
    ta.rad ~ dist.front.abs.std ,
    kappa ~ dist.front.abs.std + cross + dist.front.abs.std * cross),
    data = d, 
    family = von_mises(link = "tan_half", link_kappa = "log"), #kappa is the inverse of variance
    chains = 4, cores = 4,
    iter = 5000, warmup = 2000, thin = 10,
    control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  loo(mod.null.vm) #-842.9 ± 29.7
  loo(mod.a.lm) #-846.3 ± 29.4
  loo(mod.a.lm.distr.int) #-841.8 ± 30.7
  loo(mod.b) #-840.7 ± 30.7
  loo(mod.c) #-839.2 ± 30.9
  plot(conditional_effects(mod.b), points = TRUE)
  
  mod.d <- brm(bf(
    ta.rad ~ 1,
    kappa ~ dist.front.abs.std + cross + dist.front.abs.std * cross),
    data = d, 
    family = von_mises(link = "tan_half", link_kappa = "log"), #kappa is the inverse of variance
    chains = 4, cores = 4,
    iter = 5000, warmup = 2000, thin = 10,
    control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  loo(mod.null.vm) #-842.9 ± 29.7
  loo(mod.a.lm) #-846.3 ± 29.4
  loo(mod.a.lm.distr.int) #-841.8 ± 30.7
  loo(mod.b) #-840.7 ± 30.7
  loo(mod.c) #-839.2 ± 30.9
  loo(mod.d) #-838.0 ± 30.8
  plot(conditional_effects(mod.d), points = TRUE)
  
  mod.e <- brm(bf(
    ta.rad ~ 1,
    kappa ~ dist.front.abs.std + cross),
    data = d, 
    family = von_mises(link = "tan_half", link_kappa = "log"), #kappa is the inverse of variance
    chains = 4, cores = 4,
    iter = 5000, warmup = 2000, thin = 10,
    control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  loo(mod.null.vm) #-842.9 ± 29.7
  loo(mod.a.lm) #-846.3 ± 29.4
  loo(mod.a.lm.distr.int) #-841.8 ± 30.7
  loo(mod.b) #-840.7 ± 30.7
  loo(mod.c) #-839.2 ± 30.9
  loo(mod.d) #-838.0 ± 30.8
  loo(mod.e) #-842.4 ± 30.4
  plot(conditional_effects(mod.d), points = TRUE)
  
  mod.f <- brm(bf(
    ta.rad ~ 1,
    kappa ~ dist.front.abs.std),
    data = d, 
    family = von_mises(link = "tan_half", link_kappa = "log"), #kappa is the inverse of variance
    chains = 4, cores = 4,
    iter = 5000, warmup = 2000, thin = 10,
    control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  loo(mod.null.vm) #-842.9 ± 29.7
  loo(mod.a.lm) #-846.3 ± 29.4
  loo(mod.a.lm.distr.int) #-841.8 ± 30.7
  loo(mod.b) #-840.7 ± 30.7
  loo(mod.c) #-839.2 ± 30.9
  loo(mod.d) #-838.0 ± 30.8
  loo(mod.e) #-842.4 ± 30.4
  loo(mod.f) #-840.1  ± 30.3
  loo(mod.g) #-844.3 ± 30.0
  plot(conditional_effects(mod.d), points = TRUE)
  
  mod.g <- brm(bf(
    ta.rad ~ 1,
    kappa ~ cross),
    data = d, 
    family = von_mises(link = "tan_half", link_kappa = "log"), #kappa is the inverse of variance
    chains = 4, cores = 4,
    iter = 5000, warmup = 2000, thin = 10,
    control = list(adapt_delta = 0.99, max_treedepth = 10))
  
  bestmod <- mod.d
  summary(mod.d)
}
