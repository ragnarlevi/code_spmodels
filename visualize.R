library(tidyverse)
source("simulate_data.R")



########## Poisson --------------
Poisson_t <- plot_per_t("Poisson_ln", 1:10, 0.4)


Poisson_t$beta1_plot + ylim(c(-0.5, 1.2)) + theme_minimal(base_size = 35) +labs(title = expression(beta[1] ~ " over time by model")) + theme(
  legend.position     = c(0.60, 0.30),    # inside panel, 2% from left & bottom
  legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
  legend.background    = element_rect(fill = alpha("white", 0.5))
)

Poisson_t$beta_error_plot + ylim(c(-0.5, 1.2)) + theme_minimal(base_size = 35) +labs(title = expression(l[1] ~ " error for " ~  beta ~  " over time by model")) + theme(
  legend.position     = c(0.60, 0.70),    # inside panel, 2% from left & bottom
  legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
  legend.background    = element_rect(fill = alpha("white", 0.5))
)+  scale_y_continuous(trans = "log2")



Poisson_t$edge_plot + ylim(c(-0.5, 3)) + theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$a_{4,1}$ estimate over time by model"),
                                                                               y = latex2exp::TeX("$a_{4,1}$ estimate")) + 
  theme(
  legend.position     = c(0.60, 0.70),    # inside panel, 2% from left & bottom
  legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
  legend.background    = element_rect(fill = alpha("white", 0.5))
)

Poisson_t$edge_error_plot + ylim(c(0, 5)) + theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$l_1$ error for $A$ over time by model"),
                                                                               y = latex2exp::TeX("$l_{1}$ error")) + 
  theme(
    legend.position     = c(0.60, 0.60),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  ) +  scale_y_continuous(trans = "log2")


Poisson_t$beta2_plot + ylim(c(0, 0.15)) + theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$\\beta_2$ estimate over time by model"),
                                                                           y = latex2exp::TeX("$\\beta_2$ estimate")) + 
  theme(
    legend.position     = c(0.60, 0.60),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  )


Poisson_t$edge_plot

Poisson_lambda <- plot_per_lambda("Poisson_ln", 1:9, 0.4)
Poisson_lambda$f1_plot+ theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$\\f_1$ score over $\\lambda$ by model, time, and time points"),
                                                            y = latex2exp::TeX("$$\\f_1$ score")) + 
  theme(
    legend.position     = c(0.25, 0.75),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  )


Poisson_lambda$side_by_side+ theme_minimal(base_size = 38) +labs(x = "", y="")


Poisson_lambda <- plot_per_lambda("Poisson", 1:9, 0.8)
Poisson_lambda$f1_plot+ theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$\\f_1$ score over $\\lambda$ by model, time, and time points"),
                                                            y = latex2exp::TeX("$$\\f_1$ score")) + 
  theme(
    legend.position     = c(0.25, 0.75),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  )


# Poisson_lambda$side_by_side+ theme_minimal(base_size = 38) +labs(x = "", y="")

####
## Poisson plot error ------
###


# Chose id witih se error calculation
tmp <- new.env()
loaded_objs <- load(paste0("sim_data/", "Poisson", "_ln_per_t", 0.4, "id", 4, ".RData") , envir = tmp)
my_list <- mget(loaded_objs, envir = tmp)

# select index of a to plot error, here I choose index 1 and 8
poisson_error <- plot_error(my_list, 1, 8, is_ordinary = TRUE)
poisson_error$beta_plot


poisson_error$a_plot +
  scale_color_brewer("Parameter", palette = "Set2",
                     labels = c(latex2exp::TeX("$a_{1,1}$"), latex2exp::TeX("$a_{4,2}$"))) +
  scale_fill_brewer("Parameter", palette = "Set2",
                    labels = c(latex2exp::TeX("$a_{1,1}$"), latex2exp::TeX("$a_{4,2}$"))) +
  labs(
    x     = expression(log(t)),
    y     = "Estimate ± 2 SE",
    title = latex2exp::TeX("$a_{1,1}$ and  $a_{4,2}$ over time")
  ) +
  theme_minimal(base_size = 36) +
  theme(
    plot.title      = element_text(hjust = 0.5),
    legend.position = c(0.8, 0.6),
    legend.background = element_rect(fill = "white")
  )



######## Lambda Poisson ########

library(tidyverse)
source("simulate_data.R")


# change density to 0.2, 0.4 and 0.8
poisson_lambda <- plot_per_lambda("Poisson_ln", 1:9, 0.4, A_compare = "100-10" )
poisson_lambda$f1_plot+ theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$\\f_1$ score over $\\lambda$ by model, time, and time points"),
                                                        y = latex2exp::TeX("$$\\f_1$ score")) + 
  theme(
    legend.position     = c(0.68, 0.05),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  ) + xlim(c(0, 50)) + theme(legend.position = "none")


poisson_lambda$side_by_side+ theme_minimal(base_size = 38) +labs(x = "", y="")

poisson_lambda$true_false_graph + theme_minimal(base_size = 32) +labs(x = "", y="")

######### ZIP --------------
source("simulate_data.R")

zip_t <- plot_per_t("zip_ig", 1:10, 0.8)

zip_t$beta1_plot+ ylim(c(-0.5, 1.2)) + theme_minimal(base_size = 35) +labs(title = expression(beta[1] ~ " over time by model")) + theme(
  legend.position     = c(0.60, 0.30),    # inside panel, 2% from left & bottom
  legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
  legend.background    = element_rect(fill = alpha("white", 0.5))
)

zip_t$beta_error_plot  + ylim(c(-0.5, 15)) + theme_minimal(base_size = 35) +labs(title = expression(l[1] ~ " error for " ~  beta ~  " over time by model")) + theme(
  legend.position     = c(0.60, 0.70),    # inside panel, 2% from left & bottom
  legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
  legend.background    = element_rect(fill = alpha("white", 0.5))
)+  scale_y_continuous(trans = "log2")

zip_t$edge_plot + ylim(c(-0.5, 4)) + theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$a_{4,1}$ estimate over time by model"),
                                                                         y = latex2exp::TeX("$a_{4,1}$ estimate")) + 
  theme(
    legend.position     = c(0.60, 0.70),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  )


 
zip_t$edge_error_plot  + ylim(c(0, 5)) + theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$l_1$ error for $A$ over time by model"),
                                                                             y = latex2exp::TeX("$l_{1}$ error")) + 
  theme(
    legend.position     = c(0.60, 0.60),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  ) +  scale_y_continuous(trans = "log2")



zip_t$beta2_plot  + theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$\\beta_2$ estimate over time by model"),
                                                                              y = latex2exp::TeX("$\\beta_2$ estimate")) + 
  theme(
    legend.position     = c(0.60, 0.60),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  )






####
## ZIP plot error ------
###

tmp <- new.env()
loaded_objs <- load(paste0("sim_data/", "zip", "_ig_per_t", 0.4, "id", 1, ".RData") , envir = tmp)
my_list <- mget(loaded_objs, envir = tmp)

# Here index of a nr 2 and 7 are chosen
zip_error <- plot_error(my_list, 2, 7,is_ordinary = FALSE)
zip_error$beta_plot


zip_error$a_plot +
  scale_color_brewer("Parameter", palette = "Set2",
                     labels = c(latex2exp::TeX("$a_{2,1}$"), latex2exp::TeX("$a_{3,1}$"))) +
  scale_fill_brewer("Parameter", palette = "Set2",
                    labels = c(latex2exp::TeX("$a_{2,1}$"), latex2exp::TeX("$a_{3,1}$"))) +
  labs(
    x     = expression(log(t)),
    y     = "Estimate ± 2 SE",
    title = latex2exp::TeX("$a_{2,1}$ and  $a_{3,1}$ over time")
  ) +
  theme_minimal(base_size = 36) +
  theme(
    plot.title      = element_text(hjust = 0.5),
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(fill = "white")
  )


######## Lambda ZIP ########

# Select density 0.2, 0.4, 0.8
zip_lambda <- plot_per_lambda("2zip_ig", 1:10, 0.4, A_compare = "50-0.5" )
zip_lambda$f1_plot+ theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$\\f_1$ score over $\\lambda$ by model, time, and time points"),
                                                            y = latex2exp::TeX("$$\\f_1$ score")) + 
  theme(
    legend.position     = c(0.78, 0.45),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  ) + theme(legend.position = "none")


zip_lambda$side_by_side+ theme_minimal(base_size = 38) +labs(x = "", y="")

zip_lambda$true_false_graph + theme_minimal(base_size = 32) +labs(x = "", y="")



######### Hurdle --------------
library(tidyverse)
source("simulate_data.R")

hurdle_t <- plot_per_t("3hurdle_ig", 1:10, 0.4)

hurdle_t$beta1_plot+ ylim(c(-0.5, 1.2)) + theme_minimal(base_size = 35) +labs(title = expression(beta[1] ~ " over time by model")) + theme(
  legend.position     = c(0.60, 0.30),    # inside panel, 2% from left & bottom
  legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
  legend.background    = element_rect(fill = alpha("white", 0.5))
)

hurdle_t$beta_error_plot  + ylim(c(-0.5, 15)) + theme_minimal(base_size = 35) +labs(title = expression(l[1] ~ " error for " ~  beta ~  " over time by model")) + theme(
  legend.position     = c(0.60, 0.70),    # inside panel, 2% from left & bottom
  legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
  legend.background    = element_rect(fill = alpha("white", 0.5))
)+  scale_y_continuous(trans = "log2")

hurdle_t$edge_plot + ylim(c(-0.5, 4)) + theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$a_{4,1}$ estimate over time by model"),
                                                                         y = latex2exp::TeX("$a_{4,1}$ estimate")) + 
  theme(
    legend.position     = c(0.60, 0.70),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  )



hurdle_t$edge_error_plot  + ylim(c(0, 5)) + theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$l_1$ error for $A$ over time by model"),
                                                                             y = latex2exp::TeX("$l_{1}$ error")) + 
  theme(
    legend.position     = c(0.60, 0.60),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  ) +  scale_y_continuous(trans = "log2")



hurdle_t$beta2_plot  + theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$\\beta_2$ estimate over time by model"),
                                                        y = latex2exp::TeX("$\\beta_2$ estimate")) + 
  theme(
    legend.position     = c(0.60, 0.60),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  )


######## Lambda Hurdle ########
library(tidyverse)
source("simulate_data.R")

# select density 
hurdle_lambda <- plot_per_lambda("3hurdle_ig", 1:10, 0.2, A_compare = "200-15" )
hurdle_lambda$f1_plot+ theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$\\f_1$ score over $\\lambda$ by model, time, and time points"),
                                                            y = latex2exp::TeX("$$\\f_1$ score")) + 
  theme(
    legend.position     = c(0.78, 0.45),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  )+theme(legend.position = "none") #+ xlim(c(0, 40)) 


hurdle_lambda$side_by_side+ theme_minimal(base_size = 38) +labs(x = "", y="")
hurdle_lambda$true_false_graph  + theme_minimal(base_size = 32) +labs(x = "", y="")






########  Hurdle error plit ########

tmp <- new.env()
loaded_objs <- load(paste0("sim_data/", "3hurdle", "_ig_per_t", 0.4, "id", 1, ".RData") , envir = tmp)
my_list <- mget(loaded_objs, envir = tmp)


hurdle_error <- plot_error(my_list, 3, 9,is_ordinary = FALSE)
hurdle_error$beta_plot
hurdle_error$a_plot+scale_color_brewer("Parameter", palette = "Set2",
                     labels = c(latex2exp::TeX("$a_{2,2}$"), latex2exp::TeX("$a_{3,4}$"))) +
  scale_fill_brewer("Parameter", palette = "Set2",
                    labels = c(latex2exp::TeX("$a_{2,2}$"), latex2exp::TeX("$a_{3,4}$"))) +
  labs(
    x     = expression(log(t)),
    y     = "Estimate ± 2 SE",
    title = latex2exp::TeX("$a_{2,2}$ and  $a_{3,4}$ over time")
  ) +
  theme_minimal(base_size = 36) +
  theme(
    plot.title      = element_text(hjust = 0.5),
    legend.position = c(0.8, 0.9),
    legend.background = element_rect(fill = "white")
  )


######### ZIP PSI --------------
library(tidyverse)
source("simulate_data.R")

# Note error going slowly to zero fast as the latent effects are weak, y_latent has many zeros.

zip_t <- plot_per_t_psi("5zip_psi_ig", 1:10, 0.4)

zip_t$beta1_plot+ ylim(c(-0.5, 1.2)) + theme_minimal(base_size = 35) +labs(title = expression(beta[1] ~ " over time by model")) + theme(
  legend.position     = c(0.60, 0.30),    # inside panel, 2% from left & bottom
  legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
  legend.background    = element_rect(fill = alpha("white", 0.5))
)

zip_t$beta_error_plot  + ylim(c(-0.5, 15)) + theme_minimal(base_size = 35) +labs(title = expression(l[1] ~ " error for " ~  beta ~  " over time by model")) + theme(
  legend.position     = c(0.60, 0.70),    # inside panel, 2% from left & bottom
  legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
  legend.background    = element_rect(fill = alpha("white", 0.5))
)+  scale_y_continuous(trans = "log2")

zip_t$edge_plot + ylim(c(-0.5, 4)) + theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$a_{4,1}$ estimate over time by model"),
                                                                         y = latex2exp::TeX("$a_{4,1}$ estimate")) + 
  theme(
    legend.position     = c(0.60, 0.70),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  )



zip_t$edge_error_plot  + ylim(c(0, 5)) + theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$l_1$ error for $\\phi$ over time by model"),
                                                                             y = latex2exp::TeX("$l_{1}$ error")) + 
  theme(
    legend.position     = c(0.60, 0.60),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  ) +  scale_y_continuous(trans = "log2")



zip_t$beta2_plot  + theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$\\beta_2$ estimate over time by model"),
                                                        y = latex2exp::TeX("$\\beta_2$ estimate")) + 
  theme(
    legend.position     = c(0.60, 0.60),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  )



# ZIP PSI ERROR -------


tmp <- new.env()
loaded_objs <- load(paste0("sim_data/", "zip_psi", "_ln_per_t", 0.4, "id", 1, ".RData") , envir = tmp)
my_list <- mget(loaded_objs, envir = tmp)

zip_psi_error <- plot_error(my_list, 1, 2,is_ordinary = FALSE)
zip_psi_error$beta_plot
zip_psi_error$a_plot+   scale_color_brewer("Parameter", palette = "Set2",
                     labels = c(latex2exp::TeX("$\\phi_{1}$"), latex2exp::TeX("$\\phi_{2}$"))) +
  scale_fill_brewer("Parameter", palette = "Set2",
                    labels = c(latex2exp::TeX("$\\phi_{1}$"), latex2exp::TeX("$\\phi_{2}$"))) +
  labs(
    x     = expression(log(t)),
    y     = "Estimate ± 2 SE",
    title = latex2exp::TeX("$\\phi_{1}$ and  $\\phi_{2}$ over time")
  ) +
  theme_minimal(base_size = 36) +
  theme(
    plot.title      = element_text(hjust = 0.5),
    legend.position = c(0.8, 0.6),
    legend.background = element_rect(fill = "white")
  )



######### Speed --------------

library(tidyverse)

tmp <- new.env()
loaded_objs <- load(paste0("sim_data/2zip_speed0.4id1.RData") , envir = tmp)
my_list <- mget(loaded_objs, envir = tmp)



parts_list <- strsplit(names(my_list$time_ig), "-", fixed = TRUE)
info <- do.call(rbind, parts_list)
secs_vec <- sapply(my_list$time_ig, function(d) {
  as.numeric(d, units = "secs")
})


data_to_plot <- data.frame(nr_data_points = as.numeric(unlist(my_list$nr_data_points)), nr_param  = info[, 2], time = as.numeric(secs_vec))

library(scales)
data_to_plot %>%
  ggplot(aes(
    x = as.numeric(nr_param),
    y = time,
    color = factor(nr_data_points)
  )) +
  # Draw lines thicker and add points
  geom_line(size = 1.3) +
  geom_point(size = 3) +
  scale_color_discrete(
    name   = "Number of\nData Points",
    labels = function(x) comma(as.numeric(x))
  )+
  # Use a viridis discrete palette and give the legend a nicer title
  #scale_color_viridis_d(name = "Number of\nData Points") +
  # Force x‐axis breaks at each unique nr_param value (so ticks show up nicely)
  scale_x_continuous(
    breaks = sort(unique(as.numeric(data_to_plot$nr_param))),
    labels = scales::comma_format()  # if nr_param is large, formatting helps
  ) +
  # Add informative labels
  labs(
    title   = "Processing Time vs. Model Size",
    #subtitle = "Lines grouped by different numbers of data points",
    x       = "Number of Parameters",
    y       = "Time (seconds)"
  ) +
  # Pick a clean theme and bump up base font size
  theme_minimal(base_size = 28) +
  theme(
    plot.title       = element_text(face = "bold", size = 26, hjust = 0.5),
   # plot.subtitle    = element_text(size = 20, hjust = 0.5, color = "gray40"),
    axis.title       = element_text(face = "bold"),
    axis.text        = element_text(color = "gray20"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    legend.text      = element_text(size = 24),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )




######### EM VS direct --------


betas_mp_driect <- list()
betas2_mp_driect <- list()
a_mp_driect <- list()

time_mp_driect <- list()
dens <- 0.4
for(id in 1:10){
  
  tmp <- new.env()
  loaded_objs <- load(paste0("sim_data/", "zip_direct", "_gamma_per_t", dens, "id", id, ".RData") , envir = tmp)
  my_list <- mget(loaded_objs, envir = tmp)
  
  betas_mp_driect[[as.character(id)]] <- bind_rows(my_list$beta_est_ig)
  betas_mp_driect[[as.character(id)]]$id <- id
  betas_mp_driect[[as.character(id)]]$t <- as.numeric(names(my_list$beta_est_ig))
  
  a_mp_driect[[as.character(id)]] <- as.data.frame(t(bind_rows(my_list$A_est_ig)))
  a_mp_driect[[as.character(id)]]$id <- id
  a_mp_driect[[as.character(id)]]$t <- as.numeric(names(my_list$beta_est_ig))
  
  betas2_mp_driect[[as.character(id)]] <- as.data.frame(t(data.frame(my_list$beta2_est_ig)))
  betas2_mp_driect[[as.character(id)]]$id  <- id
  betas2_mp_driect[[as.character(id)]]$t <- as.numeric(names(my_list$beta_est_ig))
  
  
  
  time_mp_driect[[as.character(id)]] <- data.frame(time = sapply(my_list$time_ig,  function(dt) as.numeric(dt, units = "mins")),
                                                   id = id,
                                                   t = as.numeric(names(my_list$beta_est_ig)))
  
  a_true <- my_list$sim$a
  beta_true <- my_list$sim$beta1
  beta2_true <- my_list$sim$beta2
  
}


betas_mp <- list()
betas2_mp <- list()
a_mp <- list()

time_mp <- list()
dens <- 0.4
for(id in 1:5){
  
  tmp <- new.env()
  loaded_objs <- load(paste0("sim_data/", "zip", "_gamma_per_t", dens, "id", id, ".RData") , envir = tmp)
  my_list <- mget(loaded_objs, envir = tmp)
  
  betas_mp[[as.character(id)]] <- bind_rows(my_list$beta_est_ig)
  betas_mp[[as.character(id)]]$id <- id
  betas_mp[[as.character(id)]]$t <- as.numeric(names(my_list$beta_est_ig))
  a_mp[[as.character(id)]] <- as.data.frame(t(bind_rows(my_list$A_est_ig)))
  a_mp[[as.character(id)]]$id <- id
  a_mp[[as.character(id)]]$t <- as.numeric(names(my_list$beta_est_ig))
  
  betas2_mp[[as.character(id)]] <- as.data.frame(t(data.frame(my_list$beta2_est_ig)))
  betas2_mp[[as.character(id)]]$id  <- id
  betas2_mp[[as.character(id)]]$t <- as.numeric(names(my_list$beta_est_ig))
  
  
  
  time_mp[[as.character(id)]] <- data.frame(time = sapply(my_list$time_ig,  function(dt) as.numeric(dt, units = "mins")),
                                            id = id,
                                            t = as.numeric(names(my_list$beta_est_ig)))
  
  a_true <- my_list$sim$a
  beta_true <- my_list$sim$beta1
  beta2_true <- my_list$sim$beta2
  
}









betas_mp <- bind_rows(betas_mp)
a_mp <- bind_rows(a_mp)
time_mp <- bind_rows(time_mp)
betas2_mp <- bind_rows(betas2_mp)

betas_mp_driect <- bind_rows(betas_mp_driect)
a_mp_driect <- bind_rows(a_mp_driect)
time_mp_driect <- bind_rows(time_mp_driect)
betas2_mp_driect <- bind_rows(betas2_mp_driect)



# # log likelihood
# 
# additive <- TRUE
# sim <- simulate_claims(50, 5000, spatial_type = "graph",additive =  additive, area = nr_regions, 
#                        model_type = "zip", mixing = "ig", density = 0.4,  seed = 2)
# print(t)
# # subset claims
# t <- 1000
# claims <- sim$claims[sim$years <= t]
# locs <- sim$locs[sim$years <= t]
# agg_claims <- sim$agg_claims[, 1:t]
# X <- sim$X[sim$years <= t, ]
# years <- sim$years[sim$years <= t]
# exposure <- sim$exposure[sim$years <= t]
# 
# 
# tmp <- new.env()
# loaded_objs <- load(paste0("sim_data/", "zip", "_ln_per_t",  0.4, "id", 2, ".RData") , envir = tmp)
# my_list <- mget(loaded_objs, envir = tmp)
# 
# 
# # Compute spatial aggregate effects and baseline nu
# se <- get_spatial_aggregate(locs, get_W_from_array(my_list$A_est_ig[["1000"]], 10), NA, agg_claims, years, "learn_graph")
# nu <- as.numeric(exp(X %*% my_list$beta_est_ig[["1000"]]))  # baseline mu from covariates
# 
# # Get full mu by combining covariate effects (nu), spatial effects, and exposure
# mu <- get_mu(nu, se$spatial_effect, exposure, additive)
# 
# phi <- exp(my_list$beta2_est_ig[["1000"]])
# prop <- out_ig$prop
# 
# 
# sum(log_zinb(claims, mu, phi, prop)) 


# error 

beta_true_matrix_mp <- matrix(beta_true,
                              nrow = nrow(betas_mp),
                              ncol = length(beta_true),
                              byrow = TRUE)

beta_true_matrix_mp_direct <- matrix(beta_true,
                                     nrow = nrow(betas_mp_driect),
                                     ncol = length(beta_true),
                                     byrow = TRUE)

betas_mp$error <- rowSums(abs(betas_mp[, c("X11", "X12", "X13")] - beta_true_matrix_mp))
betas_mp_driect$error <- rowSums(abs(betas_mp_driect[, c("X11", "X12", "X13")] - beta_true_matrix_mp_direct))

a_true_matrix_mp <- matrix(a_true,
                           nrow = nrow(a_mp),
                           ncol = length(a_true),
                           byrow = TRUE)

a_true_matrix_direct <- matrix(a_true,
                               nrow = nrow(a_mp_driect),
                               ncol = length(a_true),
                               byrow = TRUE)


a_mp$error <- rowSums(abs(a_mp[ ,1:(10*11/2)] - a_true_matrix_mp))
a_mp_driect$error <- rowSums(abs(a_mp_driect[ ,1:(10*11/2)] - a_true_matrix_direct))


ret <- list()


df1 <- betas_mp %>% mutate(model = "Mixed ZIP")
df2 <- betas_mp_driect %>% mutate(model = "Mixed ZIP direct")

# bind into one
betas_all <- bind_rows(df1, df2)

# dodged box‐plot
ggplot(betas_all, aes(x = factor(t), y = X11, fill = model)) +
  # true parameter line + points
  geom_line(
    aes(x = factor(t), y = 1, group = 1),
    inherit.aes = FALSE,
    color       = "black",
    size        = 1,
    alpha = 0.5
  ) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", 
               position = position_dodge(width = 0.8), 
               width = 0.2) + 
  labs(
    x     = "Time (t)",
    y     = expression(beta[1] ~ " estimate"),
    fill  = "Model spec",
    title = expression(beta[1] ~ " over Time, by Model")
  )  +
  theme_minimal(base_size = 24) +
  theme(
    plot.title   = element_text(hjust = 0.5),
    axis.text.x  = element_text(angle = 45, vjust = 0.5)
  )+ ylim(c(-0.5, 1.2)) + theme_minimal(base_size = 35) +labs(title = expression(beta[1] ~ " over time by model")) + theme(
    legend.position     = c(0.60, 0.30),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  )

# dodged box‐plot
ggplot(betas_all, aes(x = factor(t), y = error, fill = model)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", 
               position = position_dodge(width = 0.8), 
               width = 0.2) + 
  labs(
    x     = "Time (t)",
    y     = expression(l[1] ~ " error"),
    fill  = "Model spec",
    title = expression(l[1] ~ " error over Time, by Model")
  ) +
  theme_minimal(base_size = 24) +
  theme(
    plot.title   = element_text(hjust = 0.5),
    axis.text.x  = element_text(angle = 45, vjust = 0.5)
  )+ ylim(c(-0.5, 15)) + theme_minimal(base_size = 35) +labs(title = expression(l[1] ~ " error for " ~  beta ~  " over time by model")) + theme(
    legend.position     = c(0.60, 0.70),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  )+  scale_y_continuous(trans = "log2")


df1 <- a_mp %>% mutate(model = "Mixed Poisson")
df2 <- a_mp_driect %>% mutate(model = "Mixed Poisson direct")

# bind into one
as_all <- bind_rows(df1, df2)

# dodged box‐plot
ggplot(as_all, aes(x = factor(t), y = V8, fill = model)) +
  # true parameter line + points
  geom_line(
    aes(x = factor(t), y = a_true[8], group = 1),
    inherit.aes = FALSE,
    color       = "black",
    size        = 1,
    alpha = 0.5
  ) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", 
               position = position_dodge(width = 0.8), 
               width = 0.2) + 
  labs(
    x     = "Time (t)",
    y     = latex2exp::TeX("$a_{1,1}$ estimate"),
    fill  = "Model spec",
    title = latex2exp::TeX("$a_{1,1}$ over time, by model")
  )  +
  theme_minimal(base_size = 24) +
  theme(
    plot.title   = element_text(hjust = 0.5),
    axis.text.x  = element_text(angle = 45, vjust = 0.5)
  )+  ylim(c(-0.5, 4)) + theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$a_{4,1}$ estimate over time by model"),
                                                             y = latex2exp::TeX("$a_{4,1}$ estimate")) + 
  theme(
    legend.position     = c(0.60, 0.70),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  )


# dodged box‐plot
ggplot(as_all, aes(x = factor(t), y = log(error), fill = model)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", 
               position = position_dodge(width = 0.8), 
               width = 0.2) + 
  labs(
    x     = "Time (t)",
    y     = latex2exp::TeX("$a_{1,1}$ estimate"),
    fill  = "Model spec",
    title = latex2exp::TeX("$a_{1,1}$ over time, by model")
  )  +
  theme_minimal(base_size = 24) +
  theme(
    plot.title   = element_text(hjust = 0.5),
    axis.text.x  = element_text(angle = 45, vjust = 0.5)
  )+ ylim(c(0, 5)) + theme_minimal(base_size = 35) +labs(title = latex2exp::TeX("$l_1$ error for $A$ over time by model"),
                                                         y = latex2exp::TeX("$l_{1}$ error")) + 
  theme(
    legend.position     = c(0.10, 0.10),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  ) +  scale_y_continuous(trans = "log2")




df1 <- time_mp %>% mutate(model = "Mixed Poisson")
df2 <- time_p %>% mutate(model = "Poisson")
times_df <- bind_rows(df1, df2)

ret$times <- ggplot(times_df, aes(x = as.factor(t), y = log(time), fill = model)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", 
               position = position_dodge(width = 0.8), 
               width = 0.2)+ 
  labs(
    x     = "Time points (t)",
    y     = latex2exp::TeX("$log$ complexity"),
    fill  = "Model spec",
    title = latex2exp::TeX("$a_{1,1}$ over time, by model")
  )  +
  theme_minimal(base_size = 24) +
  theme(
    plot.title   = element_text(hjust = 0.5),
    axis.text.x  = element_text(angle = 45, vjust = 0.5)
  )
ret$times

df1 <- betas2_mp %>% mutate(model = "Mixed ZIP")
df2 <- betas2_mp_driect %>% mutate(model = "Mixed ZIP direct")

# bind into one
beta2_all <- bind_rows(df1, df2)

# dodged box‐plot
ggplot(beta2_all, aes(x = factor(t), y = V1, fill = model)) +
  # true parameter line + points
  geom_line(
    aes(x = factor(t), y = beta2_true, group = 1),
    inherit.aes = FALSE,
    color       = "black",
    size        = 1,
    alpha = 0.5
  ) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", 
               position = position_dodge(width = 0.8), 
               width = 0.2) + 
  labs(
    x     = "Time (t)",
    y     = latex2exp::TeX("$a_{1,1}$ estimate"),
    fill  = "Model spec",
    title = latex2exp::TeX("$a_{1,1}$ over time, by model")
  )  +
  theme_minimal(base_size = 36) +
  theme(
    plot.title   = element_text(hjust = 0.5),
    axis.text.x  = element_text(angle = 45, vjust = 0.5)
  ) +labs(title = latex2exp::TeX("$\\beta_2$ estimate over time by model"),
          y = latex2exp::TeX("$\\beta_2$ estimate")) + 
  theme(
    legend.position     = c(0.60, 0.60),    # inside panel, 2% from left & bottom
    legend.justification = c(0, 0),         # anchor legend’s bottom‐left corner
    legend.background    = element_rect(fill = alpha("white", 0.5))
  ) + ylim(c(0,2.5))












