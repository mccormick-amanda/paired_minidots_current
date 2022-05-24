#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)
library(GGally)
source("model/stan_utility.R")

# stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

# select site year
type = "reyk_2018"

# read data
data_list = read_rdump(paste0("data/", type,  "/data_list.R"))

names(data_list)
data_list$obs_per_day


# set reference temperature
data_list$temp_ref = 12

# set sd
data_list$sig_b0 = 0.01
data_list$sig_r = 0.01


#==========
#========== Fit model
#==========

# initial specifications

model = "coupled_o2_test"


model_path = paste0("model/",model,".stan")
chains = 3
iter = 100

start = Sys.time()
# fit model
fit = stan(file = model_path, data = data_list, seed=1, chains = chains, 
           iter = iter, cores = chains, 
           control = list(adapt_delta = 0.8, max_treedepth = 10))

end = Sys.time()
end - start

# saveRDS(fit, paste0("model/output/", model,"/", type,"/output.rds"))

# summary of fit
fit_summary = summary(fit, probs=c(0.16, 0.5, 0.84))$summary %>% 
{as_tibble(.) %>%
    mutate(var = rownames(summary(fit)$summary))}

# check Rhat & n_eff
fit_summary %>% filter(Rhat > 1.05) %>% select(Rhat, n_eff, var) %>% arrange(-Rhat)
fit_summary %>% filter(n_eff < 0.5*(chains*iter/2)) %>% select(Rhat, n_eff, var) %>% 
  arrange(n_eff) %>%
  mutate(eff_frac = n_eff/(chains*iter/2))


# additional diagnostics
check_div(fit)
check_treedepth(fit)
check_energy(fit)




#==========
#========== Examine Chains
#==========


# fixed parameters by step

# for coupled_o2_update_sd (sig_b0 and sig_r are no longer parameters)
# fixed_par_v = c("gamma_1[1]","gamma_2[1]","gamma_1[2]","gamma_2[2]",
#                 "a[1]","a[2]","sig_proc","c0","c1","lp__")

#for coupled_o2_update_Iopt
fixed_par_v = c("gamma_1[1]","gamma_2[1]","gamma_1[2]","gamma_2[2]",
                "opt[1]","opt[2]","sig_proc","c0","c1","lp__","b0[1,1]","r[1,1]")



fixed_pars = rstan::extract(fit, pars=fixed_par_v) %>%
  lapply(as_tibble) %>%
  bind_cols() %>%
  mutate(chain = rep(1:chains, each = iter/2), step = rep(c(1:(iter/2)), chains))
names(fixed_pars) = c(fixed_par_v,"chain","step")

# examine chains for parameters
fixed_pars %>%
  gather(par, value, -chain, -step) %>%
  filter(!(par  %in% c("lp__"))) %>%
  ggplot(aes(step, value, color=factor(chain)))+
  facet_wrap(~par, scales="free_y")+
  geom_line(alpha=0.5)+
  theme_bw()



# posterior densities
fixed_pars %>%
  gather(par, value, -chain, -step) %>%
  filter(!(par  %in% c("lp__"))) %>%
  ggplot(aes(value))+
  facet_wrap(~par, scales="free")+
  stat_density(alpha=0.5, geom = "line")+
  theme_bw()



#==========
#========== Prepare Output
#==========

# extract revelvant data from summary
fit_clean = fit_summary %>%
  rename(lower = `16%`, middle = `50%`, upper = `84%`)  %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~.x[1]))


# read data
data_list = read_rdump(paste0("data/", type, "/data_list.R"))


# extract structural information 
str_data = read_csv(paste0("data/", type, "/data_check.csv")) %>%
  na.omit()



# fixed parameters
# fixed_summary = fit_clean %>%
#   filter(var %in% c(fixed_par_v, "alpha[1]", "alpha[2]","alpha")) %>%
#   mutate(layer = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
#          layer = ifelse(layer==1, "pelagic", "benthic"))

fixed_summary = fit_clean %>%
  filter(var %in% c(fixed_par_v, "omega[1]", "omega[2]","omega")) %>%
  mutate(layer = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         layer = ifelse(layer==1, "pelagic", "benthic"))



# hourly data
hourly_summary = fit_clean %>%
  filter(name %in% c("b","chi","kappa","phi","D","delta","exc","x","x_pred","beta",
                     "gpp","er","nep")) %>% 
  mutate(index = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         layer = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         layer = ifelse(layer==1, "pelagic", "benthic")) %>%
  left_join(str_data %>%
              select(year, site, yday, hour) %>%
              mutate(index = row_number())) %>%
  select(name, site, layer, year, yday, hour, lower, middle, upper)




# daily data
daily_summary = fit_clean %>%
  filter(name %in% c("log_b0","log_r","b0","r","beta0","rho","GPP","ER","NEP")) %>%
  mutate(unique_day = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         layer = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         layer = ifelse(layer==1, "pelagic", "benthic")) %>%
  full_join(str_data %>%
              group_by(unique_day) %>%
              summarize(year = unique(year),
                        yday = unique(yday),
                        site = unique(site)) %>%
              ungroup) %>%
  select(name, site, year, yday, layer, lower, middle, upper) 




# write_csv(fixed_pars, paste0("model/output/", model,"/", type,"/fixed_pars.csv"))
# write_csv(fixed_summary, paste0("model/output/", model,"/", type,"/fixed_summary.csv"))
# write_csv(hourly_summary, paste0("model/output/", model,"/", type,"/hourly_summary.csv"))
# write_csv(daily_summary, paste0("model/output/", model,"/", type,"/daily_summary.csv"))


daily_summary %>%
  filter(name == "r") %>%
  ggplot(aes(x = yday,
             y = middle,
             color = layer))+
  geom_line()+
  geom_ribbon(aes(ymin = lower,
                  ymax = upper,
                  fill = layer),
              alpha = 0.1,
              linetype = 0)

fit_summary %>% filter(str_detect(var, "opt"))


hourly_summary %>%
  filter(name == "b")