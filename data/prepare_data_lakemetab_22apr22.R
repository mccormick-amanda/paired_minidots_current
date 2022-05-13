#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(lubridate)
library(rstan)
library(LakeMetabolizer)

# desired sites
# desired years
sites = c("reyk")
years = c(2019)

# read data
data = read_csv("data/clean_coupled_o2/mini_clean_22apr22.csv") %>%
  filter(site %in% sites,
         year %in% years) %>%
  # convert hour from 0:23 to 1:24 (this makes indexing easier later)
  mutate(hour = hour + 1) %>% 
  # LakeMetabolizer requires specific names for columns and a datetime column
  mutate(wnd = wspeed,
         wtr = temp.pel,
         datetime = ymd_hms(paste(date, hour, "00:00"))) 

data %>% select(date, yday, hour, temp.ben, temp.pel) %>%
  gather(layer, temp, c(temp.ben, temp.pel)) %>%
  filter(temp<20) %>%
  ggplot(aes(yday+hour/24, temp, color=layer))+
  geom_line()+
  scale_color_manual(values=c("black", "deepskyblue3"))+
  theme(axis.title.x = element_blank())

# calculate k gas exchange and o2_equilibrium with Lake Metabolizer

data$k600 = k.cole(data.frame(datetime = data$datetime, wnd = data$wnd))$k600 #units m/d
data$k600 = data$k600/24 # we want hourly units (m/h)

#view k.cole code
#k.cole # calls on k.cole.base

#view k.cole.base code
#k.cole.base

# checking with their code for k.cole()
#data$k600_b <- 2.07 + (0.215 * (data$wnd^1.7)) # units in cm h-1
#data$k600_b <- data$k600_b*24/100 #units in m d-1
#data$k600_b <- data$k600_b/100 #units in m hr-1 (this is what we want for the model; NOT daily units)

data$k <- k600.2.kGAS.base(k600 = data$k600,
                           temperature = data$wtr,
                           gas = "O2")

data$o2_eq <- o2.at.sat(data %>% select(datetime, wtr),
                        altitude = 277)$do.sat

hist(data$wspeed)

data %>%
  select(wspeed) %>%
  mutate(log_wspeed = log(wspeed+0.1),
         sqrt_wspeed = sqrt(wspeed)) %>% 
  gather(transformation, val, c(wspeed:sqrt_wspeed)) %>%
  ggplot(aes(val))+
  facet_grid(~transformation, scales = 'free')+
  ggtitle("Reyk 2019 Wind") +
  geom_histogram() +
  theme(axis.title.x = element_blank())

names(data)
data %>% filter(!is.na(wspeed)) %>% filter(wspeed>0) %>% summarize(min(wspeed))

# omit unneccesary columns and expand to full yday x 24hrs
clean_data = data %>%
  select(year, site, yday, hour, do.pel, do.ben, temp.pel, temp.ben, par.pel, par.ben,
         wspeed, k, o2_eq) %>%
  na.omit() %>%
  full_join(data %>% 
              expand(year, site, yday, hour)) %>%
  arrange(site, year, yday, hour) %>%
  # use sqrt transformed wspeed
  mutate(wspeed = sqrt(wspeed))

hist(clean_data$wspeed)


#==========
#========== Prepare for data analysis
#==========

# prepare data
data_prep = clean_data %>%
  mutate(site_year = paste(site, year, sep="_")) %>%
  # for each year, create identifier for uninterrupted stretches of observations
  group_by(site_year) %>%
  mutate(i = ifelse(is.na(do.pel)==T, 1, 0), 
         j = c(1,abs(diff(i)))) %>% 
  filter(is.na(do.pel)==F, is.na(wspeed)==F) %>%
  mutate(series = cumsum(j)) %>% 
  ungroup() %>%
  # create unique index for each series
  # remove series with fewer than 24 observations
  mutate(unique_series = as.numeric(as.factor(site_year)) + series/length(unique(series))) %>%
  group_by(unique_series) %>%
  mutate(series_length = length(unique_series)) %>%
  filter(series_length > 23) %>%
  ungroup() %>%
  # recreate series index and make unique index for, days
  # create index for observations (for joining later)
  # replace 0 par.pel with smallest non-zero value
  mutate(unique_series = as.factor(unique_series) %>% as.numeric(),
         unique_day = paste(site_year, yday) %>% as.factor() %>% as.numeric(),
         index = 1:length(do.pel),
         par.pel = ifelse(par.pel==0, min(par.pel[which(par.pel>0)]), par.pel),
         par.ben = ifelse(par.ben==0, min(par.ben[which(par.ben>0)]), par.ben)
  ) %>%
  select(-i, -j) 

# return missing observations for check
data_check = data_prep %>% 
  expand(site_year,yday,hour) %>%
  full_join(data_prep) %>%
  arrange(site_year,yday,hour)

# check unique_series
data_check %>%
  # filter(site_year=="reyk_2018") %>%
  mutate(time = yday + hour/24) %>%
  ggplot(aes(time, do.pel, color=factor(unique_series)))+
  facet_wrap(~site_year)+
  geom_line()+
  scale_color_discrete(guide = F)+
  theme_bw()

# check unique_days
data_check %>%
  mutate(time = yday + hour/24) %>%
  filter(unique_day %in% (45 + c(1:10))) %>%
  ggplot(aes(time, do.pel, color=factor(unique_day)))+
  facet_wrap(~site_year, scale = "free", nrow=2)+
  geom_line()+
  theme_bw()

# export prepared data
data_check %>% write_csv(paste0("data/", unique(data_check$site_year), "/data_check.csv"))



#==========
#========== Package data 
#==========

# define variables in evnironment 
o2_obs = data_prep %>% with(cbind(do.pel, do.ben))
o2_eq = data_prep$o2_eq
light = data_prep %>% with(cbind(par.pel, par.ben))
temp = data_prep %>% with(cbind(temp.pel, temp.ben))
wspeed = data_prep$wspeed
k = data_prep$k
map_sites = as.numeric(as.factor(data_prep$site_year))
map_days = data_prep$unique_day
days_per_site = c({data_prep %>%
    group_by(site_year) %>%
    summarize(value = length(unique(unique_day)))}$value) 
days_per_site = array(days_per_site)
obs_per_series = c({data_prep %>%
    group_by(unique_series) %>%
    summarize(value = length(unique_series))}$value) 
obs_per_day = c({data_prep %>%
    group_by(unique_day) %>%
    summarize(value = length(unique_day))}$value) 

depth = ifelse(sites == "st33", 3.3, 4.2)
z <- 
  structure(c(depth),
            .Dim = c(1))

n_obs = length(o2_obs[,1])
n_series = length(obs_per_series) 
n_days = sum(days_per_site)
n_sites = length(days_per_site)

# export as .R
stan_rdump(c("o2_obs","o2_eq","light","temp","wspeed","k","map_days",
             "obs_per_series","days_per_site", "map_sites",
             "obs_per_day", "z","n_obs","n_series","n_days","n_sites"),
           file=paste0("data/", unique(data_check$site_year), "/data_list.R"))

