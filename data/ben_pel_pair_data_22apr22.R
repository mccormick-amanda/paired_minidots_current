#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(lubridate)


################
# OBJECTIVE:
# preparing paired benthic and pelagic data for coupled o2 model
# requires hourly observations of temp, do, do.sat, par for each habitat


# only preparing data for st33 and reyk, in 2017, 2018 and 2019
sites = c("st33", "reyk")


################
# LOAD MINIDOTS
mini = read_csv("data/raw/minidot_clean_12Dec20.csv") %>%
  filter(calibrating == 'n',
         site %in% sites) %>% 
  mutate(year = year(date_time),
         date = date(date_time),
         yday = yday(date_time),
         hour = hour(date_time)) %>% 
  select(site, layer, date_time, year, date, yday, hour, md, q:do_sat)  %>%
  filter( #restrict to the summer periods (avoid winter deployments)
    (date(date_time) >= "2017-06-29" & date(date_time) <= "2017-08-18") | 
    (date(date_time) >= "2018-06-13" & date(date_time) <= "2018-08-14") |
    (date(date_time) >= "2019-06-15" & date(date_time) <= "2019-08-14") ) %>% 
  select(-do) %>%
  rename(do = do_cor) #using the corrected do values

# more specific screening based on when minidots were deployed at each site, accounting for sonde deployment, etc.
mini = mini %>%
  filter(!(site=='reyk' & year==2018 & yday<165)) %>%
  filter(!(site=='reyk' & year==2019 & yday<168)) %>%
  filter(!(site=='st33' & year==2018 & yday<165)) %>%
  filter(!(site=='st33' & year==2019 & yday<171)) %>%
  filter(!(site=='st33' & year==2019 & yday>223)) 

mini %>% group_by(site, year) %>% summarize(min=min(yday), max=max(yday))


# plot the raw minidot data
mini %>% 
  mutate(siteyr = paste(site, year, sep='_')) %>% 
  select(-c(do_eq)) %>%
  gather(var, val, c(q:do_sat)) %>%
  ggplot(aes(yday, val, color=layer))+
  facet_grid(var~siteyr, scales = 'free_y')+
  geom_line()


 

################
# LOAD SONDE
sonde = read_csv("data/raw/sonde_clean.csv") %>%
  mutate(year = year(date_time),
         date = date(date_time),
         yday = yday(date_time),
         hour = hour(date_time),
         site = "st33",
         layer = "pel") %>%
  filter( #restrict to the periods overlapping minidots (note 2019 had some power issues before june 20)
    (date(date_time) >= "2017-06-30" & date(date_time) <= "2017-08-18") | 
    (date(date_time) >= "2018-06-14" & date(date_time) <= "2018-08-14") |
    (date(date_time) >= "2019-06-20" & date(date_time) <= "2019-08-14") ) %>%
  select(site, layer, date_time, year, date, yday, hour, temp:do_eq)

sonde %>% group_by(year) %>% summarize(min=min(yday), max=max(yday))

  
# plot the raw sonde data
sonde %>% 
  select(site:hour, temp, do_sat, do, do_sol, do_eq) %>%
  gather(var, val, c(temp:do_eq)) %>%
  ggplot(aes(yday, val))+
  facet_grid(var~year, scales = 'free_y')+
  geom_point()



################
# LOAD HOBOS
hobo = read_csv("data/raw/hobo_clean_12Oct21.csv") %>%
  mutate(year = year(datetime),
         date = date(datetime),
         yday = yday(datetime),
         hour = hour(datetime),
         par = lux*0.0185) %>% #convert lux to par
  filter(site %in% sites,
         layer %in% c("ben", "pel")) %>% 
  select(site, layer, datetime, date, year, yday, hour, par) %>%
  rename(date_time = datetime) %>%
  filter( #restrict to the summer periods
    (date(date_time) >= "2017-06-30" & date(date_time) <= "2017-08-18") | 
      (date(date_time) >= "2018-06-14" & date(date_time) <= "2018-08-14") |
      (date(date_time) >= "2019-06-20" & date(date_time) <= "2019-08-14") )

hobo %>% 
  mutate(siteyr = paste(site, year, sep='_')) %>%
  ggplot(aes(yday, par, color=layer))+
  facet_grid(layer~siteyr)+
  geom_line()

# 2017 doesn't have benthic or pelagic hobo (weren't deployed alongside the sensors this year)
# 2019 st33 pelagic hobo seems to be missing, and there isn't corresponding data in the raw file in sensor project


################
# LOAD WEATHER
weather = read_csv("data/raw/weather_2010_2019.csv") %>%
  mutate(date = date(Date_Time),
         year = year(Date_Time), 
         hour = hour(Date_Time) ) %>%
  filter( #restrict to the summer periods
    (date(Date_Time) >= "2017-06-30" & date(Date_Time) <= "2017-08-18") |
    (date(Date_Time) >= "2018-06-14" & date(Date_Time) <= "2018-08-14") |
    (date(Date_Time) >= "2019-06-20" & date(Date_Time) <= "2019-08-14") ) %>%
  rename(date_time = Date_Time)

# plot the raw weather data
weather %>% 
  select(date_time, date, year, hour, wdir:gust) %>%
  gather(var, val, c(wdir:gust)) %>%
  ggplot(aes(yday(date)+hour/24, val))+
  facet_grid(var~year, scales = 'free_y')+
  geom_line()



################
# HOBO: deal with benthic and pelagic par for 2017
# there were no layer-specific hobos, so need to estimate based on the st33 buoy hobo and light attenuation (from turb)

sonde_17 = read_csv("data/raw/sonde_clean.csv",
                    col_types = cols(turbsc = col_double())) %>% 
  filter(date(date_time) >= "2017-06-30" & date(date_time) <= "2017-08-18") %>% 
  mutate(year = year(date_time),
         date = date(date_time),
         yday = yday(date_time),
         hour = hour(date_time)) %>%
  select(date_time, date, yday, year, hour, turbsc) 

# going to remove outlier points below when aggregating to day
sonde_17 %>% 
  ggplot(aes(date_time, turbsc))+
  geom_point()

# get daily turbidity average to estimate light attenuation coef (kD) for each day
sonde_17_agg = left_join(sonde_17,
                         sonde_17 %>% 
                           group_by(date) %>%
                           summarize(turb.mean = mean(turbsc),
                                     turb.sd = sd(turbsc))) %>%
  #filter out points that are outside 2sd of the daily mean
  filter(turbsc >= turb.mean - 2*turb.sd, turbsc <= turb.mean + 2*turb.sd) %>%
  select(-c(turb.mean, turb.sd)) %>%
  group_by(date, yday, year, hour) %>% #agg to hour
  summarise(turb = mean(turbsc, na.rm=T)) %>%
  ungroup() %>%
  group_by(date, yday, year) %>% #agg to day
  summarise(turb = mean(turb, na.rm=T)) %>%
  ungroup() %>%
  mutate(kD = 0.46 + 0.05*turb) #based on empirical relationship between daily turb and LTREB light profiles (Inland Waters SuppInfo)

# not sure what is generating the large peak in mid-July (July 19)
sonde_17_agg %>%
  ggplot(aes(date, turb))+ geom_line()

# here are the days surrounding this jump in turb
read_csv("data/raw/sonde_clean.csv",
         col_types = cols(turbsc = col_double())) %>% 
  filter(date(date_time) >= "2017-07-17" & date(date_time) <= "2017-07-21") %>% 
  mutate(date = date(date_time),
         hour = hour(date_time)) %>%
  select(date_time, date, hour, temp, do_sat, do, pcyv, turbsc) %>% 
  group_by(date, hour) %>%
  summarize(temp=mean(temp), do_sat=mean(do_sat), do=mean(do), pcyv=mean(pcyv), turbsc=mean(turbsc)) %>%
  gather(var, val, temp:turbsc) %>%
  ggplot(aes(x=date+hour/24, y=val))+
  facet_wrap(~var, scales='free_y')+
  geom_line()

# possibly driven by strong wind gusts? (red line)
weather %>%
  filter(date(date_time) >= "2017-06-30" & date(date_time) <= "2017-08-18") %>%
  mutate(date = date(date_time)) %>%
  group_by(date) %>%
  summarize(wspeed=mean(wspeed), gust=mean(gust) ) %>%
  ggplot(aes(x=date, y=wspeed))+
  geom_line(color='black')+
  geom_line(aes(date, gust), color='red')

# pull out the 2017 hobo data from the sonde buoy
hobo_17 = read_csv("data/raw/hobo_clean_12Oct21.csv") %>%
  mutate(year = year(datetime),
         date = date(datetime),
         yday = yday(datetime),
         hour = hour(datetime)) %>% 
  filter(site == "st33", 
         layer == "out",
         date(datetime) >= "2017-06-30" & date(datetime) <= "2017-08-18") %>% 
  mutate(par_surf = lux*0.0084519) %>% #determined from relationship between 0m light profile readings and hobo from LTREB data (Inland Waters SuppInfo)
  select(site, layer, datetime, date, year, yday, hour, par_surf) %>%
  group_by(date, yday, year, hour) %>%
  summarize(par_surf = mean(par_surf, na.rm=T))

# determine in situ light from kD and surface par levels
par_data_17 = left_join(hobo_17, sonde_17_agg) %>%
  mutate(par.ben = par_surf*exp(-kD*3.3), #determine benthic and pelagic par based on I0*e^(-kD*z)
         par.pel = par_surf*exp(-kD*1)) %>%
  select(-c(par_surf:kD))



################
# HOBO: deal with pelagic par for 2019
# there is no pelagic hobo and I didn't see one in the raw data file in the sensors project, so not sure if there was a "pelagic" hobo

sonde_19 = read_csv("data/raw/sonde_clean.csv",
                    col_types = cols(turbsc = col_double())) %>% 
  filter(date(date_time) >= "2019-06-20" & date(date_time) <= "2019-08-14") %>% 
  mutate(year = year(date_time),
         date = date(date_time),
         yday = yday(date_time),
         hour = hour(date_time)) %>%
  select(date_time, date, yday, year, hour, turbsc) 

# remove outlier points when aggregating below; 
sonde_19 %>% filter(turbsc < 100) %>% 
  ggplot(aes(date_time, turbsc))+
  geom_point()

# get daily turbidity average to estimate light attenuation coef (kD)
sonde_19_agg = left_join(sonde_19,
                         sonde_19 %>% 
                           group_by(date) %>%
                           summarize(turb.mean = mean(turbsc),
                                     turb.sd = sd(turbsc))) %>%
  #filter out points that are outside 2sd of the daily mean
  filter(turbsc >= turb.mean - 2*turb.sd, turbsc <= turb.mean + 2*turb.sd) %>%
  select(-c(turb.mean, turb.sd)) %>%
  group_by(date, yday, year, hour) %>% #agg to hour
  summarise(turb = mean(turbsc, na.rm=T)) %>%
  ungroup() %>%
  group_by(date, yday, year) %>% #agg to day
  summarise(turb = mean(turb, na.rm=T)) %>%
  ungroup() %>%
  mutate(kD = 0.46 + 0.05*turb) #based on empirical relationship between daily turb and LTREB light profiles (Inland Waters SuppInfo)

# turbidity gets really high at end of summer, but looks like there were high winds also
sonde_19_agg %>%
  ggplot(aes(date, turb))+ geom_line()

# high wind at end of summer
weather %>%
  filter(date(date_time) >= "2019-06-20" & date(date_time) <= "2019-08-11") %>%
  mutate(date = date(date_time)) %>%
  group_by(date) %>%
  summarize(wspeed=mean(wspeed), gust=mean(gust)) %>%
  ggplot(aes(x=date, y=wspeed))+
  geom_line(color='black')+
  geom_line(aes(date, gust), color='red')

# benthic hobo seems to match the high turbidity at the end of the summer
hobo %>% filter(site=='st33', layer=='ben', year==2019) %>%
  ggplot(aes(date_time, par)) +
  geom_line()

# pull out the 2019 hobo data from the sonde buoy
hobo_19 = read_csv("data/raw/hobo_clean_12Oct21.csv") %>%
  mutate(year = year(datetime),
         date = date(datetime),
         yday = yday(datetime),
         hour = hour(datetime)) %>% 
  filter(site == "st33", 
         layer == "out",
         date(datetime) >= "2019-06-20" & date(datetime) <= "2019-08-11") %>% 
  mutate(par_surf = lux*0.0084519) %>% #determined from relationship between 0m light profile readings and hobo from LTREB data (Inland Waters Supp Info)
  select(site, layer, datetime, date, year, yday, hour, par_surf) %>%
  group_by(date, yday, year, hour) %>%
  summarize(par_surf = mean(par_surf, na.rm=T))

# determine in situ light from kD and surface par levels
par_data_19 = left_join(hobo_19, sonde_19_agg) %>%
  mutate(par.ben = par_surf*exp(-kD*3.3), #determine benthic and pelagic par based on I0*e^(-kD*z)
         par.pel = par_surf*exp(-kD*1)) %>%
  select(-c(par_surf:kD))

# compare the benthic hobo values from the in situ benthic hobo logger to the values estimated from kD and incident light
hobo %>% 
  filter(site == 'st33', layer == 'ben', year == 2019) %>%
  group_by(date, yday, year, hour) %>%
  summarize(par = mean(par, na.rm=T)) %>%
  rename(par_hobo_insitu = par) %>% 
  left_join(par_data_19) %>%
  gather(var, val, par_hobo_insitu:par.pel) %>%
  filter(var != 'par.pel') %>%
  ggplot(aes(yday+hour/24, val, color=var)) +
  geom_line()



################
# aggregate everything to hour
# also involves cleaning (e.g. removing anomalous values) 

# no 'cleaning' to the hobo data; we already checked for periods when hobos were out of the water (e.g., right before or after deployment/retrieval) when combining all hobo data
# aggregate to hour
hobo_agg = hobo %>%
  group_by(site, layer, date, yday, year, hour) %>%
  summarize(par = mean(par, na.rm=T)) %>%
  ungroup() 

hobo_agg %>%
  mutate(siteyr = paste(site, year, sep='_')) %>%
  ggplot(aes(yday+hour/24, par))+
  facet_grid(layer~siteyr)+
  geom_line()


# the sonde data appear reasonable with the exception of drastic drops and spikes in do and temp, respectively; these seem implausible for St33's water column
sonde %>% #filter(temp<20, do>5)%>%
  gather(var, val, c(do, temp)) %>%
  ggplot(aes(yday+hour/24, val))+
  facet_grid(var~year, scales="free_y")+
  geom_line()

# aggregate sonde to hour
# sonde_agg = sonde %>%
#   filter(do > 5, temp < 20) %>% # removing very low do and very high temp
#   group_by(site, layer, year, date, yday, hour) %>%
#   summarize(temp = mean(temp, na.rm=T),
#             do_sat = mean(do_sat, na.rm=T),
#             do = mean(do, na.rm=T), 
#             do_eq = mean(do_eq, na.rm=T)) %>%
#   ungroup()



#library(zoo)

head(mini)


mini2 = rbind(mini,
              sonde %>% 
                select(site:hour, temp, do_sat, do, do_eq) %>% 
                mutate(md = NA,
                       q = NA) )


mini2 = mini2 %>% mutate(site_yr_layer = paste(site, year, layer, sep="_"))

window = mini2 %>% 
  expand(nesting(site_yr_layer, date, yday)) %>% 
  mutate(do_roll_mean = 0,
         do_roll_sd = 0,
         temp_roll_mean = 0,
         temp_roll_sd = 0) %>%
  as.data.frame()

uniq = unique(window$site_yr_layer)

k = 1
for (i in 1:length(uniq)) {
  
  loop_dates = unique(unlist(window[which(window$site_yr_layer == uniq[i]),]$date))
  
  for (j in 1:length(loop_dates)){
    
    window$do_roll_mean[k] = 
      as.numeric(
        mini2 %>% 
          filter(site_yr_layer == uniq[i]) %>%
          filter(date <= loop_dates[j], date >= loop_dates[j] - 6) %>% 
          summarize(do_roll_mean = mean(do, na.rm=T)) )
    
    window$do_roll_sd[k] = 
      as.numeric(
        mini2 %>% 
          filter(site_yr_layer == uniq[i]) %>%
          filter(date <= loop_dates[j], date >= loop_dates[j] - 6) %>% 
          summarize(do_roll_sd = sd(do, na.rm=T)) )
    
    window$temp_roll_mean[k] = 
      as.numeric(
        mini2 %>% 
          filter(site_yr_layer == uniq[i]) %>%
          filter(date <= loop_dates[j], date >= loop_dates[j] - 6) %>% 
          summarize(temp_roll_mean = mean(temp, na.rm=T)) )
    
    window$temp_roll_sd[k] = 
      as.numeric(
        mini2 %>% 
          filter(site_yr_layer == uniq[i]) %>%
          filter(date <= loop_dates[j], date >= loop_dates[j] - 6) %>% 
          summarize(temp_roll_sd = sd(temp, na.rm=T)) )
    
    k = k + 1
  }
}

window

# DO running 7d ave
window %>%
  ggplot(aes(yday, do_roll_mean))+
  facet_wrap(~site_yr_layer, nrow = 5)+
  geom_ribbon(aes(ymin = do_roll_mean - 3*do_roll_sd, ymax = do_roll_mean + 3*do_roll_sd), alpha=0.5)+
  geom_line()

# temp running 7d ave
window %>%
  ggplot(aes(yday, temp_roll_mean))+
  facet_wrap(~site_yr_layer, nrow = 5)+
  geom_ribbon(aes(ymin = temp_roll_mean - 5*temp_roll_sd, ymax = temp_roll_mean + 5*temp_roll_sd), alpha=0.5)+
  geom_line()


mini3 = mini2 %>% 
  # adjoin the 7-day running ave
  left_join(window) %>%
  mutate(remove_do = ifelse((do <= do_roll_mean - 3*do_roll_sd | do >= do_roll_mean + 3*do_roll_sd), 1, 0),
         remove_temp = ifelse((temp <= temp_roll_mean - 5*temp_roll_sd | temp >= temp_roll_mean + 5*temp_roll_sd), 1, 0)) 

# DO per site-yr-layer; showing the running 7d ave with shaded sd; red points = would be filtered
mini3 %>%
  ggplot()+
  facet_wrap(~site_yr_layer, nrow=5)+
  geom_ribbon(aes(x=yday, ymin = do_roll_mean - 3*do_roll_sd, ymax = do_roll_mean + 3*do_roll_sd), alpha=0.7, fill='gray')+
  geom_point(aes(x=yday, y=do, color=factor(remove_do)), size=0.5, alpha=0.5)+
  geom_line(aes(x=yday, y=do_roll_mean), color="blue")+
  scale_color_manual(values=c("black", "red"))+
  theme_classic()+
  theme(legend.position = 'none')

# temp per site-yr-layer; showing the running 7d ave with shaded sd; red points = would be filtered
# note: the only temps being removed are clear malfunctions from the sonde
mini3 %>%
  ggplot()+
  facet_wrap(~site_yr_layer, nrow=5)+
  geom_ribbon(aes(x=yday, ymin = temp_roll_mean - 5*temp_roll_sd, ymax = temp_roll_mean + 5*temp_roll_sd), alpha=0.7, fill='gray')+
  geom_point(aes(x=yday, y=temp, color=factor(remove_temp)), size=0.5, alpha=0.5)+
  geom_line(aes(x=yday, y=temp_roll_mean), color="blue")+
  scale_color_manual(values=c("black", "red"))+
  theme_classic()+
  theme(legend.position = 'none')

mini3 %>% filter(remove_temp == 1) %>% print(n=Inf)


# COMPARE TO OLD APPROACH FOR CLEANING DO BASED ON 2SD OF THAT DAY'S MEAN
# mini_a will filter do points based on whether they are within 2sd of that day's mean
# mini_a = mini2 %>% 
#   left_join(mini2 %>% 
#               group_by(site, layer, date, year) %>%
#               summarize(do_mean = mean(do, na.rm=T),
#                         do_sd = sd(do, na.rm=T)) %>%
#               ungroup()) %>%
#   mutate(remove = ifelse((do <= do_mean - 3*do_sd | do >= do_mean + 3*do_sd), 1, 0))
# 
# mini_a %>%
#   ggplot()+
#   facet_wrap(~site_yr_layer, nrow=5)+
#   geom_ribbon(aes(x=yday, ymin = do_mean - 3*do_sd, ymax = do_mean + 3*do_sd), alpha=0.7, fill='gray')+
#   geom_point(aes(x=yday, y=do, color=factor(remove)), size=0.5, alpha=0.5)+
#   geom_line(aes(x=yday, y=do_mean), color="blue")+
#   scale_color_manual(values=c("black", "red"))+
#   theme_classic()+
#   theme(legend.position = 'none')
# 
# mini_agg_a = mini_a %>%
#   # filter out points that are not within 2SD of that day's mean 
#   filter(do <= do_mean + 2*do_sd, do >= do_mean - 2*do_sd) %>%
#   group_by(site, layer, site_yr_layer, date, yday, year, hour) %>%
#   summarize(temp = mean(temp, na.rm=T),
#             do_eq = mean(do_eq, na.rm=T),
#             do = mean(do, na.rm=T),
#             do_sat = mean(do_sat, na.rm=T)) %>%
#   ungroup() 

# #### strange st33 2019 drop in benthic do
# mini_agg_a %>% filter(site == 'st33', year==2019, yday==190) %>%
#   ggplot(aes(hour, do)) +
#   geom_point() +
#   geom_line()+
#   ylab("benthic DO") +
#   theme(axis.title.x = element_blank())
# 
# # the drop in do on day 190 really does look wonky, so I am going to hack it and filter it from 'mini_a' even though it is not actually outside of 2sd of the daily mean (because the day was so variable)
# clean_a %>% filter(site == 'st33', year==2019, yday==190) %>% print(n=Inf) 
# clean_a = clean_a %>%
#   mutate(do.ben = ifelse((site == 'st33' & year==2019 & yday==190 & hour>=10 & hour <= 13), NA, do.ben)) 


# aggregate oxygen and temp data to hour
mini_agg = mini3 %>% 
  filter(do <= do_roll_mean + 3*do_roll_sd, do >= do_roll_mean - 3*do_roll_sd) %>%
  filter(temp <= temp_roll_mean + 5*temp_roll_sd, temp >= temp_roll_mean - 5*temp_roll_sd) %>%
  group_by(site, layer, site_yr_layer, date, yday, year, hour) %>%
  summarize(temp = mean(temp, na.rm=T),
            do_eq = mean(do_eq, na.rm=T),
            do = mean(do, na.rm=T),
            do_sat = mean(do_sat, na.rm=T)) %>%
  ungroup() 

mini_agg %>%
  ggplot(aes(yday+hour/24, do)) +
  facet_wrap(~site_yr_layer, nrow=5)+
  geom_line()





# adjoin the light data (does not include light for 2017 ben and pel or 2019 pel)
mini_agg_hobo = mini_agg %>% left_join(hobo_agg) 

# create wide-format data with separate columns for ben and pel and adjoin the weather data
clean = mini_agg_hobo %>%
  select(-site_yr_layer) %>%
  pivot_wider(names_from = layer, values_from = c(temp, do_eq, do, do_sat, par), names_sep=".") %>%
  left_join(weather %>% select(date, year, hour, wspeed))


# adjoin missing 2017 and 2019 par data
full_data = clean %>% 
  left_join(par_data_17 %>% 
              rename(par.ben.est17 = par.ben, par.pel.est17 = par.pel)) %>%
  left_join(par_data_19 %>%
              rename(par.ben.est19 = par.ben, par.pel.est19 = par.pel)) %>% 
  mutate(par.ben = ifelse((year==2017 & site=="st33"), par.ben.est17, par.ben),
         par.pel = ifelse((year==2017 & site=="st33"), par.pel.est17, 
                          ifelse((year==2019 & site=="st33"), par.pel.est19, par.pel))) %>%
  select(-c(par.ben.est17:par.pel.est19))




# full plot of do, par, and temp data
full_data %>% 
  select(site, yday, year, hour, temp.ben, temp.pel, do.ben, do.pel, par.ben, par.pel) %>%
  gather(var, val, temp.ben:par.pel) %>%
  mutate(layer = ifelse(var %in% c("temp.ben", "do.ben", "par.ben"), 'ben', 'pel'),
         panel = ifelse(var %in% c("temp.ben", "temp.pel"), "temp",
                        ifelse(var %in% c("do.ben", "do.pel"), "do", "par")),
         site_yr = paste(site, year, sep="_")) %>%
  ggplot(aes(yday+hour/24, val, color=layer))+
  facet_grid(panel~site_yr, scales = 'free_y')+
  geom_line(alpha=0.7)+
  scale_color_manual(values = c("black", "deepskyblue3")) +
  theme(legend.position = 'none', axis.title = element_blank())

#write_csv(full_data, "data/clean_coupled_o2/mini_clean_22apr22.csv")


