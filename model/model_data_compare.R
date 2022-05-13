

library(tidyverse)
library(lubridate)
library(stringr)

# for figs
theme_set(theme_bw() %+replace% 
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  #legend.margin = margin(0,0,0,0),
                  strip.text = element_text(size=14),
                  legend.text = element_text(size=12),
                  legend.title = element_text(size=14),
                  axis.text = element_text(size=10, color="black"),
                  #axis.title.y = element_text(size=12,angle = 90 ,margin=margin(0,15,0,0)),
                  #axis.title.x = element_text(size=12,margin=margin(15,0,0,0)),
                  strip.text.x = element_text(margin=margin(0,0,10,0)), 
                  strip.text.y = element_text(margin=margin(0,0,0,10), angle=270),
                  axis.title = element_text(size=16),
                  plot.title = element_text(hjust = 0.5, size=16))
)

# examiine rds output
type = "reyk_2019"
gammaSD = 0.25
model = "coupled_o2_update_sd_Iopt_lmtb_c01_gamma"

output_rds = readRDS(paste0("model/output/", model, "/", type, "_gamma", gammaSD,"/output.rds"))

output2 = rstan::summary(output_rds)
#d = as.data.frame(output2$summary) %>% rownames_to_column()
#head(d)
#d$parameter = word(d$rowname,1,sep = "\\[")

as.data.frame(output2$summary) %>% rownames_to_column() %>% mutate(param = word(rowname, 1, sep = "\\[")) %>% 
  arrange(desc(Rhat)) %>% filter(Rhat > 1.05) #%>% group_by(param) %>% tally() %>% print(n=Inf)

as.data.frame(output2$summary) %>% rownames_to_column() %>% filter(is.na(Rhat)) %>% as_tibble() %>% print(n=Inf)




### MAY 2022 figures for different values of the gamma sd prior

daily_comp = bind_rows(
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2017_gamma0.05","/daily_summary.csv")) %>% mutate(gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2017_gamma0.1","/daily_summary.csv")) %>% mutate(gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2017_gamma0.25","/daily_summary.csv")) %>% mutate(gamma_sd = 0.25),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2018_gamma0.05","/daily_summary.csv")) %>% mutate(gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2018_gamma0.1","/daily_summary.csv")) %>% mutate(gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2018_gamma0.25","/daily_summary.csv")) %>% mutate(gamma_sd = 0.25),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2019_gamma0.05","/daily_summary.csv")) %>% mutate(gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2019_gamma0.1","/daily_summary.csv")) %>% mutate(gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2019_gamma0.25","/daily_summary.csv")) %>% mutate(gamma_sd = 0.25),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2018_gamma0.05","/daily_summary.csv")) %>% mutate(gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2018_gamma0.1","/daily_summary.csv")) %>% mutate(gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2018_gamma0.25","/daily_summary.csv")) %>% mutate(gamma_sd = 0.25),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2019_gamma0.05","/daily_summary.csv")) %>% mutate(gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2019_gamma0.1","/daily_summary.csv")) %>% mutate(gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2019_gamma0.25","/daily_summary.csv")) %>% mutate(gamma_sd = 0.25)
  ) 


# multi panel for each site-year (gamma sd by parameter)
type = "reyk_2018"

daily_comp %>% 
  mutate(site_yr = paste(site, year, sep = "_")) %>%
  filter(site_yr==type) %>%
  mutate(layer=factor(layer, levels=c("pelagic", "benthic"))) %>%
  filter(name %in% c("beta0", "rho", "GPP", "ER", "NEP")) %>%
  mutate(name = factor(name, levels = c("beta0", "rho", "GPP", "ER", "NEP"))) %>%
  ggplot(aes(x = yday, y = middle, color=layer))+
  facet_grid(name~gamma_sd, scales="free")+ 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=layer), 
              color = NA, alpha = 0.3)+
  geom_line(size = 0.8)+
  geom_hline(aes(yintercept=0), lty=2, alpha=0.5)+
  ggtitle(type)+
  scale_color_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  scale_fill_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  theme(axis.title = element_blank(), legend.position = 'none', #strip.text.y = element_text(size=11),
        plot.title = element_text(margin = margin(0, 0, 20, 0)))



fixed_pars = bind_rows(
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2017_gamma0.05","/fixed_pars.csv")) %>% mutate(site_yr = "st33_2017", gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2017_gamma0.1","/fixed_pars.csv")) %>% mutate(site_yr = "st33_2017", gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2017_gamma0.25","/fixed_pars.csv")) %>% mutate(site_yr = "st33_2017", gamma_sd = 0.25),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2018_gamma0.05","/fixed_pars.csv")) %>% mutate(site_yr = "st33_2018", gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2018_gamma0.1","/fixed_pars.csv")) %>% mutate(site_yr = "st33_2018", gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2018_gamma0.25","/fixed_pars.csv")) %>% mutate(site_yr = "st33_2018", gamma_sd = 0.25),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2019_gamma0.05","/fixed_pars.csv")) %>% mutate(site_yr = "st33_2019", gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2019_gamma0.1","/fixed_pars.csv")) %>% mutate(site_yr = "st33_2019", gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2019_gamma0.25","/fixed_pars.csv")) %>% mutate(site_yr = "st33_2019", gamma_sd = 0.25),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2018_gamma0.05","/fixed_pars.csv")) %>% mutate(site_yr = "reyk_2018", gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2018_gamma0.1","/fixed_pars.csv")) %>% mutate(site_yr = "reyk_2018", gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2018_gamma0.25","/fixed_pars.csv")) %>% mutate(site_yr = "reyk_2018", gamma_sd = 0.25),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2019_gamma0.05","/fixed_pars.csv")) %>% mutate(site_yr = "reyk_2019", gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2019_gamma0.1","/fixed_pars.csv")) %>% mutate(site_yr = "reyk_2019", gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2019_gamma0.25","/fixed_pars.csv")) %>% mutate(site_yr = "reyk_2019", gamma_sd = 0.25)
) 

type = "reyk_2019"

fixed_pars %>% select(site_yr, chain, step, gamma_sd, 'gamma_1[1]':'c1') %>%
  filter(site_yr == "reyk_2018") %>%
  gather(var, val, 'gamma_1[1]':'c1') %>%
  ggplot(aes(val))+
  facet_grid(var ~ gamma_sd, scales = 'free')+
  stat_density(geom = "line")+
  ggtitle(type)+
  theme(strip.text.y = element_text(size=10),
        plot.title = element_text(margin = margin(0, 0, 20, 0)))

fixed_pars %>% select(site_yr, chain, step, gamma_sd, 'gamma_1[1]':'c1') %>%
  filter(site_yr == type) %>%
  gather(var, val, 'gamma_1[1]':'c1') %>%
  ggplot(aes(step, val, color = factor(chain)))+
  facet_grid(var ~ gamma_sd, scales = 'free')+
  geom_line(alpha = 0.35)+
  ggtitle(type)+
  theme(legend.position = "none", strip.text.y = element_text(size=10),
        plot.title = element_text(margin = margin(0, 0, 20, 0)))


fixed_summary = bind_rows(
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2017_gamma0.05","/fixed_summary.csv")) %>% mutate(site_yr = "st33_2017", gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2017_gamma0.1","/fixed_summary.csv")) %>% mutate(site_yr = "st33_2017", gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2017_gamma0.25","/fixed_summary.csv")) %>% mutate(site_yr = "st33_2017", gamma_sd = 0.25),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2018_gamma0.05","/fixed_summary.csv")) %>% mutate(site_yr = "st33_2018", gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2018_gamma0.1","/fixed_summary.csv")) %>% mutate(site_yr = "st33_2018", gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2018_gamma0.25","/fixed_summary.csv")) %>% mutate(site_yr = "st33_2018", gamma_sd = 0.25),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2019_gamma0.05","/fixed_summary.csv")) %>% mutate(site_yr = "st33_2019", gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2019_gamma0.1","/fixed_summary.csv")) %>% mutate(site_yr = "st33_2019", gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2019_gamma0.25","/fixed_summary.csv")) %>% mutate(site_yr = "st33_2019", gamma_sd = 0.25),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2018_gamma0.05","/fixed_summary.csv")) %>% mutate(site_yr = "reyk_2018", gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2018_gamma0.1","/fixed_summary.csv")) %>% mutate(site_yr = "reyk_2018", gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2018_gamma0.25","/fixed_summary.csv")) %>% mutate(site_yr = "reyk_2018", gamma_sd = 0.25),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2019_gamma0.05","/fixed_summary.csv")) %>% mutate(site_yr = "reyk_2019", gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2019_gamma0.1","/fixed_summary.csv")) %>% mutate(site_yr = "reyk_2019", gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2019_gamma0.25","/fixed_summary.csv")) %>% mutate(site_yr = "reyk_2019", gamma_sd = 0.25)
) 

fixed_summary %>%
  ggplot(aes(factor(gamma_sd), middle))+
  facet_grid(var ~ site_yr, scales = 'free')+
  geom_point()+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3)+
  theme(strip.text.y = element_text(size=9), axis.title.x = element_blank(), axis.text.y = element_text(size=9))
  

hourly = bind_rows(
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2017_gamma0.05","/hourly_summary.csv")) %>% mutate(site_yr = "st33_2017", gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2017_gamma0.1","/hourly_summary.csv")) %>% mutate(site_yr = "st33_2017", gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2017_gamma0.25","/hourly_summary.csv")) %>% mutate(site_yr = "st33_2017", gamma_sd = 0.25),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2018_gamma0.05","/hourly_summary.csv")) %>% mutate(site_yr = "st33_2018", gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2018_gamma0.1","/hourly_summary.csv")) %>% mutate(site_yr = "st33_2018", gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2018_gamma0.25","/hourly_summary.csv")) %>% mutate(site_yr = "st33_2018", gamma_sd = 0.25),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2019_gamma0.05","/hourly_summary.csv")) %>% mutate(site_yr = "st33_2019", gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2019_gamma0.1","/hourly_summary.csv")) %>% mutate(site_yr = "st33_2019", gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "st33_2019_gamma0.25","/hourly_summary.csv")) %>% mutate(site_yr = "st33_2019", gamma_sd = 0.25),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2018_gamma0.05","/hourly_summary.csv")) %>% mutate(site_yr = "reyk_2018", gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2018_gamma0.1","/hourly_summary.csv")) %>% mutate(site_yr = "reyk_2018", gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2018_gamma0.25","/hourly_summary.csv")) %>% mutate(site_yr = "reyk_2018", gamma_sd = 0.25),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2019_gamma0.05","/hourly_summary.csv")) %>% mutate(site_yr = "reyk_2019", gamma_sd = 0.05),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2019_gamma0.1","/hourly_summary.csv")) %>% mutate(site_yr = "reyk_2019", gamma_sd = 0.1),
  read_csv(paste0("model/output/", "coupled_o2_update_sd_Iopt_lmtb_c01_gamma", "/", "reyk_2019_gamma0.25","/hourly_summary.csv")) %>% mutate(site_yr = "reyk_2019", gamma_sd = 0.25)
) 

hourly
unique(hourly$name)

type = "st33_2017"

hourly %>%
  filter(name == "x_pred") %>%
  mutate(site = factor(site, levels=c('st33', 'reyk'))) %>%
  mutate(layer=factor(layer, levels=c("pelagic", "benthic"))) %>%
  ggplot(aes(x = yday+hour/24, y = middle, color=layer))+
  facet_grid(site~year)+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=layer), 
              color = NA, alpha = 0.3)+
  geom_line(size = 0.8, alpha = 0.4)+
  geom_hline(aes(yintercept=0), lty=2, alpha=0.5)+
  scale_y_continuous("x_pred")+
  scale_color_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  scale_fill_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  theme(axis.title.x = element_blank(), legend.position = 'top', strip.text.y = element_text(size=11))


hourly %>% 
  filter(site_yr==type) %>%
  mutate(layer=factor(layer, levels=c("pelagic", "benthic"))) %>%
  filter(name %in% c( "gpp",  "er", "nep")) %>%
  #filter(name %in% c("b", "chi", "kappa", "phi", "D", "x_pred")) %>%
  group_by(layer, yday, name, gamma_sd) %>%
    summarize(lower=mean(lower),
              middle=mean(middle),
              upper=mean(upper)) %>%
  ggplot(aes(x = yday, y = middle, color=layer))+
  facet_grid(name~gamma_sd, scales="free")+ 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=layer), 
              color = NA, alpha = 0.3)+
  geom_line(size = 0.8)+
  geom_hline(aes(yintercept=0), lty=2, alpha=0.5)+
  ggtitle(type)+
  scale_color_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  scale_fill_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  theme(axis.title = element_blank(), legend.position = 'none', #strip.text.y = element_text(size=11),
        plot.title = element_text(margin = margin(0, 0, 20, 0)))




