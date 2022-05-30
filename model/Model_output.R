library(tidyverse)
library(lubridate)


############################
## Import output from models
############################

# type = "st33_2017_gamma0.05"
# model = "coupled_o2_update_sd_Iopt_lmtb_c01_gamma"
# 
# output_rds = readRDS(paste0("model/output/", model, "/", type,"/output.rds"))
# 
# output2 = rstan::summary(output_rds)
# as.data.frame(output2$summary) %>% rownames_to_column() %>% arrange(desc(Rhat)) %>% filter(Rhat > 1.05)
# as.data.frame(output2$summary) %>% rownames_to_column() %>% filter(is.na(Rhat)) %>% as_tibble() %>% print(n=Inf)
# 
# d = as.data.frame(output2$summary) %>% rownames_to_column()






model = "coupled_o2_test"

type = "st33_2019"

# model output (daily)
daily = read_csv(paste0("model/output/", model, "/", type,"/daily_summary.csv")) 

#daily = read_csv(paste0("model/output/", model, "/", type, "_gamma0.05","/daily_summary.csv")) 

# model output (hourly)
hourly = read_csv(paste0("model/output/", model, "/", type, "_gamma0.05", "/hourly_summary.csv")) 

fixed_pars = read_csv(paste0("model/output/", model, "/", type, "_gamma0.05", "/fixed_pars.csv")) 

fixed_summary = read_csv(paste0("model/output/", model, "/", type, "_gamma0.05", "/fixed_summary.csv")) 


# for figs
theme_set(theme_bw() %+replace% 
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  #legend.margin = margin(0,0,0,0),
                  strip.text = element_text(size=14),
                  legend.text = element_text(size=12),
                  legend.title = element_text(size=14),
                  axis.text = element_text(size=12, color="black"),
                  #axis.title.y = element_text(size=12,angle = 90 ,margin=margin(0,15,0,0)),
                  #axis.title.x = element_text(size=12,margin=margin(15,0,0,0)),
                  strip.text.x = element_text(margin=margin(0,0,10,0)), 
                  strip.text.y = element_text(margin=margin(0,0,0,10), angle=270),
                  axis.title = element_text(size=14))
)



########
# DAILY model output

## daily figure
daily %>%
  filter(name %in% c("beta0", "rho", "GPP", "ER", "NEP")) %>%
  mutate(name = factor(name, levels = c("beta0", "rho", "GPP", "ER", "NEP"))) %>%
  mutate(layer=factor(layer, levels=c("pelagic", "benthic"))) %>%
  ggplot(aes(x = yday, y = middle, color=layer))+
  facet_grid(name~., scales = "free_y") +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=layer), 
              color = NA, alpha = 0.3)+
  geom_line(size = 0.8)+
  geom_hline(aes(yintercept=0), lty=2, alpha=0.5)+
  scale_color_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  scale_fill_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  theme(axis.title.x = element_blank(), legend.position = 'top', strip.text.y = element_text(size=11))



# beta0
daily %>% 
  filter(name=='beta0') %>%
  mutate(layer=factor(layer, levels=c("pelagic", "benthic"))) %>%
  ggplot(aes(x = yday, y = middle, color=layer))+
  facet_grid(~year)+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=layer), 
              color = NA, alpha = 0.3)+
  geom_line(size = 0.8)+
  geom_hline(aes(yintercept=0), lty=2, alpha=0.5)+
  scale_y_continuous("beta0")+
  scale_color_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  scale_fill_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  theme(axis.title.x = element_blank(), legend.position = 'top', strip.text.y = element_text(size=11))

#rho
daily %>% 
  filter(name=='rho') %>%
  mutate(layer=factor(layer, levels=c("pelagic", "benthic"))) %>%
  ggplot(aes(x = yday, y = middle, color=layer))+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=layer), 
              color = NA, alpha = 0.3)+
  geom_line(size = 0.8)+
  geom_hline(aes(yintercept=0), lty=2, alpha=0.5)+
  scale_y_continuous("rho")+
  scale_color_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  scale_fill_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  theme(axis.title.x = element_blank(), legend.position = 'top', strip.text.y = element_text(size=11))

#GPP
daily %>% 
  filter(name=='GPP') %>%
  mutate(layer=factor(layer, levels=c("pelagic", "benthic"))) %>%
  ggplot(aes(x = yday, y = middle, color=layer))+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=layer), 
              color = NA, alpha = 0.3)+
  geom_line(size = 0.8)+
  geom_hline(aes(yintercept=0), lty=2, alpha=0.5)+
  scale_y_continuous("GPP")+
  scale_color_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  scale_fill_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  theme(axis.title.x = element_blank(), legend.position = 'top', strip.text.y = element_text(size=11))


#NEP
daily %>% 
  filter(name=='NEP') %>%
  mutate(layer=factor(layer, levels=c("pelagic", "benthic"))) %>%
  ggplot(aes(x = yday, y = middle, color=layer))+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=layer), 
              color = NA, alpha = 0.3)+
  geom_line(size = 0.8)+
  geom_hline(aes(yintercept=0), lty=2, alpha=0.5)+
  scale_y_continuous("NEP")+
  scale_color_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  scale_fill_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  theme(axis.title.x = element_blank(), legend.position = 'top', strip.text.y = element_text(size=11))

# ER
daily %>% 
  filter(name=='ER') %>%
  mutate(layer=factor(layer, levels=c("pelagic", "benthic"))) %>%
  ggplot(aes(x = yday, y = middle, color=layer))+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=layer), 
              color = NA, alpha = 0.3)+
  geom_line(size = 0.8)+
  geom_hline(aes(yintercept=0), lty=2, alpha=0.5)+
  scale_y_continuous("ER")+
  scale_color_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  scale_fill_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  theme(axis.title.x = element_blank(), legend.position = 'top', strip.text.y = element_text(size=11))



# x_pred; D
hourly %>%
  filter(name == "x_pred") %>%
  mutate(layer=factor(layer, levels=c("pelagic", "benthic"))) %>%
  ggplot(aes(x = yday+hour/24, y = middle, color=layer))+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=layer), 
              color = NA, alpha = 0.3)+
  geom_line(size = 0.8)+
  geom_hline(aes(yintercept=0), lty=2, alpha=0.5)+
  scale_y_continuous("x_pred")+
  scale_color_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  scale_fill_manual("",values=c("deepskyblue3","black"), labels=c("Pelagic", "Benthic"))+
  theme(axis.title.x = element_blank(), legend.position = 'top', strip.text.y = element_text(size=11))

