library(readr)
library(tidyverse)
library(DataExplorer)
library(tibble)
library(lubridate)

setwd("~/.julia/dev/Examples/datafusion/processing_in_r")
dorig <-  read_csv("daily_data.csv", 
                            col_types = cols(chl = col_double())) %>% rename(chl_water=concentration_of_chlorophyll_in_water) 

d <- dorig %>% gather(chl, chl_water,key='meastype', value=y)
#View(d)
glimpse(d)

d %>% ggplot(aes(x=time, y=y,colour=meastype)) + geom_point() + facet_wrap(~meastype,ncol=1)
d %>% filter(meastype=='chl') %>% ggplot(aes(x=time, y=y,colour=meastype)) + geom_point() + facet_wrap(~meastype)


firsttime <- "1976-01-01"
elapsed <- interval(firsttime, d$time) %>%as.numeric('years')

d <- d %>% mutate(time_elapsed = elapsed)
pl <- d %>% ggplot(aes(x=time_elapsed, y=y,colour=meastype)) + geom_point(size=0.6,alpha=0.6) +
  facet_wrap(~meastype,ncol=1,scales='free') + ylab("concentration")
pl
pdf("~/.julia/dev/Examples/datafusion/figs/vis_data.pdf",width=7,height=5)
  pl
dev.off()

elapsed2 <- interval(firsttime, dorig$time) %>%as.numeric('years')
dorig <- dorig %>% mutate(time_elapsed=elapsed2, chl_water_meas = !is.na(chl_water), chl_meas=!is.na(chl)) %>%
  unite("obsscheme", chl_water_meas, chl_meas) %>% 
  mutate(obsscheme=fct_recode(obsscheme,obs1="TRUE_FALSE", obs2="FALSE_TRUE", obs3="TRUE_TRUE")) %>%
  arrange(time_elapsed)
  

glimpse(dorig)
dorig %>% filter(obsscheme=='type3')

write_csv(dorig, paste0(getwd(),"/periodic1d/observations.csv"))

#unique(dorig$obsscheme)


# conventions
# - type 1: only chl_water measured
# - type 2: only chl measured
# - type 3: both measured

# include posterior mean in figure
pm <- read_csv("postmean_paths.csv", col_names = FALSE) 
d <- d %>% mutate(postmean=rep(pm$X1,2))
pl2 <- d %>% ggplot(aes(x=time_elapsed, y=y,colour=meastype)) + geom_point(size=0.6,alpha=0.9) +
  facet_wrap(~meastype,ncol=1,scales='free') + ylab("concentration") + 
  geom_path(mapping=aes(x=time_elapsed,y=postmean),colour='black',size=0.3,alpha=0.8)
pl2
pdf("~/.julia/dev/Examples/datafusion/figs/vis_data_fit.pdf",width=7,height=5)
pl2
dev.off()

pl3 <- d %>% ggplot(aes(x=time_elapsed, y=y,colour=meastype)) + geom_point(size=0.6,alpha=0.9) +
   ylab("concentration") + 
  geom_path(mapping=aes(x=time_elapsed,y=postmean),colour='black',size=0.3,alpha=0.8)
pl3
pdf("~/.julia/dev/Examples/datafusion/figs/vis_data_fit2.pdf",width=7,height=5)
pl3
dev.off()

d1 <- d %>% mutate(t1 = time_elapsed %% 1) 
pl4 <- d1  %>% ggplot(aes(x=t1, y=y,colour=meastype)) + geom_point(size=0.6,alpha=0.9) + geom_smooth(colour='blue') + facet_wrap(~meastype)+ 
  xlab("time elapsed") + ylab("concentration") + geom_hline(yintercept = 0)
pdf("~/.julia/dev/Examples/datafusion/figs/vis_data_periodic.pdf",width=7,height=5)
pl4
dev.off()



 }

out <- loess(y ~t1, d1)
