library(tidyverse)
library(ggplot2)
library(readr)
library(astsa)
library(lubridate)

d <- read_csv("~/.julia/dev/DataFusion/csv/fused_daily_data_log.csv")
glimpse(d)
lag2.plot(d$temp, d$fused_mean_log, max.lag=2)

d <- d %>% mutate(lograd = log(rad))

dsub <- d %>% select(time, fused_mean_log, temp, lograd) %>% gather(key=type, value=y, fused_mean_log, temp, lograd) 


p1 <- dsub %>% 
  filter(year(time)>2010) %>%
  ggplot(aes(x=time,y=y, colour=type)) +
  geom_line()
p1

p1 + facet_wrap(~type, scales="free")

plot(d$time, d$fused_mean)

dsub2 <- d #%>% filter(year(time)>2010) 
lag2.plot(dsub2$temp, dsub2$fused_mean_log, max.lag=3)
lag2.plot(diff(dsub2$temp), dsub2$fused_mean_log[-1], max.lag=3)
ccf(dsub2$temp, dsub2$fused_mean)
# differencing
lag2.plot(diff(dsub2$temp), diff(dsub2$fused_mean_log), max.lag=3)
ccf(dsub2$temp, dsub2$fused_mean)



lag2.plot(dsub2$lograd, dsub2$fused_mean_log, max.lag=3)
ccf(dsub2$lograd, dsub2$fused_mean)



