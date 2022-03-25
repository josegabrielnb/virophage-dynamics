library(tidyverse)
library(ggplot2)
library(mgcv)
library(gridExtra)

setwd("~/Desktop/Julia2.0/simulation_results/")

dynamics_summary <- read_csv("dynamics_summary_100replicates.csv")

dynamics_summary <- arrange(dynamics_summary,cell_number,treatment,run)

dynamics_summary %>% group_by(cell_number,treatment) %>%
  summarise(mean_extinction = mean(t_extinction), sd_extinction = sd(t_extinction))

p1 <- dynamics_summary %>% 
  ggplot(aes(x=factor(cell_number), y=t_extinction, fill=factor(treatment))) +
  geom_boxplot(show.legend = F) +
  xlab("Cell number") +
  ylab("Time to virus extinction") +
  theme(axis.title.y = element_text(size=12), axis.title.x = element_text(size=12)) +
  theme_bw()

p2 <- dynamics_summary %>% 
  ggplot(aes(x=cell_number, y=t_extinction, colour=factor(treatment))) +
  geom_point() +
  stat_smooth(method = "gam",formula = y~s(x,bs="fs",k=5)+x) +
  xlab("Cell number") +
  ylab("") +
  theme(axis.title.x = element_text(size=12)) +
  theme_bw() +
  scale_colour_discrete(name="Treatment",
                        breaks=c("1", "2", "3"),
                        labels=c("Neutral", "Inhibition", "PCD"))

grid.arrange(p1, p2, nrow = 1)

p3 <- dynamics_summary %>% 
  ggplot(aes(x=factor(cell_number), y=max_virus_amplitude, fill=factor(treatment))) +
  geom_boxplot(show.legend = F) +
  xlab("Cell number") +
  ylab("Maximum virus wave amplitude") +
  theme(axis.title.y = element_text(size=12), axis.title.x = element_text(size=12)) +
  theme_bw()

p4 <- dynamics_summary %>% 
  ggplot(aes(x=cell_number, y=max_virus_amplitude, colour=factor(treatment))) +
  geom_point() +
  stat_smooth(method = "gam",formula = y~s(x,bs="fs",k=5)+x) +
  xlab("Cell number") +
  ylab("") +
  theme(axis.title.x = element_text(size=12)) +
  theme_bw() +
  scale_colour_discrete(name="Treatment",
                        breaks=c("1", "2", "3"),
                        labels=c("Neutral", "Inhibition", "PCD"))

grid.arrange(p3, p4, nrow = 1)

p5 <- dynamics_summary %>% 
  ggplot(aes(x=factor(cell_number), y=surviving_cells, fill=factor(treatment))) +
  geom_boxplot(show.legend = F) +
  xlab("Cell number") +
  ylab("Surviving cells") +
  theme(axis.title.y = element_text(size=12), axis.title.x = element_text(size=12)) +
  theme_bw()

p6 <- dynamics_summary %>% 
  ggplot(aes(x=cell_number, y=surviving_cells, colour=factor(treatment))) +
  geom_point() +
  stat_smooth(method = "gam",formula = y~s(x,bs="fs",k=5)) +
  xlab("Cell number") +
  ylab("") +
  theme(axis.title.x = element_text(size=12)) +
  theme_bw() +
  scale_colour_discrete(name="Treatment",
                        breaks=c("1", "2", "3"),
                        labels=c("Neutral", "Inhibition", "PCD"))

grid.arrange(p5, p6, nrow = 1)


#General Additive Model

#Time to virus extinction

abm_gam2 <- gam(t_extinction~s(cell_number,by=factor(treatment),k=5)+factor(treatment),data=dynamics_summary)
summary(abm_gam2)
plot(abm_gam2, all.terms = TRUE, pages=1, shift=abm_gam2$coefficients[1], seWithMean=TRUE)
gam.check(abm_gam2)
concurvity(abm_gam2)
sqrt(mean(residuals.gam(abm_gam2)^2))
sqrt(mean(residuals.gam(abm_gam2)^2))/(max(dynamics_summary$t_extinction)-min(dynamics_summary$t_extinction))

#Maximum virus wave amplitude

abm_gam5 <- gam(max_virus_amplitude~s(cell_number,by=factor(treatment),k=5)+factor(treatment),data=dynamics_summary)
summary(abm_gam5)
plot(abm_gam5, all.terms = TRUE, pages=1, shift=abm_gam5$coefficients[1], seWithMean=TRUE)
gam.check(abm_gam5)
concurvity(abm_gam5)
sqrt(mean(residuals.gam(abm_gam5)^2))
sqrt(mean(residuals.gam(abm_gam5)^2))/(max(dynamics_summary$max_virus_amplitude)-min(dynamics_summary$max_virus_amplitude))

#Number of surviving cells

abm_gam8 <- gam(surviving_cells~s(cell_number,by=factor(treatment),k=5)+factor(treatment),data=dynamics_summary)
summary(abm_gam8)
plot(abm_gam8, all.terms = TRUE, pages=1, shift=abm_gam8$coefficients[1], seWithMean=TRUE)
gam.check(abm_gam8)
concurvity(abm_gam8)
sqrt(mean(residuals.gam(abm_gam8)^2))
sqrt(mean(residuals.gam(abm_gam8)^2))/(max(dynamics_summary$max_virus_amplitude)-min(dynamics_summary$max_virus_amplitude))
