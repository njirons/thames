rm (list=ls())
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

library(invgamma)
library (MASS)
attach (nlschools)
library (lme4)
library (arm)
library(bridgesampling)
library(mvtnorm)
library(tidyr)
library(dplyr)

lml_thames_lm <- read.csv('data/lml_thames_lm.csv',header=FALSE,sep='')
lml_thames_reduced <- read.csv('data/lml_thames_reduced.csv',header=FALSE,sep='')
lml_thames_full <- read.csv('data/lml_thames_full.csv',header=FALSE,sep='')

lml_bridge_reduced <- read.csv('data/lml_bridge_reduced.csv',header=FALSE,sep='')
lml_bridge_full <- read.csv('data/lml_bridge_full.csv',header=FALSE,sep='')
lml_bridge_lm <- read.csv('data/lml_bridge_lm.csv',header=FALSE,sep='')

lml_mc_reduced <- read.csv('data/lml_mc_reduced.csv',header=FALSE,sep='')
lml_mc_full <- read.csv('data/lml_mc_full.csv',header=FALSE,sep='')
lml_mc_lm <- read.csv('data/lml_mc_lm.csv',header=FALSE,sep='')

iters <- round(seq(5,5000,length.out=20)); 4*iters

# plot for paper

# set theme
library(extrafont) 
my_theme <- theme_bw() +
  theme(strip.background = element_rect(fill = "white"), 
        text = element_text(face="bold", size=16),
  )
theme_set(my_theme)
require(ggsci)

# ggplot
inds <- 1:20
results <- tibble(
  thames_full = lml_thames_full$V1[inds],
  thames_full_L = lml_thames_full$V2[inds],
  thames_full_U = lml_thames_full$V3[inds],
  thames_reduced = lml_thames_reduced$V1[inds],
  thames_reduced_L = lml_thames_reduced$V2[inds],
  thames_reduced_U = lml_thames_reduced$V3[inds],
  thames_lm = lml_thames_lm$V1[inds],
  thames_lm_L = lml_thames_lm$V2[inds],
  thames_lm_U = lml_thames_lm$V3[inds],
  bridge_full = lml_bridge_full$V1[inds],
  bridge_full_L = lml_bridge_full$V2[inds],
  bridge_full_U = lml_bridge_full$V3[inds],
  bridge_reduced = lml_bridge_reduced$V1[inds],
  bridge_reduced_L = lml_bridge_reduced$V2[inds],
  bridge_reduced_U = lml_bridge_reduced$V3[inds],
  bridge_lm = lml_bridge_lm$V1[inds],
  bridge_lm_L = lml_bridge_lm$V2[inds],
  bridge_lm_U = lml_bridge_lm$V3[inds],
  mc_full = lml_mc_full$V1[inds],
  mc_full_L = lml_mc_full$V2[inds],
  mc_full_U = lml_mc_full$V3[inds],
  mc_reduced = lml_mc_reduced$V1[inds],
  mc_reduced_L = lml_mc_reduced$V2[inds],
  mc_reduced_U = lml_mc_reduced$V3[inds],
  mc_lm = lml_mc_lm$V1[inds],
  mc_lm_L = lml_mc_lm$V2[inds],
  mc_lm_U = lml_mc_lm$V3[inds],
  iters = 4*iters[inds]
)

results_long <- tibble(
  lml = c(
    results$thames_full,
    results$bridge_full,
    results$mc_full,
    results$thames_reduced,
    results$bridge_reduced,
    results$mc_reduced,
    results$thames_lm,
    results$bridge_lm,
    results$mc_lm
  ),
  lml_lower = c(
    results$thames_full_L,
    results$bridge_full_L,
    results$mc_full_L,
    results$thames_reduced_L,
    results$bridge_reduced_L,
    results$mc_reduced_L,
    results$thames_lm_L,
    results$bridge_lm_L,
    results$mc_lm_L
  ),
  lml_upper = c(
    results$thames_full_U,
    results$bridge_full_U,
    results$mc_full_U,
    results$thames_reduced_U,
    results$bridge_reduced_U,
    results$mc_reduced_U,
    results$thames_lm_U,
    results$bridge_lm_U,
    results$mc_lm_U
  ),
  method = c(
    rep('THAMES',length(inds)),
    rep('Bridge',length(inds)),
    rep('MC',length(inds)),
    rep('THAMES',length(inds)),
    rep('Bridge',length(inds)),
    rep('MC',length(inds)),
    rep('THAMES',length(inds)),
    rep('Bridge',length(inds)),
    rep('MC',length(inds))
  ),
  model = c(
    rep('Full',3*length(inds)),
    rep('Reduced',3*length(inds)),
    rep('LM',3*length(inds))
  ),
  iters = rep(4*iters[inds],9)
)

library(ggpubr)
colors <- c("THAMES" = "black", "MC" = "red", "Bridge" = "blue")
ylim <- c(-8145,-8130)
full_plot <- results %>% 
  ggplot(aes(x=iters)) + 
  geom_point(aes(y=mc_full,color='MC'),size = 3) + 
  geom_errorbar(aes(ymax=min(ylim[1],mc_full_U), ymin = max(ylim[0],mc_full_L),color='MC'))  +
  geom_point(aes(y=bridge_full,color='Bridge'),size = 3) + 
  geom_errorbar(aes(ymax=bridge_full_U, ymin = bridge_full_L,color='Bridge'))  +
  geom_point(aes(y=thames_full,color='THAMES'),size = 3) + 
  geom_errorbar(aes(ymax=thames_full_U, ymin = thames_full_L,color='THAMES'))  +
  coord_cartesian(ylim=ylim) +
  labs(x="T", y = "Log marginal likelihood",
       title='Full LMM',color='Method') +
  scale_color_manual(values=colors)

reduced_plot <- results %>% 
  ggplot(aes(x=iters)) + 
  geom_point(aes(y=mc_reduced,color='MC'),size = 3) + 
  geom_errorbar(aes(ymax=mc_reduced_U, ymin = mc_reduced_L,color='MC'))  +
  geom_point(aes(y=bridge_reduced,color='Bridge'),size = 3) + 
  geom_errorbar(aes(ymax=bridge_reduced_U, ymin = bridge_reduced_L,color='Bridge'))  +
  geom_point(aes(y=thames_reduced,color='THAMES'),size = 3) + 
  geom_errorbar(aes(ymax=thames_reduced_U, ymin = thames_reduced_L,color='THAMES'))  +
  coord_cartesian(ylim=ylim) +
  labs(x="T", y = NULL,
       title='Reduced LMM',color='Method') +
  scale_color_manual(values=colors)

lm_plot <- results %>% 
  ggplot(aes(x=iters)) + 
  geom_point(aes(y=mc_lm,color='MC'),size = 3) + 
  geom_errorbar(aes(ymax=mc_lm_U, ymin = mc_lm_L,color='MC'))  +
  geom_point(aes(y=bridge_lm,color='Bridge'),size = 3) + 
  geom_errorbar(aes(ymax=bridge_lm_U, ymin = bridge_lm_L,color='Bridge'))  +
  geom_point(aes(y=thames_lm,color='THAMES'),size = 3) + 
  geom_errorbar(aes(ymax=thames_lm_U, ymin = thames_lm_L,color='THAMES'))  +
  coord_cartesian(ylim=c(-8284,-8276)) +
  labs(x="T", y = NULL,
       title='LM',color='Method') +
  scale_color_manual(values=colors)

ggarrange(full_plot,reduced_plot,lm_plot,nrow=1,legend ='bottom',common.legend = T)

ggsave('plots/thames-nlschools.pdf')
