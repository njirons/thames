library(ggplot2)
library(ggbeeswarm)

my_theme <- theme_bw() +
  theme(strip.background = element_rect(fill = "white"), text = element_text(face="bold", size=12),
  )
theme_set(my_theme)

res = as.matrix(read.csv('data/multinomial.csv',header=FALSE,sep=''))

d_samples = c(1,20,50,100)
y_sims = 50
res_unif_bound = res[,1:length(d_samples)]
res_unif_nobound = res[,(length(d_samples)+1):(2*length(d_samples))]
res_random_bound = res[,(2*length(d_samples)+1):(3*length(d_samples))]
res_random_nobound = res[,(3*length(d_samples)+1):(4*length(d_samples))]

res_unif_nobound = cbind(t(t(c(res_unif_nobound))),rep(d_samples,rep(y_sims,length(d_samples))))
res_unif_nobound = as.data.frame(res_unif_nobound)
names(res_unif_nobound) = c("Error","d")
plot_2 = ggplot(data=res_unif_nobound,aes(x=factor(d),y=Error)) + geom_violin() +
  geom_beeswarm() + labs(y="Error",x="d") + geom_hline(yintercept = 0, linetype =2, color ="orange") +
  theme(text=element_text(size=18))

res_random_bound = cbind(t(t(c(res_random_bound))),rep(d_samples,rep(y_sims,length(d_samples))))
res_random_bound = as.data.frame(res_random_bound)
names(res_random_bound) = c("Error","d")
plot_3 = ggplot(data=res_random_bound,aes(x=factor(d),y=Error)) + geom_violin() +
  geom_beeswarm() + labs(y="Error",x="d") + geom_hline(yintercept = 0, linetype =2, color ="orange") +
  theme(text=element_text(size=18))

res_random_nobound = cbind(t(t(c(res_random_nobound))),rep(d_samples,rep(y_sims,length(d_samples))))
res_random_nobound = as.data.frame(res_random_nobound)
names(res_random_nobound) = c("Error","d")
plot_4 = ggplot(data=res_random_nobound,aes(x=factor(d),y=Error)) + geom_violin() +
  geom_beeswarm() + labs(y="Error",x="d") + geom_hline(yintercept = 0, linetype =2, color ="orange") +
  theme(text=element_text(size=18))

grid.arrange(plot_4, plot_3, plot_2, ncol=3)
ggsave("plots/thames-simple_multinomial_3_plots.pdf")
