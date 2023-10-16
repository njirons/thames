library(ggplot2)
library(ggbeeswarm)

Res_tot = as.matrix(read.csv('data/gauss_multi.csv',header=FALSE,sep=''))
colnames(Res_tot)=c("V1",   "V2",   "d",    "diff")

Res_tot = as.data.frame(Res_tot)
p<-Res_tot %>%
  mutate(d= factor(d, level = c(1,20,50,100))) %>%
  ggplot(aes(x =d, y =diff)) + geom_violin() +
  geom_beeswarm() + labs(y="Error") +
  geom_hline(yintercept = 0, linetype =2, color ="orange") +
  theme(text=element_text(size=32))
ggsave("plots/diff_marg_log_lik_half_sig_10000_vers02.pdf")
