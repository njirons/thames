library(ggplot2)

res = as.matrix(read.csv('data/zellner.csv',header=FALSE,sep=''))

x<- c(50,100,200,500,seq(1000,10000,1000))
df = res
df = cbind(res,rep(x,8-1))
df=df[df[,dim(df)[2]]>100,]
df_1 = as.data.frame(df)
names(df_1) = c("nb_var","L","log_lik","U","bridge_L","bridge_log_lik","bridge_U",
                "naive_mc_L","naive_mc_log_lik","naive_mc_U", "real","nb_iter")
df_1_clone = df_1
names(df_1_clone) = c("nb_var","1","2","3","L","log_lik","U",
                      "naive_mc_L","naive_mc_log_lik","naive_mc_U", "real","nb_iter")
df_1_clone_2 = df_1
names(df_1_clone_2) = c("nb_var","1","2","3", "4", "5", "6",
                        "L","log_lik","U", "real","nb_iter")

colors <- c("THAMES" = "black", "MC" = "red", "Bridge" = "blue")

ggplot(data=df_1[df_1==8,],aes(x= nb_iter,y = log_lik, ymax = U, ymin = L))+
  geom_hline(aes(yintercept = real), color ="grey0", linetype=2) +
  geom_point(aes(color="THAMES")) + geom_point(data=df_1_clone[df_1==8,],aes(color="Bridge")) +
  geom_errorbar(aes(color="THAMES")) +
  geom_errorbar(data=df_1_clone[df_1==8,],aes(color="Bridge")) +
  facet_grid(nb_var~., scale ="free") +
  labs(y = "Marginal log likelihood", x= "Number of iterations") +
  theme(axis.text=element_text(size=11)) +
  theme(text=element_text(size=12))+labs(color="Method")+scale_color_manual(values=colors)
ggsave('plots/thames-zellner_thames_and_bridge.pdf')

ggplot(data=df_1,aes(x= nb_iter,y = log_lik, ymax = U, ymin = L))+
  geom_hline(aes(yintercept = real), color ="grey0", linetype=2) +
  geom_point(aes(color="THAMES")) + geom_point(data=df_1_clone,aes(color="Bridge")) + geom_point(data=df_1_clone_2,aes(color="MC"))+
  geom_errorbar(aes(color="THAMES")) +
  geom_errorbar(data=df_1_clone,aes(color="Bridge")) + geom_errorbar(data=df_1_clone_2,aes(color="MC"))+
  facet_grid(nb_var~., scale ="free") +
  labs(y = "Marginal log likelihood", x= "Number of iterations")+
  theme(axis.text=element_text(size=11)) +
  theme(text=element_text(size=12))+labs(color="Method")+scale_color_manual(values=colors)
ggsave('plots/thames-zellner_allests.pdf')
