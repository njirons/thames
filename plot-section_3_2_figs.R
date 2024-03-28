library(ggplot2)
library(mosaic)
library(tidyverse)
library(ggforce)
library(reshape2)
library(gridExtra)
library(latex2exp)
library(pracma)
library(TreeTools)  # Only need this for "LnDoubleFactorial()"
library(EnvStats)
library(mosaic)
library(reshape2)
library(latex2exp)
library(ggforce)
library(gridExtra)
source("functions/scv_figures.R")
source("functions/thames_figures_toyexample.R")

my_theme <- theme_bw() +
  theme(strip.background = element_rect(fill = "white"), text = element_text(face="bold", size=12),
  )
theme_set(my_theme)

##PLOT 1: L_d vs d
#find optimal L_d for fig_1
Q_L_plus_d = function(d,L) get_deriv_alt(d, sqrt(L+d), log_get_rec_explicit)
N=172
data_1 = numeric(N)
for (i in (1:N)){
  Q_L_plus_d_in_L = function(t) Q_L_plus_d(i,t)
  data_1[i] = findZeros(Q_L_plus_d_in_L(t)~t, xlim = c(0,2))$t
}

plot_1 = ggplot(data = data.frame(L=data_1, d=(1:length(data_1))), aes(x=d, y=L)) +
  geom_point() +
  geom_hline(yintercept=1,col="red")
plot_1 + theme(text=element_text(size=12))
ggsave("plots/section_3_2_thames_l_approx.pdf")

## PLOT 2 b: Approximate Size of the HPD region
data_2_b = pchisq(1:172+data_1,df=(1:172))
plot_2_b = ggplot(data = data.frame(HPD_region=data_2_b, d=(1:length(data_2_b)))) +
  aes(x=d,y=HPD_region) +
  geom_point() +
  geom_hline(yintercept = .5,col="red") +
  ylim(0,1) +
  labs(x="d",y="% HPD region")
plot_2_b + theme(text=element_text(size=12))
ggsave("plots/section_3_2_thames_hpd_region.pdf")

##PLOT 3: SCV(d,sqrt(d+L)) vs d for several values of L
data_3 = data.frame(d=1:100)
L_values = -4:8
for (i in L_values){
  scv_rec_explicit_test = Vectorize(function(d) scv_rec_explicit(d, sqrt(d+i)))
  data_3 = cbind(data_3, c(rep(NaN,(- min(i,0))),scv_rec_explicit_test((1 - min(i,0)):100)))
}
colnames(data_3) = c("d", L_values)
data_3.long <- melt(data_3, id = "d")
#data_3.long$variable=as.numeric(levels(data_3.long$variable))[data_3.long$variable]
data_3.long$al = ifelse(data_3.long$variable == 1,1,0.4)
plot_3 = ggplot(data_3.long, aes(d,value,colour=variable), label) +
  geom_line(alpha=data_3.long$al) +
  labs(x="d",y="SCV",color="L") +
  coord_trans(y="log") +
  scale_color_viridis_d(begin=.2,end=.8,option="magma") +
  stat_function(fun=function(x) 2*lower_bound(x),col="blue2",linetype = "dashed",alpha=1) +
  stat_function(fun=lower_bound,col="blue2",linetype = "dashed",alpha=1) +
  ylim(.5,17)
plot_3 + theme(text=element_text(size=12))
ggsave("plots/section_3_2_thames_scv_mix.pdf")

colnames(data_3)[7] = "one"
plot_3_a = ggplot(data_3, aes(x=d,y=one)) +
  geom_line(col="purple") +
  stat_function(fun=lower_bound,col="blue2",linetype="dashed")+
  coord_trans(x="sqrt") +
  labs(x="d",y="SCV")
plot_3_a + theme(text=element_text(size=12))
ggsave("plots/section_3_2_thames_sqrt_plot_thm1.pdf")

# PLOT 4 (toy example): THAMES_optimal vs THAMES_massive_c vs THAMES_tiny_c vs HARMONIC_MEAN
# The result in the final picture is chosen
# such that the peaks do not completely distort the picture,
# allthough the performance is the same each time.
set.seed(6)
d = 2
T = 10000
n = 100
est_mu = rep(0,d)
full_data = matrix(rep(0,d*n),ncol=d)
#Basic settings for a toy example

preset_random_numbers = matrix(rnorm(d*T),ncol=d)
# The estimators need to work on the same posterior sample

sigma_squ = 1
c_squ = d + 1
thames_opt = sim_mean_thames_open_simple(T, est_mu, n, c_squ, sigma_squ,
                                   type="inv_likeli",
                                   preset_random_numbers=preset_random_numbers,
                                   full_data=full_data)
values = c(50*sqrt(3), sqrt(3)*.01)#c(50,1/1000)
estimators = list()
for (i in (1:length(values))){
  c_squ = values[i]
  estimators[[i]] = sim_mean_thames_open_simple(T, est_mu, n, c_squ, sigma_squ,
                                           type="inv_likeli",
                                           full_data=full_data,
                                           preset_random_numbers)
}
estimators[[3]] = sim_mean_thames_open_simple(T, est_mu, n, c_squ, sigma_squ,
                                         type="harmonic_mean",
                                         full_data=full_data,
                                         preset_random_numbers)

sigma_n_squ = 1/(1/sigma_squ+n)

log_constant_01 = sum(apply(full_data, 1,function(x) sum((x-mean(x))^2)/2))
log_constant = log(2*pi/n)/2-n*log(2*pi)/2- log_constant_01

true_sol = exp((d/2)*log(2*pi*(sigma_squ+1/n))+
                 sum(est_mu^2)/(2*(sigma_squ+1/n))-
                 log_constant)
data_4 = thames_opt
data_4_sol = true_sol
data_4 = data.frame(matrix(c(100:T,
                             thames_opt[100:T],
                             estimators[[1]][100:T],
                             estimators[[2]][100:T],
                             estimators[[3]][100:T]),ncol=5))
#colnames(data_4) = c("t","optimal radius","small radius","large radius","Harmonic Mean Estimator")

plot_4 = ggplot(data = data_4, aes(x=X1)) +
  geom_line(aes(y=X2, col="X2")) +
  geom_line(aes(y=X3, col="X3"),alpha=.6) +
  geom_line(aes(y=X4, col="X4"),alpha=.6) +
  geom_line(aes(y=X5, col="X5"),alpha=.6) +
  geom_hline(alpha=.5,yintercept = data_4_sol,col="black", linetype="dashed") +
  scale_color_manual(name="Estimators",values = c(
    'X2' = 'blue',
    'X3' = 'green3',
    'X4' = "brown",
    'X5' = 'purple'),labels=c(TeX("\\textbf{THAMES with c=\\sqrt{3}}"),
                              TeX("\\textbf{THAMES with c=50\\sqrt{3}}"),
                              TeX("\\textbf{THAMES with c=0.1\\sqrt{3}}"),
                                  "Harmonic Mean Estimator")) +
  labs(y="Inverse marginal likelihood",x="T") + theme(text=element_text(size=18))
plot_4 = plot_4 + theme(legend.position=c(.62, .77))

# Plot A using ggplot2
# The data is just 0
test = data.frame(x=sqrt(sigma_n_squ) * preset_random_numbers[,1],
                  y=sqrt(sigma_n_squ) * preset_random_numbers[,2])
plot_4_b = ggplot(test, aes(x=x,y=y,value=1/(dnorm(x, sd=sqrt(sigma_n_squ))*dnorm(y,sd=sqrt(sigma_n_squ))))) +
  geom_point(aes(size=1/(dnorm(x/sqrt(sigma_n_squ))*dnorm(y/sqrt(sigma_n_squ)))))+
  geom_ellipse(aes(x0 = 0, y0 = 0, a = sqrt(3) * sqrt(sigma_n_squ)/2,
                   b =  sqrt(3) * sqrt(sigma_n_squ)/2, angle = 0),col="blue") +

  geom_ellipse(aes(x0 = 0, y0 = 0, a = 6*sqrt(3) * sqrt(sigma_n_squ)/2,
                   b =  6*sqrt(3) * sqrt(sigma_n_squ)/2, angle = 0),col="green3") +
  geom_ellipse(aes(x0 = 0, y0 = 0, a = .25*sqrt(3) * sqrt(sigma_n_squ)/2,
                   b =  .25*sqrt(3) * sqrt(sigma_n_squ)/2, angle = 0),col="brown") +
  labs(size="Inverse \nposterior \nlikelihood")

highest_index_1 = rev(order(test[,1]^2 + test[,2]^2))[1]
highest_index_2 = rev(order(test[,1]^2 + test[,2]^2))[2]
plot_4_b = plot_4_b + geom_text(data = test[highest_index_1,],
                     aes(label = sprintf("%s",highest_index_1)),
                     #nudge_x = 1,
                     nudge_y = -.05,
                     col="red") +
  geom_text(data = test[highest_index_2,],
            aes(label = sprintf("%s",highest_index_2)),
            #hjust = .7,
            nudge_y = -.05,
            col="red") + theme(text=element_text(size=18))
grid.arrange(plot_4, plot_4_b, ncol=2)
#grid.arrange(plot_4, plot_4_b, ncol=2)
ggsave('plots/section_3_2_thames_toyexample.jpeg', grid.arrange(plot_4, plot_4_b, ncol=2), device = "jpeg",dpi=700,height=5.6,width=11.14)

## PLOT 5 ratio of SCV of heuristic and optimal SCV v.s. d
data_5 = data.frame(d=1:100)
L_values = c(1)
for (i in L_values){
  scv_rec_explicit_test = Vectorize(function(d) scv_rec_explicit(d, sqrt(qchisq(.5,df=d))))
  data_5 = cbind(data_5, c(rep(NaN,(- min(i,0))),scv_rec_explicit_test((1 - min(i,0)):100)))
}

colnames(data_5)[2] = "one"
data_5["one"] = data_5["one"]/data_3["one"]
plot_5 = ggplot(data_5, aes(x=d,y=one)) +
  geom_point() +
  geom_hline(yintercept = 1,col="red")+
  labs(x="d",y="ratio of SCVs")
plot_5 + theme(text=element_text(size=12))
ggsave("plots/section_3_2_thames_ratio_opt_vs_heuristic.pdf")

## Other stuff ##

#rule of thumb for remark to theorem 1 statement 1
mean(abs(data_1[1:10] - 1) <= .5)
mean(abs(data_1[10:100] - 1) <= .05)
mean(abs(data_1[100:172] - 1) <= .005)

