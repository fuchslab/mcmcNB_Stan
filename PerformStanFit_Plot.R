library(rstan)
rstan_options(auto_write = TRUE)


N_anzahl <- 1000
fracPop <- c(0.2,0.8)
alpha <- c( 20,70)
prob <- c( 0.1,0.4)
beta <- prob / (1-prob)
N <- 2
Pop <- length(fracPop)

set.seed(234)
y1<-c()
y2<-c()
for(i in 1:(N_anzahl)){
  u <- runif(1,0,1)
  whichPop <- if(u<fracPop[1]) 1 else if (u<(sum(fracPop[1:2]))) 2 else 3
  y1[i]<-rnbinom(1, size=alpha[whichPop], prob=prob[whichPop]) 
  u <- runif(1,0,1)
  whichPop <- if(u<fracPop[1]) 1 else if (u<(sum(fracPop[1:2]))) 2 else 3
  y2[i]<-rnbinom(1, size=alpha[whichPop], prob=prob[whichPop]) 
}
y <- y1+y2


# Stan fit                  
iter <- 1000
chains <- 4                  
infi <- 1000                


fit_1000_1000_2_1000_2_2020_02_28_singlecell3_4 = stan(file="2020_02_28_singlecell3_third.stan", 
                                                       data=list(y = y, N_cell = 2, N_sample = N_anzahl, Pop = Pop, infi = infi),
                                                       init = list(list(theta = array(fracPop), alpha = array(alpha), beta = array(beta)),
                                                                   list(theta = array(fracPop), alpha = array(alpha), beta = array(beta)),
                                                                   list(theta = array(fracPop), alpha = array(alpha), beta = array(beta)),
                                                                   list(theta = array(fracPop), alpha = array(alpha), beta = array(beta))),
                                                       iter = iter, chains=chains, warmup = 300, seed = 345  ,save_dso = TRUE) 

saveRDS(fit_1000_1000_2_1000_2_2020_02_28_singlecell3_4, file= "fit_1000_1000_2_1000_2_2020_02_28_singlecell3.rds")



# Plot
library(colourlovers)
library(gridExtra)
library(ggplot2)
palette4 <- clpalette('694737')
fit_1000_1000_2_1000_2_2020_02_28_singlecell3_4 <- readRDS(file= "fit_1000_1000_2_1000_2_2020_02_28_singlecell3.rds")



rstan_ggtheme_options(legend.position = "none")
dens1 <- stan_dens(fit_1000_1000_2_1000_2_2020_02_28_singlecell3_4, pars = "alpha[1]", inc_warmup = FALSE, nrow = 3,fill =swatch(palette4)[[1]][2])+
  geom_segment(aes(x = alpha[1], y = 0, xend = alpha[1], yend = 0.045, linetype = "Truth"), size = 0.5)+
  geom_segment(aes(x = b["alpha[1]",1], y = 0, xend = b["alpha[1]",1], yend = 0.045, linetype = "Mean"), size = 0.5)+
  geom_segment(aes(x = b["alpha[1]","50%"], y = 0, xend = b["alpha[1]","50%"], yend = 0.045, linetype = "Median"), size = 0.5) +
  scale_linetype_manual("",values=c("Truth"=1,"Mean"=2, "Median" = 3))+ggtitle("Posterior densities")

rstan_ggtheme_options(legend.position = "top")
dens2 <-stan_dens(fit_1000_1000_2_1000_2_2020_02_28_singlecell3_4, pars = "alpha[2]", inc_warmup = FALSE, nrow = 3,fill =swatch(palette4)[[1]][2])+
  geom_segment(aes(x = alpha[2], y = 0, xend = alpha[2], yend = 0.03, linetype = "True Value"), size = 0.5)+
  geom_segment(aes(x = b["alpha[2]",1], y = 0, xend = b["alpha[2]",1], yend = 0.03, linetype = "Mean"), size = 0.5)+
  geom_segment(aes(x = b["alpha[2]","50%"], y = 0, xend = b["alpha[2]","50%"], yend = 0.03, linetype = "Median"), size = 0.5) + 
  scale_linetype_manual("",values=c("True Value"=1,"Mean"=2, "Median" = 3))

rstan_ggtheme_options(legend.position = "none")
dens3 <-stan_dens(fit_1000_1000_2_1000_2_2020_02_28_singlecell3_4, pars = "beta[1]", inc_warmup = FALSE, nrow = 3,fill =swatch(palette4)[[1]][2])+
  geom_segment(aes(x = beta[1], y = 0, xend = beta[1], yend = 10, linetype = "Truth"), size = 0.5)+
  geom_segment(aes(x = b["beta[1]",1], y = 0, xend = b["beta[1]",1], yend = 10, linetype = "Mean"), size = 0.5)+
  geom_segment(aes(x = b["beta[1]","50%"], y = 0, xend = b["beta[1]","50%"], yend = 10, linetype = "Median"), size = 0.5) +
  scale_linetype_manual("",values=c("Truth"=1,"Mean"=2, "Median" = 3))

dens4 <-stan_dens(fit_1000_1000_2_1000_2_2020_02_28_singlecell3_4, pars = "beta[2]", inc_warmup = FALSE, nrow = 3,fill =swatch(palette4)[[1]][2])+
  geom_segment(aes(x = beta[2], y = 0, xend = beta[2], yend = 4, linetype = "True Value"), size = 0.5)+
  geom_segment(aes(x = b["beta[2]",1], y = 0, xend = b["beta[2]",1], yend = 4, linetype = "Mean"), size = 0.5) +
  geom_segment(aes(x = b["beta[2]","50%"], y = 0, xend = b["beta[2]","50%"], yend = 4, linetype = "Median"), size = 0.5) +
  scale_linetype_manual("",values=c("True Value"=1,"Mean"=2, "Median" = 3))


dens5 <-stan_dens(fit_1000_1000_2_1000_2_2020_02_28_singlecell3_4, pars = "theta[1]", inc_warmup = FALSE, nrow = 3,fill =swatch(palette4)[[1]][2])+
  geom_segment(aes(x = theta[1], y = 0, xend = theta[1], yend = 10, linetype = "Truth"), size = 0.5)+
  geom_segment(aes(x = b["theta[1]",1], y = 0, xend = b["theta[1]",1], yend = 10, linetype = "Mean"), size = 0.5) +
  geom_segment(aes(x = b["theta[1]","50%"], y = 0, xend = b["theta[1]","50%"], yend = 10, linetype = "Median"), size = 0.5) +
  scale_linetype_manual("",values=c("Truth"=1,"Mean"=2, "Median" = 3))

dens6 <-stan_dens(fit_1000_1000_2_1000_2_2020_02_28_singlecell3_4, pars = "theta[2]", inc_warmup = FALSE, nrow = 3,fill =swatch(palette4)[[1]][2])+
  geom_segment(aes(x = theta[2], y = 0, xend = theta[2], yend = 10, linetype = "True Value"), size = 0.5)+
  geom_segment(aes(x = b["theta[2]",1], y = 0, xend = b["theta[2]",1], yend = 10, linetype = "Mean"), size = 0.5) +
  geom_segment(aes(x = b["theta[2]","50%"], y = 0, xend = b["theta[2]","50%"], yend = 10, linetype = "Median"), size = 0.5) +
  scale_linetype_manual("",values=c("True Value"=1,"Mean"=2, "Median" = 3))
Dens_Plot <- gridExtra::grid.arrange(dens1+ theme(plot.title = element_text(size = 13)) ,dens2,dens3,dens4,dens5,dens6, nrow =3) 


rstan_ggtheme_options(legend.position = "right")
Trace_Plot <- stan_trace(fit_1000_1000_2_1000_2_2020_02_28_singlecell3_4, inc_warmup = TRUE, nrow = 3)+scale_color_manual(values = swatch(palette4)[[1]][c(1,3:5)])+ggtitle("Traces of 4 Markov chains (300 warmup + 700 samling steps)")
Trace_Plot <- Trace_Plot + theme(plot.title = element_text(size = 13)) 


pdf(paste0("Density_Traces_fit_1000_1000_2_1000_2_2020_02_28_singlecell3_4.pdf"), width = 12, height = 7)
gridExtra::grid.arrange(Trace_Plot, Dens_Plot, ncol =2)
dev.off()