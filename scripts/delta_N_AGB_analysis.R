#########################################
# packages necessary
#########################################

library(tidyr)
library(lme4)
library(car)
library(emmeans)
library(dplyr)
library(plyr)
library(RColorBrewer)
library(multcompView)
library(ggstatsplot)
library(tidyverse)
library(ggstatsplot)
library(ggplot2)
library(lme4)
library(car)
library(emmeans)
library(nlme)
library(marginaleffects)
library(piecewiseSEM)
library(rstantools)
library(multcomp)
library(treemapify)
library(relaimpo)
library(r2glmm)
library(patchwork)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(MuMIn)
library(boot)
library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(ggfortify)
library(patchwork)
library(semPlot)
library(piecewiseSEM)
library(emmeans)
library(car)
library(tidyverse)
library(MuMIn)
library(multcomp)
library(multcompView)
library(nlme)
library(cowplot)


####################### script to explore and analyse delta Nmass and Narea and delta above ground biomass between plots
#### that received N and plots that did not : N-control, NK-K, NP-P, NPK-PK ####################################

leaf_analysis <- read.csv("../data/traits_nofence_chi_sub.csv")
names(leaf_analysis)
length(unique(leaf_analysis$Taxon)) # 196
#################################### calculate Delta Nmass per Taxon  ######################################
##### if the taxa change between treatments, deltas are calculated only for the Taxon which is kept across treatments

### subset plots which did not receive N (they can be control, K, PK, P) ################## 
traits_lowN <- subset(leaf_analysis, Ntrt_fac == '0')
nrow(traits_lowN) # 894
table(traits_lowN$trt)


### subset plots which receive N (they can be N, NP, NPK, NK, NPK)
traits_highN <- subset(leaf_analysis, Ntrt_fac == '1')
nrow(traits_highN ) # 858

#################### combine low N and high N subsets##########################
Delta <- left_join(traits_lowN, traits_highN, 
                   by = c('site_code', 
                          'block_fac','Ptrt_fac', 'Ktrt_fac', 'Taxon'))
names(Delta)
nrow(Delta) ## 894

################ calculate delta: NPK-NP, N-control, NK-K, NP-P ################
Delta$delta_nmass <- ((Delta$nmass.y - Delta$nmass.x) / Delta$nmass.x) * 100
Delta$delta_narea <- ((Delta$narea.y - Delta$narea.x) / Delta$narea.x) * 100
Delta$delta_spp_live_mass <- ((Delta$spp_live_mass.y - Delta$spp_live_mass.x) / Delta$spp_live_mass.x) * 100

#### selection of high_N columns, but we will use only site characteristics and climate + taxon characteristics, don't use the other numeric columns concerning plots since they are for high N trt only ############
names(Delta)
delta_subset<- subset(Delta, select = c(1:65, 128:130))
colnames(delta_subset) <- sub("\\.x$", "", colnames(delta_subset)) ## remove .x  
names(delta_subset)
nrow(delta_subset) # 894

hist(delta_subset$delta_nmass) #ok
hist(delta_subset$delta_narea) 
hist(delta_subset$delta_spp_live_mass) 

#############  remove instances where any delta values are 3 times higher than the MAD Leys et al., 2020 ##########################################
#################################################################################################################################
nrow(delta_subset) # 894
delta_nmass_mad <- mad(delta_subset$delta_nmass, na.rm = T)
delta_nmass_mad   # 33.68
delta_nmass_median <- median(delta_subset$delta_nmass, na.rm = T)
delta_nmass_median ### 21.62
delta_nmass_mean<- mean(delta_subset$delta_nmass, na.rm = T)
delta_nmass_mean ## 27.92

deltanmass_mad_data <- subset(delta_subset,
                              delta_nmass < 3 * delta_nmass_mad &
                                delta_nmass > 3 * -delta_nmass_mad)

nrow(deltanmass_mad_data) # 524 
hist(deltanmass_mad_data$delta_nmass) # normal dist
names(deltanmass_mad_data)
length(table(deltanmass_mad_data$Taxon)) # 87

######################## delta nmass analysis #########################################################
#### Hyp 1: Cold temperature, drought and PAR increase deltanmass, deltanmass N-fixers < fixers, deltanmass C4 < C3 #######################
#########################################################################################################################
names(deltanmass_mad_data)
deltanmass_lmer_aridity <- lmer(delta_nmass ~tmp + aridity + par2_gs + Nfix +
                                  photosynthetic_pathway + Ptrt_fac*Ktrt_fac + 
                                  (1|Taxon) + (1|Taxon:site_code) + (1|Taxon:site_code:block_fac),data = deltanmass_mad_data)
Anova(deltanmass_lmer_aridity, type = 'II')
summary(deltanmass_lmer_aridity)


plot(resid(deltanmass_lmer_aridity) ~ fitted(deltanmass_lmer_aridity))

vif(deltanmass_lmer_aridity) # no colinearity
AIC(deltanmass_lmer_aridity) ## 4983.48

plot(deltanmass_lmer_aridity)
qqnorm(residuals(deltanmass_lmer_aridity))
qqline(residuals(deltanmass_lmer_aridity))

densityPlot(residuals(deltanmass_lmer_aridity))
shapiro.test(residuals(deltanmass_lmer_aridity)) 
outlierTest(deltanmass_lmer_aridity)

residuals <- resid(deltanmass_lmer_aridity)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

plot(fitted(deltanmass_lmer_aridity), residuals, xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs. Fitted Values")  # heteroscedasticity : none

r.squaredGLMM(deltanmass_lmer_aridity) ## mariginal:0.109 , conditional: 0.280, not terrible 

delta_nmass_model_aridity <- data.frame(Var = c('Tg', 'aridity', 'PAR', 'N fixer', 
                                        'C3/C4', 'Soil P', 'Soil K', 'Soil P x Soil K'))

delta_nmass_model_aridity$df <- as.matrix(Anova(deltanmass_lmer_aridity))[1:8, 2]
delta_nmass_model_aridity$Slope <- c(
  summary(emtrends(deltanmass_lmer_aridity, ~tmp, var = "tmp"))[1, 2],
  summary(emtrends(deltanmass_lmer_aridity, ~aridity, var = "aridity"))[1, 2],
  summary(emtrends(deltanmass_lmer_aridity, ~par2_gs, var = "par2_gs"))[1, 2],
  NA, NA, NA, NA, NA)
delta_nmass_model_aridity$SE <- c(
  summary(emtrends(deltanmass_lmer_aridity, ~tmp, var = "tmp"))[1, 3],
  summary(emtrends(deltanmass_lmer_aridity, ~aridity, var = "aridity"))[1, 3],
  summary(emtrends(deltanmass_lmer_aridity, ~par2_gs, var = "par2_gs"))[1, 3],
  NA, NA, NA, NA, NA)
delta_nmass_model_aridity$p <- as.matrix(Anova(deltanmass_lmer_aridity))[1:8, 3]
delta_nmass_model_aridity$RelImp <- as.matrix(calc.relip.mm(deltanmass_lmer_aridity)$lmg)[1:8]
delta_nmass_model_aridity$RelImp <- delta_nmass_model_aridity$RelImp * 100
delta_nmass_model_aridity

###################### # percentage increase of delta nmass in plots receiving P compared to plots not receiving P (NP-P)
(summary(emmeans(deltanmass_lmer_aridity, ~Ptrt_fac))[2,2] - summary(emmeans(deltanmass_lmer_aridity, ~Ptrt_fac))[1,2])/
  summary(emmeans(deltanmass_lmer_aridity, ~Ptrt_fac))[1,2]
# - 0.444 : plots receiving P respond less to N addition than plots not receiving P 
## this effect is significant in the model 

# percentage increase of Nmass in N fixers compared to non-N fixers
(summary(emmeans(deltanmass_lmer_aridity, ~Nfix))[2,2] - summary(emmeans(deltanmass_lmer_aridity, ~Nfix))[1,2])/
  summary(emmeans(deltanmass_lmer_aridity, ~Nfix))[1,2]
# - 0.913: non fixers increase in deltaNmass is 91.36% greater than  fixers 
## highly significant in the model

# percentage increase of Nmass in C4 compared to C3
(summary(emmeans(deltanmass_lmer_aridity, ~photosynthetic_pathway))[2,2] - summary(emmeans(deltanmass_lmer_aridity, ~photosynthetic_pathway))[1,2])/
  summary(emmeans(deltanmass_lmer_aridity, ~photosynthetic_pathway))[1,2]
#  0.5323: C4 plants increase in delta Nmass is 50.91% greater than in C3 
# ns in the model

### assign treatment group labels
deltanmass_mad_data$PKgroup[deltanmass_mad_data$Ptrt_fac == '0' & deltanmass_mad_data$Ktrt_fac == '0'] <- '-P, -K+µ'
deltanmass_mad_data$PKgroup[deltanmass_mad_data$Ptrt_fac == '1' & deltanmass_mad_data$Ktrt_fac == '0'] <- '+P, -K+µ'
deltanmass_mad_data$PKgroup[deltanmass_mad_data$Ptrt_fac == '0' & deltanmass_mad_data$Ktrt_fac == '1'] <- '-P, +K+µ'
deltanmass_mad_data$PKgroup[deltanmass_mad_data$Ptrt_fac == '1' & deltanmass_mad_data$Ktrt_fac == '1'] <- '+P, +K+µ'


test = cld(emmeans(deltanmass_lmer_aridity, ~ Ptrt_fac * Ktrt_fac)) ## comp slopes : plots receiving P have less DeltaNmass than plots which don't receive P 
test$.group
deltanmass_letters <- data.frame(x = c(1, 2, 3, 4),
                                 PKgroup = c('-P, -K+µ', '+P, -K+µ', '-P, +K+µ', '+P, +K+µ'), 
                                 y = c(120, 120, 120, 120), 
                                 group = c(cld(emmeans(deltanmass_lmer_aridity, ~Ptrt_fac * Ktrt_fac))[3, 8],
                                           cld(emmeans(deltanmass_lmer_aridity, ~Ptrt_fac * Ktrt_fac))[1, 8],
                                           cld(emmeans(deltanmass_lmer_aridity, ~Ptrt_fac * Ktrt_fac))[4, 8],
                                           cld(emmeans(deltanmass_lmer_aridity, ~Ptrt_fac * Ktrt_fac))[2, 8]))

deltanmass_letters$letter[deltanmass_letters$group == " 1 "] <- "a"
deltanmass_letters$letter[deltanmass_letters$group == " 12"] <- "ab"
deltanmass_letters$letter[deltanmass_letters$group == "  2"] <- "b"

names(deltanmass_mad_data)


(deltanmass_plot_PK <- ggplot(data = deltanmass_mad_data, 
                           aes(x = PKgroup, y = delta_nmass, colour = factor(PKgroup))) +
    scale_colour_manual(values = c("darkblue", "skyblue", 'salmon4', 'orange')) +
    theme(legend.position = "none",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          axis.title.y = element_text(size = 40, colour = 'black'),
          axis.title.x = element_text(size = 30, colour = 'black'),
          axis.text.x = element_text(size = 20, colour = 'black'),
          axis.text.y = element_text(size = 30, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "white")) +
    geom_boxplot(outlier.color = NA) +
    geom_point()+
    geom_text(data = deltanmass_letters, aes(x = x, y = y, label = letter), size = 8, colour="black") +
    scale_x_discrete(limits = c('-P, -K+µ', '+P, -K+µ', '-P, +K+µ', '+P, +K+µ'))+
    labs(fill = "Soil N") +
    ylab(expression('∆'* italic('N')['mass'] * ' (%)')) +
    xlab('P x K treatment'))  
deltanmass_plot_PK

#################### deltaNmass N2 fixers ######################################################
deltanmass_mad_data$Nfix <- factor(deltanmass_mad_data$Nfix, levels = c("no", "yes"), labels = c("Non-fixers", "Fixers"))

test = cld(emmeans(deltanmass_lmer_aridity, ~ Nfix)) ## comp slopes : plots receiving P have less DeltaNmass than plots which don't receive P 
test$.group
deltanmass_letters <- data.frame(x = c(1, 2),
                                 Nfix = c('Non-fixers', 'Fixers'), 
                                 y = c(120, 120), 
                                 group = c(cld(emmeans(deltanmass_lmer_aridity, ~Nfix))[2, 7],
                                           cld(emmeans(deltanmass_lmer_aridity, ~Nfix))[1, 7]))

deltanmass_letters$letter[deltanmass_letters$group == " 1 "] <- "a"
deltanmass_letters$letter[deltanmass_letters$group == "  2"] <- "b"


names(deltanmass_mad_data)
(deltanmass_plot_Nfix <- ggplot(data = deltanmass_mad_data, 
                           aes(x = Nfix, y = delta_nmass, colour = factor(Nfix), shape = factor(Nfix))) +
    scale_colour_manual(values = c("black", "black")) +  
    scale_shape_manual(values = c(1, 8))+
                         
    theme(legend.position = "none",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          axis.title.y = element_text(size = 40, colour = 'black'),
          axis.title.x = element_text(size = 30, colour = 'black'),
          axis.text.x = element_text(size = 30, colour = 'black'),
          axis.text.y = element_text(size = 30, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "white")) +
    geom_boxplot(outlier.color = NA) +
    geom_point()+
    geom_point(size = 2) +
    geom_text(data = deltanmass_letters, aes(x = x, y = y, label = letter), size = 8, colour="black") +
    ylab(expression('∆'* italic('N')['mass'] * ' (%)')) +
    xlab('Nitrogen fixation'))

  deltanmass_plot_Nfix
  
####################### linear regression figures #######################
######### aridity######################
aridity.regline <- data.frame(emmeans(deltanmass_lmer_aridity,"aridity",
                                      at = list(aridity =seq(min(deltanmass_mad_data$aridity, na.rm = T), max(deltanmass_mad_data$aridity, na.rm = T), 0.01)),
                                      type = "response")) 

(aridity_plot <- ggplot(data = deltanmass_mad_data, aes(x = aridity, y = delta_nmass)) + 
    
    geom_point(aes(color = factor(PKgroup), size = delta_nmass, shape=factor(Nfix)))+
    scale_size_continuous(range = c(3, 1.5)) + # Set the range for the size of points
    #scale_color_gradient(low = "darkblue", high = "lightblue") +
    scale_colour_manual(values = c("darkblue", "skyblue", 'salmon4', 'orange')) +
    scale_shape_manual(values = c(1, 8))  +
    
    geom_smooth(data = aridity.regline, aes(y = emmean), 
                col = 'black', lwd = 2, alpha = 0.8) +
    
    geom_ribbon(data = aridity.regline, aes(y=emmean, ymin = lower.CL, ymax = upper.CL), fill = "gray", alpha = 0.3)+ 
    scale_x_continuous(limits = c(0, 2.2)) +
    theme(legend.position = "none", 
          axis.title.y = element_text(size = 40, colour = 'black'),
          axis.title.x = element_text(size = 40, colour = 'black'),
          axis.text.x = element_text(size = 30, colour = 'black'),
          axis.text.y = element_text(size = 30, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "white")) +
    ylab(expression('∆'* italic('N')['mass'] * ' (%)')) +
    xlab("MI"))

######################################################################################################
################ temperature regression ############

tmp.regline <- data.frame(emmeans(deltanmass_lmer_aridity,"tmp",
                                      at = list(tmp =seq(min(deltanmass_mad_data$tmp, na.rm = T), max(deltanmass_mad_data$tmp, na.rm = T), 0.1)),
                                      type = "response")) 

(tmp_plot <- ggplot(data = deltanmass_mad_data, aes(x = tmp, y = delta_nmass)) + 
    
    geom_point(aes(color = factor(PKgroup), size = delta_nmass, shape= factor(Nfix)))+
    scale_size_continuous(range = c(3, 1.5)) + # Set the range for the size of points
   # scale_color_gradient(low = "darkblue", high = "lightblue") +
    scale_colour_manual(values = c("darkblue", "skyblue", 'salmon4', 'orange'))+
    scale_shape_manual(values = c(1, 8)) +
    
    geom_smooth(data = tmp.regline, aes(y = emmean), 
                col = 'black', lwd = 2, alpha = 0.8) +
    
    geom_ribbon(data = tmp.regline, aes(y=emmean, ymin = lower.CL, ymax = upper.CL), fill = "grey", alpha = 0.3)+ 
    
    scale_x_continuous(limits = c(3, 21)) +
    theme(legend.position = "none", 
          axis.title.y = element_text(size = 40, colour = 'black'),
          axis.title.x = element_text(size = 40, colour = 'black'),
          axis.text.x = element_text(size = 30, colour = 'black'),
          axis.text.y = element_text(size = 30, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "white")) +
    ylab(expression('∆'* italic('N')['mass'] * ' (%)')) +
    xlab(expression (italic('Tg') ['°C']))) 

######################################################################################################
################ PAR regression ############

PAR.regline <- data.frame(emmeans(deltanmass_lmer_aridity,"par2_gs",
                                  at = list(par2_gs =seq(min(deltanmass_mad_data$par2_gs, na.rm = T), max(deltanmass_mad_data$par2_gs, na.rm = T), 3)),
                                  type = "response")) 

(PAR_plot <- ggplot(data = deltanmass_mad_data, aes(x = par2_gs, y = delta_nmass)) + 
    
    geom_point(aes(color = factor(PKgroup), size = delta_nmass, shape= factor(Nfix)))+
    scale_size_continuous(range = c(3, 1.5)) + # Set the range for the size of points
   # scale_color_gradient(low = "lightblue", high = "darkblue") +
    
    scale_colour_manual(values = c("darkblue", "skyblue", 'salmon4', 'orange'))+
    scale_shape_manual(values = c(1, 8)) +
    
   # geom_smooth(data = PAR.regline, aes(y = emmean), 
    #            col = 'black', lwd = 2, alpha = 0.8) +
    
   # geom_ribbon(data = PAR.regline, aes(y=emmean, ymin = lower.CL, ymax = upper.CL), fill = "gray", alpha = 0.3)+ 
    
    #scale_x_continuous(limits = c(3, 21)) +
    theme(legend.position = "none", 
          axis.title.y = element_text(size = 40, colour = 'black'),
          axis.title.x = element_text(size = 40, colour = 'black'),
          axis.text.x = element_text(size = 30, colour = 'black'),
          axis.text.y = element_text(size = 30, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "white")) +
    ylab(expression('∆'* italic('N')['mass'] * ' (%)')) +
    xlab(expression (italic('PAR')))) 


############ merge regression plots 
rel_widths <- c(2, 2, 2)
merged_plot_deltanmass<- plot_grid(aridity_plot, tmp_plot, PAR_plot, ncol = 3, 
                                   align = "vh", labels = c("(a)", "(b)", "(c)"), label_size = 25,
                                   label_x = c(0.86, 0.86, 0.86), rel_widths = rel_widths)
merged_plot_deltanmass


final_plot_deltanmass <- plot_grid(
  merged_plot_deltanmass,
  plot_grid(deltanmass_plot_PK, deltanmass_plot_Nfix, ncol = 2, labels = c("(d)", "(e)"), label_size = 25, label_x = c(0.89, 0.89), rel_widths = c(0.5, 0.5)),
  ncol = 1, align = "vh", rel_heights = c(1, 1)
)

final_plot_deltanmass <- final_plot_deltanmass + theme(text = element_text(family = "Helvetica"))


######## Structural equation model for deltanmass and delta ABG corrected by the max cover ##################
#############################################################################################################

delta_subset$PKgroup[delta_subset$Ptrt_fac == '0' & delta_subset$Ktrt_fac == '0'] <- '-P, -K+µ'
delta_subset$PKgroup[delta_subset$Ptrt_fac == '1' & delta_subset$Ktrt_fac == '0'] <- '+P, -K+µ'
delta_subset$PKgroup[delta_subset$Ptrt_fac == '0' & delta_subset$Ktrt_fac == '1'] <- '-P, +K+µ'
delta_subset$PKgroup[delta_subset$Ptrt_fac == '1' & delta_subset$Ktrt_fac == '1'] <- '+P, +K+µ'

delta_subset$Nfix_level = 0 ## creation of a column with Nfix_level= 0 
delta_subset$Nfix_level[delta_subset$Nfix == 'yes'] = 1 
delta_subset$Nfix_level[delta_subset$Nfix == 'no'] = 0 

delta_subset$ps_level = 0 ## creation of a column with photosynthetic pathway = 0 
delta_subset$ps_level[delta_subset$photosynthetic_pathway == 'C3'] = 1 
delta_subset$ps_level[delta_subset$photosynthetic_pathway == 'C4'] = 0 

delta_nmass_mad <- mad(delta_subset$delta_nmass, na.rm = T)
delta_nmass_mad   # 33.68

delta_sem_Mad <- subset(delta_subset,
                       delta_nmass < 3 * delta_nmass_mad &
                         delta_nmass > 3 * -delta_nmass_mad)
                         
nrow(delta_sem_Mad) # 524
names(delta_sem_Mad)

subset_SEM  <- subset(delta_sem_Mad, select = c(site_code,Taxon, block, tmp, par2_gs, aridity, 
                                                delta_nmass,delta_spp_live_mass,delta_narea, Ptrt_fac,Ktrt_fac,Nfix_level,ps_level))
names(subset_SEM)

subset_SEM <- na.omit(subset_SEM)
nrow(subset_SEM) # 384


############### Structural equation model #######################

########################### using psem ####################### 
delta_sem <- psem (
  
  lm1 <-  lm(aridity~tmp , data=subset_SEM), 
  
  lme2  <-  lme(delta_spp_live_mass~ aridity + tmp + par2_gs + Nfix_level + ps_level +Ptrt_fac*Ktrt_fac, 
                 random = ~ 1|Taxon,
                 data = subset_SEM),
  
  lme3 <-  lme(delta_nmass~ delta_spp_live_mass + aridity + tmp + par2_gs + Nfix_level + ps_level +Ptrt_fac*Ktrt_fac,
                   random = ~ 1|Taxon,
                   data = subset_SEM)
 
  
)

summary(delta_sem)
plot(delta_sem)

model_summary <- summary(delta_sem)

r_squared_lm1 <- summary(lm1)$r.squared # 0.28
r_squared_lme2 <- r.squaredGLMM(lme2) # 0.12
r_squared_lme3 <- r.squaredGLMM(lme3) # 0.15

################################# delta narea analysis ######################################
#############################################################################################

#############  remove instances where any delta values are 3 times higher than the MAD ##########################################
#################################################################################################################################
nrow(delta_subset) # 894
delta_narea_mad <- mad(delta_subset$delta_narea, na.rm = T)
delta_narea_mad   # 38.71
delta_narea_median <- median(delta_subset$delta_narea, na.rm = T)
delta_narea_median ### 11.42
delta_narea_mean<- mean(delta_subset$delta_narea, na.rm = T)
delta_narea_mean ## 37.07

deltanarea_mad_data <- subset(delta_subset,
                              delta_narea < 3 * delta_narea_mad &
                                delta_narea > 3 * -delta_narea_mad)

nrow(deltanarea_mad_data) # 443 
hist(deltanarea_mad_data$delta_narea) # normal dist

#### Hyp 1: Cold temperature and drought increase deltanmass, N fixation decreases deltanmass #######################
#########################################################################################################################
######### When using aridity index, tmp and aridity have significant effects but not par ##############
deltanarea_lmer_aridity <- lmer(delta_narea ~tmp + aridity + par2_gs + Nfix + 
                                  photosynthetic_pathway + Ptrt_fac*Ktrt_fac + 
                                  (1|Taxon) + (1|Taxon:site_code) + (1|Taxon:site_code:block_fac),
                                data = deltanarea_mad_data)
Anova(deltanarea_lmer_aridity, type = 'II') ## tmp, aridity, p, nfix, photo pathway effects
summary(deltanarea_lmer_aridity)


plot(resid(deltanarea_lmer_aridity) ~ fitted(deltanarea_lmer_aridity))

vif(deltanarea_lmer_aridity) # no colinearity
AIC(deltanarea_lmer_aridity) ## 4464.286

plot(deltanmass_lmer_aridity)
qqnorm(residuals(deltanmass_lmer_aridity))
qqline(residuals(deltanmass_lmer_aridity))

densityPlot(residuals(deltanarea_lmer_aridity))
shapiro.test(residuals(deltanarea_lmer_aridity)) 
outlierTest(deltanarea_lmer_aridity)

residuals <- resid(deltanarea_lmer_aridity)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

plot(fitted(deltanarea_lmer_aridity), residuals, xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs. Fitted Values")  # heteroscedasticity : none

r.squaredGLMM(deltanarea_lmer_aridity) ## Conditional: 0.20, mariginal: 0.06

############ export stat #################
delta_narea_model_aridity <- data.frame(Var = c('Tg', 'aridity', 'PAR', 'N fixer', 
                                                'C3/C4', 'Soil P', 'Soil K', 'Soil P x Soil K'))

delta_narea_model_aridity$df <- as.matrix(Anova(deltanarea_lmer_aridity))[1:8, 2]
delta_narea_model_aridity$Slope <- c(
  summary(emtrends(deltanarea_lmer_aridity, ~tmp, var = "tmp"))[1, 2],
  summary(emtrends(deltanarea_lmer_aridity, ~aridity, var = "aridity"))[1, 2],
  summary(emtrends(deltanarea_lmer_aridity, ~par2_gs, var = "par2_gs"))[1, 2],
  NA, NA, NA, NA, NA)
delta_narea_model_aridity$SE <- c(
  summary(emtrends(deltanarea_lmer_aridity, ~tmp, var = "tmp"))[1, 3],
  summary(emtrends(deltanarea_lmer_aridity, ~aridity, var = "aridity"))[1, 3],
  summary(emtrends(deltanarea_lmer_aridity, ~par2_gs, var = "par2_gs"))[1, 3],
  NA, NA, NA, NA, NA)
delta_narea_model_aridity$p <- as.matrix(Anova(deltanarea_lmer_aridity))[1:8, 3]
delta_narea_model_aridity$RelImp <- as.matrix(calc.relip.mm(deltanarea_lmer_aridity)$lmg)[1:8]
delta_narea_model_aridity$RelImp <- delta_narea_model_aridity$RelImp * 100
delta_narea_model_aridity

###################### # percentage increase of delta narea in plots receiving P compared to plots not receiving P (NP-P)
(summary(emmeans(deltanarea_lmer_aridity, ~Ptrt_fac))[2,2] - summary(emmeans(deltanarea_lmer_aridity, ~Ptrt_fac))[1,2])/
  summary(emmeans(deltanarea_lmer_aridity, ~Ptrt_fac))[1,2]
# - 0.607 : plots receiving P respond less to N addition than plots not receiving P 
## this effect is significant in the model 

# percentage increase of Nmass in N fixers compared to non-N fixers
(summary(emmeans(deltanarea_lmer_aridity, ~Nfix))[2,2] - summary(emmeans(deltanarea_lmer_aridity, ~Nfix))[1,2])/
  summary(emmeans(deltanarea_lmer_aridity, ~Nfix))[1,2]
# - 1.124: non fixers increase in deltaNarea is 112.4 % greater than  fixers 
## highly significant in the model

# percentage increase of Nmass in C4 compared to C3
(summary(emmeans(deltanarea_lmer_aridity, ~photosynthetic_pathway))[2,2] - summary(emmeans(deltanarea_lmer_aridity, ~photosynthetic_pathway))[1,2])/
  summary(emmeans(deltanarea_lmer_aridity, ~photosynthetic_pathway))[1,2]
#  3.0006: C4 plants increase in delta Nmass is 300% greater than in C3 
# ns in the model

### assign treatment group labels
deltanarea_mad_data$PKgroup[deltanarea_mad_data$Ptrt_fac == '0' & deltanarea_mad_data$Ktrt_fac == '0'] <- '-P, -K+µ'
deltanarea_mad_data$PKgroup[deltanarea_mad_data$Ptrt_fac == '1' & deltanarea_mad_data$Ktrt_fac == '0'] <- '+P, -K+µ'
deltanarea_mad_data$PKgroup[deltanarea_mad_data$Ptrt_fac == '0' & deltanarea_mad_data$Ktrt_fac == '1'] <- '-P, +K+µ'
deltanarea_mad_data$PKgroup[deltanarea_mad_data$Ptrt_fac == '1' & deltanarea_mad_data$Ktrt_fac == '1'] <- '+P, +K+µ'


test = cld(emmeans(deltanarea_lmer_aridity, ~ Ptrt_fac * Ktrt_fac)) ## comp slopes : plots receiving P have less DeltaNmass than plots which don't receive P 
test$.group
deltanarea_letters <- data.frame(x = c(1, 2, 3, 4),
                                 PKgroup = c('-P, -K+µ', '+P, -K+µ', '-P, +K+µ', '+P, +K+µ'), 
                                 y = c(120, 120, 120, 120), 
                                 group = c(cld(emmeans(deltanarea_lmer_aridity, ~Ptrt_fac * Ktrt_fac))[3, 8],
                                           cld(emmeans(deltanarea_lmer_aridity, ~Ptrt_fac * Ktrt_fac))[1, 8],
                                           cld(emmeans(deltanarea_lmer_aridity, ~Ptrt_fac * Ktrt_fac))[4, 8],
                                           cld(emmeans(deltanarea_lmer_aridity, ~Ptrt_fac * Ktrt_fac))[2, 8]))

deltanarea_letters$letter[deltanarea_letters$group == " 1"] <- "a"
### no difference in the slopes 

(deltanarea_plot <- ggplot(data = deltanarea_mad_data, 
                           aes(x = PKgroup, y = delta_narea, colour = factor(PKgroup))) +
    scale_colour_manual(values = c("darkblue", "skyblue", 'salmon4', 'orange'))+
    theme(legend.position = "none",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          axis.title.y = element_text(size = 30, colour = 'black'),
          axis.title.x = element_text(size = 30, colour = 'black'),
          axis.text.x = element_text(size = 20, colour = 'black'),
          axis.text.y = element_text(size = 20, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "white")) +
    geom_boxplot(outlier.color = NA) +
    geom_text(data = deltanarea_letters, aes(x = x, y = y, label = letter), size = 6, colour="black") +
    #scale_fill_manual(values = c("gray40", "darkgreen"), labels = c("Ambient", "Added N")) +
    geom_point()+
    scale_x_discrete(limits = c('-P, -K+µ', '+P, -K+µ', '-P, +K+µ', '+P, +K+µ'))+
    labs(fill = "Soil N") +
    ylab(expression('∆'* italic('N')['area']))  +
    xlab('P x K treatment'))

#################### deltaNarea N2 fixers

deltanarea_mad_data$Nfix <- factor(deltanarea_mad_data$Nfix, levels = c("no", "yes"), labels = c("Non-fixers", "Fixers"))

test = cld(emmeans(deltanarea_lmer_aridity, ~ Nfix)) ## comp slopes : plots receiving P have less DeltaNmass than plots which don't receive P 
test$.group
deltanarea_letters <- data.frame(x = c(1, 2),
                                 Nfix = c('Non-fixers', 'Fixers'), 
                                 y = c(150, 150), 
                                 group = c(cld(emmeans(deltanarea_lmer_aridity, ~Nfix))[2, 7],
                                           cld(emmeans(deltanarea_lmer_aridity, ~Nfix))[1, 7]))

deltanarea_letters$letter[deltanarea_letters$group == " 1 "] <- "a"
deltanarea_letters$letter[deltanarea_letters$group == "  2"] <- "b"


(deltanarea_plot_Nfix <- ggplot(data = deltanarea_mad_data, 
                                aes(x = Nfix, y = delta_narea, colour = factor(Nfix), shape = factor(Nfix))) +
    scale_colour_manual(values = c("black", "black")) +  
    scale_shape_manual(values = c(1, 8))+
    
    theme(legend.position = "none",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          axis.title.y = element_text(size = 40, colour = 'black'),
          axis.title.x = element_text(size = 30, colour = 'black'),
          axis.text.x = element_text(size = 30, colour = 'black'),
          axis.text.y = element_text(size = 30, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "white")) +
    geom_boxplot(outlier.color = NA) +
    geom_point()+
    geom_point(size = 2) +
    geom_text(data = deltanarea_letters, aes(x = x, y = y, label = letter), size = 8, colour="black") +
    ylab(expression('∆'* italic('N')['area'] * ' (%)')) +
    xlab('Nitrogen fixation'))

deltanarea_plot_Nfix

####################### linear regression figures #######################
######### aridity
aridity.regline <- data.frame(emmeans(deltanarea_lmer_aridity,"aridity",
                                      at = list(aridity =seq(min(deltanarea_mad_data$aridity, na.rm = T), max(deltanarea_mad_data$aridity, na.rm = T), 0.01)),
                                      type = "response")) 

(aridity_plot_naraea <- ggplot(data = deltanarea_mad_data, aes(x = aridity, y = delta_narea)) + 
    
    geom_point(aes(color = factor(PKgroup), size = delta_narea, shape = factor(Nfix)))+
    scale_size_continuous(range = c(3, 1.5)) + 
    #scale_color_gradient(low = "darkblue", high = "lightblue") +
  
    scale_colour_manual(values = c("darkblue", "skyblue", 'salmon4', 'orange'))+
    scale_shape_manual(values = c(1,8))  +
    
    geom_smooth(data = aridity.regline, aes(y = emmean), 
                col = 'black', lwd = 2, alpha = 0.8) +
    
    geom_ribbon(data = aridity.regline, aes(y=emmean, ymin = lower.CL, ymax = upper.CL), fill = "gray", alpha = 0.3)+ 
    scale_x_continuous(limits = c(0, 2.2)) +
    theme(legend.position = "none", 
          axis.title.y = element_text(size = 40, colour = 'black'),
          axis.title.x = element_text(size = 40, colour = 'black'),
          axis.text.x = element_text(size = 30, colour = 'black'),
          axis.text.y = element_text(size = 30, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "white")) +
    ylab(expression('∆'* italic('N')['area'] * ' (%)')) +
    xlab("MI"))

aridity_plot_naraea

######################################################################################################
################ temperature regression ############

tmp.regline <- data.frame(emmeans(deltanarea_lmer_aridity,"tmp",
                                  at = list(tmp =seq(min(deltanarea_mad_data$tmp, na.rm = T), max(deltanarea_mad_data$tmp, na.rm = T), 0.1)),
                                  type = "response")) 

(tmp_plot_narea <- ggplot(data = deltanarea_mad_data, aes(x = tmp, y = delta_narea)) + 
    
    geom_point(aes(color = factor(PKgroup), size = delta_narea, shape=factor(Nfix)))+
    scale_size_continuous(range = c(3, 1.5)) + # Set the range for the size of points
   # scale_color_gradient(low = "darkblue", high = "lightblue") +
    scale_colour_manual(values = c("darkblue", "skyblue", 'salmon4', 'orange'))+
    scale_shape_manual(values = c(1,8))  +
    geom_smooth(data = tmp.regline, aes(y = emmean), 
                col = 'black', lwd = 2, alpha = 0.8) +
    
    geom_ribbon(data = tmp.regline, aes(y=emmean, ymin = lower.CL, ymax = upper.CL), fill = "gray", alpha = 0.3)+ 
    scale_x_continuous(limits = c(3, 21)) +
    theme(legend.position = "none", 
          axis.title.y = element_text(size = 40, colour = 'black'),
          axis.title.x = element_text(size = 40, colour = 'black'),
          axis.text.x = element_text(size = 30, colour = 'black'),
          axis.text.y = element_text(size = 30, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "white")) +
    ylab(expression('∆'* italic('N')['area'] * ' (%)')) +
    xlab(expression (italic('Tg') ['°C'])))

tmp_plot_narea

##########################################################################################################

################ PAR regression ############

PAR.regline <- data.frame(emmeans(deltanarea_lmer_aridity,"par2_gs",
                                  at = list(par2_gs =seq(min(deltanarea_mad_data$par2_gs, na.rm = T), max(deltanarea_mad_data$par2_gs, na.rm = T), 3)),
                                  type = "response")) 

(PAR_plot_narea <- ggplot(data = deltanarea_mad_data, aes(x = par2_gs, y = delta_narea)) + 
    
    geom_point(aes(color = factor(PKgroup), size = delta_narea, shape=factor(Nfix)))+
    scale_size_continuous(range = c(1.5, 3)) + # Set the range for the size of points
    # scale_color_gradient(low = "darkblue", high = "lightblue") +
    scale_colour_manual(values = c("darkblue", "skyblue", 'salmon4', 'orange'))+
    scale_shape_manual(values = c(1,8))  +
    geom_smooth(data = PAR.regline, aes(y = emmean), 
                col = 'black', lwd = 2, alpha = 0.8) +
    
    geom_ribbon(data = PAR.regline, aes(y=emmean, ymin = lower.CL, ymax = upper.CL), fill = "gray", alpha = 0.3)+ 
    #scale_x_continuous(limits = c(3, 21)) +
    theme(legend.position = "none", 
          axis.title.y = element_text(size = 40, colour = 'black'),
          axis.title.x = element_text(size = 40, colour = 'black'),
          axis.text.x = element_text(size = 30, colour = 'black'),
          axis.text.y = element_text(size = 30, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "white")) +
    ylab(expression('∆'* italic('N')['area'] * ' (%)')) +
    xlab(expression (italic('PAR'))))

PAR_plot_narea


################ merge regression plots 
rel_widths <- c(2, 2, 2)
merged_plot_deltanarea<- plot_grid(aridity_plot_naraea, tmp_plot_narea, PAR_plot_narea, ncol = 3, 
                                   align = "vh", labels = c("(a)", "(b)", "(c)"), label_size = 25, 
                                   label_x = c(0.86, 0.86, 0.88), rel_widths = rel_widths)
merged_plot_deltanarea


final_plot_delta_narea <- plot_grid(
  merged_plot_deltanarea,
  plot_grid(deltanarea_plot, deltanarea_plot_Nfix, ncol = 2, labels = c("(d)", "(e)"), label_size = 25, label_x = c(0.9, 0.9), rel_widths = c(0.5, 0.5)),
  ncol = 1, align = "vh", rel_heights = c(1, 1)
)

final_plot_delta_narea

final_plot_delta_narea <- final_plot_delta_narea + theme(text = element_text(family = "Helvetica"))



