#########################################
# packages necessary
#########################################
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(emmeans)
library(RColorBrewer)
library(multcompView)
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


################################ Statistics #########################################################################
### function to calculate relative importance for mixed models 
### from https://gist.github.com/BERENZ/e9b581a4b7160357934e
calc.relip.mm <- function(model,type = 'lmg') {
  if (!isLMM(model) & !isGLMM(model)) {
    stop('Currently supports only lmer/glmer objects', call. = FALSE)
  }
  require(lme4)
  X <- getME(model,'X')
  X <- X[ , -1]
  Y <- getME(model, 'y')
  s_resid <- sigma(model)
  s_effect <- getME(model, 'theta') * s_resid
  s2 <- sum(s_resid^2, s_effect^2)
  V <- Diagonal(x = s2, n = nrow(X))
  YX <- cbind(Y, X)
  cov_XY <- solve(t(YX) %*% solve(V) %*% as.matrix(YX))
  colnames(cov_XY) <- rownames(cov_XY) <- colnames(YX)
  importances <- calc.relimp(as.matrix(cov_XY), rela = F, type = type)
  return(importances)
}

calc.relip.boot.mm <- function(model,type = 'lmg') {
  if (!isLMM(model) & !isGLMM(model)) {
    stop('Currently supports only lmer/glmer objects', call. = FALSE)
  }
  require(lme4)
  X <- getME(model,'X')
  X <- X[ , -1]
  Y <- getME(model, 'y')
  s_resid <- sigma(model)
  s_effect <- getME(model, 'theta') * s_resid
  s2 <- sum(s_resid^2, s_effect^2)
  V <- Diagonal(x = s2, n = nrow(X))
  YX <- cbind(Y, X)
  cov_XY <- solve(t(YX) %*% solve(V) %*% as.matrix(YX))
  colnames(cov_XY) <- rownames(cov_XY) <- colnames(YX)
  bootresults <- boot.relimp(as.matrix(cov_XY), b=1000, rela = F, type = type)
  importances <- booteval.relimp(bootresults, norank=T)
  return(importances)
}

multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
}

############# load data #############################################################################

leaf_analysis <- read.csv("../data/traits_nofence_chi_sub.csv")
names(leaf_analysis)
length(unique(leaf_analysis$Taxon)) # 196

################### make table summarizing climate data for each site (table S1)
selected_data <- leaf_analysis %>% 
  select(site_code, Lat, Long, tmp, aridity, par2_gs, first_nutrient_year)

calculate_mean <- function(x) if (is.numeric(x)) mean(x, na.rm = TRUE) else first(x)
summary_table <- selected_data %>%
  group_by(site_code) %>%
  summarise_all(calculate_mean)
summary_table

######################## linear mixed effects model for Nmass #######################################
#####################################################################################################

#Hypothesis: Nmass increases with N and P addition, increases with aridity and cold tmp
########## Nmass is higher in N-fixing species and C3 species #############

names(leaf_analysis)
hist(leaf_analysis$lognmass) # normal distribution

leafNmass_lmer <- lmer(lognmass~ Ntrt_fac * Ptrt_fac * Ktrt_fac + tmp + 
                         par2_gs + loglma + chi + Nfix + photosynthetic_pathway +
                         (1|Taxon) + (1|Taxon:site_code) + (1|Taxon:site_code:block_fac), 
                       data = leaf_analysis)
plot(resid(leafNmass_lmer) ~ fitted(leafNmass_lmer))
summary(leafNmass_lmer)
Anova(leafNmass_lmer)
AIC(leafNmass_lmer) # 235.448
vif(leafNmass_lmer)

plot(leafNmass_lmer)
qqnorm(residuals(leafNmass_lmer))
qqline(residuals(leafNmass_lmer))

densityPlot(residuals(leafNmass_lmer))
shapiro.test(residuals(leafNmass_lmer)) 
outlierTest(leafNmass_lmer)

residuals <- resid(leafNmass_lmer)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## good
plot(fitted(leafNmass_lmer), residuals, xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs. Fitted Values")  # heteroscedasticity : none 

r.squaredGLMM(leafNmass_lmer) ## R2 mariginal: 0.46, R2 conditional : 0.83


######### Export Nmass model ####################################

Nmass_model <- data.frame(Var = c('Soil N', 'Soil P', 'Soil K+µ', 'Tg', 
                                  'PAR', 'ln LMA', 'χ', 'N fixer', 'C3/C4',
                                  'Soil N x Soil P', 'Soil N x Soil K', 'Soil P x Soil K',
                                  'Soil N x Soil P x Soil K'))
Nmass_model$df <- as.matrix(Anova(leafNmass_lmer))[1:13, 2]
Nmass_model$Slope <- c(NA, NA, NA,
                       summary(emtrends(leafNmass_lmer, ~tmp, var = "tmp"))[1, 2],
                       summary(emtrends(leafNmass_lmer, ~par2_gs, var = "par2_gs"))[1, 2],
                       summary(emtrends(leafNmass_lmer, ~loglma, var = "loglma"))[1, 2],
                       summary(emtrends(leafNmass_lmer, ~chi, var = "chi"))[1, 2],
                       NA, NA, NA, NA, NA, NA)
Nmass_model$SE <- c(NA, NA, NA,
                    summary(emtrends(leafNmass_lmer, ~tmp, var = "tmp"))[1, 3],
                    summary(emtrends(leafNmass_lmer, ~par2_gs, var = "par2_gs"))[1, 3],
                    summary(emtrends(leafNmass_lmer, ~loglma, var = "loglma"))[1, 3],
                    summary(emtrends(leafNmass_lmer, ~chi, var = "chi"))[1, 3],
                    NA, NA, NA, NA, NA, NA)
Nmass_model$p <- as.matrix(Anova(leafNmass_lmer))[1:13, 3]
Nmass_model$RelImp <- as.matrix(calc.relip.mm(leafNmass_lmer)$lmg)[1:13]
Nmass_model$RelImp <- Nmass_model$RelImp * 100
Nmass_model

#### soil nitrogen effect for the tables in supplementary information 
# percentage increase of Nmass in plots receiving N compared to plots not receiving N
(summary(emmeans(leafNmass_lmer, ~Ntrt_fac))[2,2] - summary(emmeans(leafNmass_lmer, ~Ntrt_fac))[1,2])/
  summary(emmeans(leafNmass_lmer, ~Ntrt_fac))[1,2]
# 0.201442

# percentage increase of Nmass in plots receiving N but not P compared to plots not receiving N or P
(summary(emmeans(leafNmass_lmer, ~Ntrt_fac * Ptrt_fac))[2,3] - summary(emmeans(leafNmass_lmer, ~Ntrt_fac * Ptrt_fac))[1,3])/
  summary(emmeans(leafNmass_lmer, ~Ntrt_fac * Ptrt_fac))[1,3]
# 0.2516

# percentage increase of Nmass in N fixers compared to non-N fixers
(summary(emmeans(leafNmass_lmer, ~Nfix))[2,2] - summary(emmeans(leafNmass_lmer, ~Nfix))[1,2])/
  summary(emmeans(leafNmass_lmer, ~Nfix))[1,2]
# 0.7029515

# percentage increase of Nmass in C3s compared to C4s
(summary(emmeans(leafNmass_lmer, ~photosynthetic_pathway))[1,2] - summary(emmeans(leafNmass_lmer, ~photosynthetic_pathway))[2,2])/
  summary(emmeans(leafNmass_lmer, ~photosynthetic_pathway))[2,2]
# 0.889215

leaf_analysis$PKgroup[leaf_analysis$Ptrt_fac == '0' & leaf_analysis$Ktrt_fac == '0'] <- '-P, -K+µ'
leaf_analysis$PKgroup[leaf_analysis$Ptrt_fac == '1' & leaf_analysis$Ktrt_fac == '0'] <- '+P, -K+µ'
leaf_analysis$PKgroup[leaf_analysis$Ptrt_fac == '0' & leaf_analysis$Ktrt_fac == '1'] <- '-P, +K+µ'
leaf_analysis$PKgroup[leaf_analysis$Ptrt_fac == '1' & leaf_analysis$Ktrt_fac == '1'] <- '+P, +K+µ'

cld_test <- cld(emmeans(leafNmass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac)) ## comp slopes : plots receiving N higher than those which did not
cld_test$.group

leafnmass_letters <- data.frame(x = c(0.8, 1.2, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2),
                                NPgroup = c('-P, -K+µ', '-P, -K+µ', '+P, -K+µ', '+P, -K+µ', 
                                            '-P, +K+µ', '-P, +K+µ', '+P, +K+µ', '+P, +K+µ'),
                                Ntrt_fac = c(0, 1, 0, 1, 0, 1, 0, 1),
                                y = c(2.7, 2.7, 2.7, 2.7, 2.7, 2.7, 2.8, 2.8), 
                                group = c(cld(emmeans(leafNmass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[1, 9],
                                          cld(emmeans(leafNmass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[7, 9],
                                          cld(emmeans(leafNmass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[3, 9],
                                          cld(emmeans(leafNmass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[6, 9],
                                          cld(emmeans(leafNmass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[2, 9],
                                          cld(emmeans(leafNmass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[8, 9],
                                          cld(emmeans(leafNmass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[4, 9],
                                          cld(emmeans(leafNmass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[5, 9]))
leafnmass_letters$Ntrt_fac <- as.factor(leafnmass_letters$Ntrt_fac)
leafnmass_letters$letter[leafnmass_letters$group == " 1 "] <- "a"
leafnmass_letters$letter[leafnmass_letters$group == "  2"] <- "b"

#### Regression table emmeans + SE 
se_table <- summary(emmeans(leafNmass_lmer, ~Ntrt_fac * Ktrt_fac * Ptrt_fac))
se_table
se_table = as.data.frame(se_table)

se_table$PKgroup[se_table$Ptrt_fac == '0' & se_table$Ktrt_fac == '0'] <- '-P, -K+µ'
se_table$PKgroup[se_table$Ptrt_fac == '1' & se_table$Ktrt_fac == '0'] <- '+P, -K+µ'
se_table$PKgroup[se_table$Ptrt_fac == '0' & se_table$Ktrt_fac == '1'] <- '-P, +K+µ'
se_table$PKgroup[se_table$Ptrt_fac == '1' & se_table$Ktrt_fac == '1'] <- '+P, +K+µ'
names(se_table)[which(names(se_table) == "emmean")] <- "lognmass"


################ Nmass violin plot + emmean and SE ################################################
leaf_analysis$Ntrt_fac <- as.factor(leaf_analysis$Ntrt_fac)
se_table$Ntrt_fac <- as.factor(se_table$Ntrt_fac)

nmass_violin <- ggplot(data = leaf_analysis, 
                      aes(x = PKgroup, y = lognmass, fill = Ntrt_fac),position = position_dodge(width = 0.5), size = 0.1) +
  geom_violin(trim=FALSE)+
  geom_point(data = se_table, aes(x = PKgroup, y = lognmass, fill = Ntrt_fac), 
            position = position_dodge(width = 0.9), size = 3) +

  geom_errorbar(data = se_table, aes(x = PKgroup, y = lognmass, ymin = lognmass - SE, ymax = lognmass + SE, 
                                     fill = Ntrt_fac), position = position_dodge(width = 0.9), width = 0.3, size = 0.8) +

  geom_text(data = leafnmass_letters, aes(x = x, y = y, label = letter), size = 6) +
  
  scale_fill_manual(values = c("gray60", "darkseagreen"), labels = c("Ambient", "Added N")) +
  
    theme(legend.position = "right",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          axis.title.y = element_text(size = 30, colour = 'black'),
          axis.title.x = element_text(size = 30, colour = 'black'),
          axis.text.x = element_text(size = 20, colour = 'black'),
          axis.text.y = element_text(size = 20, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "white")) +
    labs(fill = "Soil N") +
    ylab(expression('ln ' * italic('N')['mass'])) +
    xlab('P x K treatment') 

nmass_violin
nmass_violin <- nmass_violin + theme(text = element_text(family = "Helvetica"))

################################## Tree map N mass ###################################################

calc.relip.mm(leafNmass_lmer)$lmg

relimp_leafnmass <- NULL
relimp_leafnmass$Factor <- c('Soil N', 'Soil P', 'Soil K+µ', 'Tg', 'PAR', 'LMA', 'χ',
                             'N2 fixation', 'C3/C4', 'Soil Interactions', 'Unexplained')
relimp_leafnmass$Importance <- as.numeric(as.character(c(calc.relip.mm(leafNmass_lmer)$lmg[1:9], 
                                                         sum(calc.relip.mm(leafNmass_lmer)$lmg[10:13]),
                                                         1 - sum(calc.relip.mm(leafNmass_lmer)$lmg))))

relimp_leafnmass_df <- as.data.frame(relimp_leafnmass)

tm <- treemapify(data = relimp_leafnmass_df,
                 area = "Importance", start = "topleft")
tm$x <- (tm$xmax + tm$xmin) / 2
tm$y <- (tm$ymax + tm$ymin) / 2

nmass_tm <- full_join(relimp_leafnmass_df, tm, by = "Factor")
nmass_tm$name <-c('Soil~N', 'Soil~P', 'Soil~K[+µ]', 'italic(T)[g]', 'italic(PAR)',
                  'italic(LMA)[]', 'χ', 'N[2]~fixation', 
                  'C[3]/C[4]', 'Soil~Interactions', 'Unexplained')
nmass_tm$slope <- c(1,1,1,-1, 1, 
                    -1, -1, 1, 1, 1, 1)
nmass_tm$relationship = nmass_tm$slope*nmass_tm$Importance
#nmass_tm$slope <- factor(nmass_tm$slope, levels=c("No","negative", "positive"))


nmass_tm_test <- nmass_tm %>%
  mutate(
    Importance = as.factor(Importance),
    slope = as.factor(slope)
  )

(nmass_treemap <- ggplot(nmass_tm, 
                         aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                             label = name)) +
    geom_rect(aes(fill = Importance), color = "black") +
    theme(legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "none",
          panel.background = element_rect(fill = 'white'),
          axis.title = element_text(colour = 'white'),
          axis.text = element_text(colour = 'white'),
          axis.ticks = element_line(colour = "white")) + 
    
    #scale_fill_manual(values = c(positive = "darkblue", negative ="red4", nothing = "lightblue"))+
    scale_fill_gradient(low = "lightcyan2", high = "lightcyan4") +
    
    geom_text(data = filter(nmass_tm, Factor == 'PAR'), 
              aes(x = x, y = y), parse = T, size = 14, color="darkblue") +
    
    geom_text(data = filter(nmass_tm, Factor == 'LMA'), 
              aes(x = x, y = y), parse = T, size = 10, color="red3") +
    
    geom_text(data = filter(nmass_tm, Factor == 'χ'), 
              aes(x = x, y = y), parse = TRUE, size = 20, color = "red3", label=expression(chi)) +
   
    geom_text(data = filter(nmass_tm, Factor == 'N2 fixation'), 
              aes(x = x, y = y), parse = T, size = 9, color="darkblue") +
    
    
    geom_text(data = filter(nmass_tm, Factor == 'Unexplained'), 
              aes(x = x, y = y), parse = T, size = 6) +
    
    geom_text(data = filter(nmass_tm, Factor == 'Tg'), 
              aes(x = x, y = y), parse = TRUE, size = 10, color="red3") +
    
    geom_text(data = filter(nmass_tm,  Factor == 'C3/C4'), 
              aes(x = x, y = y), parse = T, size = 5, color="darkblue") +
    
    geom_text(data = filter(nmass_tm,  Factor == 'Soil N'), 
              aes(x = x, y = y), parse = T, size = 5, color="darkblue") +
    
    geom_text(data = filter(nmass_tm,  Factor == 'Soil Interactions'), 
              aes(x = x, y = y), parse = T, size = 4) +
    
    ggrepel::geom_text_repel(data = filter(nmass_tm, Factor == 'Soil P' |
                                             Factor == 'Soil K+µ'), 
                             aes(x = x, y = y), parse = T, size = 5, 
                             direction = "y", xlim = c(1.01, NA)) +
    
    scale_x_continuous(limits = c(0, 1.2), expand = c(0, 0), 
                       name = "X (\U03C7)") +
    #scale_x_continuous(limits = c(0, 1.2), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)))

nmass_treemap
nmass_treemap <- nmass_treemap + theme(text = element_text(family = "Helvetica"))


### linear mixed effects model for Narea ####################
###############################################################################
hist(leaf_analysis$lognarea) # normal distribution
leafnarea_lmer <- lmer(lognarea ~ Ntrt_fac * Ptrt_fac * Ktrt_fac + tmp + 
                         par2_gs + loglma + chi + Nfix + photosynthetic_pathway +
                         (1|Taxon) + (1|Taxon:site_code) + (1|Taxon:site_code:block_fac), 
                       data = leaf_analysis)
plot(resid(leafnarea_lmer) ~ fitted(leafnarea_lmer))
summary(leafnarea_lmer)
Anova(leafnarea_lmer)
AIC(leafnarea_lmer) ## 235.448
vif(leafnarea_lmer)

plot(leafnarea_lmer)
qqnorm(residuals(leafnarea_lmer))
qqline(residuals(leafnarea_lmer))

densityPlot(residuals(leafnarea_lmer))
shapiro.test(residuals(leafnarea_lmer))
outlierTest(leafnarea_lmer)

residuals <- resid(leafnarea_lmer)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## good
plot(fitted(leafnarea_lmer), residuals, xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs. Fitted Values")  # heteroscedasticity : high, when fitted values are high 

r.squaredGLMM(leafnarea_lmer) ## Conditional : 0.96, marginal: 0.87


########### Export model ###########################################
Narea_model <- data.frame(Var = c('Soil N', 'Soil P', 'Soil K+µ', 'Tg', 
                                  'PAR', 'ln LMA', 'χ', 'N fixer', 'C3/C4',
                                  'Soil N x Soil P', 'Soil N x Soil K', 'Soil P x Soil K',
                                  'Soil N x Soil P x Soil K'))
Narea_model$df <- as.matrix(Anova(leafnarea_lmer))[1:13, 2]
Narea_model$Slope <- c(NA, NA, NA,
                       summary(emtrends(leafnarea_lmer, ~tmp, var = "tmp"))[1, 2],
                       summary(emtrends(leafnarea_lmer, ~par2_gs, var = "par2_gs"))[1, 2],
                       summary(emtrends(leafnarea_lmer, ~loglma, var = "loglma"))[1, 2],
                       summary(emtrends(leafnarea_lmer, ~chi, var = "chi"))[1, 2],
                       NA, NA, NA, NA, NA, NA)
Narea_model$SE <- c(NA, NA, NA,
                    summary(emtrends(leafnarea_lmer, ~tmp, var = "tmp"))[1, 3],
                    summary(emtrends(leafnarea_lmer, ~par2_gs, var = "par2_gs"))[1, 3],
                    summary(emtrends(leafnarea_lmer, ~loglma, var = "loglma"))[1, 3],
                    summary(emtrends(leafnarea_lmer, ~chi, var = "chi"))[1, 3],
                    NA, NA, NA, NA, NA, NA)
Narea_model$p <- as.matrix(Anova(leafnarea_lmer))[1:13, 3]
Narea_model$RelImp <- as.matrix(calc.relip.mm(leafnarea_lmer)$lmg)[1:13]
Narea_model$RelImp <- Narea_model$RelImp * 100
Narea_model

################################# Figure comparison between treatments ######################## 

#### soil nitrogen effect
# percentage increase of Narea in plots receiving N compared to plots not receiving N
(summary(emmeans(leafnarea_lmer, ~Ntrt_fac))[2,2] - summary(emmeans(leafnarea_lmer, ~Ntrt_fac))[1,2])/
  summary(emmeans(leafnarea_lmer, ~Ntrt_fac))[1,2]
# 0.2244

# percentage increase of Narea in plots receiving N but not P compared to plots not receiving N or P
(summary(emmeans(leafnarea_lmer, ~Ntrt_fac * Ptrt_fac))[2,3] - summary(emmeans(leafnarea_lmer, ~Ntrt_fac * Ptrt_fac))[1,3])/
  summary(emmeans(leafnarea_lmer, ~Ntrt_fac * Ptrt_fac))[1,3]
# 0.2810871

# percentage increase of Narea in plots receiving N and P compared to plots receiving P but not N
(summary(emmeans(leafnarea_lmer, ~Ntrt_fac * Ptrt_fac))[4,3] - summary(emmeans(leafnarea_lmer, ~Ntrt_fac * Ptrt_fac))[3,3])/
  summary(emmeans(leafnarea_lmer, ~Ntrt_fac * Ptrt_fac))[3,3]
# 0.1706608

# percentage increase of Narea in N fixers compared to non-N fixers
(summary(emmeans(leafnarea_lmer, ~Nfix))[2,2] - summary(emmeans(leafnarea_lmer, ~Nfix))[1,2])/
  summary(emmeans(leafnarea_lmer, ~Nfix))[1,2]
# 0.8040299

# percentage increase of Narea in C3s compared to C4s
(summary(emmeans(leafnarea_lmer, ~photosynthetic_pathway))[1,2] - summary(emmeans(leafnarea_lmer, ~photosynthetic_pathway))[2,2])/
  summary(emmeans(leafnarea_lmer, ~photosynthetic_pathway))[2,2]
# 1.027255

leaf_analysis$PKgroup[leaf_analysis$Ptrt_fac == '0' & leaf_analysis$Ktrt_fac == '0'] <- '-P, -K+µ'
leaf_analysis$PKgroup[leaf_analysis$Ptrt_fac == '1' & leaf_analysis$Ktrt_fac == '0'] <- '+P, -K+µ'
leaf_analysis$PKgroup[leaf_analysis$Ptrt_fac == '0' & leaf_analysis$Ktrt_fac == '1'] <- '-P, +K+µ'
leaf_analysis$PKgroup[leaf_analysis$Ptrt_fac == '1' & leaf_analysis$Ktrt_fac == '1'] <- '+P, +K+µ'

test = cld(emmeans(leafnarea_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac)) ## comp slopes : plots receiving N higher than those which did not
test$.group

leafnarea_letters <- data.frame(x = c(0.8, 1.2, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2),
                                NPgroup = c('-P, -K+µ', '-P, -K+µ', '+P, -K+µ', '+P, -K+µ', 
                                            '-P, +K+µ', '-P, +K+µ', '+P, +K+µ', '+P, +K+µ'),
                                Ntrt_fac = c(0, 1, 0, 1, 0, 1, 0, 1),
                                y = c(6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5), 
                                group = c(cld(emmeans(leafnarea_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[1, 9],
                                          cld(emmeans(leafnarea_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[7, 9],
                                          cld(emmeans(leafnarea_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[3, 9],
                                          cld(emmeans(leafnarea_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[6, 9],
                                          cld(emmeans(leafnarea_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[2, 9],
                                          cld(emmeans(leafnarea_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[8, 9],
                                          cld(emmeans(leafnarea_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[4, 9],
                                          cld(emmeans(leafnarea_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[5, 9]))
leafnarea_letters$Ntrt_fac <- as.factor(leafnarea_letters$Ntrt_fac)
leafnarea_letters$letter[leafnarea_letters$group == " 1 "] <- "a"
leafnarea_letters$letter[leafnarea_letters$group == "  2"] <- "b"

#### Regression table emmeans + SE 
se_table <- summary(emmeans(leafnarea_lmer, ~Ntrt_fac * Ktrt_fac * Ptrt_fac))
se_table
se_table = as.data.frame(se_table)

se_table$PKgroup[se_table$Ptrt_fac == '0' & se_table$Ktrt_fac == '0'] <- '-P, -K+µ'
se_table$PKgroup[se_table$Ptrt_fac == '1' & se_table$Ktrt_fac == '0'] <- '+P, -K+µ'
se_table$PKgroup[se_table$Ptrt_fac == '0' & se_table$Ktrt_fac == '1'] <- '-P, +K+µ'
se_table$PKgroup[se_table$Ptrt_fac == '1' & se_table$Ktrt_fac == '1'] <- '+P, +K+µ'
names(se_table)[which(names(se_table) == "emmean")] <- "lognarea"


################ Narea violin plot + emmean and SE ################################################
leaf_analysis$Ntrt_fac <- as.factor(leaf_analysis$Ntrt_fac)
se_table$Ntrt_fac <- as.factor(se_table$Ntrt_fac)

narea_violin <- ggplot(data = leaf_analysis, 
                       aes(x = PKgroup, y = lognarea, fill = Ntrt_fac),position = position_dodge(width = 0.5), size = 0.1) +
  geom_violin(trim=FALSE, outlier.color = NA)+
  geom_point(data = se_table, aes(x = PKgroup, y = lognarea, fill = Ntrt_fac), 
             position = position_dodge(width = 0.9), size = 3) +
  
  geom_errorbar(data = se_table, aes(x = PKgroup, y = lognarea, ymin = lognarea - SE, ymax = lognarea + SE, 
                                     fill = Ntrt_fac), position = position_dodge(width = 0.9), width = 0.3, size = 0.8) +
  
  geom_text(data = leafnarea_letters, aes(x = x, y = y, label = letter), size = 6) +
  
  scale_fill_manual(values = c("gray60", "darkseagreen"), labels = c("Ambient", "Added N")) +
  
  theme(legend.position = "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.y = element_text(size = 30, colour = 'black'),
        axis.title.x = element_text(size = 30, colour = 'black'),
        axis.text.x = element_text(size = 20, colour = 'black'),
        axis.text.y = element_text(size = 20, colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white")) +
  labs(fill = "Soil N") +
  ylab(expression('ln ' * italic('N')['area'])) +
  xlab('P x K treatment') 


narea_violin
narea_violin <- narea_violin + theme(text = element_text(family = "Helvetica"))

################################## Tree map Narea ###################################################

calc.relip.mm(leafnarea_lmer)$lmg

relimp_leafnarea<- NULL
relimp_leafnarea$Factor <- c('Soil N', 'Soil P', 'Soil K+µ', 'Tg', 'PAR', 'LMA', 'χ',
                             'N2 fixation', 'C3/C4', 'Soil Interactions', 'Unexplained')
relimp_leafnarea$Importance <- as.numeric(as.character(c(calc.relip.mm(leafnarea_lmer)$lmg[1:9], 
                                                         sum(calc.relip.mm(leafnarea_lmer)$lmg[10:13]),
                                                         1 - sum(calc.relip.mm(leafnarea_lmer)$lmg))))

relimp_leafnarea_df <- as.data.frame(relimp_leafnarea)

tm <- treemapify(data = relimp_leafnarea_df,
                 area = "Importance", start = "topleft")
tm$x <- (tm$xmax + tm$xmin) / 2
tm$y <- (tm$ymax + tm$ymin) / 2

narea_tm <- full_join(relimp_leafnarea_df, tm, by = "Factor")
narea_tm$name <-c('Soil~N', 'Soil~P', 'Soil~K[+µ]', 'italic(T)[g]', 'italic(PAR)',
                  'italic(LMA)[]', 'italic(χ)', 'N[2]~fixation', 
                  'C[3]/C[4]', 'Soil~Interactions', 'Unexplained')
narea_tm$slope <- c(1,1,1,-1, 1, 
                    -1, -1, 1, 1, 1, 1)
narea_tm$relationship = narea_tm$slope*narea_tm$Importance
#nmass_tm$slope <- factor(nmass_tm$slope, levels=c("No","negative", "positive"))


narea_tm_test <- narea_tm %>%
  mutate(
    Importance = as.factor(Importance),
    slope = as.factor(slope)
  )

(narea_treemap <- ggplot(narea_tm, 
                         aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                             label = name)) +
    geom_rect(aes(fill = Importance), color = "black") +
    theme(legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "none",
          panel.background = element_rect(fill = 'white'),
          axis.title = element_text(colour = 'white'),
          axis.text = element_text(colour = 'white'),
          axis.ticks = element_line(colour = "white")) + 
    
    #scale_fill_manual(values = c(positive = "darkblue", negative ="red4", nothing = "lightblue"))+
    scale_fill_gradient(low = "lightcyan2", high = "lightcyan4") +
    
    geom_text(data = filter(narea_tm, Factor == 'LMA'), 
              aes(x = x, y = y), parse = T, size = 16, color="darkblue") +
    
    geom_text(data = filter(narea_tm, Factor == 'PAR'), 
              aes(x = x, y = y), parse = T, size = 8, color="darkblue") +
    
    geom_text(data = filter(narea_tm, Factor == 'Unexplained'), 
              aes(x = x, y = y), parse = T, size = 7) +
    
    geom_text(data = filter(narea_tm, Factor == 'χ'), 
              aes(x = x, y = y), parse = T, size = 14, color="red3", label=expression(chi)) +
    
    geom_text(data = filter(narea_tm, Factor == 'Tg'), 
              aes(x = x, y = y), parse = T, size = 8, color="red3") +
    
    geom_text(data = filter(narea_tm,  Factor == 'C3/C4'), 
              aes(x = x, y = y), parse = T, size = 8, color="darkblue") +
    
    geom_text(data = filter(narea_tm, Factor == 'N2 fixation'), 
              aes(x = x, y = y), parse = T, size = 4, color="darkblue") +
    
    geom_text(data = filter(narea_tm,  Factor == 'Soil N'), 
              aes(x = x, y = y), parse = T, size = 5, color="darkblue") +
    
    geom_text(data = filter(narea_tm,  Factor == 'Soil K+µ'), 
              aes(x = x, y = y), parse = T, size = 4) +
    
    ggrepel::geom_text_repel(data = filter(narea_tm, Factor == 'Soil P' | Factor == 'Soil Interactions'), 
                             aes(x = x, y = y), parse = T, size = 4, 
                             direction = "y", xlim = c(1.03, NA)) +
    scale_x_continuous(limits = c(0, 1.2), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)))

narea_treemap <- narea_treemap + theme(text = element_text(family = "Helvetica"))
