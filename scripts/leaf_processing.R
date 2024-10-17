
################################################# Processing  ##########################################
########################################################################################################

traits <- read.csv("../data/traits.csv")
names(traits)

names(traits)
nrow(traits) # 2788

traits$Ntrt_fac <- as.factor(traits$Ntrt)
traits$Ptrt_fac <- as.factor(traits$Ptrt)
traits$Ktrt_fac <- as.factor(traits$Ktrt)
traits$block_fac <- as.factor(traits$block)
traits$trt_fac <- as.factor(traits$trt)

## add in photosynthetic pathway information
levels(as.factor(traits$Family)) # check Amaranthaceae, Asteraceae, Boraginaceae, Caryophyllaceae, Cyperaceae, Euphorbiaceae,
# Polygonaceae, Poaceae, Scrophulariaceae
# only one!
# C4
#traits$photosynthetic_pathway ='NULL'
traits$photosynthetic_pathway[traits$photosynthetic_pathway == 'NULL'
                              & traits$Family == 'Cyperaceae' & traits$Taxon == 'FIMBRISTYLIS DICHOTOMA'] <- 'C4'

traits$photosynthetic_pathway[traits$photosynthetic_pathway == 'NULL'] <- 'C3'

table(traits$photosynthetic_pathway)
## C4 : 521
## C3: 2267

### calculate SLA #############
traits$sla_m2_g = traits$SLA_v2 * (1/1000000)
hist(traits$sla_m2_g)

### calculate LMA #############
traits$lma = 1/traits$sla_m2_g
hist(traits$lma) # some extremely high values

### calculate narea #############
traits$narea = (traits$leaf_pct_N / 100) * (traits$lma)
hist(traits$narea) # some extremely high values

### calculate nmass #############
traits$nmass = traits$leaf_pct_N
hist(traits$nmass)
hist(traits$leaf_pct_N)

### calculate biomass for each species depending on the percentage max-cover for each species ###############
traits$spp_live_mass = traits$live_mass * (traits$max_cover / 100)
hist(traits$spp_live_mass)

## calculate big delta C13 from small delta 
traits$delta = ((-0.008 - traits$leaf_C13_delta_PDB * 0.001) / (1 + traits$leaf_C13_delta_PDB * 0.001)) * 1000
hist(traits$delta)

##### calculate chi for C3 plants #################### 
traits$chi[traits$photosynthetic_pathway == 'C3'] = 
  (traits$delta[traits$photosynthetic_pathway == 'C3'] * 0.001 - 0.0044) / (0.027 - 0.0044)
hist(traits$chi)

##### calculate chi for C4 plants #################### 
traits$chi[traits$photosynthetic_pathway == 'C4'] = 
  (traits$delta[traits$photosynthetic_pathway == 'C4'] * 0.001 - 0.0044) / ((-0.0057 + 0.03*0.4) - 0.0044)

hist(traits$chi)

#### create a new column for fence #################
traits$fence <- 'no'
traits$fence[traits$trt == 'Fence' | traits$trt == 'NPK+Fence'] <- 'yes'

################### transform some variables into log to get normal distribution #####################################

traits$loglma <- log(traits$lma)
hist(traits$loglma)

traits$lognarea <- log(traits$narea)
hist(traits$lognarea)

traits$lognmass <- log(traits$nmass)
hist(traits$lognmass)

names(traits)
traits$log_spp_live_mass <- log(traits$spp_live_mass)
hist(traits$log_spp_live_mass)

#################### exclude C points where chi is higher than 0.1 and lower than 0.95 ######### 
traits_chi_sub <- subset(traits, chi > 0.1 & chi < 0.95)
nrow(traits_chi_sub) ### 2106


hist(traits_chi_sub$chi) # quasi normal distribution 
hist(traits_chi_sub$aridity) 

nrow(traits_chi_sub) # 2106
min(traits_chi_sub$chi) # 0.114 
max(traits_chi_sub$chi) # 0.938

length(unique(traits_chi_sub$Taxon)) # 207 species 
length(unique(traits_chi_sub$site_code)) # 26 sites

################# exclusion of fence plots  ###########################
unique(traits_chi_sub$trt)
unique(traits_chi_sub$fence)

## exclusion of fence treatments 
traits_nofence_chi_sub = subset(traits_chi_sub, fence == "no")
traits_nofence_chi_sub$fence
nrow(traits_nofence_chi_sub) ## 1752
length(unique(traits_nofence_chi_sub$Taxon)) # 196 species 
length(unique(traits_nofence_chi_sub$site_code)) # 26 sites


view(traits_nofence_chi_sub)
write.csv(traits_nofence_chi_sub, "../data/traits_nofence_chi_sub.csv")

