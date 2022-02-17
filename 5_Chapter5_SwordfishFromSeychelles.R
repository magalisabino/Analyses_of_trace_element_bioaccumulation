##--------------------------------------------------------------------------------------------------------
## SCRIPT : Allows to reproduce all analyses of Chapter 5 dedicated to trace element bioaccumulation
##          in swordfish (Xiphias gladius) from the Seychelles.
##          This script includes the following analyses :
##            - Length-standardisation of TE concentrations for all TE for which concentrations are
##              correlated with LJFL
##            - Clustering according to FA profiles to identify FA trophic groups among all samples swordfish
##            - Correlation tests and plots to test for correlations between SI values and LJFL
##            - GAM models to identify contribution of physiological parameters, trophic ecology and season
##              in explaining TE bioaccumulation in swordfish
##            - GAM models with addition of interaction term "sex:d15N:length" to test for model improvement
##            - Test (+ associated plots) for significant difference in LJFL according to sex, FA trophic
##              group and season
##            - Test (+ associated plots) for significant difference in SI values according to sex,
##              FA trophic group and season
##            - Correlation tests (+ associated corrplot) for correlation among bioachemical tracers
##              (i.e. SI and FA)
##            - Calculation and plot of sex ratio by season
##
## As part of :
##        Magali SABINO PhD - "Bioaccumulation of trace elements in Seychelles marine food webs"
##
## Author : Magali Sabino
## First created : 2022-02-01
## Last update : 2022-02-17
##
##
####
##
## R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
## Copyright (C) 2021 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##--------------------------------------------------------------------------------------------------------


### 0 // Packages ##########################################################################################

## Open libraries
lapply(c("tidyverse", "sp", "sf", "rgdal", "ggspatial", "ggpubr",
         "FSA", "corrplot", "vegan", "ggdendro", "dendextend",
         "NbClust", "factoextra", "nicheROVER", "gam"),
       library, character.only = TRUE)

## Functions
pow.funk <- function(l, a, b, c, d) a*(l-c)^b-d
pow.funk.lin <- function(l, a, b) a*l+b



### I // Creation of database with only SWO samples from the Seychelles ##############################################################################################

SWO_data <- data_TE_SI %>% 
  filter(c_sp_fao == "SWO",
         sex %in% c("M","F"),
         !is.na(d13C),
         !is.na(length)) %>% 
  select(organism_identifier,sex,season,longitude,latitude,length,
         d13C,d15N,Cd_stat,Co_stat,Cu_stat,Fe_stat,Mn_stat,Pb_stat,
         Se_stat,TAs_stat,Zn_stat,THg_stat) %>% 
  rename(Cd = Cd_stat, Co = Co_stat, Cu = Cu_stat, Fe = Fe_stat,
         Mn = Mn_stat, Pb = Pb_stat, Se = Se_stat, As = TAs_stat,
         Zn = Zn_stat, Hg = THg_stat) %>% 
  left_join(SWO_FA_percent, by = "organism_identifier") %>% 
  select(-latitude_deg_dec, -longitude_deg_dec, -lowerjawfork_length)

SWO_data <- SWO_data %>% 
  mutate(As = (As*4.390106),
         Cd = (Cd*4.390106),
         Co = (Co*4.390106),
         Cu = (Cu*4.390106),
         Fe = (Fe*4.390106),
         Hg = (Hg*4.390106),
         Mn = (Mn*4.390106),
         Pb = (Pb*4.390106),
         Se = (Se*4.390106),
         Zn = (Zn*4.390106))



### II // Mapping sampling locations ##############################################################################################

# Projection of project
projection_UTM40S <- "+proj=utm +zone=40 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs" # Project in UTM 40S

# Projection of points
# Points: a classical data.frame 
# Determine coordinates (colums "lon", and "lat", i.e. X and Y)
Points <- SWO_data %>% 
  dplyr::select(longitude,latitude) %>% 
  rename(lon = longitude, lat = latitude) %>% 
  filter(!is.na(lon))
coordinates(Points) <- c("lon", "lat")

# Projection in WGS84
proj4string(Points) <- CRS("+proj=longlat +datum=WGS84")

# Apply your projection to the points
Points <- spTransform(Points, CRS(projection_UTM40S))

# Transform an "sp" object to an "sf"
Points_sf <- st_as_sf(Points)
st_is_valid(Points_sf)

# Import SHAPEFILES: EPSG:4326/WGS84
SEY_borders <- readOGR(dsn = "C:/Users/msabin02/Nextcloud/PhD/R_projects/Final_scripts/mapping_shp_files", layer = "gadm36_SYC_1") #dsn: your direction // layer: shapefile name WITHOUT extension
SEY_shelf <- readOGR(dsn = "C:/Users/msabin02/Nextcloud/PhD/R_projects/Final_scripts/mapping_shp_files", layer = "Sey_Shelf_200m_Poly_2March2018")
SEY_ZEE <- readOGR(dsn = "C:/Users/msabin02/Nextcloud/PhD/R_projects/Final_scripts/mapping_shp_files", layer = "eez_v10_recoupe")

# Reproject shp files
SEY_borders_UTM40S <- spTransform(SEY_borders, CRS(projection_UTM40S))
SEY_shelf_UTM40S <- spTransform(SEY_shelf, CRS(projection_UTM40S))
SEY_ZEE_UTM40S <- spTransform(SEY_ZEE, CRS(projection_UTM40S))

# Convert to sf format
SEY_borders_UTM40S <- st_as_sf(SEY_borders_UTM40S)
SEY_shelf_UTM40S <- st_as_sf(SEY_shelf_UTM40S)
SEY_ZEE_UTM40S <- st_as_sf(SEY_ZEE_UTM40S)

# Plot map
ggplot()+
  geom_sf(data = SEY_shelf_UTM40S, color = "#BBBDC0", fill = "#CAE3EF")+
  geom_sf(data = SEY_borders_UTM40S, fill = "#58585B", color = "#58585B")+
  geom_sf(data = Points_sf$geometry,
          aes(fill = Points_sf$FA_cluster),
          fill = "#FBAF3F",
          color = "black",
          shape = 21,
          size = 1.5) +
  coord_sf(crs = projection_UTM40S,
           xlim = c(min(st_coordinates(Points_sf)[,1])-50000, max(st_coordinates(Points_sf)[,1])+50000),
           ylim = c(min(st_coordinates(Points_sf)[,2])-100000, max(st_coordinates(Points_sf)[,2])+50000),
           expand = FALSE)+
  annotation_scale()+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank())

rm(projection_UTM40S,SEY_borders,SEY_borders_UTM40S,SEY_shelf,SEY_shelf_UTM40S,
   SEY_ZEE,SEY_ZEE_UTM40S,Points,Points_sf)



### III // Length-standardisation of TE concentrations ##############################################################################################

## 1 / Correlation tests between trace elements and LJFL ##############################################################################################

# Selection of data
data_tracers_x <- SWO_data %>% 
  filter(!is.na(length),
         !is.na(Cu),
         !is.na(Hg)) %>% 
  select(As,Cd,Co,Cu,Fe,Hg,Mn,Pb,Se,Zn) %>% 
  mutate(As = log10(As),
         Cd = log10(Cd),
         Co = log10(Co),
         Cu = log10(Cu),
         Fe = log10(Fe),
         Hg = log10(Hg),
         Mn = log10(Mn),
         Pb = log10(Pb),
         Se = log10(Se),
         Zn = log10(Zn))
data_tracers_x <- as.data.frame(data_tracers_x)

size <- as.data.frame(SWO_data %>% 
  filter(!is.na(length),
         !is.na(Cu),
         !is.na(Hg)) %>% 
  select(length))
size <- size$length

# Calculation of correlations
Output_ALL_cor <- NULL
for (i in 1: length(colnames(data_tracers_x))){ # metal loop
  #i = 1
  
  Metal_name <- colnames(data_tracers_x)[i]
  Metal_data <- data_tracers_x[,i]
  
  test_a <- shapiro.test(size)
  test_b <- shapiro.test(Metal_data)
  
  Output_cor <- data.frame(Metal = Metal_name,
                           N_val_used = sum(is.na(Metal_data ) == F),
                           Normality = NA,
                           Test_cor = NA,
                           P_val = NA,
                           Cor_val = NA)
  
  
  if (Output_cor$N_val_used <= 2){
    
    Output_ALL_cor <- rbind(Output_ALL_cor, Output_cor)
    
  }else{
    
    
    if(test_a$p.value >= 0.05 & test_b$p.value >= 0.05){# Pearson
      Output_cor$Normality <- T
      
      test_stat <- cor.test(size, Metal_data, method = "pearson")
      
      Output_cor$Test_cor <- "pearson"
      Output_cor$P_val <- test_stat$p.value
      Output_cor$Cor_val <- test_stat$estimate
      
    }else{
      # Kendall
      Output_cor$Normality <- F
      
      test_stat <- cor.test(size, Metal_data, method = "kendall")
      test_stat$p.value
      
      Output_cor$Test_cor <- "kendall"
      Output_cor$P_val <- test_stat$p.value
      Output_cor$Cor_val <- test_stat$estimate
      
    } # fin test correlation
    
    Output_ALL_cor <- rbind(Output_ALL_cor, Output_cor)
    rm(Output_cor, test_stat, test_a, test_b, Metal_data, Metal_name)
    
  } # fin boucle si N suffisant
}

rm(i)
rm(data_tracers_x,size)

# Dataframe creation for further plotting
corr_val <- Output_ALL_cor %>% 
  select(Metal,Cor_val) %>% 
  mutate(var = "Length") %>% 
  spread(Metal,Cor_val)
rownames(corr_val) <- corr_val$var
corr_val <- corr_val %>% 
  select(-var)
corr_val <- as.matrix(corr_val)

sig_val <- Output_ALL_cor %>% 
  select(Metal,P_val) %>% 
  mutate(var = "Length") %>% 
  spread(Metal,P_val)
rownames(sig_val) <- sig_val$var
sig_val <- sig_val %>% 
  select(-var)
sig_val <- as.matrix(sig_val)

# Correlation plot
corrplot(corr_val,
         method = "color",
         cl.pos = "n",
         addCoef.col = "black",
         tl.col="black", tl.srt = 0,
         p.mat = sig_val, sig.level = 0.05, insig = "blank",
         diag = T)

rm(Output_ALL_cor)
rm(corr_val,sig_val)


## 2 / Length-stand. for Arsenic ##############################################################################################

# Data
data_norm_As <- as.data.frame(SWO_data %>% 
  select(organism_identifier,length,As) %>% 
  filter(!is.na(length),
         !is.na(As)) %>% 
  mutate(logAs = log10(As)))

# Test linear regression
lin.reg <- lm(data_norm_As$logAs ~ data_norm_As$length)
par(mfrow = c(2,2))
plot(lin.reg)
reg.par <- summary(lin.reg)

a <- reg.par$coefficients[2,1]
b <- reg.par$coefficients[1,1]

summary(data_norm_As$length)
tbet_As <- data.frame(length=seq(80, 230, 0.2))
tbet_As[,2] <-  round(a * tbet_As[,1] + b, digits=5)

# Visual verification of adjustment curve
fit.plot.As <- data_norm_As %>% 
  mutate(x.par = max(length),
         y.par = max(logAs),
         fit.par = "Y= -0.003X + 0.96\nR2 = 0.31") %>% 
  ggplot() +
  geom_point(aes(x = length, y=logAs), size=2, alpha=0.3) +
  geom_line(data = tbet_As, aes(x = tbet_As[,1], y = tbet_As[,2]), size=0.8, col='#048C7F') +
  geom_text(aes(x = x.par, y = y.par, label = fit.par), hjust = 1, vjust =1, size = 3.5)+
  labs(x='Lower jaw-fork length (cm)', y='Log As') +
  theme_bw()+
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.99, .1),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size=9))

# Calculation of residuals and stand. concentrations
data_norm_As <- data_norm_As %>% 
  mutate(logAs_length = a*length+b,
         residuals = logAs - logAs_length,
         logAs_lmean = pow.funk.lin(150, a, b),
         logAs_stand = residuals + logAs_lmean,
         As_stand = 10^(logAs_stand))

R2 <- data_norm_As %>% 
  mutate(mean_logAs = mean(logAs),
         SSR_ind = (logAs_length - logAs)^2,
         SST_ind = (logAs - mean_logAs)^2) %>% 
  summarise(SSR = sum(SSR_ind),
            SST = sum(SST_ind),
            R2 = 1- (SSR/SST))

As_norm <- data_norm_As %>% 
  select(organism_identifier,As_stand)

rm(lin.reg,reg.par,a,b,R2)


## 3 / Length-stand. for Cadmium ##############################################################################################

# Data
data_norm_Cd <- as.data.frame(SWO_data %>% 
  filter(!is.na(length),
         !is.na(Cd)) %>% 
  mutate(logCd = log10(Cd)))

# Test linear regression
lin.reg <- lm(data_norm_Cd$logCd ~ data_norm_Cd$length)
par(mfrow = c(2,2))
plot(lin.reg)
reg.par <- summary(lin.reg)

a <- reg.par$coefficients[2,1]
b <- reg.par$coefficients[1,1]

summary(data_norm_Cd$length)
tbet_Cd <- data.frame(length=seq(80, 230, 0.2))
tbet_Cd[,2] <-  round(a * tbet_Cd[,1] + b, digits=5)

# Visual verification of adjustment curve
fit.plot.Cd <- data_norm_Cd %>% 
  mutate(x.par = min(length),
         y.par = max(logCd),
         fit.par = "Y= 0.005X - 1.49\nR2 = 0.29") %>% 
  ggplot() +
  geom_point(aes(x = length, y=logCd), size=2, alpha=0.3) +
  geom_line(data = tbet_Cd, aes(x = tbet_Cd[,1], y = tbet_Cd[,2]), size=0.8, col='#048C7F') +
  geom_text(aes(x = x.par, y = y.par, label = fit.par), hjust = 0, vjust =1, size = 3.5)+
  labs(x='Lower jaw-fork length (cm)', y='Log Cd') +
  theme_bw()+
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.99, .1),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size=9))

## Calcul des résidus et des concentrations standardisées
data_norm_Cd <- data_norm_Cd %>% 
  mutate(logCd_length = a*length+b,
         residuals = logCd - logCd_length,
         logCd_lmean = pow.funk.lin(150, a, b),
         logCd_stand = residuals + logCd_lmean,
         Cd_stand = 10^(logCd_stand))

R2 <- data_norm_Cd %>% 
  mutate(mean_logCd = mean(logCd),
         SSR_ind = (logCd_length - logCd)^2,
         SST_ind = (logCd - mean_logCd)^2) %>% 
  summarise(SSR = sum(SSR_ind),
            SST = sum(SST_ind),
            R2 = 1- (SSR/SST))

Cd_norm <- data_norm_Cd %>% 
  select(organism_identifier,Cd_stand)

rm(lin.reg,reg.par,a,b,R2)


## 4 / Length-stand. for Cobalt ##############################################################################################

# Data
data_norm_Co <- as.data.frame(SWO_data %>% 
   filter(!is.na(length),
         !is.na(Co)) %>% 
  mutate(logCo = log10(Co)))

# Test linear regression
lin.reg <- lm(data_norm_Co$logCo ~ data_norm_Co$length)
par(mfrow = c(2,2))
plot(lin.reg)
reg.par <- summary(lin.reg)

a <- reg.par$coefficients[2,1]
b <- reg.par$coefficients[1,1]

summary(data_norm_Co$length)
tbet_Co <- data.frame(length=seq(80, 230, 0.2))
tbet_Co[,2] <-  round(a * tbet_Co[,1] + b, digits=5)

# Visual verification of adjustment curve
fit.plot.Co <- data_norm_Co %>% 
  mutate(x.par = min(length),
         y.par = min(logCo),
         fit.par = "Y= 0.002X - 2.20\nR2 = 0.002") %>% 
  ggplot() +
  geom_point(aes(x = length, y=logCo), size=2, alpha=0.3) +
  geom_line(data = tbet_Co, aes(x = tbet_Co[,1], y = tbet_Co[,2]), size=0.8, col='#048C7F') +
  geom_text(aes(x = x.par, y = y.par, label = fit.par), hjust = 0, vjust =0, size = 3.5)+
  labs(x='Lower jaw-fork length (cm)', y='Log Co') +
  theme_bw()+
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.99, .1),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size=9))

# Calculation of residuals and stand. concentrations
data_norm_Co <- data_norm_Co %>% 
  mutate(logCo_length = a*length+b,
         residuals = logCo - logCo_length,
         logCo_lmean = pow.funk.lin(150, a, b),
         logCo_stand = residuals + logCo_lmean,
         Co_stand = 10^(logCo_stand))

R2 <- data_norm_Co %>% 
  mutate(mean_logCo = mean(logCo),
         SSR_ind = (logCo_length - logCo)^2,
         SST_ind = (logCo - mean_logCo)^2) %>% 
  summarise(SSR = sum(SSR_ind),
            SST = sum(SST_ind),
            R2 = 1-(SSR/SST))

Co_norm <- data_norm_Co %>% 
  select(organism_identifier,Co_stand)

rm(lin.reg,reg.par,a,b,R2)


## 5 / Length-stand. for Copper ##############################################################################################

# Data
data_norm_Cu <- as.data.frame(SWO_data %>% 
  filter(!is.na(length),
         !is.na(Cu)) %>% 
  mutate(logCu = log10(Cu)))

# Test linear regression
lin.reg <- lm(data_norm_Cu$logCu ~ data_norm_Cu$length)
par(mfrow = c(2,2))
plot(lin.reg)
reg.par <- summary(lin.reg)

a <- reg.par$coefficients[2,1]
b <- reg.par$coefficients[1,1]

summary(data_norm_Cu$length)
tbet_Cu <- data.frame(length=seq(80, 230, 0.2))
tbet_Cu[,2] <-  round(a * tbet_Cu[,1] + b, digits=5)

# Visual verification of adjustment
fit.plot.Cu <- data_norm_Cu %>% 
  mutate(x.par = min(length),
         y.par = max(logCu),
         fit.par = "Y= -0.003X + 0.69\nR2 = 0.11") %>% 
  ggplot() +
  geom_point(aes(x = length, y=logCu), size=2, alpha=0.3) +
  geom_line(data = tbet_Cu, aes(x = tbet_Cu[,1], y = tbet_Cu[,2]), size=0.8, col='#048C7F') +
  geom_text(aes(x = x.par, y = y.par, label = fit.par), hjust = 0, vjust =1, size = 3.5)+
  labs(x='Lower jaw-fork length (cm)', y='Log Cu') +
  theme_bw()+
  theme(axis.title.x = element_text(face = "italic"), 
        axis.title.y = element_text(face = "italic"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.99, .1),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size=9))

# Calculation of residuals and stand. concentrations
data_norm_Cu <- data_norm_Cu %>% 
  mutate(logCu_length = a*length+b,
         residuals = logCu - logCu_length,
         logCu_lmean = pow.funk.lin(150, a, b),
         logCu_stand = residuals + logCu_lmean,
         Cu_stand = 10^(logCu_stand))

R2 <- data_norm_Cu %>% 
  mutate(mean_logCu = mean(logCu),
         SSR_ind = (logCu_length - logCu)^2,
         SST_ind = (logCu - mean_logCu)^2) %>% 
  summarise(SSR = sum(SSR_ind),
            SST = sum(SST_ind),
            R2 = 1- (SSR/SST))

Cu_norm <- data_norm_Cu %>% 
  select(organism_identifier,Cu_stand)

rm(lin.reg,reg.par,a,b,R2)


## 6 / Length-stand. for Mercury ##############################################################################################

# Data 
data_norm_Hg <- as.data.frame(SWO_data %>% 
  filter(!is.na(length),
         !is.na(Hg)) %>% 
  mutate(logHg = log10(Hg)))

# Adjustment bioaccumulation curve logHg ~ FL
x <- data_norm_Hg$logHg
l <- data_norm_Hg$length

fit <- nls(x ~ pow.funk(l, a, b, c, d), algo='port', start=list(a=1, b=0.5, c=min(l), d=2), control=nls.control(maxiter=1000),
           lower=c(0.0001,0.001,0,-2),upper=c(1000,1,min(l),8))

a <- coef(fit)[1]
b <- coef(fit)[2]
c <- coef(fit)[3]
d <- coef(fit)[4]
summary(fit)

summary(data_norm_Hg$length)
tbet_Hg <- data.frame(length=seq(80, 230, 1))
tbet_Hg[,2] <-  round(a * (tbet_Hg[,1]-c) ^ b - d, digits=5)


# Visual verification of adjustment curve
fit.plot.Hg <- data_norm_Hg %>% 
  mutate(x.par = max(length),
         y.par = min(logHg),
         fit.par = "Y=0.29 (X - 74.85)0.30 - 0.62\nR2 = 0.35") %>% 
  ggplot() +
  geom_point(aes(x = length, y=logHg), size=2, alpha=0.3) +
  geom_line(data = tbet_Hg, aes(x = tbet_Hg[,1], y = tbet_Hg[,2]), size=0.8, col='#048C7F') +
  geom_text(aes(x = x.par, y = y.par, label = fit.par), hjust = 1, vjust =0, size = 3.5)+
  #xlim(80,235)+
  labs(x='Lower jaw-fork length (cm)', y='Log Hg') +
  theme_bw()+
  theme(axis.title.x = element_text(face = "italic"), # remove x-axis labels
        axis.title.y = element_text(face = "italic"), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.99, .1),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size=9))

# Calculation of residuals and stand. concentrations
data_norm_Hg <- data_norm_Hg %>% 
  mutate(logHg_length = a*(length-c)^b-d,
         residuals = logHg - logHg_length,
         logHg_lmean = pow.funk(150, a, b, c, d),
         logHg_stand = residuals + logHg_lmean,
         Hg_stand = 10^(logHg_stand))

R2 <- data_norm_Hg %>% 
  mutate(mean_logHg = mean(logHg),
         SSR_ind = (logHg_length - logHg)^2,
         SST_ind = (logHg - mean_logHg)^2) %>% 
  summarise(SSR = sum(SSR_ind),
            SST = sum(SST_ind),
            R2 = 1- (SSR/SST))

Hg_norm <- data_norm_Hg %>% 
  select(organism_identifier,Hg_stand)

rm(fit,a,b,c,d,l,x)


## 7 / Length-stand. for Manganese ##############################################################################################

# Data 
data_norm_Mn <- as.data.frame(SWO_data %>% 
  filter(!is.na(length),
         !is.na(Mn)) %>% 
  mutate(logMn = log10(Mn)))

# Test linear regression
lin.reg <- lm(data_norm_Mn$logMn ~ data_norm_Mn$length)
par(mfrow = c(2,2))
plot(lin.reg)
reg.par <- summary(lin.reg)

a <- reg.par$coefficients[2,1]
b <- reg.par$coefficients[1,1]

summary(data_norm_Mn$length)
tbet_Mn <- data.frame(length=seq(80, 230, 0.2))
tbet_Mn[,2] <-  round(a * tbet_Mn[,1] + b, digits=5)

# Visual verification of adjustmen curve
fit.plot.Mn <- data_norm_Mn %>% 
  mutate(x.par = min(length),
         y.par = max(logMn),
         fit.par = "Y= -0.002X - 0.55\nR2 = 0.07") %>% 
  ggplot() +
  geom_point(aes(x = length, y=logMn), size=2, alpha=0.3) +
  geom_line(data = tbet_Mn, aes(x = tbet_Mn[,1], y = tbet_Mn[,2]), size=0.8, col='#048C7F') +
  geom_text(aes(x = x.par, y = y.par, label = fit.par), hjust = 0, vjust =1, size = 3.5)+
  labs(x='Lower jaw-fork length (cm)', y='Log Mn') +
  theme_bw()+
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.99, .1),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size=9))

# Calculation of residuals and stand. concentrations
data_norm_Mn <- data_norm_Mn %>% 
  mutate(logMn_length = a*length+b,
         residuals = logMn - logMn_length,
         logMn_lmean = pow.funk.lin(150, a, b),
         logMn_stand = residuals + logMn_lmean,
         Mn_stand = 10^(logMn_stand))

R2 <- data_norm_Mn %>% 
  mutate(mean_logMn = mean(logMn),
         SSR_ind = (logMn_length - logMn)^2,
         SST_ind = (logMn - mean_logMn)^2) %>% 
  summarise(SSR = sum(SSR_ind),
            SST = sum(SST_ind),
            R2 = 1- (SSR/SST))

Mn_norm <- data_norm_Mn %>% 
  select(organism_identifier,Mn_stand)

rm(lin.reg,reg.par,a,b,R2)


## 8 / Length-stand. for Selenium ##############################################################################################

# Data
data_norm_Se <- as.data.frame(SWO_data %>% 
  filter(!is.na(length),
         !is.na(Se)) %>% 
  mutate(logSe = log10(Se)))

# Linear regression
lin.reg <- lm(data_norm_Se$logSe ~ data_norm_Se$length)
par(mfrow = c(2,2))
plot(lin.reg)
reg.par <- summary(lin.reg)

a <- reg.par$coefficients[2,1]
b <- reg.par$coefficients[1,1]

summary(data_norm_Se$length)
tbet_Se <- data.frame(length=seq(80, 230, 0.2))
tbet_Se[,2] <-  round(a * tbet_Se[,1] + b, digits=5)


# Visual verification of adjustment curve
fit.plot.Se <- data_norm_Se %>% 
  mutate(x.par = min(length),
         y.par = max(logSe),
         fit.par = "Y= -0.001X + 0.67\nR2 = 0.04") %>% 
  ggplot() +
  geom_point(aes(x = length, y=logSe), size=2, alpha=0.3) +
  geom_line(data = tbet_Se, aes(x = tbet_Se[,1], y = tbet_Se[,2]), size=0.8, col='#048C7F') +
  geom_text(aes(x = x.par, y = y.par, label = fit.par), hjust = 0, vjust =1, size = 3.5)+
  labs(x='Lower jaw-fork length (cm)', y='Log Se') +
  theme_bw()+
  theme(axis.title.x = element_text(face = "italic"), 
        axis.title.y = element_text(face = "italic"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.99, .1),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size=9))

# Calculation of residuals and stand. concentrations
data_norm_Se <- data_norm_Se %>% 
  mutate(logSe_length = a*length+b,
         residuals = logSe - logSe_length,
         logSe_lmean = pow.funk.lin(150, a, b),
         logSe_stand = residuals + logSe_lmean,
         Se_stand = 10^(logSe_stand))

R2 <- data_norm_Se %>% 
  mutate(mean_logSe = mean(logSe),
         SSR_ind = (logSe_length - logSe)^2,
         SST_ind = (logSe - mean_logSe)^2) %>% 
  summarise(SSR = sum(SSR_ind),
            SST = sum(SST_ind),
            R2 = 1- (SSR/SST))

Se_norm <- data_norm_Se %>% 
  select(organism_identifier,Se_stand)

rm(lin.reg,reg.par,a,b,R2)


## 9 / Length-stand. for Zinc ##############################################################################################

# Data
data_norm_Zn <- as.data.frame(SWO_data %>% 
  filter(!is.na(length),
         !is.na(Zn)) %>% 
  mutate(logZn = log10(Zn)))

# Linear regression
lm.fit <- lm(Zn ~ length, data = data_norm_Zn)
par(mfrow = c(2,2))
plot(lm.fit)

# Fit with NLS model
x <- data_norm_Zn$logZn
l <- data_norm_Zn$length

fit <- nls(x ~ pow.funk(l, a, b, c, d), algo='port', start=list(a=1, b=2, c=min(l), d=2), control=nls.control(maxiter=500),
           lower=c(0.0001,0,0,-10),upper=c(1000,5,min(l),20))

a <- coef(fit)[1]
b <- coef(fit)[2]
c <- coef(fit)[3]
d <- coef(fit)[4]
summary(fit)

summary(data_norm_Zn$length)
tbet_Zn <- data.frame(length=seq(80, 230, 0.2))
tbet_Zn[,2] <-  round(a * (tbet_Zn[,1]-c) ^ b - d, digits=5)

# Visual verification of adjustment curve
fit.plot.Zn <- data_norm_Zn %>% 
  mutate(x.par = min(length),
         y.par = max(logZn),
         fit.par = "Y= 0.0001 (X - 69.13)1.56 + 1.46\nR2 = 0.08") %>% 
  ggplot() +
  geom_point(aes(x = length, y=logZn), size=2, alpha=0.3) +
  geom_line(data = tbet_Zn, aes(x = tbet_Zn[,1], y = tbet_Zn[,2]), size=0.8, col='#048C7F') +
  geom_text(aes(x = x.par, y = y.par, label = fit.par), hjust = 0, vjust =1, size = 3.5)+
  labs(x='Lower jaw-fork length (cm)', y='Log Zn') +
  theme_bw()+
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.99, .1),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size=9))

# Calculation of residuals and stand. concentrations
data_norm_Zn <- data_norm_Zn %>% 
  mutate(logZn_length = a*(length-c)^b-d,
         residuals = logZn - logZn_length,
         logZn_lmean = pow.funk(150, a, b, c, d),
         logZn_stand = residuals + logZn_lmean,
         Zn_stand = 10^(logZn_stand))

R2 <- data_norm_Zn %>% 
  mutate(mean_logZn = mean(logZn),
         SSR_ind = (logZn_length - logZn)^2,
         SST_ind = (logZn - mean_logZn)^2) %>% 
  summarise(SSR = sum(SSR_ind),
            SST = sum(SST_ind),
            R2 = 1- (SSR/SST))

Zn_norm <- data_norm_Zn %>% 
  select(organism_identifier,Zn_stand)

rm(lm.fit,fit,a,b,c,d,l,x,R2)


## 10 / Plot all fits ##############################################################################################

ggarrange(fit.plot.As, fit.plot.Cd, fit.plot.Co, fit.plot.Cu,
          fit.plot.Hg,fit.plot.Mn,fit.plot.Se,fit.plot.Zn,
          labels = c("A.","B.","C.","D.","E.","F.","G.","H."),
          ncol = 2, nrow = 4, align = "hv")

rm(fit.plot.As, fit.plot.Cd, fit.plot.Co, fit.plot.Cu,
   fit.plot.Hg,fit.plot.Mn,fit.plot.Se,fit.plot.Zn)
rm(data_norm_As,data_norm_Cd,data_norm_Co,data_norm_Cu,
   data_norm_Hg,data_norm_Mn,data_norm_Se,data_norm_Zn)
rm(tbet_As,tbet_Cd,tbet_Co,tbet_Cu,tbet_Hg,tbet_Mn,tbet_Se,tbet_Zn)


## 11 / Database update with stand. concentrations ##############################################################################################

SWO_data <- SWO_data %>% 
  left_join(As_norm, by = "organism_identifier") %>% 
  left_join(Cd_norm, by = "organism_identifier") %>% 
  left_join(Co_norm, by = "organism_identifier") %>% 
  left_join(Cu_norm, by = "organism_identifier") %>% 
  left_join(Hg_norm, by = "organism_identifier") %>% 
  left_join(Mn_norm, by = "organism_identifier") %>% 
  left_join(Se_norm, by = "organism_identifier") %>% 
  left_join(Zn_norm, by = "organism_identifier")

rm(As_norm,Cd_norm,Co_norm,Cu_norm,Hg_norm,Mn_norm,Se_norm,Zn_norm)



### IV // Fatty acid trophic groups - Clustering method ##############################################################################################

## 1 / Clustering ##############################################################################################

# FA trophic marker selection
data_cluster <- SWO_data %>%
  filter(!is.na(c22_6w3_c),
         !is.na(Cu),
         !is.na(Hg),
         !length < 100) %>% 
  select(organism_identifier, season, sex, length,
         c14_c,c16_1w7_c,c16_c,c17_c,c18_1w7_c,c18_1w9_c,c18_c,
         c20_1w9_c,c20_4w6_c,c20_5w3_c,c22_5w3_c,c22_5w6_c,
         c22_6w3_c,c24_1w9_c)

# Scaling data to reduce effect of value gaps between FAs
names(data_cluster)
clustering <- as.matrix(data_cluster[,5:18])
rownames(clustering) <- data_cluster$organism_identifier
clustering <- scale(clustering)

# Need of data translation to avoid negative data (i.e. Bray Curtis matrix)
# Linear_positiv_translation <- function(matrix){
#   result <- matrix + abs(min(matrix, na.rm = TRUE)) + 0.1
# }
clustering <- Linear_positiv_translation(clustering)
clustering_mtx <- clustering
clustering <- as.data.frame(clustering)

# Clustering
cluster_Bray <- vegdist(clustering, method = "euclidean")
hc <- hclust(cluster_Bray, method = "ward.D2")

# Mise en forme data for plotting
graph.cluster <- as.dendrogram(hc) # Build dendogram object from hclust results
graph.cluster <- dendro_data(graph.cluster, type = "rectangle") # Extract the data for rectangle lines
clusters <- data_cluster %>% 
  select(organism_identifier, season, sex, length) %>% 
  rename(label = organism_identifier) # Extract infos on species
graph.cluster[["labels"]] <- merge(graph.cluster[["labels"]], clusters, by = "label") # Merge infos with labels of dendogram

# Number of clusters
res.nbclust <-
  NbClust(clustering_mtx, distance = "euclidean",
          min.nc = 2, max.nc = 10,
          method = "ward.D2", index = "all")

fviz_nbclust(res.nbclust, ggtheme = theme_minimal())

res.hc <- clustering_mtx %>%
  eclust("hclust", k = 2, graph = FALSE)

# Determination of clusters
clust_FA <- as.data.frame(cutree(hc, k = 2))
clust_FA$label <- row.names(clust_FA)
colnames(clust_FA)[1] <- c("cluster")
clust_FA$cluster <- as.character(clust_FA$cluster)
graph.cluster[["labels"]] <- merge(graph.cluster[["labels"]], clust_FA, by = "label") # Merge infos with labels of dendogram

# Plotting
dendro <- ggplot()+
  geom_segment(data = segment(graph.cluster), aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = label(graph.cluster), aes(x = x, y = y, label = label, hjust = 1,
                                             color = cluster,
                                             angle = 90), size = 3)+
  scale_color_manual(values = c("#48A4E3","#DE8E06"))+
  labs(y = "Weight", color = "Cluster")+
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "italic"),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        panel.grid = element_blank())

rm(data_cluster, clustering, hc, res.hc, cluster_Bray, clusters,
   clustering_mtx, res.nbclust)

clust_FA <- clust_FA %>% 
  mutate(cluster = as.character(cluster),
         new_cluster = "1",
         new_cluster = ifelse(cluster == "1", "2", new_cluster))


## 2 / FA niches representation ##############################################################################################

# Data selection
data_niche <- SWO_data %>% 
  rename(label = organism_identifier) %>% 
  left_join(clust_FA, by = "label") %>% 
  mutate(label = as.factor(label)) %>% 
  filter(!is.na(cluster)) %>% 
  select(label,new_cluster,c14_c,c16_1w7_c,c16_c,c17_c,
         c18_1w7_c,c18_1w9_c,c18_c,c20_1w9_c,c20_4w6_c,c20_5w3_c,
         c22_5w3_c,c22_5w6_c,c22_6w3_c,c24_1w9_c)
rownames(data_niche) <- data_niche$label

# nMDS computing
MDS_2017 <- metaMDS(data_niche[,3:16],distance = "bray",k = 2,try = 300)

# Get coordinates of individuals on the nMDS
data_niche_FA <- as.data.frame(scores(MDS_2017))
data_niche_FA$cluster <- data_niche$new_cluster
data_niche_FA$label <- data_niche$label
rm(data_niche)

# Niche ellipses computing and niche size calculation
nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_FA), data_niche_FA$cluster,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_FA[ii,1:2]))


# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})

# Extraction of points coordinates of all calculated ellipses
niche.par = fish.par
pfrac = 0.05
alpha=0.95

niso <- ncol(niche.par[[1]]$mu)
nspec <- length(niche.par)
npts <- 100
nell <- sapply(niche.par, function(x) nrow(x$mu))
species.names <- colnames(fish.size)
D <- combn(niso, 2)

Mu_all <- NULL
Sigma_all <- NULL
ell.coord_all <- NULL

for (ii in 1:nspec) {
  
  #ell.tmp <- array(NA, c(nell[ii], ncol(D), npts + 1, 2))
  for (jj in 1:nell[ii]) {
    
    for (kk in 1:ncol(D)) {
      
      Mu_prov <- as.data.frame(niche.par[[ii]]$mu[jj, D[, 
                                                        kk]])
      Mu_prov <- as.data.frame(t(Mu_prov))
      Mu_prov$SPP <- species.names[ii]
      Mu_prov$Iter <- jj
      
      Mu_all <- rbind(Mu_all, Mu_prov)
      
      
      
      Sigma_prov <- as.data.frame(niche.par[[ii]]$Sigma[D[, kk], D[, 
                                                                   kk], jj])
      
      Sigma_prov$SPP <- species.names[ii]
      Sigma_prov$Iter <- jj
      
      Sigma_all <- rbind(Sigma_all, Sigma_prov)
      
      ell.coord <- as.data.frame(ellipse(niche.par[[ii]]$mu[jj, D[, 
                                                                  kk]], V = niche.par[[ii]]$Sigma[D[, kk], D[, 
                                                                                                             kk], jj], alpha = alpha, n = npts))
      
      ell.coord_prov <- ell.coord
      ell.coord_prov$N_point <- c(1:101)
      ell.coord_prov$SPP <- species.names[ii]
      ell.coord_prov$Iter <- jj
      
      ell.coord_all <- rbind(ell.coord_all, ell.coord_prov)
    }
  }
}

rm(ii,jj,kk,alpha,niso)
rm(ell.coord)
rm(pfrac,nspec,npts,nell,species.names,D,Mu_all,niche.par,
   Sigma_all,Mu_prov,Sigma_prov,ell.coord_prov)

# Renumbering of points
Sp_names <- unique(ell.coord_all$SPP)
data_niche <- data_niche_FA %>% 
  rename(SPP = cluster)
ell.coord_all_Ntrue <- NULL

for (ii in 1:length(Sp_names)) {
  
  # data.frame avec data utilisées pour calculer les ellipses
  Ellipses_sel <- data_niche %>% 
    filter(SPP == Sp_names[ii])
  
  # data.frame avec les aires par ellipse (Area = aire, param = id de l'ellise)
  Areas_data_selec <- data.frame(Area = fish.size[,ii], param = c(1:1000))
  
  # Les points des ellipses (iso1 = valeurs des X, iso2 = valeurs des Y, group = espèce, n_point = N° du point, param = id de l'ellipse)
  Points_data_selec <- ell.coord_all %>% 
    filter(SPP == Sp_names[ii]) %>% 
    rename(iso1 = x, iso2 = y, group = SPP, n_point = N_point,
           param = Iter)
  
  # Moyennes des parametres Mu_1 et Mu_2
  # Ici on pourrait aussi faire les moyennes des X et Y des points d'ellipses
  MU1 <- mean(Ellipses_sel$NMDS1)
  MU2 <- mean(Ellipses_sel$NMDS2)
  Coords_MU <- data.frame(iso1 = MU1, iso2 = MU2)
  rm(MU1,MU2)
  
  # Calcul de la différence entre valeurs des X des points d'ellipse, et Mu1
  Points_data_selec0 <- Points_data_selec
  Points_data_selec0$DIFF <-  Points_data_selec0$iso1 - Coords_MU$iso1
  
  # Ne prendre que les points inférieurs à Mu2 en axe des Y
  Points_data_selec1 <- Points_data_selec0[Points_data_selec0$iso2 < Coords_MU$iso2,]
  Points_data_selec1 <- Points_data_selec1[Points_data_selec1$DIFF >= 0,]  # Différence positive, par défaut
  
  # Créer objet avec le min_DIFF associé au 1er point souhaité
  Mins_param <- as.data.frame(Points_data_selec1 %>% group_by(param) %>%
                                summarise(min_DIFF = min(DIFF)))
  
  Points_data_selec2 <- left_join(Points_data_selec0, Mins_param, by = "param")
  
  # Détection des premiers points (avec 2 filtres)
  First_points <- Points_data_selec2[Points_data_selec2$DIFF == Points_data_selec2$min_DIFF,]
  First_points <- First_points[First_points$iso2 < Coords_MU$iso2,]
  First_points <- First_points %>% dplyr::select(n_point, param) %>% dplyr::rename(First_point = n_point)
  # Ici nrow(First_points) doit correspondre au nombre d'ellipses
  
  # Renumérotation des points
  Points_data_selec <- left_join(Points_data_selec, First_points, by = "param")
  Points_data_selec$n_point_new <- ifelse(Points_data_selec$n_point == Points_data_selec$First_point, 1,
                                          ifelse(Points_data_selec$n_point > Points_data_selec$First_point,
                                                 Points_data_selec$n_point - Points_data_selec$First_point + 1,
                                                 ifelse(Points_data_selec$n_point < Points_data_selec$First_point,
                                                        Points_data_selec$n_point + max(Points_data_selec$n_point) - Points_data_selec$First_point + 1, NA)))
  
  ell.coord_all_Ntrue <- rbind(ell.coord_all_Ntrue,Points_data_selec)
  
  rm(Ellipses_sel,Areas_data_selec,Points_data_selec)
  rm(Coords_MU,Points_data_selec0,Points_data_selec1,Mins_param,
     Points_data_selec2,First_points)
}

rm(Sp_names,data_niche)

ell.coord_all_Ntrue <- ell.coord_all_Ntrue %>% 
  rename(x = iso1, y = iso2, N_point = n_point, SPP = group, Iter = param,
         N_true = n_point_new)

# Recovery of point coordinates for CI95% ellipses
ell.coord_inf_all <- NULL
ell.coord_sup_all <- NULL
ell.coord_mean95_all <- NULL

for (i in 1:ncol(fish.size)){
  
  # Selection of ellipse sizes for only one species
  fish_selec <- fish.size[,i]
  fish_area <- data.frame(area = fish_selec, Iter = c(1:1000))
  
  
  # Selection of point coordinates of all ellipses for one species
  fish_name <- colnames(fish.size)[i]
  ell.coord_fish <- ell.coord_all_Ntrue[ell.coord_all_Ntrue$SPP == fish_name,]
  
  
  # Merging ellipse sizes with data related to ellipses
  ell.coord_fish <- ell.coord_fish %>%
    left_join(fish_area, by = "Iter")
  
  
  # Identification of 2.5% smallest and 2.5% largest ellipses to get point coordinates
  Middle_area <- fish_area[order(fish_area$area, decreasing = F),] # Ranking ellipse sizes from smallest to largest
  QUANT_inf_area <- Middle_area[0:250,] # Selection of size of 2.5% smallest ellipses
  QUANT_sup_area <- Middle_area[750:1000,] # Selection of size of 2.5% largest ellipses
  MEAN_95_area <- Middle_area[25:975,] # Selection of size 95% middle sizes
  
  ell.coord_fish$CI_inf <- ifelse(ell.coord_fish$area %in% QUANT_inf_area[,1], 1, 0)
  ell.coord_fish$CI_sup <- ifelse(ell.coord_fish$area %in% QUANT_sup_area[,1], 1, 0)
  ell.coord_fish$MEAN_95 <- ifelse(ell.coord_fish$area %in% MEAN_95_area[,1], 1, 0)
  
  ell.coord_inf <- ell.coord_fish[ell.coord_fish$CI_inf == 1,] # Selection of point coordinates or 2.5% smallest ellipses
  ell.coord_sup <- ell.coord_fish[ell.coord_fish$CI_sup == 1,] # Selection of point coordinates or 2.5% smallest ellipses
  ell.coord_mean95 <- ell.coord_fish[ell.coord_fish$MEAN_95 == 1,]
  
  ell.coord_inf_all <- rbind(ell.coord_inf_all,ell.coord_inf)
  ell.coord_sup_all <- rbind(ell.coord_sup_all,ell.coord_sup)
  ell.coord_mean95_all <- rbind(ell.coord_mean95_all,ell.coord_mean95)
}

rm(i)
rm(fish_selec, fish_area, fish_name, ell.coord_fish,
   Middle_area, QUANT_inf_area, QUANT_sup_area,
   ell.coord_inf, ell.coord_sup)
rm(ell.coord_mean95, MEAN_95_area)
rm(fish.par, fish.size)

# Calculation of mean and CI95% ellipses for each cluster
ell.coord_MEAN_FA <- as.data.frame(
  ell.coord_mean95_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("1","2")))

ell.coord_CIinf_FA <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>%
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("1","2")))

ell.coord_CIsup_FA <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>%
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("1","2")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)

# Get nMDS data for further nMDS plotting
FA.scores <- as.data.frame(scores(MDS_2017, "species"))
FA.scores$FA <- rownames(FA.scores)

rm(MDS_2017)

# Plot nMDS and Bayes ellipses
clust_FAniches <- data_niche_FA %>% 
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(fill = cluster), size = 2.5, shape = 21, color = "grey50") + # add all individual points
  geom_polygon(data = ell.coord_MEAN_FA, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_FA, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_FA, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_text(data = FA.scores, aes(x = NMDS1, y = NMDS2,label = FA), size = 3, fontface = "bold", color = "grey23") + # add the FA labels
  scale_color_manual(values = c("#DE8E06","#48A4E3"))+
  scale_fill_manual(values = c("#DE8E06","#48A4E3"))+
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", color = "Cluster")+
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.title = element_blank())


## 3 / Plot dendrogram and FA niches ##############################################################################################

ggarrange(dendro, clust_FAniches, ncol = 2, labels = c("A.","B."),
          align = "hv")

rm(dendro, graph.cluster,clust_FAniches)
rm(data_niche_FA, ell.coord_MEAN_FA, ell.coord_CIinf_FA,
   ell.coord_CIsup_FA, FA.scores)


## 4 / Update database with cluster numbers ##############################################################################################

clust_FA <- clust_FA %>% 
  rename(organism_identifier = label)

SWO_data <- SWO_data %>% 
  left_join(clust_FA, by = "organism_identifier")

rm(clust_FA)



### V // Relationship between SI values and LJFL ##############################################################################################

## Test for correlations between SI values and LJFL
data_test <- SWO_data %>% 
  filter(!is.na(d13C),
         !is.na(length),
         !is.na(d13C),
         sex %in% c("M","F"))

shapiro.test(data_test$d13C)
shapiro.test(data_test$d15N)
shapiro.test(data_test$length)

cor.test(data_test$d13C, data_test$length, method = "kendall")
cor.test(data_test$d15N, data_test$length, method = "kendall")

data_test_M <- SWO_data %>% 
  filter(!is.na(d13C),
         !is.na(length),
         !is.na(d13C),
         sex == "M")

shapiro.test(data_test_M$d13C)
shapiro.test(data_test_M$d15N)
shapiro.test(data_test_M$length)

cor.test(data_test_M$d13C, data_test_M$length, method = "pearson")
cor.test(data_test_M$d15N, data_test_M$length, method = "kendall")

data_test_F <- SWO_data %>% 
  filter(!is.na(d13C),
         !is.na(length),
         !is.na(d13C),
         sex == "F")

shapiro.test(data_test_F$d13C)
shapiro.test(data_test_F$d15N)
shapiro.test(data_test_F$length)

cor.test(data_test_F$d13C, data_test_F$length, method = "kendall")
cor.test(data_test_F$d15N, data_test_F$length, method = "pearson")

## Test for linearity for each relationship
glm_d13C_length <- glm(d13C ~ length, data = data_test)
par(mfrow=c(2,2))
plot(glm_d13C_length)

glm_d15N_length <- glm(d15N ~ length, data = data_test)
par(mfrow=c(2,2))
plot(glm_d15N_length)

glm_d13C_length <- glm(d13C ~ length, data = data_test_M)
par(mfrow=c(2,2))
plot(glm_d13C_length)

glm_d15N_length <- glm(d15N ~ length, data = data_test_M)
par(mfrow=c(2,2))
plot(glm_d15N_length)

glm_d13C_length <- glm(d13C ~ length, data = data_test_F)
par(mfrow=c(2,2))
plot(glm_d13C_length)

glm_d15N_length <- glm(d15N ~ length, data = data_test_F)
par(mfrow=c(2,2))
plot(glm_d15N_length)

rm(data_test,data_test_M,data_test_F)
rm(glm_d13C_length,glm_d15N_length)

# Plot all relationships
d13C_length <- SWO_data %>% 
  filter(!is.na(length),
         !is.na(d13C),
         sex %in% c("M","F")) %>% 
   ggplot()+
  geom_point(aes(x = length, y = d13C), size = 2, alpha = 0.3)+
  geom_smooth(aes(x = length, y = d13C), color = "#DE8E06", method = "gam", se = TRUE)+
  geom_text(aes(x = max(length), y = max(d13C), label = "Kendall correlation\np > 0.05"),
            hjust = 1, vjust = 1, size = 3.5)+
  labs(x = "Fish length (cm)", y = "d13C")+
  theme_bw()+
  theme(legend.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(face = "italic"),
    axis.title.y = element_text(face = "italic"))

d13C_length_sex <- SWO_data %>% 
  filter(!is.na(d15N),
         !is.na(length),
         !is.na(d13C),
         sex %in% c("M","F")) %>% 
  mutate(sex = factor(sex, levels = c("M","F"))) %>% 
  ggplot()+
  geom_point(aes(x = length, y = d13C, color = sex), size = 2, alpha = 0.3)+
  geom_smooth(aes(x = length, y = d13C, color = sex), method = "gam", se = TRUE)+
  geom_text(aes(x = min(length), y = max(d13C), label = "Pearson and Kendall correlations\np > 0.05"),
            hjust = 0, vjust = 1, size = 3.5)+
  scale_color_manual(values = c("#257898","#E05F25"))+
  labs(x = "Fish length (cm)", y = "d13C")+
  theme_bw()+
  theme(legend.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(face = "italic"),
    axis.title.y = element_text(face = "italic"))

d15N_length <- SWO_data %>% 
  filter(!is.na(d15N),
         !is.na(length),
         !is.na(d13C),
         sex %in% c("M","F")) %>% 
  ggplot()+
  geom_point(aes(x = length, y = d15N), size = 2, alpha = 0.3)+
  geom_smooth(aes(x = length, y = d15N), color = "#DE8E06", method = "gam", se = TRUE)+ #
  geom_text(aes(x = max(length), y = min(d15N), label = "Kendall correlation\np < 0.001, cor = 0.24"),
            hjust = 1, vjust = 0, size = 3.5)+
  labs(x = "Lower jaw-fork length (cm)", y = "d15N")+
  theme_bw()+
  theme(legend.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(face = "italic"),
    axis.title.y = element_text(face = "italic"))

d15N_length_sex <- SWO_data %>% 
  filter(!is.na(d15N),
         !is.na(length),
         !is.na(d13C),
         sex %in% c("M","F")) %>% 
  mutate(sex = factor(sex, levels = c("M","F"))) %>% 
  ggplot()+
  geom_point(aes(x = length, y = d15N, color = sex), size = 2, alpha = 0.3)+
  geom_smooth(aes(x = length, y = d15N, color = sex, group = sex), method = "gam", se = TRUE)+
  geom_text(aes(x = max(length), y = min(d15N), label = "Kendall and Pearson correlations\np < 0.001, cor = 0.31\np = 0.003, cor = 0.52"),
            hjust = 1, vjust = 0, size = 3.5)+
  scale_color_manual(values = c("#257898","#E05F25"))+
  labs(x = "Lower jaw-fork length (cm)", y = "d15N")+
  theme_bw()+
  theme(legend.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(face = "italic"),
    axis.title.y = element_text(face = "italic"))

ggarrange(d13C_length, d15N_length, d13C_length_sex, d15N_length_sex,
          ncol = 2, nrow = 2, align = "hv",
          labels = c("A.","B.","C.","D."))

rm(d13C_length, d15N_length, d13C_length_sex, d15N_length_sex)



### VI // Contribution of physiological parameters, trophic ecology and season in explaining trace element bioaccumulation ##############################################################################################

## 1 / Arsenic ##############################################################################################

database_As <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,As)

## Step 1
mod_gam_As0 <- mgcv::gam(log10(As) ~ 1,
                         data = database_As, method = "REML")

mod_gam_As1 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_As, method = "REML")
summary(mod_gam_As1) # Complete model

mod_gam_As2 <- mgcv::gam(log10(As) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_As, method = "REML")
summary(mod_gam_As2) # length

mod_gam_As3 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_As, method = "REML")
summary(mod_gam_As3) # sex

mod_gam_As4 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                           factor(sex)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_As, method = "REML")
summary(mod_gam_As4) # new_cluster

mod_gam_As5 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_As, method = "REML")
summary(mod_gam_As5) # d13C

mod_gam_As6 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_As, method = "REML")
summary(mod_gam_As6) # d15N

mod_gam_As7 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_As, method = "REML")
summary(mod_gam_As7) # season

mod_gam_As8 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_As, method = "REML")
summary(mod_gam_As8) # Sex:length

mod_gam_As9 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_As, method = "REML")
summary(mod_gam_As9) # new_cluster:length

mod_gam_As10 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_As, method = "REML")
summary(mod_gam_As10) # season:length

mod_gam_As11 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_As, method = "REML")
summary(mod_gam_As11) # season:d15N

mod_gam_As12 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            d13C:length+
                            d15N:length,
                          data = database_As, method = "REML")
summary(mod_gam_As12) # season:d13C

mod_gam_As13 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_As, method = "REML")
summary(mod_gam_As13) # d13C:length

mod_gam_As14 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_As, method = "REML")
summary(mod_gam_As14) # d15N:length

AIC(mod_gam_As0,mod_gam_As1,mod_gam_As2,mod_gam_As3,mod_gam_As4,
    mod_gam_As5,mod_gam_As6,mod_gam_As7,mod_gam_As8,mod_gam_As9,
    mod_gam_As10,mod_gam_As11,mod_gam_As12,mod_gam_As13,mod_gam_As14)

## Step 2

# Begin again from n°11, remove d13C (mod 4), d15N (mod 6),
# new_cluster:length (mod 9), season:d13C (mod 12)

mod_gam_As11bis <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                               factor(sex)+
                               factor(new_cluster)+
                               factor(season)+
                               factor(sex):length+
                               factor(season):length+
                               d13C:length+
                               d15N:length,
                             data = database_As, method = "REML")
summary(mod_gam_As11bis)

AIC(mod_gam_As0,mod_gam_As1,mod_gam_As11,mod_gam_As11bis)

mod_gam_As15 <- mgcv::gam(log10(As) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            d13C:length+
                            d15N:length,
                          data = database_As, method = "REML")
summary(mod_gam_As15) # length

mod_gam_As16 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            d13C:length+
                            d15N:length,
                          data = database_As, method = "REML")
summary(mod_gam_As16) # sex

mod_gam_As17 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            d13C:length+
                            d15N:length,
                          data = database_As, method = "REML")
summary(mod_gam_As17) # clustFA

mod_gam_As18 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):length+
                            d13C:length+
                            d15N:length,
                          data = database_As, method = "REML")
summary(mod_gam_As18) # season

mod_gam_As19 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(season):length+
                            d13C:length+
                            d15N:length,
                          data = database_As, method = "REML")
summary(mod_gam_As19) # sex:length

mod_gam_As20 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            d13C:length+
                            d15N:length,
                          data = database_As, method = "REML")
summary(mod_gam_As20) # season:length

mod_gam_As21 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            d15N:length,
                          data = database_As, method = "REML")
summary(mod_gam_As21) # d13C:length

mod_gam_As22 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            d13C:length,
                          data = database_As, method = "REML")
summary(mod_gam_As22) # d15N:length

AIC(mod_gam_As0,mod_gam_As1,mod_gam_As11bis,mod_gam_As15,mod_gam_As16,
    mod_gam_As17,mod_gam_As18,mod_gam_As19,mod_gam_As20,mod_gam_As21,
    mod_gam_As22)

## Step 3

# Begin again from n°21

mod_gam_As23 <- mgcv::gam(log10(As) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            d15N:length,
                          data = database_As, method = "REML")
summary(mod_gam_As23) # length

mod_gam_As24 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            d15N:length,
                          data = database_As, method = "REML")
summary(mod_gam_As24) # sex

mod_gam_As25 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            d15N:length,
                          data = database_As, method = "REML")
summary(mod_gam_As25) # clustFA

mod_gam_As26 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):length+
                            d15N:length,
                          data = database_As, method = "REML")
summary(mod_gam_As26) # season

mod_gam_As27 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(season):length+
                            d15N:length,
                          data = database_As, method = "REML")
summary(mod_gam_As27) # sex:length

mod_gam_As28 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            d15N:length,
                          data = database_As, method = "REML")
summary(mod_gam_As28) # season:length

mod_gam_As29 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length,
                          data = database_As, method = "REML")
summary(mod_gam_As29) # d15N:length

AIC(mod_gam_As0,mod_gam_As1,mod_gam_As11bis,mod_gam_As21,mod_gam_As23,
    mod_gam_As24,mod_gam_As25,mod_gam_As26,mod_gam_As27,mod_gam_As28,
    mod_gam_As29)

## Step 4

# Begin again from n°29

mod_gam_As30 <- mgcv::gam(log10(As) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length,
                          data = database_As, method = "REML")
summary(mod_gam_As30) # length

mod_gam_As31 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length,
                          data = database_As, method = "REML")
summary(mod_gam_As31) # sex

mod_gam_As32 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length,
                          data = database_As, method = "REML")
summary(mod_gam_As32) # clustFA

mod_gam_As33 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):length,
                          data = database_As, method = "REML")
summary(mod_gam_As33) # season

mod_gam_As34 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(season):length,
                          data = database_As, method = "REML")
summary(mod_gam_As34) # sex:length

mod_gam_As35 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length,
                          data = database_As, method = "REML")
summary(mod_gam_As35) # season:length

AIC(mod_gam_As0,mod_gam_As1,mod_gam_As11bis,mod_gam_As21,mod_gam_As29,
    mod_gam_As30,mod_gam_As31,mod_gam_As32,mod_gam_As33,mod_gam_As34,
    mod_gam_As35)

## Step 5

# BEgin from 31 and 24, and remove from each

mod_gam_As36 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(season):length,
                          data = database_As, method = "REML")
summary(mod_gam_As36)

AIC(mod_gam_As0,mod_gam_As1,mod_gam_As11bis,mod_gam_As21,mod_gam_As29,
    mod_gam_As36)

# Best model is 31

rm(mod_gam_As0,mod_gam_As1,mod_gam_As2,mod_gam_As3,mod_gam_As4,
   mod_gam_As5,mod_gam_As6,mod_gam_As7,mod_gam_As8,mod_gam_As9,
   mod_gam_As10,mod_gam_As11,mod_gam_As12,mod_gam_As13,mod_gam_As14,
   mod_gam_As11bis,mod_gam_As15,mod_gam_As16,
   mod_gam_As17,mod_gam_As18,mod_gam_As19,mod_gam_As20,mod_gam_As21,
   mod_gam_As22,mod_gam_As23,
   mod_gam_As24,mod_gam_As25,mod_gam_As26,mod_gam_As27,mod_gam_As28,
   mod_gam_As29,mod_gam_As30,mod_gam_As31,mod_gam_As32,mod_gam_As33,mod_gam_As34,
   mod_gam_As35,mod_gam_As36)
rm(database_As)

## 2 / Length stand. As ##############################################################################################

database_As <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster)) %>% 
  select(length,sex,new_cluster,d13C,d15N,
         season,As_stand) %>% 
  rename(As = As_stand)

## Step 1

mod_gam_As0 <- mgcv::gam(log10(As) ~ 1,
                         data = database_As, method = "REML")

mod_gam_As1 <- mgcv::gam(log10(As) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_As, method = "REML")
summary(mod_gam_As1) # Complete model

mod_gam_As2 <- mgcv::gam(log10(As) ~ factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_As, method = "REML")
summary(mod_gam_As2) # sex

mod_gam_As3 <- mgcv::gam(log10(As) ~ factor(sex)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_As, method = "REML")
summary(mod_gam_As3) # new_cluster

mod_gam_As4 <- mgcv::gam(log10(As) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_As, method = "REML")
summary(mod_gam_As4) # d13C

mod_gam_As5 <- mgcv::gam(log10(As) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_As, method = "REML")
summary(mod_gam_As5) # d15N

mod_gam_As6 <- mgcv::gam(log10(As) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_As, method = "REML")
summary(mod_gam_As6) # season

mod_gam_As7 <- mgcv::gam(log10(As) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d13C,
                         data = database_As, method = "REML")
summary(mod_gam_As7) # season:d15N

mod_gam_As8 <- mgcv::gam(log10(As) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N,
                         data = database_As, method = "REML")
summary(mod_gam_As8) # season:d13C

AIC(mod_gam_As0,mod_gam_As1,mod_gam_As2,mod_gam_As3,mod_gam_As4,
    mod_gam_As5,mod_gam_As6,mod_gam_As7,mod_gam_As8)

## Step 2

# Begin from n°7, remove sex (mod 2, d13C (mod 4), d15N (mod 5)

mod_gam_As7bis <- mgcv::gam(log10(As) ~ factor(new_cluster)+
                              factor(season)+
                              factor(season):d13C,
                            data = database_As, method = "REML")
summary(mod_gam_As7bis) # season:d15N

AIC(mod_gam_As0,mod_gam_As1,mod_gam_As7,mod_gam_As7bis)

# Begin from n°7

mod_gam_As9 <- mgcv::gam(log10(As) ~ factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d13C,
                         data = database_As, method = "REML")
summary(mod_gam_As9) # sex

mod_gam_As10 <- mgcv::gam(log10(As) ~ factor(sex)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(season):d13C,
                          data = database_As, method = "REML")
summary(mod_gam_As10) # clustFA

mod_gam_As11 <- mgcv::gam(log10(As) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d15N)+
                            factor(season)+
                            factor(season):d13C,
                          data = database_As, method = "REML")
summary(mod_gam_As11) # d13C

mod_gam_As12 <- mgcv::gam(log10(As) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            factor(season)+
                            factor(season):d13C,
                          data = database_As, method = "REML")
summary(mod_gam_As12) # d15N

mod_gam_As13 <- mgcv::gam(log10(As) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season):d13C,
                          data = database_As, method = "REML")
summary(mod_gam_As13) # season

mod_gam_As14 <- mgcv::gam(log10(As) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season),
                          data = database_As, method = "REML")
summary(mod_gam_As14) # season:d13C

AIC(mod_gam_As0,mod_gam_As1,mod_gam_As7,mod_gam_As9,mod_gam_As10,mod_gam_As11,
    mod_gam_As12,mod_gam_As13,mod_gam_As14)

## Step 3

# Begin from n°12, remove d13C (mod 11)

mod_gam_As12bis <- mgcv::gam(log10(As) ~ factor(sex)+
                               factor(new_cluster)+
                               factor(season)+
                               factor(season):d13C,
                             data = database_As, method = "REML")
summary(mod_gam_As12bis)

AIC(mod_gam_As0,mod_gam_As1,mod_gam_As7,mod_gam_As12,mod_gam_As12bis)

mod_gam_As15 <- mgcv::gam(log10(As) ~ factor(new_cluster)+
                            factor(season)+
                            factor(season):d13C,
                          data = database_As, method = "REML")
summary(mod_gam_As15) # sex

mod_gam_As16 <- mgcv::gam(log10(As) ~ factor(sex)+
                            factor(season)+
                            factor(season):d13C,
                          data = database_As, method = "REML")
summary(mod_gam_As16) # clustFA

mod_gam_As17 <- mgcv::gam(log10(As) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(season):d13C,
                          data = database_As, method = "REML")
summary(mod_gam_As17) # season

mod_gam_As18 <- mgcv::gam(log10(As) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(season),
                          data = database_As, method = "REML")
summary(mod_gam_As18) # season:d13C

AIC(mod_gam_As0,mod_gam_As1,mod_gam_As7,mod_gam_As12bis,mod_gam_As15,
    mod_gam_As16,mod_gam_As17,mod_gam_As18)

## Step 4

# Begin from n°18
mod_gam_As19 <- mgcv::gam(log10(As) ~ factor(new_cluster)+
                            factor(season),
                          data = database_As, method = "REML")
summary(mod_gam_As19) # sex

mod_gam_As20 <- mgcv::gam(log10(As) ~ factor(sex)+
                            factor(season),
                          data = database_As, method = "REML")
summary(mod_gam_As20) # clustFA

mod_gam_As21 <- mgcv::gam(log10(As) ~ factor(sex)+
                            factor(new_cluster),
                          data = database_As, method = "REML")
summary(mod_gam_As21) # season

AIC(mod_gam_As0,mod_gam_As1,mod_gam_As7,mod_gam_As12bis,mod_gam_As18,
    mod_gam_As19,mod_gam_As20,mod_gam_As21)

## Step 5

# Begin from n°19

mod_gam_As22 <- mgcv::gam(log10(As) ~ factor(season),
                          data = database_As, method = "REML")
summary(mod_gam_As22) # clustFA

mod_gam_As23 <- mgcv::gam(log10(As) ~ factor(new_cluster),
                          data = database_As, method = "REML")
summary(mod_gam_As23) # season

AIC(mod_gam_As0,mod_gam_As1,mod_gam_As7,mod_gam_As12bis,mod_gam_As18,
    mod_gam_As19,mod_gam_As22,mod_gam_As23)

# Best model is n°19

rm(mod_gam_As0,mod_gam_As1,mod_gam_As2,mod_gam_As3,mod_gam_As4,
   mod_gam_As5,mod_gam_As6,mod_gam_As7,mod_gam_As8,mod_gam_As9,
   mod_gam_As10,mod_gam_As11,mod_gam_As12,mod_gam_As13,mod_gam_As14,
   mod_gam_As15,mod_gam_As16,mod_gam_As17,mod_gam_As18,mod_gam_As19,
   mod_gam_As20,mod_gam_As21,mod_gam_As22,mod_gam_As23,
   mod_gam_As7bis,mod_gam_As12bis)
rm(database_As)

## 3 / As & length-stand. As - Plot response variable vs explaining variables ##############################################################################################

## LOG-TRANSFORMED As
database_As_ls <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,As)

# Best model is n°31

mod_gam_As31 <- gam::gam(log10(As) ~ length +
                           factor(new_cluster)+
                           factor(season)+
                           factor(sex):length+
                           factor(season):length,
                         data = database_As_ls)
summary(mod_gam_As31)

mod_gam_As31bis <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                               factor(new_cluster)+
                               factor(season)+
                               factor(sex):length+
                               factor(season):length,
                             data = database_As_ls, method = "REML")
summary(mod_gam_As31bis)

## LOG-TRANSFORMED LENGTH-STAND. As
database_As <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster)) %>% 
  select(length,sex,new_cluster,d13C,d15N,
         season,As_stand) %>% 
  rename(As = As_stand)

# Best model is n°19
mod_gam_As19 <- gam::gam(log10(As) ~ factor(new_cluster)+
                           factor(season),
                         data = database_As)
summary(mod_gam_As19)

mod_gam_As19bis <- mgcv::gam(log10(As) ~ factor(new_cluster)+
                               factor(season),
                             data = database_As, method = "REML")
summary(mod_gam_As19bis)

par(mfrow=c(2,3))
plot(mod_gam_As31,se = TRUE, ylab = "", xlab = "", las = 1)
plot(mod_gam_As19,se = TRUE, ylab = "", xlab = "", las = 1)

rm(database_As_ls,mod_gam_As31,mod_gam_As31bis)
rm(database_As,mod_gam_As19,mod_gam_As19bis)


## 4 / Cadmium ##############################################################################################

database_Cd <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Cd)

## Step 1
mod_gam_Cd0 <- mgcv::gam(log10(Cd) ~ 1,
                         data = database_Cd, method = "REML")

mod_gam_Cd1 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd1) # Complete model

mod_gam_Cd2 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd2) # length

mod_gam_Cd3 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd3) # sex

mod_gam_Cd4 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                           factor(sex)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd4) # new_cluster

mod_gam_Cd5 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd5) # d13C

mod_gam_Cd6 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd6) # d15N

mod_gam_Cd7 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd7) # season

mod_gam_Cd8 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd8) # Sex:length

mod_gam_Cd9 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd9) # new_cluster:length

mod_gam_Cd10 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd10) # season:length

mod_gam_Cd11 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd11) # season:d15N

mod_gam_Cd12 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            d13C:length+
                            d15N:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd12) # season:d13C

mod_gam_Cd13 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd13) # d13C:length

mod_gam_Cd14 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd14) # d15N:length

AIC(mod_gam_Cd0,mod_gam_Cd1,mod_gam_Cd2,mod_gam_Cd3,mod_gam_Cd4,
    mod_gam_Cd5,mod_gam_Cd6,mod_gam_Cd7,mod_gam_Cd8,mod_gam_Cd9,
    mod_gam_Cd10,mod_gam_Cd11,mod_gam_Cd12,mod_gam_Cd13,mod_gam_Cd14)

## Step 2

# Begin from n°7, remove d13C (mod 5), d15N (mod 6), clustFA:length (mod 9)

mod_gam_Cd7bis <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                              factor(sex)+
                              factor(new_cluster)+
                              factor(sex):length+
                              factor(season):length+
                              factor(season):d15N+
                              factor(season):d13C+
                              d13C:length+
                              d15N:length,
                            data = database_Cd, method = "REML")
summary(mod_gam_Cd7bis)

AIC(mod_gam_Cd0,mod_gam_Cd1,mod_gam_Cd7,mod_gam_Cd7bis)

mod_gam_Cd15 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd15) # length

mod_gam_Cd16 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd13) # sex

mod_gam_Cd17 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd17) # clustFA

mod_gam_Cd18 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd18) # sex:length

mod_gam_Cd19 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd19) # season:length

mod_gam_Cd20 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd20) # season:d15N

mod_gam_Cd21 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            d13C:length+
                            d15N:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd21) # season:d13C

mod_gam_Cd22 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd22) # d13C:length

mod_gam_Cd23 <- mgcv::gam(log10(Cd) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd23) # d15N:length

AIC(mod_gam_Cd0,mod_gam_Cd1,mod_gam_Cd7bis,mod_gam_Cd15,mod_gam_Cd16,
    mod_gam_Cd17,mod_gam_Cd18,mod_gam_Cd19,mod_gam_Cd20,mod_gam_Cd21,
    mod_gam_Cd22,mod_gam_Cd23)

## Step 3

# Begin again from n°23, remove length (mod 15)

mod_gam_Cd23bis <- mgcv::gam(log10(Cd) ~ factor(sex)+
                               factor(new_cluster)+
                               factor(sex):length+
                               factor(season):length+
                               factor(season):d15N+
                               factor(season):d13C+
                               d13C:length,
                             data = database_Cd, method = "REML")
summary(mod_gam_Cd23bis)

AIC(mod_gam_Cd0,mod_gam_Cd1,mod_gam_Cd7bis,mod_gam_Cd23,mod_gam_Cd23bis)

mod_gam_Cd24 <- mgcv::gam(log10(Cd) ~ factor(new_cluster)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd24) # sex

mod_gam_Cd25 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd25) # clustFA

mod_gam_Cd26 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd26) # sex:length

mod_gam_Cd27 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd27) # season:length

mod_gam_Cd28 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd28) # season:d15N

mod_gam_Cd29 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            d13C:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd29) # season:d13C

mod_gam_Cd30 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd30) # d13C:length

AIC(mod_gam_Cd0,mod_gam_Cd1,mod_gam_Cd7bis,mod_gam_Cd23bis,mod_gam_Cd24,mod_gam_Cd25,
    mod_gam_Cd26,mod_gam_Cd27,mod_gam_Cd28,mod_gam_Cd29,mod_gam_Cd30)

## Step 4

# Begin from n°29

mod_gam_Cd31 <- mgcv::gam(log10(Cd) ~ factor(new_cluster)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            d13C:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd31) # sex

mod_gam_Cd32 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            d13C:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd32) # clustFA

mod_gam_Cd33 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(season):length+
                            factor(season):d15N+
                            d13C:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd33) # sex:length

mod_gam_Cd34 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):d15N+
                            d13C:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd34) # season:length

mod_gam_Cd35 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):length+
                            d13C:length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd35) # season:d15N

mod_gam_Cd36 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd36) # d13C:length

AIC(mod_gam_Cd0,mod_gam_Cd1,mod_gam_Cd7bis,mod_gam_Cd23bis,mod_gam_Cd29,mod_gam_Cd31,
    mod_gam_Cd32,mod_gam_Cd33,mod_gam_Cd34,mod_gam_Cd35,mod_gam_Cd36)

## Step 5

# Begin from n°36

mod_gam_Cd37 <- mgcv::gam(log10(Cd) ~ factor(new_cluster)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd37) # sex

mod_gam_Cd38 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd38) # clust FA

mod_gam_Cd39 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(season):length+
                            factor(season):d15N,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd39) # sex:length

mod_gam_Cd40 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):d15N,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd40) # season:length

mod_gam_Cd41 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):length,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd41) # season:d15N

AIC(mod_gam_Cd0,mod_gam_Cd1,mod_gam_Cd7bis,mod_gam_Cd23bis,mod_gam_Cd29,
    mod_gam_Cd36,mod_gam_Cd37,mod_gam_Cd38,mod_gam_Cd39,mod_gam_Cd40,
    mod_gam_Cd41)

# Best model is 36

rm(mod_gam_Cd0,mod_gam_Cd1,mod_gam_Cd2,mod_gam_Cd3,mod_gam_Cd4,
   mod_gam_Cd5,mod_gam_Cd6,mod_gam_Cd7,mod_gam_Cd8,mod_gam_Cd9,
   mod_gam_Cd10,mod_gam_Cd11,mod_gam_Cd12,mod_gam_Cd13,mod_gam_Cd14,
   mod_gam_Cd7bis,mod_gam_Cd15,mod_gam_Cd16,
   mod_gam_Cd17,mod_gam_Cd18,mod_gam_Cd19,mod_gam_Cd20,mod_gam_Cd21,
   mod_gam_Cd22,mod_gam_Cd23,mod_gam_Cd24,mod_gam_Cd23bis,mod_gam_Cd25,
   mod_gam_Cd26,mod_gam_Cd27,mod_gam_Cd28,mod_gam_Cd29,mod_gam_Cd30,
   mod_gam_Cd31,mod_gam_Cd32,mod_gam_Cd33,mod_gam_Cd34,mod_gam_Cd35,
   mod_gam_Cd36,mod_gam_Cd37,mod_gam_Cd38,mod_gam_Cd39,mod_gam_Cd40,
   mod_gam_Cd41)
rm(database_Cd)


## 5 / Length-stand. Cd ##############################################################################################

database_Cd <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Cd_stand) %>% 
  rename(Cd = Cd_stand)

## Step 1

mod_gam_Cd0 <- mgcv::gam(log10(Cd) ~ 1,
                         data = database_Cd, method = "REML")

mod_gam_Cd1 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd1) # Complete model

mod_gam_Cd2 <- mgcv::gam(log10(Cd) ~ factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd2) # sex

mod_gam_Cd3 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd3) # new_cluster

mod_gam_Cd4 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd4) # d13C

mod_gam_Cd5 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd5) # d15N

mod_gam_Cd6 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd6) # season

mod_gam_Cd7 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d13C,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd7) # season:d15N

mod_gam_Cd8 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd8) # season:d13C

AIC(mod_gam_Cd0,mod_gam_Cd1,mod_gam_Cd2,mod_gam_Cd3,mod_gam_Cd4,
    mod_gam_Cd5,mod_gam_Cd6,mod_gam_Cd7,mod_gam_Cd8)

## Step 2

# Begin from n°8, remove d13C (mod 4), d15N (mod 5), season:d15N (mod 7)

mod_gam_Cd8bis <- mgcv::gam(log10(Cd) ~ factor(sex)+
                              factor(new_cluster)+
                              factor(season)+
                              factor(season):d15N,
                            data = database_Cd, method = "REML")
summary(mod_gam_Cd8bis)

AIC(mod_gam_Cd0,mod_gam_Cd1,mod_gam_Cd8,mod_gam_Cd8bis)

mod_gam_Cd9 <- mgcv::gam(log10(Cd) ~ factor(new_cluster)+
                           factor(season)+
                           factor(season):d15N,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd9) # sex

mod_gam_Cd10 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(season)+
                            factor(season):d15N,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd10) # clustFA

mod_gam_Cd11 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(season):d15N,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd11) # season

mod_gam_Cd12 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(season),
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd12) # season:d15N

AIC(mod_gam_Cd0,mod_gam_Cd1,mod_gam_Cd8bis,mod_gam_Cd9,mod_gam_Cd10,
    mod_gam_Cd11,mod_gam_Cd12)

## Step 3

# BEgin from n°9

mod_gam_Cd13 <- mgcv::gam(log10(Cd) ~ factor(season)+
                            factor(season):d15N,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd13) # clustFA

mod_gam_Cd14 <- mgcv::gam(log10(Cd) ~ factor(new_cluster)+
                            factor(season):d15N,
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd14) # season

mod_gam_Cd15 <- mgcv::gam(log10(Cd) ~ factor(new_cluster)+
                            factor(season),
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd15) # season:d15N

AIC(mod_gam_Cd0,mod_gam_Cd1,mod_gam_Cd8bis,mod_gam_Cd9,mod_gam_Cd13,
    mod_gam_Cd14,mod_gam_Cd15)

# Best model is n°9

rm(mod_gam_Cd0,mod_gam_Cd1,mod_gam_Cd2,mod_gam_Cd3,mod_gam_Cd4,
   mod_gam_Cd5,mod_gam_Cd6,mod_gam_Cd7,mod_gam_Cd8,mod_gam_Cd9,
   mod_gam_Cd10,mod_gam_Cd11,mod_gam_Cd12,mod_gam_Cd13,mod_gam_Cd14,
   mod_gam_Cd15,mod_gam_Cd8bis)
rm(database_Cd)

## 6 / Cd & length-stand. Cd - Plot response variable vs explaining variables ##############################################################################################

## LOG-TRANSFORMED Cd
database_Cd_ls <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Cd)

# Best model is n°36

mod_gam_Cd36 <- gam::gam(log10(Cd) ~ factor(sex)+
                           factor(new_cluster)+
                           factor(sex):length+
                           factor(season):length+
                           factor(season):d15N,
                         data = database_Cd_ls)
summary(mod_gam_Cd36)

mod_gam_Cd36bis <- mgcv::gam(log10(Cd) ~ factor(sex)+
                               factor(new_cluster)+
                               factor(sex):length+
                               factor(season):length+
                               factor(season):d15N,
                             data = database_Cd_ls, method = "REML")
summary(mod_gam_Cd36bis)

## LOG-TRANSFORMED LENGTH-STAND. Cd
database_Cd <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Cd_stand) %>% 
  rename(Cd = Cd_stand)

# Best model is n°9
mod_gam_Cd9 <- gam::gam(log10(Cd) ~ factor(new_cluster)+
                          factor(season)+
                          factor(season):d15N,
                        data = database_Cd)
summary(mod_gam_Cd9)

mod_gam_Cd9bis <- mgcv::gam(log10(Cd) ~ factor(new_cluster)+
                              factor(season)+
                              factor(season):d15N,
                            data = database_Cd, method = "REML")
summary(mod_gam_Cd9bis)

par(mfrow=c(2,2))
plot(mod_gam_Cd36,se = TRUE, ylab = "", xlab = "",las = 1)
plot(mod_gam_Cd9,se = TRUE, ylab = "", xlab = "", las = 1)

rm(database_Cd_ls,mod_gam_Cd36,mod_gam_Cd36bis)
rm(database_Cd,mod_gam_Cd9,mod_gam_Cd9bis)


## 7 / Cobalt ##############################################################################################

database_Co <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  dplyr::select(length,sex,new_cluster,d13C,d15N,season,Co)

## Step 1
mod_gam_Co0 <- mgcv::gam(log10(Co) ~ 1,
                         data = database_Co, method = "REML")

mod_gam_Co1 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Co, method = "REML")
summary(mod_gam_Co1) # Complete model

mod_gam_Co2 <- mgcv::gam(log10(Co) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Co, method = "REML")
summary(mod_gam_Co2) # length

mod_gam_Co3 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Co, method = "REML")
summary(mod_gam_Co3) # sex

mod_gam_Co4 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                           factor(sex)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Co, method = "REML")
summary(mod_gam_Co4) # new_cluster

mod_gam_Co5 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Co, method = "REML")
summary(mod_gam_Co5) # d13C

mod_gam_Co6 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Co, method = "REML")
summary(mod_gam_Co6) # d15N

mod_gam_Co7 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Co, method = "REML")
summary(mod_gam_Co7) # season

mod_gam_Co8 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Co, method = "REML")
summary(mod_gam_Co8) # Sex:length

mod_gam_Co9 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Co, method = "REML")
summary(mod_gam_Co9) # new_cluster:length

mod_gam_Co10 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co10) # season:length

mod_gam_Co11 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co11) # season:d15N

mod_gam_Co12 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            d13C:length+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co12) # season:d13C

mod_gam_Co13 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co13) # d13C:length

mod_gam_Co14 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co14) # d15N:length

AIC(mod_gam_Co0,mod_gam_Co1,mod_gam_Co2,mod_gam_Co3,mod_gam_Co4,
    mod_gam_Co5,mod_gam_Co6,mod_gam_Co7,mod_gam_Co8,mod_gam_Co9,
    mod_gam_Co10,mod_gam_Co11,mod_gam_Co12,mod_gam_Co13,mod_gam_Co14)

## Step 2

# Begin from n°13, remove sex (mod 3), clustFA:length (mod 9)

mod_gam_Co13bis <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                               factor(new_cluster)+
                               s(d13C)+
                               s(d15N)+
                               factor(season)+
                               factor(sex):length+
                               factor(season):length+
                               factor(season):d15N+
                               factor(season):d13C+
                               d15N:length,
                             data = database_Co, method = "REML")
summary(mod_gam_Co13bis)

AIC(mod_gam_Co0,mod_gam_Co1,mod_gam_Co13,mod_gam_Co13bis)

mod_gam_Co15 <- mgcv::gam(log10(Co) ~ factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co15) # length

mod_gam_Co16 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co16) # clustFA

mod_gam_Co17 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co17) # d13C

mod_gam_Co18 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co18) # d15N

mod_gam_Co19 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co19) # season

mod_gam_Co20 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co20) # sex:length

mod_gam_Co21 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co21) # season:length

mod_gam_Co22 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co22) # season:d15N

mod_gam_Co23 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co23) # season:d13C

mod_gam_Co24 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co24) # d15N:length

AIC(mod_gam_Co0,mod_gam_Co1,mod_gam_Co13bis,mod_gam_Co15,mod_gam_Co16,
    mod_gam_Co17,mod_gam_Co18,mod_gam_Co19,mod_gam_Co20,mod_gam_Co21,
    mod_gam_Co22,mod_gam_Co23,mod_gam_Co24)

## Step 3

# Begin from n°17, remove clustFA

mod_gam_Co17bis <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                               s(d15N)+
                               factor(season)+
                               factor(sex):length+
                               factor(season):length+
                               factor(season):d15N+
                               factor(season):d13C+
                               d15N:length,
                             data = database_Co, method = "REML")
summary(mod_gam_Co17bis)

AIC(mod_gam_Co0,mod_gam_Co1,mod_gam_Co13bis,mod_gam_Co17,mod_gam_Co17bis)



mod_gam_Co25 <- mgcv::gam(log10(Co) ~ s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co25) # length

mod_gam_Co26 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co26) # d15N

mod_gam_Co27 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            s(d15N)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co27) # season

mod_gam_Co28 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            s(d15N)+
                            factor(season)+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co28) # sex:length

mod_gam_Co29 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co29) # season:length

mod_gam_Co30 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co30) # season:d15N

mod_gam_Co31 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            d15N:length,
                          data = database_Co, method = "REML")
summary(mod_gam_Co31) # season:d13C

mod_gam_Co31 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co31) # d15N:length

AIC(mod_gam_Co0,mod_gam_Co1,mod_gam_Co13bis,mod_gam_Co17bis,mod_gam_Co25,
    mod_gam_Co26,mod_gam_Co27,mod_gam_Co28,mod_gam_Co29,mod_gam_Co30,
    mod_gam_Co31)

## Step 4

# Begin from 31

mod_gam_Co32 <- mgcv::gam(log10(Co) ~ s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co32) # length

mod_gam_Co33 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co33) # d15N

mod_gam_Co34 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            s(d15N)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co34) # season

mod_gam_Co35 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            s(d15N)+
                            factor(season)+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co35) # sex:length

mod_gam_Co36 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co36) # season:length

mod_gam_Co37 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co37) # season:d15N

mod_gam_Co38 <- mgcv::gam(log10(Co) ~ s(length, k = 3)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N,
                          data = database_Co, method = "REML")
summary(mod_gam_Co38) # season:d13C

AIC(mod_gam_Co0,mod_gam_Co1,mod_gam_Co13bis,mod_gam_Co17bis,mod_gam_Co31,
    mod_gam_Co32,mod_gam_Co33,mod_gam_Co34,mod_gam_Co35,mod_gam_Co36,
    mod_gam_Co37,mod_gam_Co38)

## Step 5

# Begin from n°32

mod_gam_Co39 <- mgcv::gam(log10(Co) ~ factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co39) # d15N

mod_gam_Co40 <- mgcv::gam(log10(Co) ~ s(d15N)+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co40) # sex:length

mod_gam_Co41 <- mgcv::gam(log10(Co) ~ s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co41) # season:length

mod_gam_Co42 <- mgcv::gam(log10(Co) ~ s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co42) # season:d15N

mod_gam_Co43 <- mgcv::gam(log10(Co) ~ s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N,
                          data = database_Co, method = "REML")
summary(mod_gam_Co43) # season:d13C

AIC(mod_gam_Co0,mod_gam_Co1,mod_gam_Co13bis,mod_gam_Co17bis,mod_gam_Co31,
    mod_gam_Co32,mod_gam_Co39,mod_gam_Co40,mod_gam_Co41,mod_gam_Co42,
    mod_gam_Co43)

## Best model is 32

rm(mod_gam_Co0,mod_gam_Co1,mod_gam_Co2,mod_gam_Co3,mod_gam_Co4,
   mod_gam_Co5,mod_gam_Co6,mod_gam_Co7,mod_gam_Co8,mod_gam_Co9,
   mod_gam_Co10,mod_gam_Co11,mod_gam_Co12,mod_gam_Co13,mod_gam_Co14,
   mod_gam_Co13bis,mod_gam_Co15,mod_gam_Co16,
   mod_gam_Co17,mod_gam_Co18,mod_gam_Co19,mod_gam_Co20,mod_gam_Co21,
   mod_gam_Co22,mod_gam_Co23,mod_gam_Co24,
   mod_gam_Co25,
   mod_gam_Co26,mod_gam_Co27,mod_gam_Co28,mod_gam_Co29,mod_gam_Co30,
   mod_gam_Co31,mod_gam_Co32,
   mod_gam_Co33,mod_gam_Co34,mod_gam_Co35,mod_gam_Co36,mod_gam_Co37,
   mod_gam_Co38,mod_gam_Co39,
   mod_gam_Co40,mod_gam_Co41,mod_gam_Co42,mod_gam_Co43,
   mod_gam_Co17bis)
rm(database_Co)


## 8 / Length-stand. Co ##############################################################################################

database_Co <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Co_stand) %>% 
  rename(Co = Co_stand)

## Step 1

mod_gam_Co0 <- mgcv::gam(log10(Co) ~ 1,
                         data = database_Co, method = "REML")

mod_gam_Co1 <- mgcv::gam(log10(Co) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Co, method = "REML")
summary(mod_gam_Co1) # Complete model

mod_gam_Co2 <- mgcv::gam(log10(Co) ~ factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Co, method = "REML")
summary(mod_gam_Co2) # sex

mod_gam_Co3 <- mgcv::gam(log10(Co) ~ factor(sex)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Co, method = "REML")
summary(mod_gam_Co3) # new_cluster

mod_gam_Co4 <- mgcv::gam(log10(Co) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Co, method = "REML")
summary(mod_gam_Co4) # d13C

mod_gam_Co5 <- mgcv::gam(log10(Co) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Co, method = "REML")
summary(mod_gam_Co5) # d15N

mod_gam_Co6 <- mgcv::gam(log10(Co) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Co, method = "REML")
summary(mod_gam_Co6) # season

mod_gam_Co7 <- mgcv::gam(log10(Co) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d13C,
                         data = database_Co, method = "REML")
summary(mod_gam_Co7) # season:d15N

mod_gam_Co8 <- mgcv::gam(log10(Co) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N,
                         data = database_Co, method = "REML")
summary(mod_gam_Co8) # season:d13C

AIC(mod_gam_Co0,mod_gam_Co1,mod_gam_Co2,mod_gam_Co3,mod_gam_Co4,
    mod_gam_Co5,mod_gam_Co6,mod_gam_Co7,mod_gam_Co8)

## Step 2

# Begin from n°7, remove season (mod 6)

mod_gam_Co7bis <- mgcv::gam(log10(Co) ~ factor(sex)+
                              factor(new_cluster)+
                              s(d13C)+
                              s(d15N)+
                              factor(season):d13C,
                            data = database_Co, method = "REML")
summary(mod_gam_Co7bis)

AIC(mod_gam_Co0,mod_gam_Co1,mod_gam_Co7,mod_gam_Co7bis)

mod_gam_Co9 <- mgcv::gam(log10(Co) ~ factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season):d13C,
                         data = database_Co, method = "REML")
summary(mod_gam_Co9) # sex

mod_gam_Co10 <- mgcv::gam(log10(Co) ~ factor(sex)+
                            s(d13C)+
                            s(d15N)+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co10) # clustFA

mod_gam_Co11 <- mgcv::gam(log10(Co) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d15N)+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co11) # d13C

mod_gam_Co12 <- mgcv::gam(log10(Co) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co12) # d15N

mod_gam_Co13 <- mgcv::gam(log10(Co) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N),
                          data = database_Co, method = "REML")
summary(mod_gam_Co13) # season:d13C

AIC(mod_gam_Co0,mod_gam_Co1,mod_gam_Co7bis,mod_gam_Co9,mod_gam_Co10,
    mod_gam_Co11,mod_gam_Co12,mod_gam_Co13)

## Step 3

# Begin from n°10

mod_gam_Co14 <- mgcv::gam(log10(Co) ~ s(d13C)+
                            s(d15N)+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co14) # sex

mod_gam_Co15 <- mgcv::gam(log10(Co) ~ factor(sex)+
                            s(d15N)+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co15) # d13C

mod_gam_Co16 <- mgcv::gam(log10(Co) ~ factor(sex)+
                            s(d13C)+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co16) # d15N

mod_gam_Co17 <- mgcv::gam(log10(Co) ~ factor(sex)+
                            s(d13C)+
                            s(d15N),
                          data = database_Co, method = "REML")
summary(mod_gam_Co17) # season:d13C

AIC(mod_gam_Co0,mod_gam_Co1,mod_gam_Co7bis,mod_gam_Co10,mod_gam_Co14,
    mod_gam_Co15,mod_gam_Co16,mod_gam_Co17)

## Step 4

# Begin from n°14

mod_gam_Co18 <- mgcv::gam(log10(Co) ~ s(d15N)+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co18) # d13C

mod_gam_Co19 <- mgcv::gam(log10(Co) ~ s(d13C)+
                            factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co19) # d15N

mod_gam_Co20 <- mgcv::gam(log10(Co) ~ s(d13C)+
                            s(d15N),
                          data = database_Co, method = "REML")
summary(mod_gam_Co20) # season:d13C

AIC(mod_gam_Co0,mod_gam_Co1,mod_gam_Co7bis,mod_gam_Co10,mod_gam_Co14,
    mod_gam_Co18,mod_gam_Co19,mod_gam_Co20)

## Step 5

# Begin from n°18

mod_gam_Co21 <- mgcv::gam(log10(Co) ~ factor(season):d13C,
                          data = database_Co, method = "REML")
summary(mod_gam_Co21)

mod_gam_Co22 <- mgcv::gam(log10(Co) ~ s(d15N),
                          data = database_Co, method = "REML")
summary(mod_gam_Co22)

AIC(mod_gam_Co0,mod_gam_Co1,mod_gam_Co7bis,mod_gam_Co10,mod_gam_Co14,
    mod_gam_Co18,mod_gam_Co21,mod_gam_Co22)

## Best model is n°18

rm(mod_gam_Co0,mod_gam_Co1,mod_gam_Co2,mod_gam_Co3,mod_gam_Co4,
   mod_gam_Co5,mod_gam_Co6,mod_gam_Co7,mod_gam_Co8,mod_gam_Co9,
   mod_gam_Co10,mod_gam_Co11,mod_gam_Co12,mod_gam_Co13,mod_gam_Co14,
   mod_gam_Co15,mod_gam_Co16,mod_gam_Co17,mod_gam_Co18,mod_gam_Co19,
   mod_gam_Co20,mod_gam_Co21,mod_gam_Co22,mod_gam_Co7bis)
rm(database_Co)

## 9 / Co & length-stand. Co - Plot response variable vs explaining variables ##############################################################################################

## LOG-TRANSFORMED Co
database_Co <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  dplyr::select(length,sex,new_cluster,d13C,d15N,season,Co)

# Best model is n°32

mod_gam_Co32 <- gam::gam(log10(Co) ~ d15N+
                           factor(season)+
                           factor(sex):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Co)
summary(mod_gam_Co32)

mod_gam_Co32bis <- mgcv::gam(log10(Co) ~ s(d15N)+
                               factor(season)+
                               factor(sex):length+
                               factor(season):length+
                               factor(season):d15N+
                               factor(season):d13C,
                             data = database_Co, method = "REML")
summary(mod_gam_Co32bis)

## LOG-TRANSFORMED LENGTH-STAND. Co
database_Co_ls <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Co_stand) %>% 
  rename(Co = Co_stand)

# Best model is n°18
mod_gam_Co18 <- gam::gam(log10(Co) ~ d15N+
                           #factor(season)+
                           factor(season):d13C,
                         data = database_Co_ls)
summary(mod_gam_Co18)
AIC(mod_gam_Co18)

mod_gam_Co18bis <- mgcv::gam(log10(Co) ~ s(d15N)+
                               factor(season):d13C,
                             data = database_Co_ls, method = "REML")
summary(mod_gam_Co18bis)

par(mfrow=c(2,2))
plot(mod_gam_Co32,se = TRUE, ylab = "", xlab = "", las = 1)
plot(mod_gam_Co18,se = TRUE, ylab = "", xlab = "", las = 1)

rm(database_Co,mod_gam_Co32,mod_gam_Co32bis)
rm(database_Co_ls,mod_gam_Co18,mod_gam_Co18bis)


## 10 / Copper ##############################################################################################

database_Cu <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Cu)

## Step 1
mod_gam_Cu0 <- mgcv::gam(log10(Cu) ~ 1,
                         data = database_Cu, method = "REML")

mod_gam_Cu1 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu1) # Complete model

mod_gam_Cu2 <- mgcv::gam(log10(Cu) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu2) # length

mod_gam_Cu3 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu3) # sex

mod_gam_Cu4 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                           factor(sex)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu4) # new_cluster

mod_gam_Cu5 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu5) # d13C

mod_gam_Cu6 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu6) # d15N

mod_gam_Cu7 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu7) # season

mod_gam_Cu8 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu8) # Sex:length

mod_gam_Cu9 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu9) # new_cluster:length

mod_gam_Cu10 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu10) # season:length

mod_gam_Cu11 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu11) # season:d15N

mod_gam_Cu12 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            d13C:length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu12) # season:d13C

mod_gam_Cu13 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu13) # d13C:length

mod_gam_Cu14 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu14) # d15N:length

AIC(mod_gam_Cu0,mod_gam_Cu1,mod_gam_Cu2,mod_gam_Cu3,mod_gam_Cu4,
    mod_gam_Cu5,mod_gam_Cu6,mod_gam_Cu7,mod_gam_Cu8,mod_gam_Cu9,
    mod_gam_Cu10,mod_gam_Cu11,mod_gam_Cu12,mod_gam_Cu13,mod_gam_Cu14)

## Step 2

# Begin from n°11, remove d15N (mod 5), d13C (mod 6), season:length
# (mod 10)

mod_gam_Cu11bis <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                               factor(sex)+
                               factor(new_cluster)+
                               factor(season)+
                               factor(sex):length+
                               factor(new_cluster):length+
                               factor(season):d13C+
                               d13C:length+
                               d15N:length,
                             data = database_Cu, method = "REML")
summary(mod_gam_Cu11bis)

AIC(mod_gam_Cu0,mod_gam_Cu1,mod_gam_Cu11bis)

mod_gam_Cu15 <- mgcv::gam(log10(Cu) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu15) # length

mod_gam_Cu16 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu16) # sex

mod_gam_Cu17 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu17) # clustFA

mod_gam_Cu18 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu18) # season

mod_gam_Cu19 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu19) # sex:length

mod_gam_Cu20 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu20) # clustFA:length

mod_gam_Cu21 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            d13C:length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu21) # season:d13C

mod_gam_Cu22 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu22) # d13C:length

mod_gam_Cu23 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu23) # d15N:length

AIC(mod_gam_Cu0,mod_gam_Cu1,mod_gam_Cu11bis,mod_gam_Cu15,mod_gam_Cu16,
    mod_gam_Cu17,mod_gam_Cu18,mod_gam_Cu19,mod_gam_Cu20,mod_gam_Cu21,
    mod_gam_Cu22,mod_gam_Cu23)

## Step 3

# Begin from n°21

mod_gam_Cu24 <- mgcv::gam(log10(Cu) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            d13C:length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu24) # length

mod_gam_Cu25 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            d13C:length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu25) # sex

mod_gam_Cu26 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            d13C:length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu26) # clustFA

mod_gam_Cu27 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            d13C:length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu27) # season

mod_gam_Cu28 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(new_cluster):length+
                            d13C:length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu28) # sex:length

mod_gam_Cu29 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            d13C:length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu29) # clustFA:length

mod_gam_Cu30 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu30) # d13C:length

mod_gam_Cu31 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            d13C:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu31) # d15N:length

AIC(mod_gam_Cu0,mod_gam_Cu1,mod_gam_Cu11bis,mod_gam_Cu21,mod_gam_Cu24,
    mod_gam_Cu25,mod_gam_Cu26,mod_gam_Cu27,mod_gam_Cu28,mod_gam_Cu29,
    mod_gam_Cu30,mod_gam_Cu31)

## Step 4

# Begin from n°30, remove clustFA (mod 26), clustFA:length (mod 29)

mod_gam_Cu30bis <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                               factor(sex)+
                               factor(season)+
                               factor(sex):length+
                               d15N:length,
                             data = database_Cu, method = "REML")
summary(mod_gam_Cu30bis)

AIC(mod_gam_Cu0,mod_gam_Cu1,mod_gam_Cu11bis,mod_gam_Cu21,mod_gam_Cu30,
    mod_gam_Cu30bis)

mod_gam_Cu32 <- mgcv::gam(log10(Cu) ~ factor(sex)+
                            factor(season)+
                            factor(sex):length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu32) # length

mod_gam_Cu33 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(season)+
                            factor(sex):length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu33) # sex

mod_gam_Cu34 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(sex):length+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu34) # season

mod_gam_Cu35 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(season)+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu35) # sex:length

mod_gam_Cu36 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(season)+
                            factor(sex):length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu36) # d15N:length

AIC(mod_gam_Cu0,mod_gam_Cu1,mod_gam_Cu11bis,mod_gam_Cu21,mod_gam_Cu30bis,
    mod_gam_Cu32,mod_gam_Cu33,mod_gam_Cu34,mod_gam_Cu35,mod_gam_Cu36)

## Step 5

# Begin from n°35, remove sex (mod 33)

mod_gam_Cu35bis <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                               factor(season)+
                               d15N:length,
                             data = database_Cu, method = "REML")
summary(mod_gam_Cu35bis)

AIC(mod_gam_Cu0,mod_gam_Cu1,mod_gam_Cu11bis,mod_gam_Cu21,mod_gam_Cu30bis,
    mod_gam_Cu35,mod_gam_Cu35bis)

mod_gam_Cu37 <- mgcv::gam(log10(Cu) ~ factor(season)+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu37) # length

mod_gam_Cu38 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            d15N:length,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu38) # season

mod_gam_Cu39 <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                            factor(season),
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu39) # d15N:length

AIC(mod_gam_Cu0,mod_gam_Cu1,mod_gam_Cu11bis,mod_gam_Cu21,mod_gam_Cu30bis,
    mod_gam_Cu35bis,mod_gam_Cu37,mod_gam_Cu38,mod_gam_Cu39)

# Best model is 35bis

rm(mod_gam_Cu0,mod_gam_Cu1,mod_gam_Cu2,mod_gam_Cu3,mod_gam_Cu4,
   mod_gam_Cu5,mod_gam_Cu6,mod_gam_Cu7,mod_gam_Cu8,mod_gam_Cu9,
   mod_gam_Cu10,mod_gam_Cu11,mod_gam_Cu12,mod_gam_Cu13,mod_gam_Cu14,
   mod_gam_Cu11bis,mod_gam_Cu15,mod_gam_Cu16,
   mod_gam_Cu17,mod_gam_Cu18,mod_gam_Cu19,mod_gam_Cu20,mod_gam_Cu21,
   mod_gam_Cu22,mod_gam_Cu23,
   mod_gam_Cu24,
   mod_gam_Cu25,mod_gam_Cu26,mod_gam_Cu27,mod_gam_Cu28,mod_gam_Cu29,
   mod_gam_Cu30,mod_gam_Cu31,
   mod_gam_Cu30bis,
   mod_gam_Cu32,mod_gam_Cu33,mod_gam_Cu34,mod_gam_Cu35,mod_gam_Cu36,
   mod_gam_Cu35bis,mod_gam_Cu37,mod_gam_Cu38,mod_gam_Cu39)
rm(database_Cu)


## 11 / Length-stand. Cu ##############################################################################################

database_Cu <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Cu_stand) %>% 
  rename(Cu = Cu_stand)

## Step 1

mod_gam_Cu0 <- mgcv::gam(log10(Cu) ~ 1,
                         data = database_Cu, method = "REML")

mod_gam_Cu1 <- mgcv::gam(log10(Cu) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu1) # Complete model

mod_gam_Cu2 <- mgcv::gam(log10(Cu) ~ factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu2) # sex

mod_gam_Cu3 <- mgcv::gam(log10(Cu) ~ factor(sex)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu3) # new_cluster

mod_gam_Cu4 <- mgcv::gam(log10(Cu) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu4) # d13C

mod_gam_Cu5 <- mgcv::gam(log10(Cu) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu5) # d15N

mod_gam_Cu6 <- mgcv::gam(log10(Cu) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu6) # season

mod_gam_Cu7 <- mgcv::gam(log10(Cu) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d13C,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu7) # season:d15N

mod_gam_Cu8 <- mgcv::gam(log10(Cu) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu8) # season:d13C

AIC(mod_gam_Cu0,mod_gam_Cu1,mod_gam_Cu2,mod_gam_Cu3,mod_gam_Cu4,
    mod_gam_Cu5,mod_gam_Cu6,mod_gam_Cu7,mod_gam_Cu8)

## Step 2

# BEgin from n°7, remove d13C (mod 4)

mod_gam_Cu7bis <- mgcv::gam(log10(Cu) ~ factor(sex)+
                              s(d15N)+
                              factor(new_cluster)+
                              factor(season)+
                              factor(season):d13C,
                            data = database_Cu, method = "REML")
summary(mod_gam_Cu7bis)

AIC(mod_gam_Cu0,mod_gam_Cu1,mod_gam_Cu7,mod_gam_Cu7bis)

mod_gam_Cu9 <- mgcv::gam(log10(Cu) ~ s(d15N)+
                           factor(new_cluster)+
                           factor(season)+
                           factor(season):d13C,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu9) # sex

mod_gam_Cu10 <- mgcv::gam(log10(Cu) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(season):d13C,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu10) # d15N

mod_gam_Cu11 <- mgcv::gam(log10(Cu) ~ factor(sex)+
                            s(d15N)+
                            factor(season)+
                            factor(season):d13C,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu11) # clustFA

mod_gam_Cu12 <- mgcv::gam(log10(Cu) ~ factor(sex)+
                            s(d15N)+
                            factor(new_cluster)+
                            factor(season):d13C,
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu12) # season

mod_gam_Cu13 <- mgcv::gam(log10(Cu) ~ factor(sex)+
                            s(d15N)+
                            factor(new_cluster)+
                            factor(season),
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu13) # season:d13C

AIC(mod_gam_Cu0,mod_gam_Cu1,mod_gam_Cu7bis,mod_gam_Cu9,mod_gam_Cu10,
    mod_gam_Cu11,mod_gam_Cu12,mod_gam_Cu13)

## Step 3

# Begin from n°13

mod_gam_Cu14 <- mgcv::gam(log10(Cu) ~ s(d15N)+
                            factor(new_cluster)+
                            factor(season),
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu14) # sex

mod_gam_Cu15 <- mgcv::gam(log10(Cu) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(season),
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu15) # d15N

mod_gam_Cu16 <- mgcv::gam(log10(Cu) ~ factor(sex)+
                            s(d15N)+
                            factor(season),
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu16) # clustFA

mod_gam_Cu17 <- mgcv::gam(log10(Cu) ~ factor(sex)+
                            s(d15N)+
                            factor(new_cluster),
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu17) # season

AIC(mod_gam_Cu0,mod_gam_Cu1,mod_gam_Cu7bis,mod_gam_Cu13,mod_gam_Cu14,
    mod_gam_Cu15,mod_gam_Cu16,mod_gam_Cu17)

## Step 4

# Begin from n°14, remove clustFA (mod 16)

mod_gam_Cu14bis <- mgcv::gam(log10(Cu) ~ s(d15N)+
                               factor(season),
                             data = database_Cu, method = "REML")
summary(mod_gam_Cu14bis)

AIC(mod_gam_Cu0,mod_gam_Cu1,mod_gam_Cu7bis,mod_gam_Cu13,mod_gam_Cu14,
    mod_gam_Cu14bis)

mod_gam_Cu18 <- mgcv::gam(log10(Cu) ~ factor(season),
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu18) # d15N

mod_gam_Cu19 <- mgcv::gam(log10(Cu) ~ s(d15N),
                          data = database_Cu, method = "REML")
summary(mod_gam_Cu19) # season

AIC(mod_gam_Cu0,mod_gam_Cu1,mod_gam_Cu7bis,mod_gam_Cu13,mod_gam_Cu14bis,
    mod_gam_Cu18,mod_gam_Cu19)

# Best model is 14bis

rm(mod_gam_Cu0,mod_gam_Cu1,mod_gam_Cu2,mod_gam_Cu3,mod_gam_Cu4,
   mod_gam_Cu5,mod_gam_Cu6,mod_gam_Cu7,mod_gam_Cu8,mod_gam_Cu9,
   mod_gam_Cu10,mod_gam_Cu11,mod_gam_Cu12,mod_gam_Cu13,mod_gam_Cu14,
   mod_gam_Cu7bis,mod_gam_Cu15,mod_gam_Cu16,
   mod_gam_Cu17,mod_gam_Cu18,mod_gam_Cu19,mod_gam_Cu14bis)
rm(database_Cu)


## 12 / Cu & length-stand. Cu - Plot response variable vs explaining variables ##############################################################################################

## LOG-TRANSFORMED Cu
database_Cu_ls <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Cu)

# Best model is n°35bis

mod_gam_Cu35bis <- gam::gam(log10(Cu) ~ s(length)+
                              factor(season)+
                              d15N:length,
                            data = database_Cu_ls)
summary(mod_gam_Cu35bis)

mod_gam_Cu35ter <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                               factor(season)+
                               d15N:length,
                             data = database_Cu_ls, method = "REML")
summary(mod_gam_Cu35ter)

## LOG-TRANSFORMED LENGTH-STAND. Cu
database_Cu <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Cu_stand) %>% 
  rename(Cu = Cu_stand)

# Best model is n°14bis
mod_gam_Cu14bis <- gam::gam(log10(Cu) ~ s(d15N)+
                              factor(season),
                            data = database_Cu)
summary(mod_gam_Cu14bis)

mod_gam_Cu14ter <- mgcv::gam(log10(Cu) ~ s(d15N)+
                               factor(season),
                             data = database_Cu, method = "REML")
summary(mod_gam_Cu14ter)

par(mfrow=c(2,2))
plot(mod_gam_Cu35bis,se = TRUE, ylab = "", xlab = "", las = 1)
plot(mod_gam_Cu14bis,se = TRUE, ylab = "", xlab = "", las = 1)

rm(database_Cu_ls,mod_gam_Cu35bis,mod_gam_Cu35ter)
rm(database_Cu,mod_gam_Cu14bis,mod_gam_Cu14ter)


## 13 / Iron ##############################################################################################

database_Fe <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Fe)

## Step 1
mod_gam_Fe0 <- mgcv::gam(log10(Fe) ~ 1,
                         data = database_Fe, method = "REML")

mod_gam_Fe1 <- mgcv::gam(log10(Fe) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Fe, method = "REML")
summary(mod_gam_Fe1) # Complete model

mod_gam_Fe2 <- mgcv::gam(log10(Fe) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Fe, method = "REML")
summary(mod_gam_Fe2) # length

mod_gam_Fe3 <- mgcv::gam(log10(Fe) ~ s(length, k = 3)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Fe, method = "REML")
summary(mod_gam_Fe3) # sex

mod_gam_Fe4 <- mgcv::gam(log10(Fe) ~ s(length, k = 3)+
                           factor(sex)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Fe, method = "REML")
summary(mod_gam_Fe4) # new_cluster

mod_gam_Fe5 <- mgcv::gam(log10(Fe) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Fe, method = "REML")
summary(mod_gam_Fe5) # d13C

mod_gam_Fe6 <- mgcv::gam(log10(Fe) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Fe, method = "REML")
summary(mod_gam_Fe6) # d15N

mod_gam_Fe7 <- mgcv::gam(log10(Fe) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Fe, method = "REML")
summary(mod_gam_Fe7) # season

mod_gam_Fe8 <- mgcv::gam(log10(Fe) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Fe, method = "REML")
summary(mod_gam_Fe8) # Sex:length

mod_gam_Fe9 <- mgcv::gam(log10(Fe) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Fe, method = "REML")
summary(mod_gam_Fe9) # new_cluster:length

mod_gam_Fe10 <- mgcv::gam(log10(Fe) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe10) # season:length

mod_gam_Fe11 <- mgcv::gam(log10(Fe) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe11) # season:d15N

mod_gam_Fe12 <- mgcv::gam(log10(Fe) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            d13C:length+
                            d15N:length,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe12) # season:d13C

mod_gam_Fe13 <- mgcv::gam(log10(Fe) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe13) # d13C:length

mod_gam_Fe14 <- mgcv::gam(log10(Fe) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe14) # d15N:length

AIC(mod_gam_Fe0,mod_gam_Fe1,mod_gam_Fe2,mod_gam_Fe3,mod_gam_Fe4,
    mod_gam_Fe5,mod_gam_Fe6,mod_gam_Fe7,mod_gam_Fe8,mod_gam_Fe9,
    mod_gam_Fe10,mod_gam_Fe11,mod_gam_Fe12,mod_gam_Fe13,mod_gam_Fe14)

## Step 2

# Begin from n°14, remove length (mod 2), sex (mod 3), d13C (mod 5),
# d15N (mod 6)

mod_gam_Fe14bis <- mgcv::gam(log10(Fe) ~ factor(new_cluster)+
                               factor(season)+
                               factor(sex):length+
                               factor(new_cluster):length+
                               factor(season):length+
                               factor(season):d15N+
                               factor(season):d13C+
                               d13C:length,
                             data = database_Fe, method = "REML")
summary(mod_gam_Fe14bis)

AIC(mod_gam_Fe0,mod_gam_Fe1,mod_gam_Fe14,mod_gam_Fe14bis)

mod_gam_Fe15 <- mgcv::gam(log10(Fe) ~ factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe15) # clustFA

mod_gam_Fe16 <- mgcv::gam(log10(Fe) ~ factor(new_cluster)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe16) # season

mod_gam_Fe17 <- mgcv::gam(log10(Fe) ~ factor(new_cluster)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe17) # sex:length

mod_gam_Fe18 <- mgcv::gam(log10(Fe) ~ factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe18) # clustFA:length

mod_gam_Fe19 <- mgcv::gam(log10(Fe) ~ factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe19) # season:length

mod_gam_Fe20 <- mgcv::gam(log10(Fe) ~ factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe20) # season:d15N

mod_gam_Fe21 <- mgcv::gam(log10(Fe) ~ factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            d13C:length,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe21) # season:d13C

mod_gam_Fe22 <- mgcv::gam(log10(Fe) ~ factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe22) # d13C:length

AIC(mod_gam_Fe0,mod_gam_Fe1,mod_gam_Fe14bis,mod_gam_Fe15,mod_gam_Fe16,
    mod_gam_Fe17,mod_gam_Fe18,mod_gam_Fe19,mod_gam_Fe20,mod_gam_Fe21,
    mod_gam_Fe22)

## Step 3

# Begin from n°22, remove length:season (mod 19)

mod_gam_Fe22bis <- mgcv::gam(log10(Fe) ~ factor(new_cluster)+
                               factor(season)+
                               factor(sex):length+
                               factor(new_cluster):length+
                               factor(season):d15N+
                               factor(season):d13C,
                             data = database_Fe, method = "REML")
summary(mod_gam_Fe22bis)

AIC(mod_gam_Fe0,mod_gam_Fe1,mod_gam_Fe14bis,mod_gam_Fe22,mod_gam_Fe22bis)

mod_gam_Fe23 <- mgcv::gam(log10(Fe) ~ factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe23) # clustFA

mod_gam_Fe24 <- mgcv::gam(log10(Fe) ~ factor(new_cluster)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe24) # season

mod_gam_Fe25 <- mgcv::gam(log10(Fe) ~ factor(new_cluster)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe25) # sex:length

mod_gam_Fe26 <- mgcv::gam(log10(Fe) ~ factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe26) # clustFA:length

mod_gam_Fe27 <- mgcv::gam(log10(Fe) ~ factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d13C,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe27) # season:d15N

mod_gam_Fe28 <- mgcv::gam(log10(Fe) ~ factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe28) # season:d13C

mod_gam_Fe29 <- mgcv::gam(log10(Fe) ~ factor(new_cluster)+
                            factor(season)+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Fe, method = "REML")
summary(mod_gam_Fe29) # all with length

AIC(mod_gam_Fe0,mod_gam_Fe1,mod_gam_Fe14bis,mod_gam_Fe22bis,mod_gam_Fe23,
    mod_gam_Fe24,mod_gam_Fe25,mod_gam_Fe26,mod_gam_Fe27,mod_gam_Fe28,
    mod_gam_Fe29)

## Best model is 22bis

rm(mod_gam_Fe0,mod_gam_Fe1,mod_gam_Fe2,mod_gam_Fe3,mod_gam_Fe4,
   mod_gam_Fe5,mod_gam_Fe6,mod_gam_Fe7,mod_gam_Fe8,mod_gam_Fe9,
   mod_gam_Fe10,mod_gam_Fe11,mod_gam_Fe12,mod_gam_Fe13,mod_gam_Fe14,
   mod_gam_Fe14bis,mod_gam_Fe15,mod_gam_Fe16,
   mod_gam_Fe17,mod_gam_Fe18,mod_gam_Fe19,mod_gam_Fe20,mod_gam_Fe21,
   mod_gam_Fe22,
   mod_gam_Fe22bis,mod_gam_Fe23,
   mod_gam_Fe24,mod_gam_Fe25,mod_gam_Fe26,mod_gam_Fe27,mod_gam_Fe28,
   mod_gam_Fe29)
rm(database_Fe)


## 14 / Fe - Plot response variable vs explaining variables ##############################################################################################

## LOG-TRANSFORMED Fe
database_Fe <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Fe)

# Best model is n°35bis

mod_gam_Fe22bis <- gam::gam(log10(Fe) ~ factor(new_cluster)+
                              factor(season)+
                              factor(sex):length+
                              factor(new_cluster):length+
                              factor(season):d15N+
                              factor(season):d13C,
                            data = database_Fe)
summary(mod_gam_Fe22bis)

mod_gam_Fe22ter <- mgcv::gam(log10(Fe) ~ factor(new_cluster)+
                               factor(season)+
                               factor(sex):length+
                               factor(new_cluster):length+
                               factor(season):d15N+
                               factor(season):d13C,
                             data = database_Fe, method = "REML")
summary(mod_gam_Fe22ter)

par(mfrow=c(1,2))
plot(mod_gam_Fe22bis,se = TRUE, ylab = "", xlab = "", las = 1)

rm(database_Fe,mod_gam_Fe22bis,mod_gam_Fe22ter)


## 15 / Mercury ##############################################################################################

database_Hg <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Hg)

## Step 1
mod_gam_Hg0 <- mgcv::gam(log10(Hg) ~ 1,
                          data = database_Hg, method = "REML")

mod_gam_Hg1 <- mgcv::gam(log10(Hg) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg1) # Complete model

mod_gam_Hg2 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg2) # length

mod_gam_Hg3 <- mgcv::gam(log10(Hg) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg3) # sex

mod_gam_Hg4 <- mgcv::gam(log10(Hg) ~ s(length, k = 3)+
                            factor(sex)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg4) # new_cluster

mod_gam_Hg5 <- mgcv::gam(log10(Hg) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg5) # d13C

mod_gam_Hg6 <- mgcv::gam(log10(Hg) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg6) # d15N

mod_gam_Hg7 <- mgcv::gam(log10(Hg) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg7) # season

mod_gam_Hg8 <- mgcv::gam(log10(Hg) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg8) # Sex:length

mod_gam_Hg9 <- mgcv::gam(log10(Hg) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg9) # new_cluster:length

mod_gam_Hg10 <- mgcv::gam(log10(Hg) ~ s(length, k = 3)+
                             factor(sex)+
                             factor(new_cluster)+
                             s(d13C)+
                             s(d15N)+
                             factor(season)+
                             factor(sex):length+
                             factor(new_cluster):length+
                             factor(season):d15N+
                             factor(season):d13C+
                             d13C:length+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg10) # season:length

mod_gam_Hg11 <- mgcv::gam(log10(Hg) ~ s(length, k = 3)+
                             factor(sex)+
                             factor(new_cluster)+
                             s(d13C)+
                             s(d15N)+
                             factor(season)+
                             factor(sex):length+
                             factor(new_cluster):length+
                             factor(season):length+
                             factor(season):d13C+
                             d13C:length+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg11) # season:d15N

mod_gam_Hg12 <- mgcv::gam(log10(Hg) ~ s(length, k = 3)+
                             factor(sex)+
                             factor(new_cluster)+
                             s(d13C)+
                             s(d15N)+
                             factor(season)+
                             factor(sex):length+
                             factor(new_cluster):length+
                             factor(season):length+
                             factor(season):d15N+
                             d13C:length+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg12) # season:d13C

mod_gam_Hg13 <- mgcv::gam(log10(Hg) ~ s(length, k = 3)+
                             factor(sex)+
                             factor(new_cluster)+
                             s(d13C)+
                             s(d15N)+
                             factor(season)+
                             factor(sex):length+
                             factor(new_cluster):length+
                             factor(season):length+
                             factor(season):d15N+
                             factor(season):d13C+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg13) # d13C:length

mod_gam_Hg14 <- mgcv::gam(log10(Hg) ~ s(length, k = 3)+
                             factor(sex)+
                             factor(new_cluster)+
                             s(d13C)+
                             s(d15N)+
                             factor(season)+
                             factor(sex):length+
                             factor(new_cluster):length+
                             factor(season):length+
                             factor(season):d15N+
                             factor(season):d13C+
                             d13C:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg14) # d15N:length

AIC(mod_gam_Hg0,mod_gam_Hg1,mod_gam_Hg2,mod_gam_Hg3,mod_gam_Hg4,
    mod_gam_Hg5,mod_gam_Hg6,mod_gam_Hg7,mod_gam_Hg8,mod_gam_Hg9,
    mod_gam_Hg10,mod_gam_Hg11,mod_gam_Hg12,mod_gam_Hg13,mod_gam_Hg14)

## Step 2

# Begin from n°10, remove length (mod 2), d13C (mod 5), d15N (mod 6),
# season:d13C (mod 12), d13C:length (mod 13)

mod_gam_Hg10bis <- mgcv::gam(log10(Hg) ~ factor(sex)+
                                factor(new_cluster)+
                                factor(season)+
                                factor(sex):length+
                                factor(new_cluster):length+
                                factor(season):d15N+
                                d15N:length,
                              data = database_Hg, method = "REML")
summary(mod_gam_Hg10bis) # season:length

AIC(mod_gam_Hg1,mod_gam_Hg10,mod_gam_Hg10bis)

mod_gam_Hg15 <- mgcv::gam(log10(Hg) ~ factor(new_cluster)+
                             factor(season)+
                             factor(sex):length+
                             factor(new_cluster):length+
                             factor(season):d15N+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg15) # sex

mod_gam_Hg16 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(season)+
                             factor(sex):length+
                             factor(new_cluster):length+
                             factor(season):d15N+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg16) # clustFA

mod_gam_Hg17 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(new_cluster)+
                             factor(sex):length+
                             factor(new_cluster):length+
                             factor(season):d15N+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg17) # season

mod_gam_Hg18 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(new_cluster)+
                             factor(season)+
                             factor(new_cluster):length+
                             factor(season):d15N+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg18) # sex:length

mod_gam_Hg19 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(new_cluster)+
                             factor(season)+
                             factor(sex):length+
                             factor(season):d15N+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg19) # clustFA:length

mod_gam_Hg20 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(new_cluster)+
                             factor(season)+
                             factor(sex):length+
                             factor(new_cluster):length+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg20) # season:d15N

mod_gam_Hg21 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(new_cluster)+
                             factor(season)+
                             factor(sex):length+
                             factor(new_cluster):length+
                             factor(season):d15N,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg21) # d15N:length

AIC(mod_gam_Hg1,mod_gam_Hg10bis,mod_gam_Hg15,mod_gam_Hg16,mod_gam_Hg17,
    mod_gam_Hg18,mod_gam_Hg19,mod_gam_Hg20,mod_gam_Hg21)

## Step 3

# Begin from n°18
mod_gam_Hg22 <- mgcv::gam(log10(Hg) ~ factor(new_cluster)+
                             factor(season)+
                             factor(new_cluster):length+
                             factor(season):d15N+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg22) # sex

mod_gam_Hg23 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(season)+
                             factor(new_cluster):length+
                             factor(season):d15N+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg23) # clustFA

mod_gam_Hg24 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(new_cluster)+
                             factor(new_cluster):length+
                             factor(season):d15N+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg24) # season

mod_gam_Hg25 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(new_cluster)+
                             factor(season)+
                             factor(season):d15N+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg25) # clustFA:length

mod_gam_Hg26 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(new_cluster)+
                             factor(season)+
                             factor(new_cluster):length+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg26) # season:d15N

mod_gam_Hg27 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(new_cluster)+
                             factor(season)+
                             factor(new_cluster):length+
                             factor(season):d15N,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg27) # d15N:length

AIC(mod_gam_Hg1,mod_gam_Hg10bis,mod_gam_Hg18,
    mod_gam_Hg23,mod_gam_Hg24,mod_gam_Hg25,mod_gam_Hg26,mod_gam_Hg27)

## Step 4

# Bengin from n°13

mod_gam_Hg28 <- mgcv::gam(log10(Hg) ~ factor(season)+
                             factor(new_cluster):length+
                             factor(season):d15N+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg28) # sex

mod_gam_Hg29 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(new_cluster):length+
                             factor(season):d15N+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg29) # season

mod_gam_Hg30 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(season)+
                             factor(season):d15N+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg30) # clustFA:length

mod_gam_Hg31 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(season)+
                             factor(new_cluster):length+
                             d15N:length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg31) # season:d15N

mod_gam_Hg32 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(season)+
                             factor(new_cluster):length+
                             factor(season):d15N,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg32) # d15N:length


AIC(mod_gam_Hg1,mod_gam_Hg10bis,mod_gam_Hg18bis,mod_gam_Hg23,
    mod_gam_Hg28,mod_gam_Hg29,mod_gam_Hg30,mod_gam_Hg31,mod_gam_Hg32)

## Step 5

# Begin from n°32, remove season:d15N (mod 31)

mod_gam_Hg32bis <- mgcv::gam(log10(Hg) ~ factor(sex)+
                                factor(season)+
                                factor(new_cluster):length,
                              data = database_Hg, method = "REML")
summary(mod_gam_Hg32bis)

AIC(mod_gam_Hg1,mod_gam_Hg10bis,mod_gam_Hg18bis,mod_gam_Hg23,
    mod_gam_Hg32,mod_gam_Hg32bis)

mod_gam_Hg33 <- mgcv::gam(log10(Hg) ~ factor(season)+
                             factor(new_cluster):length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg33) # sex

mod_gam_Hg34 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(new_cluster):length,
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg34) # season

mod_gam_Hg35 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(season),
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg35) # clustFA:length

AIC(mod_gam_Hg1,mod_gam_Hg10bis,mod_gam_Hg18bis,mod_gam_Hg23,
    mod_gam_Hg32bis,mod_gam_Hg33,mod_gam_Hg34,mod_gam_Hg35)

## Best model is n°34

rm(mod_gam_Hg0,mod_gam_Hg1,mod_gam_Hg2,mod_gam_Hg3,mod_gam_Hg4,
   mod_gam_Hg5,mod_gam_Hg6,mod_gam_Hg7,mod_gam_Hg8,mod_gam_Hg9,
   mod_gam_Hg10,mod_gam_Hg11,mod_gam_Hg12,mod_gam_Hg13,mod_gam_Hg14,
   mod_gam_Hg15,mod_gam_Hg16,mod_gam_Hg17,mod_gam_Hg18,mod_gam_Hg19,
   mod_gam_Hg20,mod_gam_Hg21,mod_gam_Hg22,mod_gam_Hg23,mod_gam_Hg24,
   mod_gam_Hg25,mod_gam_Hg26,mod_gam_Hg27,mod_gam_Hg28,mod_gam_Hg29,
   mod_gam_Hg30,mod_gam_Hg31,mod_gam_Hg32,mod_gam_Hg33,mod_gam_Hg34,
   mod_gam_Hg35,
   mod_gam_Hg10bis,mod_gam_Hg18bis,mod_gam_Hg32bis)
rm(database_Hg)


## 16 / Length-stand. Hg ##############################################################################################

database_Hg <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Hg_stand) %>% 
  rename(Hg = Hg_stand)

## Step 1

mod_gam_Hg0 <- mgcv::gam(log10(Hg) ~ 1,
                          data = database_Hg, method = "REML")

mod_gam_Hg1 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg1) # Complete model

mod_gam_Hg2 <- mgcv::gam(log10(Hg) ~ factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg2) # sex

mod_gam_Hg3 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg3) # new_cluster

mod_gam_Hg4 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d15N)+
                            factor(season)+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg4) # d13C

mod_gam_Hg5 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            factor(season)+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg5) # d15N

mod_gam_Hg6 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg6) # season

mod_gam_Hg7 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(season):d13C,
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg7) # season:d15N

mod_gam_Hg8 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(season):d15N,
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg8) # season:d13C

AIC(mod_gam_Hg0,mod_gam_Hg1,mod_gam_Hg2,mod_gam_Hg3,mod_gam_Hg4,
    mod_gam_Hg5,mod_gam_Hg6,mod_gam_Hg7,mod_gam_Hg8)

## Step 2

# Begin from n°8, remove d13C (mod 54), d15N (mod 5), season:d15N (mod 7)

mod_gam_Hg8bis <- mgcv::gam(log10(Hg) ~ factor(sex)+
                               factor(new_cluster)+
                               factor(season),
                             data = database_Hg, method = "REML")
summary(mod_gam_Hg8bis)

AIC(mod_gam_Hg0,mod_gam_Hg1,mod_gam_Hg8,mod_gam_Hg8bis)

mod_gam_Hg9 <- mgcv::gam(log10(Hg) ~ factor(new_cluster)+
                            factor(season),
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg9) # sex

mod_gam_Hg10 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(season),
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg10) # clustFA

mod_gam_Hg11 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                             factor(new_cluster),
                           data = database_Hg, method = "REML")
summary(mod_gam_Hg11) # season

AIC(mod_gam_Hg0,mod_gam_Hg1,mod_gam_Hg8bis,mod_gam_Hg9,mod_gam_Hg10,
    mod_gam_Hg11)

## Best model is 11

rm(mod_gam_Hg0,mod_gam_Hg1,mod_gam_Hg2,mod_gam_Hg3,mod_gam_Hg4,
   mod_gam_Hg5,mod_gam_Hg6,mod_gam_Hg7,mod_gam_Hg8,mod_gam_Hg8bis,
   mod_gam_Hg9,mod_gam_Hg10,mod_gam_Hg11)

rm(database_Hg)


## 17 / Hg & length-stand. Hg - Plot response variable vs explaining variables ##############################################################################################

## LOG-TRANSFORMED Hg
database_Hg_ls <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Hg)

# Best model is n°34

mod_gam_Hg34 <- gam::gam(log10(Hg) ~ factor(sex)+
                            factor(new_cluster):length,
                          data = database_Hg_ls)
summary(mod_gam_Hg34)

mod_gam_Hg34bis <- mgcv::gam(log10(Hg) ~ factor(sex)+
                                factor(new_cluster):length,
                              data = database_Hg_ls, method = "REML")
summary(mod_gam_Hg34bis)

## LOG-TRANSFORMED LENGTH-STAND. Hg
database_Hg <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Hg_stand) %>% 
  rename(Hg = Hg_stand)

# Best model is n°11
mod_gam_Hg11 <- gam::gam(log10(Hg) ~ factor(sex)+
                            factor(new_cluster),
                          data = database_Hg)
summary(mod_gam_Hg11)

mod_gam_Hg11bis <- mgcv::gam(log10(Hg) ~ factor(sex)+
                                factor(new_cluster),
                              data = database_Hg, method = "REML")
summary(mod_gam_Hg11bis)

par(mfrow=c(2,2))
plot(mod_gam_Hg34,se = TRUE, ylab = "", xlab = "", las = 1)
plot(mod_gam_Hg34,se = TRUE, ylab = "", xlab = "", las = 1)
plot(mod_gam_Hg11,se = TRUE, ylab = "", xlab = "", las = 1)

rm(database_Hg_ls,mod_gam_Hg34,mod_gam_Hg34bis)
rm(database_Hg,mod_gam_Hg11,mod_gam_Hg11bis)


## 18 / Manganese ##############################################################################################

database_Mn <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Mn)

## Step 1
mod_gam_Mn0 <- mgcv::gam(log10(Mn) ~ 1,
                         data = database_Mn, method = "REML")

mod_gam_Mn1 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           #d13C:length+
                           d15N:length,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn1) # Complete model

mod_gam_Mn2 <- mgcv::gam(log10(Mn) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           #d13C:length+
                           d15N:length,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn2) # length

mod_gam_Mn3 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           #d13C:length+
                           d15N:length,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn3) # sex

mod_gam_Mn4 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                           factor(sex)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           #d13C:length+
                           d15N:length,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn4) # new_cluster

mod_gam_Mn5 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           #d13C:length+
                           d15N:length,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn5) # d13C

mod_gam_Mn6 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           #d13C:length+
                           d15N:length,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn6) # d15N

mod_gam_Mn7 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           #d13C:length+
                           d15N:length,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn7) # season

mod_gam_Mn8 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           #d13C:length+
                           d15N:length,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn8) # Sex:length

mod_gam_Mn9 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           #d13C:length+
                           d15N:length,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn9) # new_cluster:length

mod_gam_Mn10 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            #d13C:length+
                            d15N:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn10) # season:length

mod_gam_Mn11 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d13C+
                            #d13C:length+
                            d15N:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn11) # season:d15N

mod_gam_Mn12 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            #d13C:length+
                            d15N:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn12) # season:d13C

mod_gam_Mn14 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          #d13C:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn14) # d15N:length

AIC(mod_gam_Mn0,mod_gam_Mn1,mod_gam_Mn2,mod_gam_Mn3,mod_gam_Mn4,
    mod_gam_Mn5,mod_gam_Mn6,mod_gam_Mn7,mod_gam_Mn8,mod_gam_Mn9,
    mod_gam_Mn10,mod_gam_Mn11,mod_gam_Mn12,mod_gam_Mn14) #mod_gam_Mn13,

## Step 2

# BEgin from n°10, remove d13C (mod 5), d15N (mod 6)

mod_gam_Mn10bis <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                               factor(sex)+
                               factor(new_cluster)+
                               factor(season)+
                               factor(sex):length+
                               factor(new_cluster):length+
                               factor(season):d15N+
                               factor(season):d13C+
                               #d13C:length+
                               d15N:length,
                             data = database_Mn, method = "REML")
summary(mod_gam_Mn10bis)

AIC(mod_gam_Mn0,mod_gam_Mn1,mod_gam_Mn10,mod_gam_Mn10bis)

mod_gam_Mn15 <- mgcv::gam(log10(Mn) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            #d13C:length+
                            d15N:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn15) # length

mod_gam_Mn16 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            #d13C:length+
                            d15N:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn16) # sex

mod_gam_Mn17 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            #d13C:length+
                            d15N:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn17) # clustFA

mod_gam_Mn18 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            #d13C:length+
                            d15N:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn18) # season

mod_gam_Mn19 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            #d13C:length+
                            d15N:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn19) # sex:length

mod_gam_Mn20 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            #d13C:length+
                            d15N:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn20) # clustFA:length

mod_gam_Mn21 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d13C+
                            #d13C:length+
                            d15N:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn21) # season:d15N

mod_gam_Mn22 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            #d13C:length+
                            d15N:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn22) # season:d13C

mod_gam_Mn24 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn24) # d15N:length

AIC(mod_gam_Mn0,mod_gam_Mn1,mod_gam_Mn10bis,mod_gam_Mn15,mod_gam_Mn16,
    mod_gam_Mn17,mod_gam_Mn18,mod_gam_Mn19,mod_gam_Mn20,mod_gam_Mn21,
    mod_gam_Mn22,mod_gam_Mn24)

## Step 3

# Begin from n°20

mod_gam_Mn25 <- mgcv::gam(log10(Mn) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            #d13C:length+
                            d15N:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn25) # length

mod_gam_Mn26 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            #d13C:length+
                            d15N:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn26) # sex

mod_gam_Mn27 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            #d13C:length+
                            d15N:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn27) # clustFA

mod_gam_Mn28 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(sex):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            #d13C:length+
                            d15N:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn28) # season

mod_gam_Mn29 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(season):d15N+
                            factor(season):d13C+
                            #d13C:length+
                            d15N:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn29) # sex:length

mod_gam_Mn30 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):d13C+
                            #d13C:length+
                            d15N:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn30) # season:d15N

mod_gam_Mn31 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):d15N+
                            #d13C:length+
                            d15N:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn31) # season:d13C

mod_gam_Mn33 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          #d13C:length,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn33) # d15N:length

AIC(mod_gam_Mn0,mod_gam_Mn1,mod_gam_Mn10bis,mod_gam_Mn20,mod_gam_Mn25,
    mod_gam_Mn26,mod_gam_Mn27,mod_gam_Mn28,mod_gam_Mn29,mod_gam_Mn30,
    mod_gam_Mn31,mod_gam_Mn33)

## Best model is n°20

rm(mod_gam_Mn0,mod_gam_Mn1,mod_gam_Mn2,mod_gam_Mn3,mod_gam_Mn4,
   mod_gam_Mn5,mod_gam_Mn6,mod_gam_Mn7,mod_gam_Mn8,mod_gam_Mn9,
   mod_gam_Mn10,mod_gam_Mn11,mod_gam_Mn12,mod_gam_Mn13,mod_gam_Mn14,
   mod_gam_Mn10bis,mod_gam_Mn15,mod_gam_Mn16,
   mod_gam_Mn17,mod_gam_Mn18,mod_gam_Mn19,mod_gam_Mn20,mod_gam_Mn21,
   mod_gam_Mn22,mod_gam_Mn23,mod_gam_Mn24,
   mod_gam_Mn25,
   mod_gam_Mn26,mod_gam_Mn27,mod_gam_Mn28,mod_gam_Mn29,mod_gam_Mn30,
   mod_gam_Mn31,mod_gam_Mn32,mod_gam_Mn33)
rm(database_Mn)


## 19 / Length-stand. Mn ##############################################################################################

database_Mn <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Mn_stand) %>% 
  rename(Mn = Mn_stand)

## Step 1

mod_gam_Mn0 <- mgcv::gam(log10(Mn) ~ 1,
                         data = database_Mn, method = "REML")

mod_gam_Mn1 <- mgcv::gam(log10(Mn) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn1) # Complete model

mod_gam_Mn2 <- mgcv::gam(log10(Mn) ~ factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn2) # sex

mod_gam_Mn3 <- mgcv::gam(log10(Mn) ~ factor(sex)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn3) # new_cluster

mod_gam_Mn4 <- mgcv::gam(log10(Mn) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn4) # d13C

mod_gam_Mn5 <- mgcv::gam(log10(Mn) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn5) # d15N

mod_gam_Mn6 <- mgcv::gam(log10(Mn) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn6) # season

mod_gam_Mn7 <- mgcv::gam(log10(Mn) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d13C,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn7) # season:d15N

mod_gam_Mn8 <- mgcv::gam(log10(Mn) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn8) # season:d13C

AIC(mod_gam_Mn0,mod_gam_Mn1,mod_gam_Mn2,mod_gam_Mn3,mod_gam_Mn4,
    mod_gam_Mn5,mod_gam_Mn6,mod_gam_Mn7,mod_gam_Mn8)

## Step 2

# Begin from n°2, remove d13C (mod 4)

mod_gam_Mn2bis <- mgcv::gam(log10(Mn) ~ factor(new_cluster)+
                              s(d15N)+
                              factor(season)+
                              factor(season):d15N+
                              factor(season):d13C,
                            data = database_Mn, method = "REML")
summary(mod_gam_Mn2bis)

AIC(mod_gam_Mn0,mod_gam_Mn1,mod_gam_Mn2,mod_gam_Mn2bis)

mod_gam_Mn9 <- mgcv::gam(log10(Mn) ~ s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn9) # clustFA

mod_gam_Mn10 <- mgcv::gam(log10(Mn) ~ factor(new_cluster)+
                            factor(season)+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn10) # d15N

mod_gam_Mn11 <- mgcv::gam(log10(Mn) ~ factor(new_cluster)+
                            s(d15N)+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn11) # season

mod_gam_Mn12 <- mgcv::gam(log10(Mn) ~ factor(new_cluster)+
                            s(d15N)+
                            factor(season)+
                            factor(season):d13C,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn12) # season:d15N

mod_gam_Mn13 <- mgcv::gam(log10(Mn) ~ factor(new_cluster)+
                            s(d15N)+
                            factor(season)+
                            factor(season):d15N,
                          data = database_Mn, method = "REML")
summary(mod_gam_Mn13) # season:d13C

AIC(mod_gam_Mn0,mod_gam_Mn1,mod_gam_Mn2bis,mod_gam_Mn9,mod_gam_Mn10,
    mod_gam_Mn11,mod_gam_Mn12,mod_gam_Mn13)

## Best model is 2bis

rm(mod_gam_Mn0,mod_gam_Mn1,mod_gam_Mn2,mod_gam_Mn3,mod_gam_Mn4,
   mod_gam_Mn5,mod_gam_Mn6,mod_gam_Mn7,mod_gam_Mn8,mod_gam_Mn9,
   mod_gam_Mn10,mod_gam_Mn11,mod_gam_Mn12,mod_gam_Mn13,mod_gam_Mn2bis)
rm(database_Mn)


## 20 / Mn & length-stand. Mn - Plot response variable vs explaining variables ##############################################################################################

## LOG-TRANSFORMED Mn
database_Mn_ls <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Mn)

# Best model is n°20

mod_gam_Mn20 <- gam::gam(log10(Mn) ~ s(length, df = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           factor(season)+
                           factor(sex):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Mn_ls)
summary(mod_gam_Mn20)

mod_gam_Mn20bis <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                               factor(sex)+
                               factor(new_cluster)+
                               factor(season)+
                               factor(sex):length+
                               factor(season):d15N+
                               factor(season):d13C+
                               d13C:length+
                               d15N:length,
                             data = database_Mn_ls, method = "REML")
summary(mod_gam_Mn20bis)

## LOG-TRANSFORMED LENGTH-STAND. Mn
database_Mn <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Mn_stand) %>% 
  rename(Mn = Mn_stand)

# Best model is n°2bis
mod_gam_Mn2bis <- gam::gam(log10(Mn) ~ s(d15N)+
                             factor(new_cluster)+
                             factor(season)+
                             factor(season):d15N+
                             factor(season):d13C,
                           data = database_Mn)
summary(mod_gam_Mn2bis)

mod_gam_Mn2ter <- mgcv::gam(log10(Mn) ~ factor(new_cluster)+
                              s(d15N)+
                              factor(season)+
                              factor(season):d15N+
                              factor(season):d13C,
                            data = database_Mn, method = "REML")
summary(mod_gam_Mn2ter)

par(mfrow=c(4,2))
plot(mod_gam_Mn20,se = TRUE, ylab = "", xlab = "", las = 1)
plot(mod_gam_Mn2bis,se = TRUE, ylab = "", xlab = "", las = 1)

rm(database_Mn_ls,mod_gam_Mn20,mod_gam_Mn20bis)
rm(database_Mn,mod_gam_Mn2bis,mod_gam_Mn2ter)


## 21 / Lead ##############################################################################################

database_Pb <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Pb)

## Step 1
mod_gam_Pb0 <- mgcv::gam(log10(Pb) ~ 1,
                         data = database_Pb, method = "REML")

mod_gam_Pb1 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Pb, method = "REML")
summary(mod_gam_Pb1) # Complete model

mod_gam_Pb2 <- mgcv::gam(log10(Pb) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Pb, method = "REML")
summary(mod_gam_Pb2) # length

mod_gam_Pb3 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Pb, method = "REML")
summary(mod_gam_Pb3) # sex

mod_gam_Pb4 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                           factor(sex)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Pb, method = "REML")
summary(mod_gam_Pb4) # new_cluster

mod_gam_Pb5 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Pb, method = "REML")
summary(mod_gam_Pb5) # d13C

mod_gam_Pb6 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Pb, method = "REML")
summary(mod_gam_Pb6) # d15N

mod_gam_Pb7 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Pb, method = "REML")
summary(mod_gam_Pb7) # season

mod_gam_Pb8 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Pb, method = "REML")
summary(mod_gam_Pb8) # Sex:length

mod_gam_Pb9 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Pb, method = "REML")
summary(mod_gam_Pb9) # new_cluster:length

mod_gam_Pb10 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb10) # season:length

mod_gam_Pb11 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb11) # season:d15N

mod_gam_Pb12 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            d13C:length+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb12) # season:d13C

mod_gam_Pb13 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb13) # d13C:length

mod_gam_Pb14 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb14) # d15N:length

AIC(mod_gam_Pb0,mod_gam_Pb1,mod_gam_Pb2,mod_gam_Pb3,mod_gam_Pb4,
    mod_gam_Pb5,mod_gam_Pb6,mod_gam_Pb7,mod_gam_Pb8,mod_gam_Pb9,
    mod_gam_Pb10,mod_gam_Pb11,mod_gam_Pb12,mod_gam_Pb13,mod_gam_Pb14)

## Step 2

# Begin from n°8, remove sex (mod 3)

mod_gam_Pb8bis <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                              factor(new_cluster)+
                              s(d13C)+
                              s(d15N)+
                              factor(season)+
                              factor(new_cluster):length+
                              factor(season):length+
                              factor(season):d15N+
                              factor(season):d13C+
                              d13C:length+
                              d15N:length,
                            data = database_Pb, method = "REML")
summary(mod_gam_Pb8bis)

AIC(mod_gam_Pb0,mod_gam_Pb1,mod_gam_Pb8,mod_gam_Pb8bis)

mod_gam_Pb15 <- mgcv::gam(log10(Pb) ~ factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb15) # length

mod_gam_Pb16 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb16) # clustFA

mod_gam_Pb17 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d15N)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb17) # d13C

mod_gam_Pb18 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb18) # d15N

mod_gam_Pb19 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb19) # season

mod_gam_Pb20 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb20) # clustFA:length

mod_gam_Pb21 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb21) # season:length

mod_gam_Pb22 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb22) # season:d15N

mod_gam_Pb23 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            d13C:length+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb23) # season:d13C

mod_gam_Pb24 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb24) # d13C:length

mod_gam_Pb25 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb25) # d15N:length

AIC(mod_gam_Pb0,mod_gam_Pb1,mod_gam_Pb8bis,mod_gam_Pb15,mod_gam_Pb16,
    mod_gam_Pb17,mod_gam_Pb18,mod_gam_Pb19,mod_gam_Pb20,mod_gam_Pb21,
    mod_gam_Pb22,mod_gam_Pb23,mod_gam_Pb24,mod_gam_Pb25)

## Step 3

# Begin from n°16, remove d13C (mod 17), d15N (mod 18),
# clustFA:length (mod 20)

mod_gam_Pb16bis <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                               factor(season)+
                               factor(new_cluster):length+
                               factor(season):length+
                               factor(season):d15N+
                               factor(season):d13C+
                               d13C:length+
                               d15N:length,
                             data = database_Pb, method = "REML")
summary(mod_gam_Pb16bis)

AIC(mod_gam_Pb0,mod_gam_Pb1,mod_gam_Pb8bis,mod_gam_Pb16,mod_gam_Pb16bis)

mod_gam_Pb26 <- mgcv::gam(log10(Pb) ~ factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb26) # length

mod_gam_Pb27 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb27) # season

mod_gam_Pb28 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(season)+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb28) # clustFA:length

mod_gam_Pb29 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb29) # season:length

mod_gam_Pb30 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb30) # season:d15N

mod_gam_Pb31 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            d13C:length+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb31) # season:d13C

mod_gam_Pb32 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb32) # d13C:length

mod_gam_Pb33 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb33) # d15N:length

AIC(mod_gam_Pb0,mod_gam_Pb1,mod_gam_Pb8bis,mod_gam_Pb16bis,mod_gam_Pb26,
    mod_gam_Pb27,mod_gam_Pb28,mod_gam_Pb29,mod_gam_Pb30,mod_gam_Pb31,
    mod_gam_Pb32,mod_gam_Pb33)

## Step 4

# Begin from n°32

mod_gam_Pb34 <- mgcv::gam(log10(Pb) ~ factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb34) # length

mod_gam_Pb35 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb35) # season

mod_gam_Pb36 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(season)+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb36) # clustFA:length

mod_gam_Pb37 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb37) # season:length

mod_gam_Pb38 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb38) # season:d15N

mod_gam_Pb39 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb39) # season:d13C

mod_gam_Pb40 <- mgcv::gam(log10(Pb) ~ s(length, k = 3)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb40) # d15N:length

AIC(mod_gam_Pb0,mod_gam_Pb1,mod_gam_Pb8bis,mod_gam_Pb16bis,mod_gam_Pb32,
    mod_gam_Pb34,mod_gam_Pb35,mod_gam_Pb36,mod_gam_Pb37,mod_gam_Pb38,
    mod_gam_Pb39,mod_gam_Pb40)

## Step 5

# Begin from n°35, remove length (mod 34)

mod_gam_Pb35bis <- mgcv::gam(log10(Pb) ~ factor(new_cluster):length+
                               factor(season):length+
                               factor(season):d15N+
                               factor(season):d13C+
                               d15N:length,
                             data = database_Pb, method = "REML")
summary(mod_gam_Pb35bis)

AIC(mod_gam_Pb0,mod_gam_Pb1,mod_gam_Pb8bis,mod_gam_Pb16bis,mod_gam_Pb32,
    mod_gam_Pb35,mod_gam_Pb35bis)

mod_gam_Pb41 <- mgcv::gam(log10(Pb) ~ factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb41) # clustFA:length

mod_gam_Pb42 <- mgcv::gam(log10(Pb) ~ factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb42) #season:length

mod_gam_Pb43 <- mgcv::gam(log10(Pb) ~ factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb43) # season:d15N

mod_gam_Pb44 <- mgcv::gam(log10(Pb) ~ factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb44) # season:d13C

mod_gam_Pb45 <- mgcv::gam(log10(Pb) ~ factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb45) # d15N:length

AIC(mod_gam_Pb0,mod_gam_Pb1,mod_gam_Pb8bis,mod_gam_Pb16bis,mod_gam_Pb32,
    mod_gam_Pb35bis,mod_gam_Pb41,mod_gam_Pb42,mod_gam_Pb43,mod_gam_Pb44,
    mod_gam_Pb45)

## Step 6

# Begin from n°42

mod_gam_Pb46 <- mgcv::gam(log10(Pb) ~ factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb46) # clustFA:length

mod_gam_Pb47 <- mgcv::gam(log10(Pb) ~ factor(new_cluster):length+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb47) # season:d15N

mod_gam_Pb48 <- mgcv::gam(log10(Pb) ~ factor(new_cluster):length+
                            factor(season):d15N+
                            d15N:length,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb48) # season:d13C

mod_gam_Pb49 <- mgcv::gam(log10(Pb) ~ factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb49) # d15N:length

AIC(mod_gam_Pb0,mod_gam_Pb1,mod_gam_Pb8bis,mod_gam_Pb16bis,mod_gam_Pb32,
    mod_gam_Pb35bis,mod_gam_Pb42,mod_gam_Pb46,mod_gam_Pb47,mod_gam_Pb48,
    mod_gam_Pb49)

## Step 7

# Begin from n°49

mod_gam_Pb50 <- mgcv::gam(log10(Pb) ~ factor(season):d15N+
                            factor(season):d13C,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb50) # clustFA:length

mod_gam_Pb51 <- mgcv::gam(log10(Pb) ~ factor(new_cluster):length+
                            factor(season):d13C,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb51)  # season:d15N

mod_gam_Pb52 <- mgcv::gam(log10(Pb) ~ factor(new_cluster):length+
                            factor(season):d15N,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb52) # season:d13C

AIC(mod_gam_Pb0,mod_gam_Pb1,mod_gam_Pb8bis,mod_gam_Pb16bis,mod_gam_Pb32,
    mod_gam_Pb35bis,mod_gam_Pb42,mod_gam_Pb49,mod_gam_Pb50,mod_gam_Pb51,
    mod_gam_Pb52)

## Best model is 35bis

rm(mod_gam_Pb0,mod_gam_Pb1,mod_gam_Pb2,mod_gam_Pb3,mod_gam_Pb4,
   mod_gam_Pb5,mod_gam_Pb6,mod_gam_Pb7,mod_gam_Pb8,mod_gam_Pb9,
   mod_gam_Pb10,mod_gam_Pb11,mod_gam_Pb12,mod_gam_Pb13,mod_gam_Pb14,
   mod_gam_Pb8bis,mod_gam_Pb15,mod_gam_Pb16,
   mod_gam_Pb17,mod_gam_Pb18,mod_gam_Pb19,mod_gam_Pb20,mod_gam_Pb21,
   mod_gam_Pb22,mod_gam_Pb23,mod_gam_Pb24,mod_gam_Pb25,
   mod_gam_Pb16bis,mod_gam_Pb26,
   mod_gam_Pb27,mod_gam_Pb28,mod_gam_Pb29,mod_gam_Pb30,mod_gam_Pb31,
   mod_gam_Pb32,mod_gam_Pb33,
   mod_gam_Pb34,mod_gam_Pb35,mod_gam_Pb36,mod_gam_Pb37,mod_gam_Pb38,
   mod_gam_Pb39,mod_gam_Pb40,
   mod_gam_Pb35bis,mod_gam_Pb41,mod_gam_Pb42,mod_gam_Pb43,mod_gam_Pb44,
   mod_gam_Pb45,
   mod_gam_Pb46,mod_gam_Pb47,mod_gam_Pb48,mod_gam_Pb49,
   mod_gam_Pb50,mod_gam_Pb51,mod_gam_Pb52)
rm(database_Pb)


## 22 / Length-stand. Pb ##############################################################################################

database_Pb <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Pb_stand) %>% 
  rename(Pb = Pb_stand)

## Step 1

mod_gam_Pb0 <- mgcv::gam(log10(Pb) ~ 1,
                         data = database_Pb, method = "REML")

mod_gam_Pb1 <- mgcv::gam(log10(Pb) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Pb, method = "REML")
summary(mod_gam_Pb1) # Complete model

mod_gam_Pb2 <- mgcv::gam(log10(Pb) ~ factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Pb, method = "REML")
summary(mod_gam_Pb2) # sex

mod_gam_Pb3 <- mgcv::gam(log10(Pb) ~ factor(sex)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Pb, method = "REML")
summary(mod_gam_Pb3) # new_cluster

mod_gam_Pb4 <- mgcv::gam(log10(Pb) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Pb, method = "REML")
summary(mod_gam_Pb4) # d13C

mod_gam_Pb5 <- mgcv::gam(log10(Pb) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Pb, method = "REML")
summary(mod_gam_Pb5) # d15N

mod_gam_Pb6 <- mgcv::gam(log10(Pb) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Pb, method = "REML")
summary(mod_gam_Pb6) # season

mod_gam_Pb7 <- mgcv::gam(log10(Pb) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d13C,
                         data = database_Pb, method = "REML")
summary(mod_gam_Pb7) # season:d15N

mod_gam_Pb8 <- mgcv::gam(log10(Pb) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N,
                         data = database_Pb, method = "REML")
summary(mod_gam_Pb8) # season:d13C

AIC(mod_gam_Pb0,mod_gam_Pb1,mod_gam_Pb2,mod_gam_Pb3,mod_gam_Pb4,
    mod_gam_Pb5,mod_gam_Pb6,mod_gam_Pb7,mod_gam_Pb8)

## Step 2

# Begin from n°6, remove sex (mod 2)

mod_gam_Pb6bis <- mgcv::gam(log10(Pb) ~ factor(new_cluster)+
                              s(d13C)+
                              s(d15N)+
                              factor(season):d15N+
                              factor(season):d13C,
                            data = database_Pb, method = "REML")
summary(mod_gam_Pb6bis)

AIC(mod_gam_Pb0,mod_gam_Pb1,mod_gam_Pb6,mod_gam_Pb6bis)

mod_gam_Pb9 <- mgcv::gam(log10(Pb) ~ s(d13C)+
                           s(d15N)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Pb, method = "REML")
summary(mod_gam_Pb9) # clustFA

mod_gam_Pb10 <- mgcv::gam(log10(Pb) ~ factor(new_cluster)+
                            s(d15N)+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb10) # d13C

mod_gam_Pb11 <- mgcv::gam(log10(Pb) ~ factor(new_cluster)+
                            s(d13C)+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb11) # d15N

mod_gam_Pb12 <- mgcv::gam(log10(Pb) ~ factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season):d13C,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb12) # season:d15N

mod_gam_Pb13 <- mgcv::gam(log10(Pb) ~ factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season):d15N,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb13) # season:d13C

AIC(mod_gam_Pb0,mod_gam_Pb1,mod_gam_Pb6bis,mod_gam_Pb9,mod_gam_Pb10,
    mod_gam_Pb11,mod_gam_Pb12,mod_gam_Pb13)

## Step 3

# Begin from 6bis, remove d13C (mod 10), d15N (mod 11)

mod_gam_Pb6ter <- mgcv::gam(log10(Pb) ~ factor(new_cluster)+
                              factor(season):d15N+
                              factor(season):d13C,
                            data = database_Pb, method = "REML")
summary(mod_gam_Pb6ter)

AIC(mod_gam_Pb0,mod_gam_Pb1,mod_gam_Pb6bis,mod_gam_Pb6ter)


mod_gam_Pb14 <- mgcv::gam(log10(Pb) ~ factor(season):d15N+
                            factor(season):d13C,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb14) # clustFA

mod_gam_Pb15 <- mgcv::gam(log10(Pb) ~ factor(new_cluster)+
                            factor(season):d13C,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb15) # season:d15N

mod_gam_Pb16 <- mgcv::gam(log10(Pb) ~ factor(new_cluster)+
                            factor(season):d15N,
                          data = database_Pb, method = "REML")
summary(mod_gam_Pb16) # season:d13C

AIC(mod_gam_Pb0,mod_gam_Pb1,mod_gam_Pb6ter,mod_gam_Pb14,mod_gam_Pb15,
    mod_gam_Pb16)

## Best model is 6ter

rm(mod_gam_Pb0,mod_gam_Pb1,mod_gam_Pb2,mod_gam_Pb3,mod_gam_Pb4,
   mod_gam_Pb5,mod_gam_Pb6,mod_gam_Pb7,mod_gam_Pb8,
   mod_gam_Pb6bis,mod_gam_Pb9,mod_gam_Pb10,
   mod_gam_Pb11,mod_gam_Pb12,mod_gam_Pb13,
   mod_gam_Pb14,mod_gam_Pb15,
   mod_gam_Pb16)
rm(database_Pb)


## 23 / Selenium ##############################################################################################

database_Se <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Se)

## Step 1
mod_gam_Se0 <- mgcv::gam(log10(Se) ~ 1,
                         data = database_Se, method = "REML")

mod_gam_Se1 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           #d13C:length+
                           d15N:length,
                         data = database_Se, method = "REML")
summary(mod_gam_Se1) # Complete model

mod_gam_Se2 <- mgcv::gam(log10(Se) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           #d13C:length+
                           d15N:length,
                         data = database_Se, method = "REML")
summary(mod_gam_Se2) # length

mod_gam_Se3 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           #d13C:length+
                           d15N:length,
                         data = database_Se, method = "REML")
summary(mod_gam_Se3) # sex

mod_gam_Se4 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                           factor(sex)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           #d13C:length+
                           d15N:length,
                         data = database_Se, method = "REML")
summary(mod_gam_Se4) # new_cluster

mod_gam_Se5 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           #d13C:length+
                           d15N:length,
                         data = database_Se, method = "REML")
summary(mod_gam_Se5) # d13C

mod_gam_Se6 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           #d13C:length+
                           d15N:length,
                         data = database_Se, method = "REML")
summary(mod_gam_Se6) # d15N

mod_gam_Se7 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           #d13C:length+
                           d15N:length,
                         data = database_Se, method = "REML")
summary(mod_gam_Se7) # season

mod_gam_Se8 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           #d13C:length+
                           d15N:length,
                         data = database_Se, method = "REML")
summary(mod_gam_Se8) # Sex:length

mod_gam_Se9 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           #d13C:length+
                           d15N:length,
                         data = database_Se, method = "REML")
summary(mod_gam_Se9) # new_cluster:length

mod_gam_Se10 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            #d13C:length+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se10) # season:length

mod_gam_Se11 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d13C+
                            #d13C:length+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se11) # season:d15N

mod_gam_Se12 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            #d13C:length+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se12) # season:d13C

mod_gam_Se14 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          #d13C:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se14) # d15N:length

AIC(mod_gam_Se0,mod_gam_Se1,mod_gam_Se2,mod_gam_Se3,mod_gam_Se4,
    mod_gam_Se5,mod_gam_Se6,mod_gam_Se7,mod_gam_Se8,mod_gam_Se9,
    mod_gam_Se10,mod_gam_Se11,mod_gam_Se12,mod_gam_Se14)

## Step 2

# Begin from n°10

mod_gam_Se15 <- mgcv::gam(log10(Se) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se15) # length

mod_gam_Se16 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se16) # sex

mod_gam_Se17 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                            factor(sex)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se17) # clustFA

mod_gam_Se18 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se18) # d13C

mod_gam_Se19 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se19) # d15N

mod_gam_Se20 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se20) # season

mod_gam_Se21 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se21) # sex:length

mod_gam_Se22 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se22) # clustFA:length

mod_gam_Se23 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se23) # season:d15N

mod_gam_Se24 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se24) # season:d13C

mod_gam_Se26 <- mgcv::gam(log10(Se) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Se, method = "REML")
summary(mod_gam_Se26) # d15N:length

AIC(mod_gam_Se0,mod_gam_Se1,mod_gam_Se10,mod_gam_Se15,mod_gam_Se16,
    mod_gam_Se17,mod_gam_Se18,mod_gam_Se19,mod_gam_Se20,mod_gam_Se21,
    mod_gam_Se22,mod_gam_Se23,mod_gam_Se24,mod_gam_Se26) #mod_gam_Se25,

## Step 3

# BEgin from n°23, remove length (mod 15), d15N (mod 19)

mod_gam_Se23bis <- mgcv::gam(log10(Se) ~ factor(sex)+
                               factor(new_cluster)+
                               s(d13C)+
                               factor(season)+
                               factor(sex):length+
                               factor(new_cluster):length+
                               factor(season):d13C+
                               d15N:length,
                             data = database_Se, method = "REML")
summary(mod_gam_Se23bis)

AIC(mod_gam_Se0,mod_gam_Se1,mod_gam_Se10,mod_gam_Se23,mod_gam_Se23bis)

mod_gam_Se27 <- mgcv::gam(log10(Se) ~ factor(new_cluster)+
                            s(d13C)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se27) # sex

mod_gam_Se28 <- mgcv::gam(log10(Se) ~ factor(sex)+
                            s(d13C)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se28) # clustFA

mod_gam_Se29 <- mgcv::gam(log10(Se) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se29) # d13C

mod_gam_Se30 <- mgcv::gam(log10(Se) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se30) # season

mod_gam_Se31 <- mgcv::gam(log10(Se) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se31) # sex:length

mod_gam_Se32 <- mgcv::gam(log10(Se) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se32) # clustFA:length

mod_gam_Se33 <- mgcv::gam(log10(Se) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            d15N:length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se33) # season:d13C

mod_gam_Se35 <- mgcv::gam(log10(Se) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d13C,
                          data = database_Se, method = "REML")
summary(mod_gam_Se35) # d15N:length

AIC(mod_gam_Se0,mod_gam_Se1,mod_gam_Se10,mod_gam_Se23bis,mod_gam_Se27,
    mod_gam_Se28,mod_gam_Se29,mod_gam_Se30,mod_gam_Se31,mod_gam_Se32,
    mod_gam_Se33,mod_gam_Se35)

## Step 4

# Begin from n°27, remove d15N:length (mod 35)

mod_gam_Se27bis <- mgcv::gam(log10(Se) ~ factor(new_cluster)+
                               s(d13C)+
                               factor(season)+
                               factor(sex):length+
                               factor(new_cluster):length+
                               factor(season):d13C,
                             data = database_Se, method = "REML")
summary(mod_gam_Se27bis)

AIC(mod_gam_Se0,mod_gam_Se1,mod_gam_Se10,mod_gam_Se23bis,mod_gam_Se27bis)

# Begin from 27bis
mod_gam_Se36 <- mgcv::gam(log10(Se) ~ s(d13C)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d13C,
                          data = database_Se, method = "REML")
summary(mod_gam_Se36) # clustFA

mod_gam_Se37 <- mgcv::gam(log10(Se) ~ factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d13C,
                          data = database_Se, method = "REML")
summary(mod_gam_Se37) # d13C

mod_gam_Se38 <- mgcv::gam(log10(Se) ~ factor(new_cluster)+
                            s(d13C)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d13C,
                          data = database_Se, method = "REML")
summary(mod_gam_Se38) # season

mod_gam_Se39 <- mgcv::gam(log10(Se) ~ factor(new_cluster)+
                            s(d13C)+
                            factor(season)+
                            factor(new_cluster):length+
                            factor(season):d13C,
                          data = database_Se, method = "REML")
summary(mod_gam_Se39) # sex:length

mod_gam_Se40 <- mgcv::gam(log10(Se) ~ factor(new_cluster)+
                            s(d13C)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):d13C,
                          data = database_Se, method = "REML")
summary(mod_gam_Se40) # clustFA:length

mod_gam_Se41 <- mgcv::gam(log10(Se) ~ factor(new_cluster)+
                            s(d13C)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length,
                          data = database_Se, method = "REML")
summary(mod_gam_Se41) # season:d13C

AIC(mod_gam_Se0,mod_gam_Se1,mod_gam_Se10,mod_gam_Se23bis,mod_gam_Se27,
    mod_gam_Se36,mod_gam_Se37,mod_gam_Se38,mod_gam_Se39,mod_gam_Se40,
    mod_gam_Se41)

## Best model is 27bis

rm(mod_gam_Se0,mod_gam_Se1,mod_gam_Se2,mod_gam_Se3,mod_gam_Se4,
   mod_gam_Se5,mod_gam_Se6,mod_gam_Se7,mod_gam_Se8,mod_gam_Se9,
   mod_gam_Se10,mod_gam_Se11,mod_gam_Se12,mod_gam_Se13,mod_gam_Se14,
   mod_gam_Se15,mod_gam_Se16,
   mod_gam_Se17,mod_gam_Se18,mod_gam_Se19,mod_gam_Se20,mod_gam_Se21,
   mod_gam_Se22,mod_gam_Se23,mod_gam_Se24,mod_gam_Se25,mod_gam_Se26,
   mod_gam_Se23bis,mod_gam_Se27,
   mod_gam_Se28,mod_gam_Se29,mod_gam_Se30,mod_gam_Se31,mod_gam_Se32,
   mod_gam_Se33,mod_gam_Se34,mod_gam_Se35,
   mod_gam_Se36,mod_gam_Se37,mod_gam_Se38,mod_gam_Se39,mod_gam_Se40)
rm(database_Se)


## 24 / Length-stand. Se ##############################################################################################

database_Se <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Se_stand) %>% 
  rename(Se = Se_stand)

## Step 1

mod_gam_Se0 <- mgcv::gam(log10(Se) ~ 1,
                         data = database_Se, method = "REML")

mod_gam_Se1 <- mgcv::gam(log10(Se) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Se, method = "REML")
summary(mod_gam_Se1) # Complete model

mod_gam_Se2 <- mgcv::gam(log10(Se) ~ factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Se, method = "REML")
summary(mod_gam_Se2) # sex

mod_gam_Se3 <- mgcv::gam(log10(Se) ~ factor(sex)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Se, method = "REML")
summary(mod_gam_Se3) # new_cluster

mod_gam_Se4 <- mgcv::gam(log10(Se) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Se, method = "REML")
summary(mod_gam_Se4) # d13C

mod_gam_Se5 <- mgcv::gam(log10(Se) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Se, method = "REML")
summary(mod_gam_Se5) # d15N

mod_gam_Se6 <- mgcv::gam(log10(Se) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Se, method = "REML")
summary(mod_gam_Se6) # season

mod_gam_Se7 <- mgcv::gam(log10(Se) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d13C,
                         data = database_Se, method = "REML")
summary(mod_gam_Se7) # season:d15N

mod_gam_Se8 <- mgcv::gam(log10(Se) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N,
                         data = database_Se, method = "REML")
summary(mod_gam_Se8) # season:d13C

AIC(mod_gam_Se0,mod_gam_Se1,mod_gam_Se2,mod_gam_Se3,mod_gam_Se4,
    mod_gam_Se5,mod_gam_Se6,mod_gam_Se7,mod_gam_Se8)

## Step 2

# Begin from n°6, remove sex (mod 2), d13C (mod 4), d15N (mod 5)

mod_gam_Se6bis <- mgcv::gam(log10(Se) ~ factor(new_cluster)+
                              factor(season):d15N+
                              factor(season):d13C,
                            data = database_Se, method = "REML")
summary(mod_gam_Se6bis)

AIC(mod_gam_Se0,mod_gam_Se1,mod_gam_Se6,mod_gam_Se6bis)

mod_gam_Se9 <- mgcv::gam(log10(Se) ~ factor(season):d15N+
                           factor(season):d13C,
                         data = database_Se, method = "REML")
summary(mod_gam_Se9) # clustFA

mod_gam_Se10 <- mgcv::gam(log10(Se) ~ factor(new_cluster)+
                            factor(season):d13C,
                          data = database_Se, method = "REML")
summary(mod_gam_Se10) # season:d15N

mod_gam_Se11 <- mgcv::gam(log10(Se) ~ factor(new_cluster)+
                            factor(season):d15N,
                          data = database_Se, method = "REML")
summary(mod_gam_Se11) # season:d13C

AIC(mod_gam_Se0,mod_gam_Se1,mod_gam_Se6bis,mod_gam_Se9,mod_gam_Se10,
    mod_gam_Se11)

## Best model is 6bis

rm(mod_gam_Se0,mod_gam_Se1,mod_gam_Se2,mod_gam_Se3,mod_gam_Se4,
   mod_gam_Se5,mod_gam_Se6,mod_gam_Se7,mod_gam_Se8,
   mod_gam_Se6bis,mod_gam_Se9,mod_gam_Se10,
   mod_gam_Se11)
rm(database_Se)


## 24 / Se & length-stand. Se - Plot response variable vs explaining variables ##############################################################################################

## LOG-TRANSFORMED Se
database_Se_ls <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Se)

# Best model is n°43

mod_gam_Se43 <- gam::gam(log10(Se) ~ s(d13C)+
                           factor(new_cluster)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):d13C,
                         #d13C:length,
                         data = database_Se_ls)
summary(mod_gam_Se43)

mod_gam_Se43bis <- mgcv::gam(log10(Se) ~ factor(new_cluster)+
                               s(d13C)+
                               factor(season)+
                               factor(sex):length+
                               factor(new_cluster):length+
                               factor(season):d13C+
                               d13C:length,
                             data = database_Se_ls, method = "REML")
summary(mod_gam_Se43bis)

## LOG-TRANSFORMED LENGTH-STAND. Se
database_Se <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Se_stand) %>% 
  rename(Se = Se_stand)

# Best model is n°6bis
mod_gam_Se6bis <- gam::gam(log10(Se) ~ factor(new_cluster)+
                             factor(season):d15N+
                             factor(season):d13C,
                           data = database_Se)
summary(mod_gam_Se6bis)

mod_gam_Se6ter <- mgcv::gam(log10(Se) ~ factor(new_cluster)+
                              factor(season):d15N+
                              factor(season):d13C,
                            data = database_Se, method = "REML")
summary(mod_gam_Se6ter)

par(mfrow=c(2,3))
plot(mod_gam_Se43,se = TRUE, ylab = "", xlab = "", las = 1)
plot(mod_gam_Se6bis,se = TRUE, ylab = "", xlab = "", las = 1)

rm(database_Se_ls,mod_gam_Se43,mod_gam_Se43bis)
rm(database_Se,mod_gam_Se6bis,mod_gam_Se6ter)


## 25 / Zinc ##############################################################################################

database_Zn <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Zn)

## Step 1
mod_gam_Zn0 <- mgcv::gam(log10(Zn) ~ 1,
                         data = database_Zn, method = "REML")

mod_gam_Zn1 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn1) # Complete model

mod_gam_Zn2 <- mgcv::gam(log10(Zn) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn2) # length

mod_gam_Zn3 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn3) # sex

mod_gam_Zn4 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                           factor(sex)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn4) # new_cluster

mod_gam_Zn5 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn5) # d13C

mod_gam_Zn6 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           factor(season)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn6) # d15N

mod_gam_Zn7 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(sex):length+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn7) # season

mod_gam_Zn8 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(new_cluster):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn8) # Sex:length

mod_gam_Zn9 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                           factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(sex):length+
                           factor(season):length+
                           factor(season):d15N+
                           factor(season):d13C+
                           d13C:length+
                           d15N:length,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn9) # new_cluster:length

mod_gam_Zn10 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn10) # season:length

mod_gam_Zn11 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn11) # season:d15N

mod_gam_Zn12 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            d13C:length+
                            d15N:length,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn12) # season:d13C

mod_gam_Zn13 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn13) # d13C:length

mod_gam_Zn14 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn14) # d15N:length

AIC(mod_gam_Zn0,mod_gam_Zn1,mod_gam_Zn2,mod_gam_Zn3,mod_gam_Zn4,
    mod_gam_Zn5,mod_gam_Zn6,mod_gam_Zn7,mod_gam_Zn8,mod_gam_Zn9,
    mod_gam_Zn10,mod_gam_Zn11,mod_gam_Zn12,mod_gam_Zn13,mod_gam_Zn14)

## Step 2

# Begin from n°7

mod_gam_Zn15 <- mgcv::gam(log10(Zn) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn15) # length

mod_gam_Zn16 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn16) # sex

mod_gam_Zn17 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn17) # clustFA

mod_gam_Zn18 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d15N)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn18) # d13C

mod_gam_Zn19 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn19) # d15N

mod_gam_Zn20 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn20) # sex:length

mod_gam_Zn21 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn21) # clustFA:length

mod_gam_Zn22 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn22) # season:length

mod_gam_Zn23 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d13C+
                            d13C:length+
                            d15N:length,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn23) # season:d15N

mod_gam_Zn24 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            d13C:length+
                            d15N:length,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn24) # season:d13C

mod_gam_Zn25 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d15N:length,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn25) # d13C:length

mod_gam_Zn26 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C+
                            d13C:length,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn26) # d15N:length

AIC(mod_gam_Zn0,mod_gam_Zn1,mod_gam_Zn7,mod_gam_Zn15,mod_gam_Zn16,
    mod_gam_Zn17,mod_gam_Zn18,mod_gam_Zn19,mod_gam_Zn20,mod_gam_Zn21,
    mod_gam_Zn22,mod_gam_Zn23,mod_gam_Zn24,mod_gam_Zn25,mod_gam_Zn26)

## Step 3

# Begin from n°26, remove d13C:length (mod 25)

mod_gam_Zn26bis <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                               factor(sex)+
                               factor(new_cluster)+
                               s(d13C)+
                               s(d15N)+
                               factor(sex):length+
                               factor(new_cluster):length+
                               factor(season):length+
                               factor(season):d15N+
                               factor(season):d13C,
                             data = database_Zn, method = "REML")
summary(mod_gam_Zn26bis)

AIC(mod_gam_Zn0,mod_gam_Zn1,mod_gam_Zn7,mod_gam_Zn26,mod_gam_Zn26bis)

mod_gam_Zn27 <- mgcv::gam(log10(Zn) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn27) # length

mod_gam_Zn28 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn28) # sex

mod_gam_Zn29 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn29) # clustFA

mod_gam_Zn30 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d15N)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn30) # d13C

mod_gam_Zn31 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn31) # d15N

mod_gam_Zn32 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn32) # sex:length

mod_gam_Zn33 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn33) # clustFA:length

mod_gam_Zn34 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn34) # season:length

mod_gam_Zn35 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn35) # season:d15N

mod_gam_Zn36 <- mgcv::gam(log10(Zn) ~ s(length, k = 3)+
                            factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(sex):length+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn36) # season:d13C

AIC(mod_gam_Zn0,mod_gam_Zn1,mod_gam_Zn7,mod_gam_Zn26bis,mod_gam_Zn27,
    mod_gam_Zn28,mod_gam_Zn29,mod_gam_Zn30,mod_gam_Zn31,mod_gam_Zn32,
    mod_gam_Zn33,mod_gam_Zn34,mod_gam_Zn35,mod_gam_Zn36)

## Step 4

# Begin from 26bis, remove length (mod 27), sex (mod 28), sex:length (mod 32)

mod_gam_Zn26ter <- mgcv::gam(log10(Zn) ~ factor(new_cluster)+
                               s(d13C)+
                               s(d15N)+
                               factor(new_cluster):length+
                               factor(season):length+
                               factor(season):d15N+
                               factor(season):d13C,
                             data = database_Zn, method = "REML")
summary(mod_gam_Zn26ter)

AIC(mod_gam_Zn0,mod_gam_Zn1,mod_gam_Zn7,mod_gam_Zn26bis,mod_gam_Zn26ter)

# Begin from 26ter
mod_gam_Zn37 <- mgcv::gam(log10(Zn) ~ s(d13C)+
                            s(d15N)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn37) # clustFA

mod_gam_Zn38 <- mgcv::gam(log10(Zn) ~ factor(new_cluster)+
                            s(d15N)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn38) # d13C

mod_gam_Zn39 <- mgcv::gam(log10(Zn) ~ factor(new_cluster)+
                            s(d13C)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn39) # d15N

mod_gam_Zn40 <- mgcv::gam(log10(Zn) ~ factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(season):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn40) # clustFA:length

mod_gam_Zn41 <- mgcv::gam(log10(Zn) ~ factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(new_cluster):length+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn41) # season:length

mod_gam_Zn42 <- mgcv::gam(log10(Zn) ~ factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn42) # season:d15N

mod_gam_Zn43 <- mgcv::gam(log10(Zn) ~ factor(new_cluster)+
                            s(d13C)+
                            s(d15N)+
                            factor(new_cluster):length+
                            factor(season):length+
                            factor(season):d15N,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn43) # season:d13C

AIC(mod_gam_Zn0,mod_gam_Zn1,mod_gam_Zn7,mod_gam_Zn26bis,mod_gam_Zn26ter,
    mod_gam_Zn37,mod_gam_Zn38,mod_gam_Zn39,mod_gam_Zn40,mod_gam_Zn41,
    mod_gam_Zn42,mod_gam_Zn43)

## Best model is 26ter

rm(mod_gam_Zn0,mod_gam_Zn1,mod_gam_Zn2,mod_gam_Zn3,mod_gam_Zn4,
   mod_gam_Zn5,mod_gam_Zn6,mod_gam_Zn7,mod_gam_Zn8,mod_gam_Zn9,
   mod_gam_Zn10,mod_gam_Zn11,mod_gam_Zn12,mod_gam_Zn13,mod_gam_Zn14,
   mod_gam_Zn15,mod_gam_Zn16,
   mod_gam_Zn17,mod_gam_Zn18,mod_gam_Zn19,mod_gam_Zn20,mod_gam_Zn21,
   mod_gam_Zn22,mod_gam_Zn23,mod_gam_Zn24,mod_gam_Zn25,mod_gam_Zn26,
   mod_gam_Zn26bis,mod_gam_Zn27,
   mod_gam_Zn28,mod_gam_Zn29,mod_gam_Zn30,mod_gam_Zn31,mod_gam_Zn32,
   mod_gam_Zn33,mod_gam_Zn34,mod_gam_Zn35,mod_gam_Zn36,
   mod_gam_Zn26ter,
   mod_gam_Zn37,mod_gam_Zn38,mod_gam_Zn39,mod_gam_Zn40,mod_gam_Zn41,
   mod_gam_Zn42,mod_gam_Zn43)
rm(database_Zn)


## 26 / Length-stand. Zn ##############################################################################################

database_Zn <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Zn_stand) %>% 
  rename(Zn = Zn_stand)

## Step 1

mod_gam_Zn0 <- mgcv::gam(log10(Zn) ~ 1,
                         data = database_Zn, method = "REML")

mod_gam_Zn1 <- mgcv::gam(log10(Zn) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn1) # Complete model

mod_gam_Zn2 <- mgcv::gam(log10(Zn) ~ factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn2) # sex

mod_gam_Zn3 <- mgcv::gam(log10(Zn) ~ factor(sex)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn3) # new_cluster

mod_gam_Zn4 <- mgcv::gam(log10(Zn) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn4) # d13C

mod_gam_Zn5 <- mgcv::gam(log10(Zn) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           factor(season)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn5) # d15N

mod_gam_Zn6 <- mgcv::gam(log10(Zn) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn6) # season

mod_gam_Zn7 <- mgcv::gam(log10(Zn) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d13C,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn7) # season:d15N

mod_gam_Zn8 <- mgcv::gam(log10(Zn) ~ factor(sex)+
                           factor(new_cluster)+
                           s(d13C)+
                           s(d15N)+
                           factor(season)+
                           factor(season):d15N,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn8) # season:d13C

AIC(mod_gam_Zn0,mod_gam_Zn1,mod_gam_Zn2,mod_gam_Zn3,mod_gam_Zn4,
    mod_gam_Zn5,mod_gam_Zn6,mod_gam_Zn7,mod_gam_Zn8)

## Step 2

# BEgin from n°6, remove d15N (mod 5)

mod_gam_Zn6bis <- mgcv::gam(log10(Zn) ~ factor(sex)+
                              factor(new_cluster)+
                              s(d13C)+
                              factor(season):d15N+
                              factor(season):d13C,
                            data = database_Zn, method = "REML")
summary(mod_gam_Zn6bis)

AIC(mod_gam_Zn0,mod_gam_Zn1,mod_gam_Zn6,mod_gam_Zn6bis)

mod_gam_Zn9 <- mgcv::gam(log10(Zn) ~ factor(new_cluster)+
                           s(d13C)+
                           factor(season):d15N+
                           factor(season):d13C,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn9) # sex

mod_gam_Zn10 <- mgcv::gam(log10(Zn) ~ factor(sex)+
                            s(d13C)+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn10) # clustFA

mod_gam_Zn11 <- mgcv::gam(log10(Zn) ~ factor(sex)+
                            factor(new_cluster)+
                            factor(season):d15N+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn11) # d13C

mod_gam_Zn12 <- mgcv::gam(log10(Zn) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn12) # season:d15N

mod_gam_Zn13 <- mgcv::gam(log10(Zn) ~ factor(sex)+
                            factor(new_cluster)+
                            s(d13C)+
                            factor(season):d15N,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn13) # season:d13C

AIC(mod_gam_Zn0,mod_gam_Zn1,mod_gam_Zn6bis,mod_gam_Zn9,mod_gam_Zn10,
    mod_gam_Zn11,mod_gam_Zn12,mod_gam_Zn13)

## Step 3

mod_gam_Zn10bis <- mgcv::gam(log10(Zn) ~ s(d13C)+
                               factor(season):d15N+
                               factor(season):d13C,
                             data = database_Zn, method = "REML")
summary(mod_gam_Zn10bis)

AIC(mod_gam_Zn0,mod_gam_Zn1,mod_gam_Zn6bis,mod_gam_Zn10,mod_gam_Zn10bis)

mod_gam_Zn14 <- mgcv::gam(log10(Zn) ~ factor(season):d15N+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn14) # d13C

mod_gam_Zn15 <- mgcv::gam(log10(Zn) ~ s(d13C)+
                            factor(season):d13C,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn15) # season:d15N

mod_gam_Zn16 <- mgcv::gam(log10(Zn) ~ s(d13C)+
                            factor(season):d15N,
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn16) # season:d13C

AIC(mod_gam_Zn0,mod_gam_Zn1,mod_gam_Zn6bis,mod_gam_Zn10bis,mod_gam_Zn14,
    mod_gam_Zn15,mod_gam_Zn16)

# Best model is 15

rm(mod_gam_Zn0,mod_gam_Zn1,mod_gam_Zn2,mod_gam_Zn3,mod_gam_Zn4,
   mod_gam_Zn5,mod_gam_Zn6,mod_gam_Zn7,mod_gam_Zn8,
   mod_gam_Zn6bis,mod_gam_Zn9,mod_gam_Zn10,
   mod_gam_Zn11,mod_gam_Zn12,mod_gam_Zn13,
   mod_gam_Zn10bis,mod_gam_Zn14,
   mod_gam_Zn15,mod_gam_Zn16)
rm(database_Zn)


## 27 / Zn & length-stand. Zn - Plot response variable vs explaining variables ##############################################################################################

## LOG-TRANSFORMED Zn
database_Zn_ls <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Zn)

# Best model is n°26ter

mod_gam_Zn26ter <- gam::gam(log10(Zn) ~ s(d13C)+
                              s(d15N)+
                              factor(new_cluster)+
                              factor(new_cluster):length+
                              factor(season):length+
                              factor(season):d15N+
                              factor(season):d13C,
                            data = database_Zn_ls)
summary(mod_gam_Zn26ter)

mod_gam_Zn26terbis <- mgcv::gam(log10(Zn) ~ factor(new_cluster)+
                                  s(d13C)+
                                  s(d15N)+
                                  factor(new_cluster):length+
                                  factor(season):length+
                                  factor(season):d15N+
                                  factor(season):d13C,
                                data = database_Zn_ls, method = "REML")
summary(mod_gam_Zn26terbis)

## LOG-TRANSFORMED LENGTH-STAND. Zn
database_Zn <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Zn_stand) %>% 
  rename(Zn = Zn_stand)

# Best model is n°15
mod_gam_Zn15 <- gam::gam(log10(Zn) ~ s(d13C)+
                           factor(season):d13C,
                         data = database_Zn)
summary(mod_gam_Zn15)

mod_gam_Zn15bis <- mgcv::gam(log10(Zn) ~ s(d13C)+
                               factor(season):d13C,
                             data = database_Zn, method = "REML")
summary(mod_gam_Zn15bis)

par(mfrow=c(2,3))
plot(mod_gam_Zn26ter,se = TRUE, ylab = "", xlab = "", las = 1)
plot(mod_gam_Zn15,se = TRUE, ylab = "", xlab = "", las = 1)

rm(database_Zn_ls,mod_gam_Zn26ter,mod_gam_Zn26terbis)
rm(database_Zn,mod_gam_Zn15,mod_gam_Zn15bis)



### VII // GAM - Adding interaction term d15N:sex:LJFL in best selected model ##############################################################################################

## 1 / Arsenic ##############################################################################################

database_As <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,As)

mod_gam_As31 <- mgcv::gam(log10(As) ~ s(length, k = 3)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length,
                          data = database_As, method = "REML")
summary(mod_gam_As31)

mod_gam_As31bis <- gam::gam(log10(As) ~ s(length)+
                            factor(new_cluster)+
                            factor(season)+
                            factor(sex):length+
                            factor(season):length+
                            factor(sex):length:d13C,
                          data = database_As)
summary(mod_gam_As31bis)

AIC(mod_gam_As31,mod_gam_As31bis)

rm(database_As, mod_gam_As31, mod_gam_As31bis)


## 2 / Cadmium ##############################################################################################

database_Cd <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Cd)

mod_gam_Cd36 <- mgcv::gam(log10(Cd) ~ factor(sex)+
                               factor(new_cluster)+
                               factor(sex):length+
                               factor(season):length+
                               factor(season):d15N,
                             data = database_Cd, method = "REML")
summary(mod_gam_Cd36)

mod_gam_Cd36bis <- mgcv::gam(log10(Cd) ~ factor(sex)+
                               factor(new_cluster)+
                               factor(sex):length+
                               factor(season):length+
                               factor(season):d15N+
                               factor(sex):length:d13C,
                             data = database_Cd, method = "REML")
summary(mod_gam_Cd36bis)

AIC(mod_gam_Cd36,mod_gam_Cd36bis)

rm(database_Cd, mod_gam_Cd36, mod_gam_Cd36bis)


## 3 / Cobalt ##############################################################################################

database_Co <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  dplyr::select(length,sex,new_cluster,d13C,d15N,season,Co)

mod_gam_Co32 <- mgcv::gam(log10(Co) ~ s(d15N)+
                               factor(season)+
                               factor(sex):length+
                               factor(season):length+
                               factor(season):d15N+
                               factor(season):d13C,
                             data = database_Co, method = "REML")
summary(mod_gam_Co32)

mod_gam_Co32bis <- mgcv::gam(log10(Co) ~ s(d15N)+
                               factor(season)+
                               factor(sex):length+
                               factor(season):length+
                               factor(season):d15N+
                               factor(season):d13C+
                               factor(sex):length:d15N,
                             data = database_Co, method = "REML")
summary(mod_gam_Co32bis)

AIC(mod_gam_Co32,mod_gam_Co32bis)

rm(database_Co, mod_gam_Co32,mod_gam_Co32bis)


## 4 / Copper ##############################################################################################

database_Cu <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Cu)

mod_gam_Cu35bis <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                               factor(season)+
                               d15N:length,
                             data = database_Cu_ls, method = "REML")
summary(mod_gam_Cu35bis)

mod_gam_Cu35ter <- mgcv::gam(log10(Cu) ~ s(length, k = 3)+
                               factor(season)+
                               d15N:length+
                               factor(sex):length:d15N,
                             data = database_Cu, method = "REML")
summary(mod_gam_Cu35ter)

AIC(mod_gam_Cu35bis,mod_gam_Cu35ter)

rm(database_Cu, mod_gam_Cu35bis, mod_gam_Cu35ter)


## 5 / Iron ##############################################################################################

database_Fe <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Fe)

mod_gam_Fe22bis <- mgcv::gam(log10(Fe) ~ factor(new_cluster)+
                               factor(season)+
                               factor(sex):length+
                               factor(new_cluster):length+
                               factor(season):d15N+
                               factor(season):d13C,
                             data = database_Fe, method = "REML")
summary(mod_gam_Fe22bis)

mod_gam_Fe22ter <- mgcv::gam(log10(Fe) ~ factor(new_cluster)+
                               factor(season)+
                               factor(sex):length+
                               factor(new_cluster):length+
                               factor(season):d15N+
                               factor(season):d13C+
                               factor(sex):length:d15N,
                             data = database_Fe, method = "REML")
summary(mod_gam_Fe22ter)

AIC(mod_gam_Fe22bis,mod_gam_Fe22ter)

rm(database_Fe, mod_gam_Fe22bis, mod_gam_Fe22ter)


## 6 / Mercury ##############################################################################################

database_Hg <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Hg)

mod_gam_Hg34 <- mgcv::gam(log10(Hg) ~ factor(sex)+
                               factor(new_cluster):length,
                             data = database_Hg, method = "REML")
summary(mod_gam_Hg34)

mod_gam_Hg34bis <- mgcv::gam(log10(Hg) ~ factor(sex)+
                                factor(new_cluster):length+
                                factor(sex):length:d15N,
                             data = database_Hg, method = "REML")
summary(mod_gam_Hg34bis)

AIC(mod_gam_Hg34,mod_gam_Hg34bis)

rm(database_Hg, mod_gam_Hg34, mod_gam_Hg34bis)


## 7 / Manganese ##############################################################################################

database_Mn <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Mn)

mod_gam_Mn20 <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                               factor(sex)+
                               factor(new_cluster)+
                               factor(season)+
                               factor(sex):length+
                               factor(season):d15N+
                               factor(season):d13C+
                               d13C:length+
                               d15N:length,
                             data = database_Mn_ls, method = "REML")
summary(mod_gam_Mn20)

mod_gam_Mn20bis <- mgcv::gam(log10(Mn) ~ s(length, k = 3)+
                               factor(sex)+
                               factor(new_cluster)+
                               factor(season)+
                               factor(sex):length+
                               factor(season):d15N+
                               factor(season):d13C+
                               factor(sex):length:d15N,
                             data = database_Mn, method = "REML")
summary(mod_gam_Mn20bis) 

AIC(mod_gam_Mn20,mod_gam_Mn20bis)

rm(database_Mn, mod_gam_Mn20, mod_gam_Mn20bis)


## 8 / Lead ##############################################################################################

database_Pb <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Pb)

mod_gam_Pb35bis <- mgcv::gam(log10(Pb) ~ factor(new_cluster):length+
                               factor(season):length+
                               factor(season):d15N+
                               factor(season):d13C+
                               d15N:length,
                             data = database_Pb, method = "REML")
summary(mod_gam_Pb35bis)

mod_gam_Pb35ter <- mgcv::gam(log10(Pb) ~ factor(new_cluster):length+
                                 factor(season):length+
                                 factor(season):d15N+
                                 factor(season):d13C+
                                 d15N:length+
                                 factor(sex):length:d15N,
                             data = database_Pb, method = "REML")
summary(mod_gam_Pb35ter)

AIC(mod_gam_Pb35bis,mod_gam_Pb35ter)

rm(database_Pb, mod_gam_Pb35bis, mod_gam_Pb35ter)


## 9 / Selenium ##############################################################################################

database_Se <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Se)

mod_gam_Se43 <- mgcv::gam(log10(Se) ~ factor(new_cluster)+
                               s(d13C)+
                               factor(season)+
                               factor(sex):length+
                               factor(new_cluster):length+
                               factor(season):d13C+
                               d13C:length,
                             data = database_Se_ls, method = "REML")
summary(mod_gam_Se43)

mod_gam_Se43bis <- mgcv::gam(log10(Se) ~ factor(new_cluster)+
                               s(d13C)+
                               factor(season)+
                               factor(sex):length+
                               factor(new_cluster):length+
                               factor(season):d13C+
                               d13C:length+
                               factor(sex):length:d15N,
                               data = database_Se, method = "REML")
summary(mod_gam_Se43bis)

AIC(mod_gam_Se43,mod_gam_Se43bis)

rm(database_Se, mod_gam_Se43, mod_gam_Se43bis)


## 10 / Zinc ##############################################################################################

database_Zn <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2")),
         season = factor(season, levels = c("IMS","NWM","IMA","SEM")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  filter(!is.na(new_cluster),
         !is.na(d13C),
         !is.na(sex)) %>% 
  select(length,sex,new_cluster,d13C,d15N,season,Zn)

mod_gam_Zn26bis <- mgcv::gam(log10(Zn) ~ factor(new_cluster)+
                                  s(d13C)+
                                  s(d15N)+
                                  factor(new_cluster):length+
                                  factor(season):length+
                                  factor(season):d15N+
                                  factor(season):d13C,
                                data = database_Zn_ls, method = "REML")
summary(mod_gam_Zn26bis)

mod_gam_Zn26ter <- mgcv::gam(log10(Zn) ~ factor(new_cluster)+
                                  s(d13C)+
                                  s(d15N)+
                                  factor(new_cluster):length+
                                  factor(season):length+
                                  factor(season):d15N+
                                  factor(season):d13C+
                                  factor(sex):length:d15N,
                             data = database_Zn, method = "REML")
summary(mod_gam_Zn26ter)

AIC(mod_gam_Zn26bis,mod_gam_Zn26ter)

rm(database_Zn, mod_gam_Zn26bis, mod_gam_Zn26ter)



### VIII // Appendices ##############################################################################################

## 1 / LJFL according to sex, FA trophic group and season ##############################################################################################

## Statistical tests for SEX
data_test_sex <- SWO_data %>% 
  mutate(sex = factor(sex, levels = c("M","F")))

shapiro.test(data_test_sex$length) # Normality
fligner.test(data_test_sex$length ~ data_test_sex$sex) # Homoscedasticité

wilcox.test(data_test_sex$length ~ data_test_sex$sex)

## Statistical tests for FA TROPHIC GROUP
data_test_FA <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2"))) %>% 
  filter(!is.na(new_cluster))

shapiro.test(data_test_FA$length) # Normality
fligner.test(data_test_FA$length ~ data_test_FA$new_cluster) # Homoscedasticité

t.test(data_test_FA$length ~ data_test_FA$new_cluster)

## Statistical tests for SEASON
data_test_season <- SWO_data %>% 
  mutate(season = factor(season, levels = c("IMS","NWM","IMA","SEM")))

shapiro.test(data_test_season$length) # Normality
fligner.test(data_test_season$length ~ data_test_season$season) # Homoscedasticité

kruskal.test(data_test_season$length ~ data_test_season$season)
dunnTest(data_test_season$length ~ data_test_season$season, method = "bh")


## Plots
sex <- data_test_sex %>% 
  mutate(sex = ifelse(sex == "M", "Males", "Females")) %>% 
  group_by(sex) %>% 
  mutate(max_value = max(length,na.rm = TRUE),
         mean = mean(length, na.rm = TRUE),
         letter = ifelse(sex == "Females", "a", "b")) %>% 
  ungroup() %>% 
  mutate(sex = factor(sex, levels = c("Males","Females"))) %>% 
  ggplot()+
  geom_boxplot(aes(x = sex, y = length), fill = "#A6A4A4")+
  geom_text(aes(x = sex, y = max_value+10, label = letter), color = "#58585D")+
  labs(x = "Sex", y = "Lower jaw-fork length (cm)")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "italic"),
        legend.title = element_blank())

FAcluster <- data_test_FA %>% 
  mutate(letter = "a",
         letter = ifelse(new_cluster == "2", "b", letter),
         new_cluster = as.character(new_cluster),
         new_cluster = ifelse(new_cluster == "1", "FA trophic group 1", "FA trophic group 2")) %>%
  group_by(new_cluster) %>% 
  mutate(position = max(length)) %>% 
  ungroup() %>% 
  ggplot()+
  geom_boxplot(aes(x = new_cluster, y = length), fill = "#A6A4A4")+
  geom_text(aes(x = new_cluster, y = position+10, label = letter))+
  labs(y = "Lower jaw-fork length (cm)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(face = "italic"))

season <- data_test_season %>% 
  group_by(season) %>% 
  mutate(max_value = max(length,na.rm = TRUE),
         mean = mean(length, na.rm = TRUE),
         letter = NA,
         letter = ifelse(season == "NWM", "a", letter),
         letter = ifelse(season == "SEM" | season == "IMS", "ab", letter),
         letter = ifelse(season == "IMA", "b", letter)) %>%
  ungroup() %>% 
  mutate(season = as.character(season),
         season = ifelse(season == "IMS", "Pre-NWM", season),
         season = ifelse(season == "IMA", "Pre-SEM", season),
         season = factor(season, levels = c("Pre-NWM","NWM","Pre-SEM","SEM"))) %>% 
  ggplot()+
  geom_boxplot(aes(x = season, y = length), fill = "#A6A4A4")+
  geom_text(aes(x = season, y = max_value+10, label = letter), color = "#58585B")+
  scale_fill_manual(values = c("#3D4C53","#E64A45","#4DB3B3","#F2C249"))+
  labs(x = "Season", y = "Lower jaw-fork length (cm)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
    axis.title.y = element_text(face = "italic"),
    legend.title = element_blank(),
    legend.position = "none")

empty <- as_ggplot(NULL)

ggarrange(sex, FAcluster, season, empty,
          nrow = 2, ncol = 2, labels = c("A.", "B.", "C.", ""), align = "hv")

rm(data_test_sex, data_test_FA, data_test_season)
rm(sex, FAcluster, season, empty)


## 2 / d13C and d15N values according to sex, FA trophic group and season ##############################################################################################

## Statistical tests for SEX
data_test_sex <- SWO_data %>% 
  mutate(sex = factor(sex, levels = c("M","F")))

shapiro.test(data_test_sex$d13C) # Normality
shapiro.test(data_test_sex$d15N) # Normality

fligner.test(data_test_sex$d13C ~ data_test_sex$sex) # Homoscedasticité
fligner.test(data_test_sex$d15N ~ data_test_sex$sex) # Homoscedasticité

wilcox.test(data_test_sex$d13C ~ data_test_sex$sex) #d13C
wilcox.test(data_test_sex$d15N ~ data_test_sex$sex) #d15N

## Statistical tests for FA TROPHIC GROUP
data_test_FA <- SWO_data %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2"))) %>% 
  filter(!is.na(new_cluster))

shapiro.test(data_test_FA$d13C) # Normality
fligner.test(data_test_FA$d13C ~ data_test_FA$new_cluster) # Homoscedasticité
wilcox.test(data_test_FA$d13C ~ data_test_FA$new_cluster)

shapiro.test(data_test_FA$d15N) # Normality
fligner.test(data_test_FA$d15N ~ data_test_FA$new_cluster) # Homoscedasticité
t.test(data_test_FA$d15N ~ data_test_FA$new_cluster)

## Statistical tests for SEASON
data_test_season <- SWO_data %>% 
  mutate(season = factor(season, levels = c("IMS","NWM","IMA","SEM")))

shapiro.test(data_test_season$d13C) # Normality
shapiro.test(data_test_season$d15N) # Normality

fligner.test(data_test_season$d13C ~ data_test_season$season) # Homoscedasticité
fligner.test(data_test_season$d15N ~ data_test_season$season) # Homoscedasticité

kruskal.test(SWO_data$d13C ~ SWO_data$season) #d13C
dunnTest(SWO_data$d13C ~ SWO_data$season, method = "bh")

kruskal.test(SWO_data$d15N ~ SWO_data$season) #d15N
dunnTest(SWO_data$d15N ~ SWO_data$season, method = "bh")

## Plots
sex_SIvalues <- data_test_sex %>% 
  select(sex,d13C,d15N) %>% 
  gather(SI,value,-sex) %>% 
  mutate(sex = ifelse(sex == "M", "Males", "Females")) %>% 
  group_by(sex,SI) %>% 
  mutate(max_value = max(value,na.rm = TRUE),
         mean = mean(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(sex = factor(sex, levels = c("Males","Females"))) %>% 
  ggplot()+
  geom_boxplot(aes(x = sex, y = value), fill = "#A6A4A4")+
  labs(x = "Sex", y = "Values (permil)")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "italic"),
        legend.title = element_blank())+
  facet_wrap(.~SI, scales = "free")+
  theme(strip.text.x = element_text(face = "bold"))

clustFA_SIvalues <- data_test_FA %>%
  select(new_cluster,d13C,d15N) %>% 
  gather(SI, value, -new_cluster) %>% 
  mutate(new_cluster = as.character(new_cluster),
         new_cluster = ifelse(new_cluster == "1", "FA trophic group 1", new_cluster),
         new_cluster = ifelse(new_cluster == "2", "FA trophic group 2", new_cluster),
         new_cluster = factor(new_cluster, levels = c("FA trophic group 1","FA trophic group 2"))) %>% 
  ggplot()+
  geom_boxplot(aes(x = new_cluster, y = value), fill = "#A6A4A4")+
  labs(y = "Values (permil)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(face = "italic"))+
  facet_wrap(.~SI, scales = "free")+
  theme(strip.text.x = element_text(face = "bold"))

season_SIvalues <- data_test_season %>% 
  select(season,d13C,d15N) %>% 
  gather(SI,value,-season) %>% 
  group_by(season,SI) %>% 
  mutate(max_value = max(value,na.rm = TRUE),
         mean = mean(value, na.rm = TRUE),
         letter = NA,
         letter = ifelse(season == "SEM", "ab", letter),
         letter = ifelse(season == "IMS", "c", letter),
         letter = ifelse(season == "NWM", "a", letter),
         letter = ifelse(season == "IMA", "b", letter)) %>% 
  ungroup() %>% 
  mutate(season = as.character(season),
         season = ifelse(season == "IMS", "Pre-NWM", season),
         season = ifelse(season == "IMA", "Pre-SEM", season),
         season = factor(season, levels = c("Pre-NWM","NWM","Pre-SEM","SEM"))) %>% 
  ggplot()+
  geom_boxplot(aes(x = season, y = value), fill = "#A6A4A4")+
  geom_text(aes(x = season, y = max_value+0.5, label = letter))+
  labs(x = "Season", y = "Values (permil)")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "italic"),
        legend.title = element_blank())+
  facet_wrap(.~SI, scales = "free")+
  theme(strip.text.y = element_text(face = "bold"))

ggarrange(ggarrange(sex_SIvalues,clustFA_SIvalues, ncol = 2,
                    labels = c("A.","B."), align = "hv"),
          season_SIvalues, nrow = 2,
          labels= c("","C."), align = "hv")

rm(data_test_sex, data_test_FA,data_test_season)
rm(sex_SIvalues,clustFA_SIvalues,season_SIvalues)


## 3 / Correlations among bioachemical tracers (SI + FAs) ##############################################################################################

## Calcul correlation coeff and p-values
data_tracers_x <- SWO_data %>% 
  filter(!is.na(c22_6w3_c),
         !is.na(d15N)) %>% 
  select(d13C,d15N,c14_c,c16_1w7_c,c18_1w7_c,c18_1w9_c,c18_c,
         c20_1w9_c,c20_4w6_c,c20_5w3_c,c22_5w3_c,c22_5w6_c,
         c22_6w3_c,c24_1w9_c)
data_tracers_x <- as.data.frame(data_tracers_x)

data_tracers_y <- SWO_data %>% 
  filter(!is.na(c22_6w3_c),
         !is.na(d15N)) %>% 
  select(d13C,d15N,c14_c,c16_1w7_c,c18_1w7_c,c18_1w9_c,c18_c,
         c20_1w9_c,c20_4w6_c,c20_5w3_c,c22_5w3_c,c22_5w6_c,
         c22_6w3_c,c24_1w9_c)
data_tracers_y <- as.data.frame(data_tracers_y)

Output_ALL_cor <- NULL # création objet final
for (j in 1 : length(colnames(data_tracers_y))){ # tracer loop
  #j = 1
  
  # Sélection du traceur trophique
  Trophic_tracer <- data_tracers_y[,j]
  Tracer_name <- colnames(data_tracers_y)[j]
  
  
  for (i in 1: length(colnames(data_tracers_x))){ # metal loop
    #i = 1
    
    Metal_name <- colnames(data_tracers_x)[i]
    Metal_data <- data_tracers_x[,i]
    
    test_a <- shapiro.test(Trophic_tracer)
    test_b <- shapiro.test(Metal_data)
    
    Output_cor <- data.frame(Trophic_tracer = Tracer_name,
                             Metal = Metal_name,
                             N_val_used = sum((is.na(Trophic_tracer) | is.na(Metal_data )) == F),
                             Normality = NA,
                             Test_cor = NA,
                             P_val = NA,
                             Cor_val = NA)
    
    
    if (Output_cor$N_val_used <= 2){
      
      Output_ALL_cor <- rbind(Output_ALL_cor, Output_cor)
      
    }else{
      
      
      if(test_a$p.value >= 0.05 & test_b$p.value >= 0.05){# Pearson
        Output_cor$Normality <- T
        
        test_stat <- cor.test(Trophic_tracer, Metal_data, method = "pearson")
        
        Output_cor$Test_cor <- "pearson"
        Output_cor$P_val <- test_stat$p.value
        Output_cor$Cor_val <- test_stat$estimate
        
      }else{
        # Kendall
        Output_cor$Normality <- F
        
        test_stat <- cor.test(Trophic_tracer, Metal_data, method = "kendall")
        test_stat$p.value
        
        Output_cor$Test_cor <- "kendall"
        Output_cor$P_val <- test_stat$p.value
        Output_cor$Cor_val <- test_stat$estimate
        
      } # fin test correlation
      
      Output_ALL_cor <- rbind(Output_ALL_cor, Output_cor)
      rm(Output_cor, test_stat, test_a, test_b, Metal_data, Metal_name)
      
    } # fin boucle si N suffisant
  } # fin boucle metal
  rm(Tracer_name,Trophic_tracer)
} # fin boucle 

rm(i,j)
rm(data_tracers_x,data_tracers_y)

## Création matrices de correlation et p_values
cor_tracers <- Output_ALL_cor %>% 
  select(Metal,Trophic_tracer,Cor_val) %>% 
  spread(Metal,Cor_val)
rownames(cor_tracers) <- cor_tracers$Trophic_tracer
cor_tracers <- cor_tracers %>% 
  select(-Trophic_tracer)
cor_tracers <- as.matrix(cor_tracers)

sig_tracers <- Output_ALL_cor %>% 
  select(Metal,Trophic_tracer,P_val) %>% 
  spread(Metal,P_val)
rownames(sig_tracers) <- sig_tracers$Trophic_tracer
sig_tracers <- sig_tracers %>% 
  select(-Trophic_tracer)
sig_tracers <- as.matrix(sig_tracers)

## Correlation plot
corrplot(cor_tracers, method = "color",
         type = "lower",
         addCoef.col = "black",
         tl.col="black", tl.srt = 45,
         p.mat = sig_tracers, sig.level = 0.05, insig = "blank",
         diag = F)

rm(cor_tracers,sig_tracers)
rm(Output_ALL_cor)


## 4 / SWO sex ratio by season ##############################################################################################

SWO_data %>% 
  select(season,sex) %>%
  group_by(season,sex) %>% 
  mutate(n_sex = n()) %>% 
  ungroup() %>% 
  distinct(season,sex, .keep_all = TRUE) %>% 
  mutate(sex = ifelse(sex == "M", "Males", "Females"),
         season = as.character(season),
         season = ifelse(season == "IMS", "Pre-NWM", season),
         season = ifelse(season == "IMA", "Pre-SEM", season),
         season = factor(season, levels = c("Pre-NWM","NWM","Pre-SEM","SEM"))) %>% 
  spread(sex, n_sex) %>% 
  mutate(SexRatio_M_F = Females/Males) %>% 
  ggplot(aes(x = season, y = SexRatio_M_F))+
  geom_point()+
  geom_line(group = 1)+
  geom_text(aes(y = SexRatio_M_F+0.2, label = round(SexRatio_M_F, 1)))+
  theme_bw()+
  labs(y = "Sex ratio (Females/Males)")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "italic"),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())

