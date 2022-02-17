##--------------------------------------------------------------------------------------------------------
## SCRIPT : Allows to reproduce all analyses of Chapter 4 dedicated to Hg and Se bioaccumulation
##          in swordfish (Xiphias gladius) from the Indian Ocean. Includes complementary analysis
##          used only in the related paper (Sabino et al. "Regional patterns in mercury (Hg) and
##          selenium (Se) concentration of swordfish in the Indian Ocean").
##          This script includes the following analyses :
##            - Plotting the sampling locations across the Indian Ocean
##            - Calculation of MHg:MSe and HBVSe ratios, as well as theoretically available Se
##            - Statistical tests for significant difference in Hg, Se, length-stand. Hg
##              and theoretically available Se concentrations, d13C and d15N values, lower
##              jaw-fork length (LJFL), MHg:MSe and HBVSe among areas in the Indian Ocean
##            - Correlation test for relationship between log(Hg) and log(Se)
##            - Length-standardisation method for Hg concentrations
##            - Plot of Hg, Se, length-stand. Hg and theoretically available Se concentrations,
##              d13C and d15N values, lower jaw-fork length (LJFL), MHG:MSe and HBVSe ratios in
##              each sampling area
##            - Plot correlation between log(Hg) and log(Se)
##            - GAM to test relationship between Hg or Se concentrations and LJFL, d13C/d15N values
##              and longitude and latitude
##            - Calculation of number of servings to reach PTI and/or percentage of covered RDI
##              for children, young aduls and adults according to the sampling area
##
## As part of :
##        Magali SABINO PhD - "Bioaccumulation of trace elements in Seychelles marine food webs"
##
## Author : Magali Sabino
## First created : 2022-01-26
## Last update : 2022-02-10
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
lapply(c("tidyverse", "openxlsx", "sp", "sf", "rnaturalearth",
         "rnaturalearthdata", "ggspatial", "rgeos", "ggpubr",
         "FSA", "corrplot", "scatterpie", "mgcv"),
       library, character.only = TRUE)



### I // Creation of database with SWO samples from the Indian Ocean ##############################################################################################

## 1 / Creation of database with data SI and TE ##############################################################################################

SWOIO_data <- read.xlsx("C:/Users/msabin02/Nextcloud/PhD/R_projects/Final_scripts/SWOIO_TE_SI.xlsx", sheet = "Data", colNames = TRUE)
SWOIO_data <- SWOIO_data %>% 
  rename(longitude = long_mean, latitude = lat_mean) %>% 
  select(organism_identifier,project,date,lowerjawfork_length,longitude,
         latitude,area,d13C,d15N,THg,Seww,TP) %>% 
  mutate(longitude = round(longitude, 2),
         latitude = round(latitude, 2)) %>% 
  filter(!is.na(d13C),
         !is.na(Seww),
         !project %in% c("Paco","CONSWO_SEY")) %>% 
  mutate(area = factor(area, levels = c("BENG","ISLU","WTIO","MOZ","SSG","MCSA")))


## 2 / Add trophic level ##############################################################################################

baseline <- read.xlsx("C:/Users/msabin02/Nextcloud/PhD/R_projects/Final_scripts/MOBI_Output_IO_Swordfish_2004-2018_06Jun2021.xlsx", colNames = TRUE)%>% 
  select(organism_identifier,d15N_phyt)

SWOIO_data <- SWOIO_data %>% 
  left_join(baseline, by = "organism_identifier") %>% 
  mutate(trophic_level = 1 + ((d15N-d15N_phyt)/3.4))

rm(baseline)


## 3 / Add MHg:MSe, HBVSe and Se avail ##############################################################################################

SWOIO_data <- SWOIO_data %>%
  mutate(MHg = THg/200590000,
         MSe = Seww/78960000,
         Se_avail = (MSe-MHg)*78960000,
         MHg_MSe = MHg/MSe,
         HBVSe = ((MSe-MHg)/MSe)*(MSe+MHg))



### II // Mapping sampling locations in the Indian Ocean ##############################################################################################

world <- ne_countries(scale = "medium", returnclass = "sf", type = "map_units")
class(world)

Coords_lim <- matrix(data = c(min(SWOIO_data$latitude), max(SWOIO_data$latitude),
                              min(SWOIO_data$longitude), max(SWOIO_data$longitude)),
                     nrow = 2, ncol = 2)
colnames(Coords_lim) <- c("latitude", "longitude")
row.names(Coords_lim) <- c("min", "max")


theme_set(theme_bw(base_size = 14))
ggplot(data = world) +
  geom_sf() +
  geom_point(data = SWOIO_data, aes(x=longitude, y=latitude, fill=area), color = "black", shape = 21, size=2)+
  coord_sf(xlim = c(floor(Coords_lim[1,2]), ceiling(Coords_lim[2,2])),
           ylim = c(floor(Coords_lim[1,1]), ceiling(Coords_lim[2,1]))) +
  scale_fill_manual(values = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"))+
  labs(x = 'Longitude', y = "Latitude", fill = 'Sampling area') +
  theme(axis.title = element_blank(),
        axis.text = element_text(size=12, colour="black"),
        legend.title = element_blank(),
        legend.text = element_text(size=14, face="bold", colour="black"))

rm(Coords_lim, world)



### III // Plot Hg ans Se concentrations, d13C and d15N values, trophic level and relationship between log(Hg) and log(Se) ##############################################################################################

## 1 / Statistical tests - Hg concentrations among areas ##############################################################################################

a1 <- aov(SWOIO_data$THg ~ SWOIO_data$area)
shapiro.test(resid(a1)) # Normality of residuals
fligner.test(data = SWOIO_data, THg ~ area) # Homoscedasticity

kruskal.test(data = SWOIO_data, THg ~ area)
dunnTest(data = SWOIO_data, THg ~ area, method = "bh")

rm(a1)

## 2 / Statistical tests - Se concentrations among areas ##############################################################################################

a1 <- aov(SWOIO_data$Seww ~ SWOIO_data$area)
shapiro.test(resid(a1)) # Normality of residuals
fligner.test(data = SWOIO_data, Seww ~ area) # Homoscedasticity

kruskal.test(data = SWOIO_data, Seww ~ area)
dunnTest(data = SWOIO_data, Seww ~ area, method = "bh")

rm(a1)


## 3 / Statistical tests - d13C values among areas ##############################################################################################

a1 <- aov(SWOIO_data$d13C ~ SWOIO_data$area)
shapiro.test(resid(a1)) # Normality of residuals
fligner.test(data = SWOIO_data, d13C ~ area) # Homoscedasticity

kruskal.test(data = SWOIO_data, d13C ~ area)
dunnTest(data = SWOIO_data, d13C ~ area, method = "bh")

rm(a1)


## 4 / Statistical tests - d15N values among areas ##############################################################################################

a1 <- aov(SWOIO_data$d15N ~ SWOIO_data$area)
shapiro.test(resid(a1)) # Normality of residuals
fligner.test(data = SWOIO_data, d15N ~ area) # Homoscedasticity

kruskal.test(data = SWOIO_data, d15N ~ area)
dunnTest(data = SWOIO_data, d15N ~ area, method = "bh")

lmFA <- lm(SWOIO_data$d15N ~ SWOIO_data$area)
anova(lmFA)
TukeyHSD(a1, 'SWOIO_data$area', conf.level = 0.95)

rm(a1,lmFA)


## 6 / Statistical tests - Trophic positions among areas ##############################################################################################

a1 <- aov(SWOIO_data$trophic_level ~ SWOIO_data$area)
shapiro.test(resid(a1)) # Normality of residuals
fligner.test(data = SWOIO_data, trophic_level ~ area) # Homoscedasticity

kruskal.test(data = SWOIO_data, trophic_level ~ area)
dunnTest(data = SWOIO_data, trophic_level ~ area, method = "bh")

rm(a1)


## 7 / Correlation between log(Hg) and log(Se) ##############################################################################################

shapiro.test(log10(SWOIO_data$THg))
shapiro.test(SWOIO_data$Seww)
cor.test(log10(SWOIO_data$THg), log10(SWOIO_data$Seww), method = "pearson")


## 8 / Plot all on same graph ##############################################################################################

# THg concentrations
THg_concentrations <- SWOIO_data %>%
  group_by(area) %>% 
  mutate(mean = mean(THg, na.rm = TRUE),
         se = std(THg),
         y_letter = max(THg, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(letter = "a",
         letter = ifelse(area == "ISLU", "ab", letter),
         letter = ifelse(area %in% c("BENG","WTIO","MCSA"), "b", letter),
         letter = ifelse(area == "MOZ", "c", letter)) %>%
  ggplot(aes(x = area))+
  geom_boxplot(aes(y = THg, fill = area))+
  geom_point(aes(y = mean), shape = 21, fill = "white", color = "black", size = 2.5)+
  geom_text(aes(y = y_letter+0.3, vjust = 0, label = letter), color = "#58585B")+
  scale_fill_manual(values = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"))+
  labs(y = "Hg (µg.g-1 ww)")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title.y = element_text(face = "italic"),
        axis.text.x = element_text(color = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"),
                                   angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank())

# Se concentrations
Se_concentrations <- SWOIO_data %>%
  group_by(area) %>% 
  mutate(mean = mean(Seww, na.rm = TRUE),
         se = std(Seww),
         y_letter = max(Seww, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(letter = "a",
         letter = ifelse(area %in% c("ISLU","WTIO","SSG"), "b", letter),
         letter = ifelse(area %in% c("MOZ","MCSA"), "c",letter)) %>%
  ggplot(aes(x = area))+
  geom_boxplot(aes(y = Seww, fill = area))+
  geom_point(aes(y = mean), shape = 21, fill = "white", color = "black", size = 2.5)+
  geom_text(aes(y = y_letter+0.3, vjust = 0, label = letter), color = "#58585B")+
  scale_fill_manual(values = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"))+
  labs(y = "Se (µg.g-1 ww)")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title.y = element_text(face = "italic"),
        axis.text.x = element_text(color = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"),
                                   angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank())

# d13C values
d13C_values <- SWOIO_data %>% 
  group_by(area) %>% 
  mutate(mean = mean(d13C, na.rm = TRUE),
         y_letter = max(d13C)) %>% 
  ungroup() %>% 
  mutate(letter = "a",
         letter = ifelse(area == "MOZ", "b", letter),
         letter = ifelse(area == "SSG", "bc", letter),
         letter = ifelse(area == "WTIO", "c", letter),
         letter = ifelse(area == "BENG", "cd", letter),
         letter = ifelse(area == "MCSA", "d", letter)) %>%
  ggplot(aes(x = area)) +
  geom_boxplot(aes(y = d13C, fill = area)) +
  geom_point(aes(y = mean), shape = 21, fill = "white", color = "black", size = 2)+
  geom_text(aes(y = y_letter+0.2, vjust = 0, label = letter), color = "#58585B")+
  scale_fill_manual(values = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"))+
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) +
  labs(y = "d13C (permil)")+
  theme(legend.position = "none",
        axis.title.y = element_text(face = "italic"),
        axis.text.x = element_text(color = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"),
                                   angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank())

# d15N values
d15N_values <- SWOIO_data %>% 
  filter(!is.na(d15N)) %>%
  group_by(area) %>% 
  mutate(mean = mean(d15N, na.rm = TRUE),
         y_letter = max(d15N)) %>% 
  ungroup() %>% 
  mutate(letter = "a",
         letter = ifelse(area == "MCSA", "ab", letter),
         letter = ifelse(area == "ISLU", "abc", letter),
         letter = ifelse(area == "MOZ", "bc", letter),
         letter = ifelse(area == "BENG", "c", letter)) %>%
  ggplot(aes(x = area)) +
  geom_boxplot(aes(y = d15N, fill = area)) +
  geom_point(aes(y = mean), shape = 21, fill = "white", color = "black", size = 2)+
  geom_text(aes(y = y_letter+0.2, vjust = 0, label = letter), color = "#58585B")+
  scale_fill_manual(values = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"))+
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) +
  labs(y = "d15N (permil)")+
  theme(legend.position = "none",
        axis.title.y = element_text(face = "italic"),
        axis.text.x = element_text(color = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"),
                                   angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank())

# Trophic levels
TP <- SWOIO_data %>% 
  filter(!is.na(trophic_level)) %>%
  group_by(area) %>% 
  mutate(mean = mean(trophic_level, na.rm = TRUE),
         y_letter = max(trophic_level)) %>% 
  ungroup() %>% 
  mutate(letter = "a",
         letter = ifelse(area == "WTIO", "b", letter),
         letter = ifelse(area == "MCSA", "bc", letter),
         letter = ifelse(area == "MOZ", "c", letter),
         letter = ifelse(area %in% c("BENG","ISLU"), "d", letter)) %>%
  ggplot(aes(x = area)) +
  geom_boxplot(aes(y = trophic_level, fill = area)) +
  geom_point(aes(y = mean), shape = 21, fill = "white", color = "black", size = 2)+
  geom_text(aes(y = y_letter+0.1, vjust = 0, label = letter), color = "#58585B")+
  scale_fill_manual(values = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"))+
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) +
  labs(y = "Trophic level")+
  theme(legend.position = "none",
        axis.title.y = element_text(face = "italic"),
        axis.text.x = element_text(color = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"),
                                   angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank())

# Relationship log(Hg) and log(Se)
logHg_logSe <- SWOIO_data %>% 
  filter(!is.na(Seww)) %>% 
  ggplot()+
  geom_point(aes(x = log10(THg), y = log10(Seww), fill = area),  shape = 21, size = 3, color = "black")+
  geom_smooth(aes(x = log10(THg), y = log10(Seww)), color = "red", method = "lm", se = FALSE)+
  labs(x = "log(Hg)", y = "log(Se)")+
  scale_fill_manual(values = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"))+
  theme_bw()+
  theme(axis.title = element_text(face = "italic"),
        legend.position = "none",
        legend.title = element_blank())

# All on same graph
ggarrange(
  ggarrange(
    ggarrange(THg_concentrations, Se_concentrations,nrow = 2, labels = c("A.","B.")),
              logHg_logSe, ncol = 2, labels = c("","C."), widths = c(1,2)),
  ggarrange(d13C_values,d15N_values,TP, ncol = 3, labels = c("D.","E.","F.")),
  nrow = 2, heights = c(2,1))

rm(THg_concentrations, Se_concentrations, d13C_values,d15N_values,TP, logHg_logSe)



### IV // Plot Relationship between LJFL and Hg concentrations, and LJFL and length-stand. Hg by area ##############################################################################################

## 1 / Length-standardisation of Hg concentrations ##############################################################################################

# Database
data_norm <- SWOIO_data %>% 
  select(organism_identifier,area,lowerjawfork_length,THg) %>% 
  rename(Hg = THg) %>% 
  filter(!is.na(lowerjawfork_length),
         !is.na(Hg),
         !organism_identifier %in% c("RUNSWO22")) %>% 
  mutate(logHg = log10(Hg))
data_norm <- as.data.frame(data_norm)

# Adjusting bioaccumulation curve - logHg ~ LJFL
pow.funk <- function(l, a, b, c, d) a*(l-c)^b-d

x <- data_norm$logHg
l <- data_norm$lowerjawfork_length

fit <- nls(x ~ pow.funk(l, a, b, c, d), algo='port', start=list(a=1, b=0.5, c=min(l), d=2), control=nls.control(maxiter=500),
           lower=c(0.0001,0.001,0,-2),upper=c(1000,1,min(l),8))

a <- coef(fit)[1]
b <- coef(fit)[2]
c <- coef(fit)[3]
d <- coef(fit)[4]

summary(fit)

fit_coeff_all <- data.frame(formula = "a*(l-c)^b-d",
                            a = a,
                            b = b,
                            c = c,
                            d = d)

summary(data_norm$lowerjawfork_length)
tbet <- data.frame(length=seq(63, 300, 1))
tbet[,2] <-  round(a * (tbet[,1]-c) ^ b - d, digits=5)

# Visual verification of adjusted curve 
summary(data_norm$logHg)
data_norm %>% 
  ggplot() +
  geom_point(aes(x = lowerjawfork_length, y=logHg, color = area), size=2) +
  geom_line(data = tbet, aes(x = tbet[,1], y = tbet[,2]), size=0.8, col='black') +
  labs(x='Lower jaw-fork length (cm)', y='log(Hg)') +
  theme_classic(base_size=14) +
  theme(legend.title=element_blank(),
        legend.position = c(1, .1),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size=7))

# Calculation of residuals and of length-stand. Hg concentrations
data_norm <- data_norm %>% 
  mutate(logHg_length = a*(lowerjawfork_length-c)^b-d,
         residuals = logHg - logHg_length,
         logHg_lmean = pow.funk(150, a, b, c, d),
         logHg_stand = residuals + logHg_lmean,
         Hg_stand = 10^(logHg_stand))

R2 <- data_norm %>% 
  mutate(mean_logHg = mean(logHg),
         SSR_ind = (logHg_length - logHg)^2,
         SST_ind = (logHg - mean_logHg)^2) %>% 
  summarise(SSR = sum(SSR_ind),
            SST = sum(SST_ind),
            R2 = 1- (SSR/SST))

# Database update with length-stand. Hg concentrations
Hg_norm_ALL <- data_norm %>% 
  select(organism_identifier,Hg_stand)

SWOIO_data <- SWOIO_data %>%
  left_join(Hg_norm_ALL, by = "organism_identifier")

rm(fit,a,b,c,d,l,x,fit_coeff_all,R2,Hg_norm_ALL,pow.funk)


## 2 / Plot relationship between LJFL and log(Hg) (with adjusted curve) ##############################################################################################

fit.plot <- data_norm %>% 
  mutate(x.par = min(lowerjawfork_length),
         y.par = max(logHg),
         fit.par = "Y=5.70 (X - 48.54)0.07 - 8\nR2 = 0.46") %>% 
  ggplot() +
  geom_point(aes(x = lowerjawfork_length, y=logHg, color = area), size=2) +
  geom_line(data = tbet, aes(x = tbet[,1], y = tbet[,2]), size=0.8, col='black') +
  geom_text(aes(x = x.par, y = y.par, label = fit.par), hjust = 0, vjust = 1, size = 3.5)+
  labs(x='Lower jaw-fork length (cm)', y='Log Hg') +
  scale_color_manual(values = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"))+
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


## 3 / Plot LJFL by area ##############################################################################################

# Test for significant difference between areas
a1 <- aov(SWOIO_data$lowerjawfork_length ~ SWOIO_data$area)
shapiro.test(resid(a1)) # Normality of residuals
fligner.test(data = SWOIO_data, lowerjawfork_length ~ area) # Homoscedasticity

kruskal.test(data = SWOIO_data, lowerjawfork_length ~ area)
dunnTest(data = SWOIO_data, lowerjawfork_length ~ area, method = "bh")

rm(a1)

# Plot LJFL by area
LJFL <- SWOIO_data %>% 
  filter(!is.na(lowerjawfork_length)) %>%
  group_by(area) %>% 
  mutate(mean = mean(lowerjawfork_length, na.rm = TRUE),
         y_letter = max(lowerjawfork_length)) %>% 
  ungroup() %>% 
  mutate(letter = "a",
         letter = ifelse(area %in%c("WTIO","SSG","MCSA"), "b", letter),
         letter = ifelse(area %in% c("BENG","MOZ"), "c", letter),
         area = factor(area, levels = c("MCSA","SSG","MOZ","WTIO","ISLU","BENG"))) %>%
  ggplot(aes(x = area)) +
  geom_boxplot(aes(y = lowerjawfork_length, fill = area)) +
  geom_point(aes(y = mean), shape = 21, fill = "white", color = "black", size = 2.5)+
  geom_text(aes(y = y_letter+10, vjust = 0, label = letter), color = "#404041")+
  scale_fill_manual(values = c("#354B5E","#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340"))+
  coord_flip()+
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) +
  labs(y = "Lower jaw-fork length (cm)")+
  theme(legend.position = "none",
        axis.title.x = element_text(face = "italic"),
        axis.text.y = element_text(color = c("#354B5E","#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340")),
        axis.title.y = element_blank())


## 4 / Plot Length-stand. Hg concentrations by area ##############################################################################################

# Test for significant difference among areas
a1 <- aov(SWOIO_data$Hg_stand ~ SWOIO_data$area)
shapiro.test(resid(a1)) # Normality of residuals
fligner.test(data = SWOIO_data, Hg_stand ~ area) # Homoscedasticity

kruskal.test(data = SWOIO_data, Hg_stand ~ area)
dunnTest(data = SWOIO_data, Hg_stand ~ area, method = "bh")

rm(a1)

# Plot Length-stand. Hg concentrations
Hg_norm <- SWOIO_data %>% 
  mutate(letter = "a",
         letter = ifelse(area == "BENG", "ab", letter),
         letter = ifelse(area == "WTIO", "bc", letter),
         letter = ifelse(area %in% c("ISLU","MCSA"), "cd", letter),
         letter = ifelse(area == "MOZ", "d", letter)) %>%
  group_by(area) %>% 
  mutate(mean = mean(Hg_stand),
         max = max(Hg_stand)) %>% 
  ungroup() %>% 
  ggplot()+
  geom_boxplot(aes(x = area, y = Hg_stand, fill = area))+
  geom_point(aes(x = area, y = mean), shape = 21, fill = "white", color = "black", size = 2.5)+
  geom_text(aes(x = area, y = max+0.5, label = letter), color = "#404041")+
  scale_fill_manual(values = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"))+
  labs(y = "Length-standardised Hg (µg.g-1 ww)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "italic"), 
        axis.text.x = element_text(color = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"),
                                   angle = 45, hjust = 1, vjust = 1),
        plot.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")


## 5 / Plot all on same graph ##############################################################################################

empty <- as_ggplot(NULL)

ggarrange(ggarrange(fit.plot, Hg_norm,
                    ncol = 2, widths = c(2,1), labels = c("A.","C.")),
          ggarrange(LJFL, empty,
                    ncol = 2, widths = c(2,1), labels = c("B.","")),
          nrow = 2, heights = c(2,1))

rm(empty, fit.plot, Hg_norm, LJFL, data_norm, tbet)



### V // GAM - Relationships between Hg or Se concentrations and LJFL, d15N and d13C values, longitude and latitude, TP ##############################################################################################

## 1 / Test correlations between long/lat and LJFL, d13C and trophic level ##############################################################################################

## Calcul correlation coeff and p-values
data_tracers_x <- SWOIO_data %>% 
  filter(!is.na(THg)) %>% 
  select(lowerjawfork_length,d13C,trophic_level)
data_tracers_x <- as.data.frame(data_tracers_x)

data_tracers_y <- SWOIO_data %>% 
  filter(!is.na(THg)) %>% 
  select(latitude,longitude)
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
         addCoef.col = "black",
         tl.col="black", tl.srt = 45,
         p.mat = sig_tracers, sig.level = 0.05, insig = "blank")

rm(cor_tracers,sig_tracers)
rm(Output_ALL_cor)


## 2 / Mercury ##############################################################################################

database_THg <- SWOIO_data %>% 
  filter(!is.na(THg))

## Test concurvity
mod_gam_THg_conc <- mgcv::gam(log10(THg) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(d15N) +
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            d13C:longitude +
                            d13C:latitude +
                            d13C:latitude +
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_THg, method = "REML")
summary(mod_gam_THg_conc) # Complete model

concurvity(mod_gam_THg_conc)

rm(mod_gam_THg_conc)

## Step 1

mod_gam_THg0 <- mgcv::gam(log10(THg) ~ 1,
                          data = database_THg, method = "REML")
summary(mod_gam_THg0) # Null  model

mod_gam_THg1 <- mgcv::gam(log10(THg) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_THg, method = "REML")
summary(mod_gam_THg1) # Complete model

mod_gam_THg2 <- mgcv::gam(log10(THg) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_THg, method = "REML")
summary(mod_gam_THg2) # length

mod_gam_THg3 <- mgcv::gam(log10(THg) ~ s(lowerjawfork_length)+
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_THg, method = "REML")
summary(mod_gam_THg3) # d13C

mod_gam_THg4 <- mgcv::gam(log10(THg) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(latitude) +
                            s(longitude)+
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_THg, method = "REML")
summary(mod_gam_THg4) # TP

mod_gam_THg5 <- mgcv::gam(log10(THg) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(longitude)+
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_THg, method = "REML")
summary(mod_gam_THg5) # lat

mod_gam_THg6 <- mgcv::gam(log10(THg) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_THg, method = "REML")
summary(mod_gam_THg6) # long

mod_gam_THg7 <- mgcv::gam(log10(THg) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            d13C:latitude +
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_THg, method = "REML")
summary(mod_gam_THg7) # d13C:long

mod_gam_THg8 <- mgcv::gam(log10(THg) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            d13C:longitude +
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_THg, method = "REML")
summary(mod_gam_THg8) # d13C:lat

mod_gam_THg9 <- mgcv::gam(log10(THg) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            d13C:longitude +
                            d13C:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_THg, method = "REML")
summary(mod_gam_THg9) # TP:lat

mod_gam_THg10 <- mgcv::gam(log10(THg) ~ s(lowerjawfork_length)+
                             s(d13C) + 
                             s(trophic_level) + 
                             s(latitude) +
                             s(longitude)+
                             d13C:longitude +
                             d13C:latitude +
                             trophic_level:latitude,
                           data = database_THg, method = "REML")
summary(mod_gam_THg10) # length:TP

AIC(mod_gam_THg0,mod_gam_THg1,mod_gam_THg2,mod_gam_THg3,mod_gam_THg4,mod_gam_THg5,mod_gam_THg6,
    mod_gam_THg7,mod_gam_THg8,mod_gam_THg9,mod_gam_THg10)

## Step 2

# Begin from n°10
mod_gam_THg11 <- mgcv::gam(log10(THg) ~ s(d13C) + 
                             s(trophic_level) + 
                             s(latitude) +
                             s(longitude)+
                             d13C:longitude +
                             d13C:latitude +
                             trophic_level:latitude,
                           data = database_THg, method = "REML")
summary(mod_gam_THg11) # length

mod_gam_THg12 <- mgcv::gam(log10(THg) ~ s(lowerjawfork_length)+
                             s(trophic_level) + 
                             s(latitude) +
                             s(longitude)+
                             d13C:longitude +
                             d13C:latitude +
                             trophic_level:latitude,
                           data = database_THg, method = "REML")
summary(mod_gam_THg12) # d13C

mod_gam_THg13 <- mgcv::gam(log10(THg) ~ s(lowerjawfork_length)+
                             s(d13C) + 
                             s(latitude) +
                             s(longitude)+
                             d13C:longitude +
                             d13C:latitude +
                             trophic_level:latitude,
                           data = database_THg, method = "REML")
summary(mod_gam_THg13) # TP

mod_gam_THg14 <- mgcv::gam(log10(THg) ~ s(lowerjawfork_length)+
                             s(d13C) + 
                             s(trophic_level) + 
                             s(longitude)+
                             d13C:longitude +
                             d13C:latitude +
                             trophic_level:latitude,
                           data = database_THg, method = "REML")
summary(mod_gam_THg14) # Latitude

mod_gam_THg15 <- mgcv::gam(log10(THg) ~ s(lowerjawfork_length)+
                             s(d13C) + 
                             s(trophic_level) + 
                             s(latitude) +
                             d13C:longitude +
                             d13C:latitude +
                             trophic_level:latitude,
                           data = database_THg, method = "REML")
summary(mod_gam_THg15) # Longitude

mod_gam_THg16 <- mgcv::gam(log10(THg) ~ s(lowerjawfork_length)+
                             s(d13C) + 
                             s(trophic_level) + 
                             s(latitude) +
                             s(longitude)+
                             d13C:latitude +
                             trophic_level:latitude,
                           data = database_THg, method = "REML")
summary(mod_gam_THg16) # d13C:long

mod_gam_THg17 <- mgcv::gam(log10(THg) ~ s(lowerjawfork_length)+
                             s(d13C) + 
                             s(trophic_level) + 
                             s(latitude) +
                             s(longitude)+
                             d13C:longitude +
                             trophic_level:latitude,
                           data = database_THg, method = "REML")
summary(mod_gam_THg17) # d13C:lat

mod_gam_THg18 <- mgcv::gam(log10(THg) ~ s(lowerjawfork_length)+
                             s(d13C) + 
                             s(trophic_level) + 
                             s(latitude) +
                             s(longitude)+
                             d13C:longitude +
                             d13C:latitude,
                           data = database_THg, method = "REML")
summary(mod_gam_THg18) # TP:lat

AIC(mod_gam_THg0,mod_gam_THg10,mod_gam_THg11,mod_gam_THg12,
    mod_gam_THg13,mod_gam_THg14,mod_gam_THg15,mod_gam_THg16,
    mod_gam_THg17,mod_gam_THg18)

# Best model is n°10

rm(mod_gam_THg0,mod_gam_THg1,mod_gam_THg2,mod_gam_THg3,mod_gam_THg4,mod_gam_THg5,mod_gam_THg6,
   mod_gam_THg7,mod_gam_THg8,mod_gam_THg9,mod_gam_THg10,mod_gam_THg11,mod_gam_THg12,mod_gam_THg13,
   mod_gam_THg14,mod_gam_THg15,mod_gam_THg16,mod_gam_THg17,
   mod_gam_THg18)

## Plotting GAM results
mod_gam_THg10 <- mgcv::gam(log10(THg) ~ s(lowerjawfork_length)+
                             s(d13C) + 
                             s(trophic_level) + 
                             s(latitude) +
                             s(longitude)+
                             d13C:longitude +
                             d13C:latitude +
                             trophic_level:latitude,
                           data = database_THg, method = "REML")
summary(mod_gam_THg10)

par(mfrow=c(2,3))
plot(mod_gam_THg10, se = TRUE, ylab = "", xlab = "", las = 1)

rm(database_THg,mod_gam_THg10)

## 4 / Length-stand. Hg ##############################################################################################

database_THgstand <- SWOIO_data %>% 
  filter(!is.na(Hg_stand))

## Test concurvity
mod_gam_THgstand_conc <- mgcv::gam(log10(Hg_stand) ~ s(lowerjawfork_length)+
                                s(d13C) + 
                                s(d15N) +
                                s(trophic_level) + 
                                s(latitude) +
                                s(longitude)+
                                d13C:longitude +
                                d13C:latitude +
                                d13C:latitude +
                                trophic_level:latitude +
                                lowerjawfork_length:trophic_level,
                              data = database_THgstand, method = "REML")
summary(mod_gam_THgstand_conc) # Complete model

concurvity(mod_gam_THgstand_conc)

rm(mod_gam_THgstand_conc)

## Step 1

mod_gam_THg0 <- mgcv::gam(log10(Hg_stand) ~ 1,
                          data = database_THgstand, method = "REML")
summary(mod_gam_THg0) # Null  model

mod_gam_THg1 <- mgcv::gam(log10(Hg_stand) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude,
                          data = database_THgstand, method = "REML")
summary(mod_gam_THg1) # Complete model

mod_gam_THg2 <- mgcv::gam(log10(Hg_stand) ~ s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude,
                          data = database_THgstand, method = "REML")
summary(mod_gam_THg2) # d13C

mod_gam_THg3 <- mgcv::gam(log10(Hg_stand) ~ s(d13C) + 
                            s(latitude) +
                            s(longitude)+
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude,
                          data = database_THgstand, method = "REML")
summary(mod_gam_THg3) # TP

mod_gam_THg4 <- mgcv::gam(log10(Hg_stand) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(longitude)+
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude,
                          data = database_THgstand, method = "REML")
summary(mod_gam_THg4) # lat

mod_gam_THg5 <- mgcv::gam(log10(Hg_stand) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude,
                          data = database_THgstand, method = "REML")
summary(mod_gam_THg5) # long

mod_gam_THg6 <- mgcv::gam(log10(Hg_stand) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            d13C:latitude +
                            trophic_level:latitude,
                          data = database_THgstand, method = "REML")
summary(mod_gam_THg6) # d13C:long

mod_gam_THg7 <- mgcv::gam(log10(Hg_stand) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            d13C:longitude +
                            trophic_level:latitude,
                          data = database_THgstand, method = "REML")
summary(mod_gam_THg7) # d13C:lat

mod_gam_THg8 <- mgcv::gam(log10(Hg_stand) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            d13C:longitude +
                            d13C:latitude,
                          data = database_THgstand, method = "REML")
summary(mod_gam_THg8) # TP:lat

AIC(mod_gam_THg0,mod_gam_THg1,mod_gam_THg2,mod_gam_THg3,mod_gam_THg4,mod_gam_THg5,mod_gam_THg6,
    mod_gam_THg7,mod_gam_THg8)

# Best model is complete model

rm(mod_gam_THg0,mod_gam_THg1,mod_gam_THg2,mod_gam_THg3,mod_gam_THg4,mod_gam_THg5,mod_gam_THg6,
   mod_gam_THg7,mod_gam_THg8)

## Plotting GAM results
mod_gam_THg1 <- mgcv::gam(log10(Hg_stand) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude,
                          data = database_THgstand, method = "REML")
summary(mod_gam_THg1)

par(mfrow=c(2,2))
plot(mod_gam_THg1, se = TRUE, ylab = "", xlab = "", las = 1)

rm(database_THgstand,mod_gam_THg1)


## 5 / Selenium ##############################################################################################

database_Se <- SWOIO_data %>% 
  filter(!is.na(Seww))

## Test concurvity
mod_gam_Se_conc <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                                     s(d13C) + 
                                     s(d15N) +
                                     s(trophic_level) + 
                                     s(latitude) +
                                     s(longitude)+
                                     d13C:longitude +
                                     d13C:latitude +
                                     d13C:latitude +
                                     trophic_level:latitude +
                                     lowerjawfork_length:trophic_level,
                                   data = database_Se, method = "REML")
summary(mod_gam_Se_conc) # Complete model

concurvity(mod_gam_Se_conc)

rm(mod_gam_Se_conc)

## Step 1

mod_gam_Se0 <- mgcv::gam(log10(Seww) ~ 1,
                         data = database_Se, method = "REML")
summary(mod_gam_Se0) # Null  model

mod_gam_Se1 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                           s(d13C) + 
                           s(trophic_level) + 
                           s(latitude) +
                           s(longitude)+
                           d13C:longitude +
                           d13C:latitude +
                           trophic_level:latitude +
                           lowerjawfork_length:trophic_level,
                         data = database_Se, method = "REML")
summary(mod_gam_Se1) # Complete model

mod_gam_Se2 <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                           s(trophic_level) + 
                           s(latitude) +
                           s(longitude)+
                           d13C:longitude +
                           d13C:latitude +
                           trophic_level:latitude +
                           lowerjawfork_length:trophic_level,
                         data = database_Se, method = "REML")
summary(mod_gam_Se2) # length

mod_gam_Se3 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                           s(trophic_level) + 
                           s(latitude) +
                           s(longitude)+
                           d13C:longitude +
                           d13C:latitude +
                           trophic_level:latitude +
                           lowerjawfork_length:trophic_level,
                         data = database_Se, method = "REML")
summary(mod_gam_Se3) # d13C

mod_gam_Se4 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                           s(d13C) + 
                           s(latitude) +
                           s(longitude)+
                           d13C:longitude +
                           d13C:latitude +
                           trophic_level:latitude +
                           lowerjawfork_length:trophic_level,
                         data = database_Se, method = "REML")
summary(mod_gam_Se4) # TP

mod_gam_Se5 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                           s(d13C) + 
                           s(trophic_level) + 
                           s(longitude)+
                           d13C:longitude +
                           d13C:latitude +
                           trophic_level:latitude +
                           lowerjawfork_length:trophic_level,
                         data = database_Se, method = "REML")
summary(mod_gam_Se5) # lat

mod_gam_Se6 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                           s(d13C) + 
                           s(trophic_level) + 
                           s(latitude) +
                           d13C:longitude +
                           d13C:latitude +
                           trophic_level:latitude +
                           lowerjawfork_length:trophic_level,
                         data = database_Se, method = "REML")
summary(mod_gam_Se6) # long

mod_gam_Se7 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                           s(d13C) + 
                           s(trophic_level) + 
                           s(latitude) +
                           s(longitude)+
                           d13C:latitude +
                           trophic_level:latitude +
                           lowerjawfork_length:trophic_level,
                         data = database_Se, method = "REML")
summary(mod_gam_Se7) # d13C:long

mod_gam_Se8 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                           s(d13C) + 
                           s(trophic_level) + 
                           s(latitude) +
                           s(longitude)+
                           d13C:longitude +
                           trophic_level:latitude +
                           lowerjawfork_length:trophic_level,
                         data = database_Se, method = "REML")
summary(mod_gam_Se8) # d13C:lat

mod_gam_Se9 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                           s(d13C) + 
                           s(trophic_level) + 
                           s(latitude) +
                           s(longitude)+
                           d13C:longitude +
                           d13C:latitude +
                           lowerjawfork_length:trophic_level,
                         data = database_Se, method = "REML")
summary(mod_gam_Se9) # TP:lat

mod_gam_Se10 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude,
                          data = database_Se, method = "REML")
summary(mod_gam_Se10) # length:TP

AIC(mod_gam_Se0,mod_gam_Se1,mod_gam_Se2,mod_gam_Se3,mod_gam_Se4,mod_gam_Se5,mod_gam_Se6,
    mod_gam_Se7,mod_gam_Se8,mod_gam_Se9,mod_gam_Se10)

## Step 2

# Begin from n°7, remove d13C:lat (mod 8)
mod_gam_Se7bis <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                              s(d13C) + 
                              s(trophic_level) + 
                              s(latitude) +
                              s(longitude)+
                              trophic_level:latitude +
                              lowerjawfork_length:trophic_level,
                            data = database_Se, method = "REML")
summary(mod_gam_Se7bis)

AIC(mod_gam_Se0,mod_gam_Se7,mod_gam_Se7bis)

# Begin from 7bis
mod_gam_Se11 <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_Se, method = "REML")
summary(mod_gam_Se11) # Length

mod_gam_Se12 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_Se, method = "REML")
summary(mod_gam_Se12) # d13C

mod_gam_Se13 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(latitude) +
                            s(longitude)+
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_Se, method = "REML")
summary(mod_gam_Se13) # TP

mod_gam_Se14 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(longitude)+
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_Se, method = "REML")
summary(mod_gam_Se14) # lat

mod_gam_Se15 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_Se, method = "REML")
summary(mod_gam_Se15) # long

mod_gam_Se16 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            lowerjawfork_length:trophic_level,
                          data = database_Se, method = "REML")
summary(mod_gam_Se16) # TP:lat

mod_gam_Se17 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            trophic_level:latitude,
                          data = database_Se, method = "REML")
summary(mod_gam_Se17) # Length:TP

AIC(mod_gam_Se0,mod_gam_Se7bis,mod_gam_Se11,mod_gam_Se12,mod_gam_Se13,
    mod_gam_Se14,mod_gam_Se15,mod_gam_Se16,mod_gam_Se17)

## Step 3

# Begin from 11, remove length:TP (mod 17)
mod_gam_Se11bis <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                               s(trophic_level) + 
                               s(latitude) +
                               s(longitude)+
                               trophic_level:latitude,
                             data = database_Se, method = "REML")
summary(mod_gam_Se11bis)

AIC(mod_gam_Se0,mod_gam_Se7bis,mod_gam_Se11,mod_gam_Se11bis)

# Begin from 11bis
mod_gam_Se18 <- mgcv::gam(log10(Seww) ~ s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            trophic_level:latitude,
                          data = database_Se, method = "REML")
summary(mod_gam_Se18) # d13C

mod_gam_Se19 <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                            s(latitude) +
                            s(longitude)+
                            trophic_level:latitude,
                          data = database_Se, method = "REML")
summary(mod_gam_Se19) # TP

mod_gam_Se20 <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(longitude)+
                            trophic_level:latitude,
                          data = database_Se, method = "REML")
summary(mod_gam_Se20) # lat

mod_gam_Se21 <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            trophic_level:latitude,
                          data = database_Se, method = "REML")
summary(mod_gam_Se21) # long

mod_gam_Se22 <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude),
                          data = database_Se, method = "REML")
summary(mod_gam_Se22) #TP:lat

AIC(mod_gam_Se0,mod_gam_Se7bis,mod_gam_Se11bis,mod_gam_Se18,mod_gam_Se19,
    mod_gam_Se20,mod_gam_Se21,mod_gam_Se22)

## Step 4

# Begin from n°22
mod_gam_Se23 <- mgcv::gam(log10(Seww) ~ s(trophic_level) + 
                            s(latitude) +
                            s(longitude),
                          data = database_Se, method = "REML")
summary(mod_gam_Se23) # d13C

mod_gam_Se24 <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                            s(latitude) +
                            s(longitude),
                          data = database_Se, method = "REML")
summary(mod_gam_Se24) # TP

mod_gam_Se25 <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(longitude),
                          data = database_Se, method = "REML")
summary(mod_gam_Se25) # lat

mod_gam_Se26 <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(latitude),
                          data = database_Se, method = "REML")
summary(mod_gam_Se26) # long

AIC(mod_gam_Se0,mod_gam_Se7bis,mod_gam_Se11bis,mod_gam_Se22,mod_gam_Se23,
    mod_gam_Se24,mod_gam_Se25,mod_gam_Se26)

# Best model is n°22

rm(mod_gam_Se0,mod_gam_Se1,mod_gam_Se2,mod_gam_Se3,mod_gam_Se4,mod_gam_Se5,mod_gam_Se6,
   mod_gam_Se7,mod_gam_Se8,mod_gam_Se9,mod_gam_Se10,
   mod_gam_Se7bis,mod_gam_Se11,mod_gam_Se12,mod_gam_Se13,
   mod_gam_Se14,mod_gam_Se15,mod_gam_Se16,mod_gam_Se17,
   mod_gam_Se11bis,mod_gam_Se18,mod_gam_Se19,
   mod_gam_Se20,mod_gam_Se21,mod_gam_Se22,
   mod_gam_Se23,
   mod_gam_Se24,mod_gam_Se25,mod_gam_Se26)

## Plotting GAM results
mod_gam_Se22 <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude),
                          data = database_Se, method = "REML")
summary(mod_gam_Se22)

par(mfrow=c(2,2))
plot(mod_gam_Se22, se = TRUE, ylab = "", xlab = "", las = 1)

rm(database_Se,mod_gam_Se22)


## 6 / Selenium - Adding log(Hg) as explaining variable ##############################################################################################

database_Se <- SWOIO_data %>% 
  filter(!is.na(Seww))

## Test concurvity
mod_gam_Se_conc <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                               s(d13C) + 
                               s(d15N) +
                               s(trophic_level) + 
                               s(latitude) +
                               s(longitude)+
                               s(log10(THg))+
                               d13C:longitude +
                               d13C:latitude +
                               d13C:latitude +
                               trophic_level:latitude +
                               lowerjawfork_length:trophic_level,
                             data = database_Se, method = "REML")
summary(mod_gam_Se_conc) # Complete model

concurvity(mod_gam_Se_conc)

rm(mod_gam_Se_conc)

## Step 1

mod_gam_Se0 <- mgcv::gam(log10(Seww) ~ 1,
                         data = database_Se, method = "REML")
summary(mod_gam_Se0) # Null  model

mod_gam_Se1 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                           s(d13C) + 
                           s(trophic_level) + 
                           s(latitude) +
                           s(longitude)+
                           s(log10(THg))+
                           d13C:longitude +
                           d13C:latitude +
                           trophic_level:latitude +
                           lowerjawfork_length:trophic_level,
                         data = database_Se, method = "REML")
summary(mod_gam_Se1) # Complete model

mod_gam_Se2 <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                           s(trophic_level) + 
                           s(latitude) +
                           s(longitude)+
                           s(log10(THg))+
                           d13C:longitude +
                           d13C:latitude +
                           trophic_level:latitude +
                           lowerjawfork_length:trophic_level,
                         data = database_Se, method = "REML")
summary(mod_gam_Se2) # length

mod_gam_Se3 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                           s(trophic_level) + 
                           s(latitude) +
                           s(longitude)+
                           s(log10(THg))+
                           d13C:longitude +
                           d13C:latitude +
                           trophic_level:latitude +
                           lowerjawfork_length:trophic_level,
                         data = database_Se, method = "REML")
summary(mod_gam_Se3) # d13C

mod_gam_Se4 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                           s(d13C) + 
                           s(latitude) +
                           s(longitude)+
                           s(log10(THg))+
                           d13C:longitude +
                           d13C:latitude +
                           trophic_level:latitude +
                           lowerjawfork_length:trophic_level,
                         data = database_Se, method = "REML")
summary(mod_gam_Se4) # TP

mod_gam_Se5 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                           s(d13C) + 
                           s(trophic_level) + 
                           s(longitude)+
                           s(log10(THg))+
                           d13C:longitude +
                           d13C:latitude +
                           trophic_level:latitude +
                           lowerjawfork_length:trophic_level,
                         data = database_Se, method = "REML")
summary(mod_gam_Se5) # lat

mod_gam_Se6 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                           s(d13C) + 
                           s(trophic_level) + 
                           s(latitude) +
                           s(log10(THg))+
                           d13C:longitude +
                           d13C:latitude +
                           trophic_level:latitude +
                           lowerjawfork_length:trophic_level,
                         data = database_Se, method = "REML")
summary(mod_gam_Se6) # long

mod_gam_Se7 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                           s(d13C) + 
                           s(trophic_level) + 
                           s(latitude) +
                           s(longitude)+
                           s(log10(THg))+
                           d13C:latitude +
                           trophic_level:latitude +
                           lowerjawfork_length:trophic_level,
                         data = database_Se, method = "REML")
summary(mod_gam_Se7) # d13C:long

mod_gam_Se8 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                           s(d13C) + 
                           s(trophic_level) + 
                           s(latitude) +
                           s(longitude)+
                           s(log10(THg))+
                           d13C:longitude +
                           trophic_level:latitude +
                           lowerjawfork_length:trophic_level,
                         data = database_Se, method = "REML")
summary(mod_gam_Se8) # d13C:lat

mod_gam_Se9 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                           s(d13C) + 
                           s(trophic_level) + 
                           s(latitude) +
                           s(longitude)+
                           s(log10(THg))+
                           d13C:longitude +
                           d13C:latitude +
                           lowerjawfork_length:trophic_level,
                         data = database_Se, method = "REML")
summary(mod_gam_Se9) # TP:lat

mod_gam_Se10 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            s(log10(THg))+
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude,
                          data = database_Se, method = "REML")
summary(mod_gam_Se10) # length:TP

mod_gam_Se11 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(longitude)+
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude,
                          data = database_Se, method = "REML")
summary(mod_gam_Se11) # log10(Hg)

AIC(mod_gam_Se0,mod_gam_Se1,mod_gam_Se2,mod_gam_Se3,mod_gam_Se4,mod_gam_Se5,mod_gam_Se6,
    mod_gam_Se7,mod_gam_Se8,mod_gam_Se9,mod_gam_Se10,mod_gam_Se11)

## Step 2

# Begin from n°6
mod_gam_Se12 <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(log10(THg))+
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_Se, method = "REML")
summary(mod_gam_Se12) # Length

mod_gam_Se13 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                            s(trophic_level) + 
                            s(latitude) +
                            s(log10(THg))+
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_Se, method = "REML")
summary(mod_gam_Se13) # d13C

mod_gam_Se14 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(latitude) +
                            s(log10(THg))+
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_Se, method = "REML")
summary(mod_gam_Se14) # TP

mod_gam_Se15 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(log10(THg))+
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_Se, method = "REML")
summary(mod_gam_Se15) # lat

mod_gam_Se16 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_Se, method = "REML")
summary(mod_gam_Se16) # log10(THg)

mod_gam_Se17 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(log10(THg))+
                            d13C:latitude +
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_Se, method = "REML")
summary(mod_gam_Se17) # d13C:long

mod_gam_Se18 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(log10(THg))+
                            d13C:longitude +
                            trophic_level:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_Se, method = "REML")
summary(mod_gam_Se18) # d13C:lat

mod_gam_Se19 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(log10(THg))+
                            d13C:longitude +
                            d13C:latitude +
                            lowerjawfork_length:trophic_level,
                          data = database_Se, method = "REML")
summary(mod_gam_Se19) # TP:lat

mod_gam_Se20 <- mgcv::gam(log10(Seww) ~ s(lowerjawfork_length)+
                            s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(log10(THg))+
                            d13C:longitude +
                            d13C:latitude +
                            trophic_level:latitude,
                          data = database_Se, method = "REML")
summary(mod_gam_Se20) # length:TP

AIC(mod_gam_Se0,mod_gam_Se6,mod_gam_Se12,mod_gam_Se13,mod_gam_Se14,
    mod_gam_Se15,mod_gam_Se16,mod_gam_Se17,mod_gam_Se18,mod_gam_Se19,
    mod_gam_Se20)

## Step 3

# Begin from 18, remove length (mod 12), TP:lat (mod 19), TP:length (mod 20)
mod_gam_Se18bis <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                               s(trophic_level) + 
                               s(latitude) +
                               s(log10(THg))+
                               d13C:longitude,
                             data = database_Se, method = "REML")
summary(mod_gam_Se18bis)

AIC(mod_gam_Se0,mod_gam_Se6,mod_gam_Se18,mod_gam_Se18bis)

# Begin from 18bis
mod_gam_Se21 <- mgcv::gam(log10(Seww) ~ s(trophic_level) + 
                            s(latitude) +
                            s(log10(THg))+
                            d13C:longitude,
                          data = database_Se, method = "REML")
summary(mod_gam_Se21) # d13C

mod_gam_Se22 <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                            s(latitude) +
                            s(log10(THg))+
                            d13C:longitude,
                          data = database_Se, method = "REML")
summary(mod_gam_Se22) # TP

mod_gam_Se23 <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(log10(THg))+
                            d13C:longitude,
                          data = database_Se, method = "REML")
summary(mod_gam_Se23) # latitude

mod_gam_Se24 <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            d13C:longitude,
                          data = database_Se, method = "REML")
summary(mod_gam_Se24) # log10(THg)

mod_gam_Se25 <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                            s(trophic_level) + 
                            s(latitude) +
                            s(log10(THg)),
                          data = database_Se, method = "REML")
summary(mod_gam_Se25) # d13C:longitude

AIC(mod_gam_Se0,mod_gam_Se6,mod_gam_Se18,mod_gam_Se18bis,mod_gam_Se21,
    mod_gam_Se22,mod_gam_Se23,mod_gam_Se24,mod_gam_Se25)

# Best model is 18bis

rm(mod_gam_Se0,mod_gam_Se1,mod_gam_Se2,mod_gam_Se3,mod_gam_Se4,mod_gam_Se5,mod_gam_Se6,
   mod_gam_Se7,mod_gam_Se8,mod_gam_Se9,mod_gam_Se10,mod_gam_Se11,
   mod_gam_Se12,mod_gam_Se13,mod_gam_Se14,
   mod_gam_Se15,mod_gam_Se16,mod_gam_Se17,mod_gam_Se18,mod_gam_Se19,
   mod_gam_Se20,
   mod_gam_Se18bis,mod_gam_Se21,
   mod_gam_Se22,mod_gam_Se23,mod_gam_Se24,mod_gam_Se25)

## Plotting GAM results
mod_gam_Se18bis <- mgcv::gam(log10(Seww) ~ s(d13C) + 
                               s(trophic_level) + 
                               s(latitude) +
                               s(log10(THg))+
                               d13C:longitude,
                             data = database_Se, method = "REML")
summary(mod_gam_Se18bis)

par(mfrow=c(2,2))
plot(mod_gam_Se18bis, se = TRUE, ylab = "", xlab = "", las = 1)

rm(database_Se,mod_gam_Se18bis)



### VI // Plot molar Hg and Se concentrations, MHg:MSe, HBVSe and theoretically available Se for each area ##############################################################################################

## 1 / Mapping molar Hg and Se concentrations by area ##############################################################################################

world <- ne_countries(scale = "medium", returnclass = "sf", type = "map_units")
class(world)

Coords_lim <- matrix(data = c(min(SWOIO_data$latitude), max(SWOIO_data$latitude),
                              min(SWOIO_data$longitude), max(SWOIO_data$longitude)),
                     nrow = 2, ncol = 2)
colnames(Coords_lim) <- c("latitude", "longitude")
row.names(Coords_lim) <- c("min", "max")

data_pie <- SWOIO_data %>% 
  select(area, latitude, longitude, MHg, MSe) %>% 
  mutate(total = MHg+MSe,
         percent_Hg = MHg*100/total,
         percent_Se = MSe*100/total) %>% 
  group_by(area) %>% 
  summarise(tot = mean(total, na.rm = TRUE),
            lat = mean(latitude, na.rm = TRUE),
            lon = mean(longitude, na.rm = TRUE),
            Se = mean(percent_Se, na.rm = TRUE),
            Hg = mean(percent_Hg, na.rm = TRUE)) %>% 
  mutate(tot = tot*100000000,
         tot = round(tot, 0))

map <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(floor(Coords_lim[1,2]), ceiling(Coords_lim[2,2])),
           ylim = c(floor(Coords_lim[1,1]), ceiling(Coords_lim[2,1]))) +
  geom_scatterpie(data = data_pie, aes(x = lon, y = lat, group = area,
                                       r = tot*2), cols = c("Se","Hg"))+
  geom_scatterpie_legend(data_pie$tot*2, x=80, y=-25)+
  theme_bw()+
  theme(axis.title = element_blank(),
        legend.position = "none") +
  labs(x = 'Longitude', y = "Latitude")+
  scale_fill_manual(values = c("#D1D2D4","#404041"))


## 2 / Plot MHg:MSe ##############################################################################################

MHg_MSe <- SWOIO_data %>%
  group_by(area) %>%
  mutate(mean = mean(MHg_MSe, na.rm = TRUE)) %>% 
  ungroup() %>% 
  ggplot(aes(x = area))+
  geom_boxplot(aes(y = MHg_MSe, fill = area))+
  geom_point(aes(y = mean), shape = 21, fill = "white", color = "black", size = 2.5)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey")+
  scale_fill_manual(values = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"))+
  labs(y = "MHg:MSe")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title.y = element_text(face = "italic"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


## 3 / Plot HBVSe ##############################################################################################

HBVSe <- SWOIO_data %>%
  group_by(area) %>%
  mutate(mean = mean(HBVSe, na.rm = TRUE)) %>% 
  ungroup() %>% 
  ggplot(aes(x = area))+
  geom_boxplot(aes(y = HBVSe, fill = area))+
  geom_point(aes(y = mean), shape = 21, fill = "white", color = "black", size = 2.5)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey")+
  scale_fill_manual(values = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"))+
  labs(y = "HBVSe")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title.y = element_text(face = "italic"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


## 4 / Plot theoretically available Se ##############################################################################################

# Statistical tests
a1 <- aov(SWOIO_data$Se_avail ~ SWOIO_data$area)
shapiro.test(resid(a1)) # Normality of residuals
fligner.test(data = SWOIO_data, Se_avail ~ area) # Homoscedasticity

kruskal.test(data = SWOIO_data, Se_avail ~ area)
dunnTest(data = SWOIO_data, Se_avail ~ area, method = "bh")

rm(a1)

# Plot theoretically available Se
Se_avail <- SWOIO_data %>%
  mutate(letter = "bc",
         letter = ifelse(area == "BENG", "a",letter),
         letter = ifelse(area == "WTIO", "b",letter),
         letter = ifelse(area == "MCSA", "c",letter)) %>%
  group_by(area) %>%
  mutate(mean = mean(Se_avail, na.rm = TRUE),
         y_letter = max(Se_avail, na.rm = TRUE)) %>% 
  ungroup() %>% 
  ggplot(aes(x = area))+
  geom_boxplot(aes(y = Se_avail, fill = area))+
  geom_point(aes(y = mean), shape = 21, fill = "white", color = "black", size = 2.5)+
  geom_text(aes(y = y_letter+0.2, vjust = 0, label = letter), color = "#58585B")+
  scale_fill_manual(values = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"))+
  labs(y = "Theoretically available Se (µg.g-1 ww)")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title.y = element_text(face = "italic"),
        axis.text.x = element_text(color = c("#D2062E","#048C7F","#88BB9A","#FBAB19","#ED7340","#354B5E"),
                                   angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank())


## 5 / Plot all on same graph ##############################################################################################

ggarrange(map, ggarrange(MHg_MSe, HBVSe, Se_avail,
                         ncol = 1, nrow = 3, align = "hv",
                         labels = c("B.","C.","D.")),
          ncol = 2, align = "hv", labels = c("A.",""),
          widths = c(2,1))

rm(world,Coords_lim,data_pie,map)
rm(MHg_MSe, HBVSe, Se_avail)



### VII // Calculation of max number of servings and/or contribution to RDI ##############################################################################################

## 1 / For Hg ##############################################################################################

servings <- SWOIO_data %>% 
  select(area, THg) %>% 
  mutate(# THg in one serving
    THg_YC = THg*25,
    THg_C = THg*50,
    THg_ADO = THg*75,
    THg_YA_A = THg*100,
    # Number of servings before reaching PTI
    serv_THg_YC = 48/THg_YC,
    serv_THg_C = 88/THg_C,
    serv_THg_ADO = 144/THg_ADO,
    serv_THg_YA = 208/THg_YA_A,
    serv_THg_A = 280/THg_YA_A)

# Calculating for Children
servings_YC <- servings %>% select(area,serv_THg_YC) %>% rename(serv_THg_C = serv_THg_YC)
servings_C <- servings %>% select(area,serv_THg_C)
servings_ADO <- servings %>% select(area,serv_THg_ADO) %>% rename(serv_THg_C = serv_THg_ADO)

servings_C <- rbind(servings_YC,servings_C,servings_ADO) %>% 
  group_by(area) %>% 
  summarise(mean_serv_THg_C = mean(serv_THg_C),
            se_serv_THg_C = std(serv_THg_C)) %>% 
  mutate(CIinf_serv_THg_C = mean_serv_THg_C-(1.96*se_serv_THg_C),
         CIsup_serv_THg_C = mean_serv_THg_C+(1.96*se_serv_THg_C),
         serv_THg_C = paste0(round(mean_serv_THg_C,0)," (",round(CIinf_serv_THg_C,0),"-",round(CIsup_serv_THg_C,0),")")) %>% 
  select(area,serv_THg_C)

rm(servings_YC,servings_ADO)

# Calculating for Young adults and Adults
servings_THg <- servings %>% 
  group_by(area) %>% 
  summarise(mean_serv_THg_YA = mean(serv_THg_YA),
            se_serv_THg_YA = std(serv_THg_YA),
            mean_serv_THg_A = mean(serv_THg_A),
            se_serv_THg_A = std(serv_THg_A)) %>% 
  mutate(CIinf_serv_THg_YA = mean_serv_THg_YA-(1.96*se_serv_THg_YA),
         CIsup_serv_THg_YA = mean_serv_THg_YA+(1.96*se_serv_THg_YA),
         CIinf_serv_THg_A = mean_serv_THg_A-(1.96*se_serv_THg_A),
         CIsup_serv_THg_A = mean_serv_THg_A+(1.96*se_serv_THg_A),
         serv_THg_YA = paste0(round(mean_serv_THg_YA,0)," (",round(CIinf_serv_THg_YA,0),"-",round(CIsup_serv_THg_YA,0),")"),
         serv_THg_A = paste0(round(mean_serv_THg_A,0)," (",round(CIinf_serv_THg_A,0),"-",round(CIsup_serv_THg_A,0),")")) %>% 
  left_join(servings_C, by = "area") %>% 
  select(area,serv_THg_C,serv_THg_YA,serv_THg_A)

rm(servings_C,servings)


## 2 / For Se ##############################################################################################

# %RDI
percent_RDI <- SWOIO_data %>% 
  select(area, Seww) %>% 
  mutate(# Se in one serving
    Se_YC = Seww*25,
    Se_C = Seww*50,
    Se_ADO = Seww*75,
    Se_YA_A = Seww*100,
    # Number of servings before reaching PTI
    RDI_Se_YC = Se_YC*100/20,
    RDI_Se_C = Se_C*100/30,
    RDI_Se_ADO = Se_ADO*100/40,
    RDI_Se_YA = Se_YA_A*100/55,
    RDI_Se_A = Se_YA_A*100/55,
    RDI_Se_Apre = Se_YA_A*100/60,
    RDI_Se_Alact = Se_YA_A*100/70)

# Calculating for Children
percent_RDI_YC <- percent_RDI %>% select(area,RDI_Se_YC) %>% rename(RDI_Se_C = RDI_Se_YC)
percent_RDI_C <- percent_RDI %>% select(area,RDI_Se_C)
percent_RDI_ADO <- percent_RDI %>% select(area,RDI_Se_ADO) %>% rename(RDI_Se_C = RDI_Se_ADO)

percent_RDI_C <- rbind(percent_RDI_YC,percent_RDI_C,percent_RDI_ADO) %>% 
  group_by(area) %>% 
  summarise(mean_RDI_Se_C = mean(RDI_Se_C),
            se_RDI_Se_C = std(RDI_Se_C)) %>% 
  mutate(CIinf_RDI_Se_C = mean_RDI_Se_C-(1.96*se_RDI_Se_C),
         CIsup_RDI_Se_C = mean_RDI_Se_C+(1.96*se_RDI_Se_C),
         RDI_Se_C = paste0(round(mean_RDI_Se_C,0)," (",round(CIinf_RDI_Se_C,0),"-",round(CIsup_RDI_Se_C,0),")")) %>% 
  select(area,RDI_Se_C)

rm(percent_RDI_YC,percent_RDI_ADO)

# Calculation for Adult women
percent_RDI_Apre <- percent_RDI %>% select(area,RDI_Se_Apre) %>% rename(RDI_Se_A = RDI_Se_Apre)
percent_RDI_A <- percent_RDI %>% select(area,RDI_Se_A)
percent_RDI_Alact <- percent_RDI %>% select(area,RDI_Se_Alact) %>% rename(RDI_Se_A = RDI_Se_Alact)

percent_RDI_A <- rbind(percent_RDI_Apre,percent_RDI_A,percent_RDI_Alact) %>% 
  group_by(area) %>% 
  summarise(mean_RDI_Se_A = mean(RDI_Se_A),
            se_RDI_Se_A = std(RDI_Se_A)) %>% 
  mutate(CIinf_RDI_Se_A = mean_RDI_Se_A-(1.96*se_RDI_Se_A),
         CIsup_RDI_Se_A = mean_RDI_Se_A+(1.96*se_RDI_Se_A),
         RDI_Se_A = paste0(round(mean_RDI_Se_A,0)," (",round(CIinf_RDI_Se_A,0),"-",round(CIsup_RDI_Se_A,0),")")) %>% 
  select(area,RDI_Se_A)

rm(percent_RDI_Apre,percent_RDI_Alact)

# Calculating for Young adults
percent_RDI <- percent_RDI %>% 
  group_by(area) %>% 
  summarise(mean_RDI_Se_YA = mean(RDI_Se_YA),
            se_RDI_Se_YA = std(RDI_Se_YA)) %>% 
  mutate(CIinf_RDI_Se_YA = mean_RDI_Se_YA-(1.96*se_RDI_Se_YA),
         CIsup_RDI_Se_YA = mean_RDI_Se_YA+(1.96*se_RDI_Se_YA),
         RDI_Se_YA = paste0(round(mean_RDI_Se_YA,0)," (",round(CIinf_RDI_Se_YA,0),"-",round(CIsup_RDI_Se_YA,0),")")) %>% 
  left_join(percent_RDI_C, by = "area") %>% 
  left_join(percent_RDI_A, by = "area") %>% 
  select(area,RDI_Se_C,RDI_Se_YA,RDI_Se_A)

rm(percent_RDI_C,percent_RDI_A)

## Servings
servings <- SWOIO_data %>% 
  select(area, Seww) %>% 
  mutate(# Se in one serving
    Se_YC = Seww*25,
    Se_C = Seww*50,
    Se_ADO = Seww*75,
    Se_YA_A = Seww*100,
    # Number of servings before reaching PTI
    serv_Se_YC = 90/Se_YC,
    serv_Se_C = 150/Se_C,
    serv_Se_ADO = 280/Se_ADO,
    serv_Se_YA = 400/Se_YA_A,
    serv_Se_A = 400/Se_YA_A)

# Calculating for Children
servings_YC <- servings %>% select(area,serv_Se_YC) %>% rename(serv_Se_C = serv_Se_YC)
servings_C <- servings %>% select(area,serv_Se_C)
servings_ADO <- servings %>% select(area,serv_Se_ADO) %>% rename(serv_Se_C = serv_Se_ADO)

servings_C <- rbind(servings_YC,servings_C,servings_ADO) %>% 
  group_by(area) %>% 
  summarise(mean_serv_Se_C = mean(serv_Se_C),
            se_serv_Se_C = std(serv_Se_C)) %>% 
  mutate(CIinf_serv_Se_C = mean_serv_Se_C-(1.96*se_serv_Se_C),
         CIsup_serv_Se_C = mean_serv_Se_C+(1.96*se_serv_Se_C),
         serv_Se_C = paste0(round(mean_serv_Se_C,0)," (",round(CIinf_serv_Se_C,0),"-",round(CIsup_serv_Se_C,0),")")) %>% 
  select(area,serv_Se_C)

rm(servings_YC,servings_ADO)

# Calculating for Young adults and Adults
servings <- servings %>% 
  group_by(area) %>% 
  summarise(mean_serv_Se_YA = mean(serv_Se_YA),
            se_serv_Se_YA = std(serv_Se_YA),
            mean_serv_Se_A = mean(serv_Se_A),
            se_serv_Se_A = std(serv_Se_A)) %>% 
  mutate(CIinf_serv_Se_YA = mean_serv_Se_YA-(1.96*se_serv_Se_YA),
         CIsup_serv_Se_YA = mean_serv_Se_YA+(1.96*se_serv_Se_YA),
         CIinf_serv_Se_A = mean_serv_Se_A-(1.96*se_serv_Se_A),
         CIsup_serv_Se_A = mean_serv_Se_A+(1.96*se_serv_Se_A),
         serv_Se_YA = paste0(round(mean_serv_Se_YA,0)," (",round(CIinf_serv_Se_YA,0),"-",round(CIsup_serv_Se_YA,0),")"),
         serv_Se_A = paste0(round(mean_serv_Se_A,0)," (",round(CIinf_serv_Se_A,0),"-",round(CIsup_serv_Se_A,0),")")) %>% 
  left_join(servings_C, by = "area") %>% 
  select(area,serv_Se_C,serv_Se_YA,serv_Se_A)

rm(servings_C)

## Both %RDI and servings
Se_percentRDI_servings <- merge(percent_RDI,servings,by = "area") %>% 
  mutate(area = factor(area, levels = c("BENG","ISLU","WTIO","MOZ","SSG","MCSA")))

rm(percent_RDI,servings)


## 3 / For theoretically available Se ##############################################################################################

# %RDI
percent_RDI <- SWOIO_data %>% 
  select(area, Se_avail) %>% 
  mutate(# Se in one serving
    Se_YC = Se_avail*25,
    Se_C = Se_avail*50,
    Se_ADO = Se_avail*75,
    Se_YA_A = Se_avail*100,
    # Number of servings before reaching PTI
    RDI_Se_YC = Se_YC*100/20,
    RDI_Se_C = Se_C*100/30,
    RDI_Se_ADO = Se_ADO*100/40,
    RDI_Se_YA = Se_YA_A*100/55,
    RDI_Se_A = Se_YA_A*100/55,
    RDI_Se_Apre = Se_YA_A*100/60,
    RDI_Se_Alact = Se_YA_A*100/70)

# Calculating for Children
percent_RDI_YC <- percent_RDI %>% select(area,RDI_Se_YC) %>% rename(RDI_Se_C = RDI_Se_YC)
percent_RDI_C <- percent_RDI %>% select(area,RDI_Se_C)
percent_RDI_ADO <- percent_RDI %>% select(area,RDI_Se_ADO) %>% rename(RDI_Se_C = RDI_Se_ADO)

percent_RDI_C <- rbind(percent_RDI_YC,percent_RDI_C,percent_RDI_ADO) %>% 
  group_by(area) %>% 
  summarise(mean_RDI_Se_C = mean(RDI_Se_C),
            se_RDI_Se_C = std(RDI_Se_C)) %>% 
  mutate(CIinf_RDI_Se_C = mean_RDI_Se_C-(1.96*se_RDI_Se_C),
         CIsup_RDI_Se_C = mean_RDI_Se_C+(1.96*se_RDI_Se_C),
         RDI_Se_C = paste0(round(mean_RDI_Se_C,0)," (",round(CIinf_RDI_Se_C,0),"-",round(CIsup_RDI_Se_C,0),")")) %>% 
  select(area,RDI_Se_C)

rm(percent_RDI_YC,percent_RDI_ADO)

# Calculation for Adult women
percent_RDI_Apre <- percent_RDI %>% select(area,RDI_Se_Apre) %>% rename(RDI_Se_A = RDI_Se_Apre)
percent_RDI_A <- percent_RDI %>% select(area,RDI_Se_A)
percent_RDI_Alact <- percent_RDI %>% select(area,RDI_Se_Alact) %>% rename(RDI_Se_A = RDI_Se_Alact)

percent_RDI_A <- rbind(percent_RDI_Apre,percent_RDI_A,percent_RDI_Alact) %>% 
  group_by(area) %>% 
  summarise(mean_RDI_Se_A = mean(RDI_Se_A),
            se_RDI_Se_A = std(RDI_Se_A)) %>% 
  mutate(CIinf_RDI_Se_A = mean_RDI_Se_A-(1.96*se_RDI_Se_A),
         CIsup_RDI_Se_A = mean_RDI_Se_A+(1.96*se_RDI_Se_A),
         RDI_Se_A = paste0(round(mean_RDI_Se_A,0)," (",round(CIinf_RDI_Se_A,0),"-",round(CIsup_RDI_Se_A,0),")")) %>% 
  select(area,RDI_Se_A)

rm(percent_RDI_Apre,percent_RDI_Alact)

# Calculating for Young adults
percent_RDI <- percent_RDI %>% 
  group_by(area) %>% 
  summarise(mean_RDI_Se_YA = mean(RDI_Se_YA),
            se_RDI_Se_YA = std(RDI_Se_YA)) %>% 
  mutate(CIinf_RDI_Se_YA = mean_RDI_Se_YA-(1.96*se_RDI_Se_YA),
         CIsup_RDI_Se_YA = mean_RDI_Se_YA+(1.96*se_RDI_Se_YA),
         RDI_Se_YA = paste0(round(mean_RDI_Se_YA,0)," (",round(CIinf_RDI_Se_YA,0),"-",round(CIsup_RDI_Se_YA,0),")")) %>% 
  left_join(percent_RDI_C, by = "area") %>% 
  left_join(percent_RDI_A, by = "area") %>% 
  select(area,RDI_Se_C,RDI_Se_YA,RDI_Se_A)

rm(percent_RDI_C,percent_RDI_A)

## Servings
servings <- SWOIO_data %>% 
  select(area, Se_avail) %>% 
  mutate(# Se in one serving
    Se_YC = Se_avail*25,
    Se_C = Se_avail*50,
    Se_ADO = Se_avail*75,
    Se_YA_A = Se_avail*100,
    # Number of servings before reaching PTI
    serv_Se_YC = 90/Se_YC,
    serv_Se_C = 150/Se_C,
    serv_Se_ADO = 280/Se_ADO,
    serv_Se_YA = 400/Se_YA_A,
    serv_Se_A = 400/Se_YA_A)

# Calculating for Children
servings_YC <- servings %>% select(area,serv_Se_YC) %>% rename(serv_Se_C = serv_Se_YC)
servings_C <- servings %>% select(area,serv_Se_C)
servings_ADO <- servings %>% select(area,serv_Se_ADO) %>% rename(serv_Se_C = serv_Se_ADO)

servings_C <- rbind(servings_YC,servings_C,servings_ADO) %>% 
  filter(!serv_Se_C < 0,
         !serv_Se_C > 100) %>% 
  group_by(area) %>% 
  summarise(mean_serv_Se_C = mean(serv_Se_C),
            se_serv_Se_C = std(serv_Se_C)) %>% 
  mutate(CIinf_serv_Se_C = mean_serv_Se_C-(1.96*se_serv_Se_C),
         CIsup_serv_Se_C = mean_serv_Se_C+(1.96*se_serv_Se_C),
         serv_Se_C = paste0(round(mean_serv_Se_C,0)," (",round(CIinf_serv_Se_C,0),"-",round(CIsup_serv_Se_C,0),")")) %>% 
  select(area,serv_Se_C)

rm(servings_YC,servings_ADO)

# Calculating for Young adults and Adults
servings <- servings %>% 
  filter(!serv_Se_YA < 0,
         !serv_Se_YA > 100) %>% 
  group_by(area) %>% 
  summarise(mean_serv_Se_YA = mean(serv_Se_YA),
            se_serv_Se_YA = std(serv_Se_YA),
            mean_serv_Se_A = mean(serv_Se_A),
            se_serv_Se_A = std(serv_Se_A)) %>% 
  mutate(CIinf_serv_Se_YA = mean_serv_Se_YA-(1.96*se_serv_Se_YA),
         CIsup_serv_Se_YA = mean_serv_Se_YA+(1.96*se_serv_Se_YA),
         CIinf_serv_Se_A = mean_serv_Se_A-(1.96*se_serv_Se_A),
         CIsup_serv_Se_A = mean_serv_Se_A+(1.96*se_serv_Se_A),
         serv_Se_YA = paste0(round(mean_serv_Se_YA,0)," (",round(CIinf_serv_Se_YA,0),"-",round(CIsup_serv_Se_YA,0),")"),
         serv_Se_A = paste0(round(mean_serv_Se_A,0)," (",round(CIinf_serv_Se_A,0),"-",round(CIsup_serv_Se_A,0),")")) %>% 
  left_join(servings_C, by = "area") %>% 
  select(area,serv_Se_C,serv_Se_YA,serv_Se_A)

rm(servings_C)

## Both %RDI and servings
Se_avail_percentRDI_servings <- merge(percent_RDI,servings,by = "area") %>% 
  mutate(area = factor(area, levels = c("BENG","ISLU","WTIO","MOZ","SSG","MCSA"))) %>% 
  rename(RDI_Seavail_C = RDI_Se_C,
         RDI_Seavail_YA = RDI_Se_YA,
         RDI_Seavail_A = RDI_Se_A,
         serv_Seavail_C = serv_Se_C,
         serv_Seavail_YA = serv_Se_YA,
         serv_Seavail_A = serv_Se_A)

rm(percent_RDI,servings)


## 4 / Table with all data ##############################################################################################

PTI_RDI <- servings_THg %>% 
  left_join(Se_percentRDI_servings, by = "area") %>% 
  left_join(Se_avail_percentRDI_servings, by = "area")

write.csv(PTI_RDI,file = "C:/Users/msabin02/Nextcloud/PhD/R_projects/Final_scripts/SWOIO_Hg_Se_RDI_PTI.csv")

rm(servings_THg,Se_percentRDI_servings,Se_avail_percentRDI_servings)
rm(PTI_RDI)

rm(SWOIO_data)