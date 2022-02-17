##--------------------------------------------------------------------------------------------------------
## SCRIPT : Allows to reproduce all analyses of Chapter 3 dedicated to trace element bioaccumulation
##          in Emperor red snapper (Lutjanus sebae)
##          This script includes the following analyses :
##            - Correlation tests to test for relationship between size and SI, and among SI
##            - Statistical tests to test for significant difference in d13C or d15N between sexes
##              and among seasons
##            - nMDS to visualise TE profiles for each sex and during each season
##            - Plot TE concentrations (mean SD) in each sex and during each season
##            - Statistical tests to test for significant difference in TE concentrations between
##              sexes and among seasons
##            - GLM and GAM models to test for relationships between TE concentrations and fork
##              length, d13C and d15N values + plotting GAM results
##            - Plotting mean fork length for each sex and each season
##            - Determining size at shift in trend in d13C/d15N relationship
##
## As part of :
##        Magali SABINO PhD - "Bioaccumulation of trace elements in Seychelles marine food webs"
##
## Author : Magali Sabino
## First created : 2022-01-25
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
lapply(c("tidyverse", "mgcv", "vegan", "nicheROVER", "ggpubr","FSA"),
       library, character.only = TRUE)



### I // Creation of database with only emperor red snapper samples ##############################################################################################

LUB_data <- data_TE_SI %>% 
  filter(c_sp_fao == "LUB") %>% 
  select(organism_identifier,length,sex,season,d13C,d15N,Ag_stat,Cd_stat,Co_stat,
         Cr_stat,Cu_stat,Fe_stat,Mn_stat,Ni_stat,Pb_stat,Se_stat,TAs_stat,
         Zn_stat,THg_stat) %>% 
  rename(Ag = Ag_stat, Cd = Cd_stat, Co = Co_stat, Cr = Cr_stat,
         Cu = Cu_stat, Fe = Fe_stat, Mn = Mn_stat, Ni = Ni_stat,
         Pb = Pb_stat, Se = Se_stat, As = TAs_stat, Zn = Zn_stat,
         Hg = THg_stat)



### II // Stable isotopes: intercorrelation, relationship with size and effect of sex and season ##############################################################################################

## 1 / Intercorrelation and relationship with size ##############################################################################################

# Calculation of correlations - For males only
database_test_M <- LUB_data %>% 
  filter(sex %in% c("M"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N)

shapiro.test(database_test_M$d13C)
shapiro.test(database_test_M$d15N)
shapiro.test(database_test_M$length)

cor.test(database_test_M$d13C, database_test_M$length, method = "pearson")
cor.test(database_test_M$d15N, database_test_M$length, method = "kendall")
cor.test(database_test_M$d15N, database_test_M$d13C, method = "pearson")

rm(database_test_M)

# Calculation of correlations - For females only
database_test_F <- LUB_data %>% 
  filter(sex %in% c("F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N)

shapiro.test(database_test_M$d13C)
shapiro.test(database_test_M$d15N)
shapiro.test(database_test_M$length)

cor.test(database_test_M$d13C, database_test_M$length, method = "pearson")
cor.test(database_test_M$d15N, database_test_M$length, method = "kendall")
cor.test(database_test_M$d15N, database_test_M$d13C, method = "pearson")

rm(database_test_F)

# Plot all relationships
length_d13C <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(sex = factor(sex, levels = c("M","F"))) %>% 
  select(length,d13C,d15N, sex) %>% 
  ggplot(aes(x = length, y = d13C, color = sex))+
  geom_point(alpha = 0.3, size = 3)+
  geom_smooth(method = "glm", se = TRUE, aes(color = sex))+
  scale_color_manual(values = c("#257898","#E05F25"))+
  labs(x = "Fork length (cm)")+
  theme_bw()+
  theme(axis.title = element_text(face = "italic"),
        legend.position = "none",
        panel.grid = element_blank())

length_d15N <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(sex = factor(sex, levels = c("M","F"))) %>% 
  select(length,d13C,d15N, sex) %>% 
  ggplot(aes(x = length, y = d15N, color = sex))+
  geom_point(alpha = 0.3, size = 3)+
  geom_smooth(method = "glm", se = TRUE, aes(color = sex))+
  scale_color_manual(values = c("#257898","#E05F25"))+
  labs(x = "Fork length (cm)")+
  theme_bw()+
  theme(axis.title = element_text(face = "italic"),
        legend.position = "none",
        panel.grid = element_blank())


d13C_d15N <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(sex = factor(sex, levels = c("M","F"))) %>% 
  select(sex,length,d13C,d15N) %>% 
  ggplot(aes(x = d13C, y = d15N, color = sex))+
  geom_point(alpha = 0.3, size = 3)+
  geom_smooth(method = "gam", se = TRUE, aes(color = sex))+
  scale_color_manual(values = c("#257898","#E05F25"))+
  theme_bw()+
  theme(axis.title = element_text(face = "italic"),
        legend.position = "none",
        panel.grid = element_blank())

ggarrange(length_d13C, length_d15N, d13C_d15N,
          ncol = 3, align = "hv", labels = c("A.","B.","C."))

rm(length_d13C,length_d15N,d13C_d15N)


## 2 / Difference in d13C and d15N mean value between sexes and between seasons ##############################################################################################

# Statistical tests - BETWEEN MALES AND FEMALES
Tracer_data <- as.data.frame(LUB_data %>% 
  select(sex,d13C,d15N) %>% 
  filter(sex %in% c("M","F"),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(sex = factor(sex, levels = c("M","F"))))

Output_ALL_test <- NULL

for (i in 2: length(colnames(Tracer_data))){ # metal loop
  #i = 2
  
  Tracer_name <- colnames(Tracer_data)[i]
  Tracer_data_F1 <- subset(Tracer_data, sex == "M")
  Tracer_data_F1 <- Tracer_data_F1[,i]
  Tracer_data_F1 <- Tracer_data_F1[!is.na(Tracer_data_F1)]
  Tracer_data_F2 <- subset(Tracer_data, sex == "F")
  Tracer_data_F2 <- Tracer_data_F2[,i]
  Tracer_data_F2 <- Tracer_data_F2[!is.na(Tracer_data_F2)]
  
  test_a <- shapiro.test(Tracer_data_F1) # Normality
  test_b <- shapiro.test(Tracer_data_F2) # Normality
  test_c <- fligner.test(Tracer_data[,i] ~ Tracer_data[,1]) # Homoscedasticity
  
  Output_test <- data.frame(Trace_element = Tracer_name,
                            N1 = length(Tracer_data_F1),
                            N2 = length(Tracer_data_F2),
                            Test = NA,
                            P_val = NA,
                            stat_val = NA)
  
  
  if (Output_test$N1 <= 2 | Output_test$N2 <= 2){ # pas de test pour 5 ech ou moins
    
    Output_ALL_test <- rbind(Output_ALL_test, Output_test)
    
  }else{
    
    if (Output_test$N1 <= 10 | Output_test$N2 <= 10 | test_a$p.value < 0.05 | test_b$p.value < 0.05 | test_c$p.value < 0.05){ # non paramétrique
      # Wilcox
      
      test_stat <- wilcox.test(Tracer_data_F1, Tracer_data_F2, paired = F)
      
      Output_test$Test <- "Wilcoxon"
      Output_test$P_val <- test_stat$p.value
      Output_test$stat_val <- test_stat$statistic
      
    }else{
      # t.test
      test_stat <- t.test(Tracer_data_F1, Tracer_data_F2, paired = F)
      
      Output_test$Test <- "t test"
      Output_test$P_val <- test_stat$p.value
      Output_test$stat_val <- test_stat$statistic
      
    }
    
    Output_ALL_test <- rbind(Output_ALL_test, Output_test)
    
  } # fin test
  
  rm(Output_test, test_stat, Tracer_name, Tracer_data_F2,Tracer_data_F1,test_a, test_b, test_c)
  
}

Output_ALL_test[,5] <- round(Output_ALL_test[,5], 3)
rm(Tracer_data,Output_ALL_test,i)


# Statistical tests - AMONG SEASONS, MALES ONLY
Tracer_data <- as.data.frame(LUB_data %>% 
  filter(sex %in% c("M"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(season = factor(season, levels = c("IMS","NWM","IMA"))) %>% 
  select(season,d13C,d15N))

Output_ALL_test <- NULL

for (i in 2: length(colnames(Tracer_data))){ # metal loop
  #i = 3
  Tracer_name <- colnames(Tracer_data)[i]
  IMS <- subset(Tracer_data, season == "IMS")
  IMS <- IMS[,i]
  NWM <- subset(Tracer_data, season == "NWM")
  NWM <- NWM[,i]
  IMA <- subset(Tracer_data, season == "IMA")
  IMA <- IMA[,i]
  # SEM <- subset(Tracer_data, season == "SEM")
  # SEM <- SEM[,i]
  
  test_a <- fligner.test(Tracer_data[,i] ~ Tracer_data[,1]) # Homoscedasticity
  a1 <- aov(Tracer_data[,i] ~ Tracer_data[,1])
  test_b <- shapiro.test(resid(a1)) # Normality
  
  Output_test <- data.frame(Tracer = Tracer_name,
                            IMS = length(IMS[!is.na(IMS)]),
                            NWM = length(NWM[!is.na(NWM)]),
                            IMA = length(IMA[!is.na(IMA)]),
                            Test = NA,
                            P_val = NA,
                            stat_val = NA,
                            IMA_IMS = NA,
                            IMA_NWM = NA,
                            IMS_NWM = NA)
  
  rm(IMS,NWM,IMA)
  
  
  if (Output_test$IMS <= 2 | Output_test$NWM <= 2 | Output_test$IMA <= 2){ # pas de test pour 5 ech ou moins
    
    Output_ALL_test <- rbind(Output_ALL_test, Output_test)
    
  }else{
    
    if (Output_test$IMS <= 5 | Output_test$NWM <= 5 | Output_test$IMA <= 5 | test_a$p.value < 0.05 | test_b$p.value < 0.05) { # non paramétrique
      # Kruskal
      
      test_stat <- kruskal.test(Tracer_data[,i] ~ Tracer_data[,1])
      
      Output_test$Test <- "Kruskal"
      Output_test$P_val <- test_stat$p.value
      Output_test$stat_val <- test_stat$statistic
      
      if (Output_test$P_val < 0.05) { # post_hoc
        post_hoc <- dunnTest(Tracer_data[,i] ~ Tracer_data[,1], method = "bh")
        Output_test$IMA_IMS <- post_hoc[[2]][1,4]
        Output_test$IMA_NWM <- post_hoc[[2]][2,4]
        Output_test$IMS_NWM <- post_hoc[[2]][3,4]
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
        rm(post_hoc)
      }else{
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
      }
      
    }else{
      # test normalité
      
      lmFA <- lm(Tracer_data[,i] ~ Tracer_data$season)
      test_stat <- anova(lmFA)
      rm(lmFA)
      
      Output_test$Test <- "ANOVA"
      Output_test$P_val <- test_stat[1,5]
      Output_test$stat_val <- test_stat[1,4]
      
      if (Output_test$P_val < 0.05) {
        post_hoc <- TukeyHSD(a1, 'Tracer_data$season', conf.level = 0.95)
        Output_test$IMA_IMS <- post_hoc[[2]][1,4]
        Output_test$IMA_NWM <- post_hoc[[2]][2,4]
        Output_test$IMS_NWM <- post_hoc[[2]][3,4]
        Output_test$IMA_SEM <- post_hoc[[2]][4,4]
        Output_test$IMS_SEM <- post_hoc[[2]][5,4]
        Output_test$NWM_SEM <- post_hoc[[2]][6,4]
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
        rm(post_hoc)
      } else {
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
      }
      
    } # fin tests stat
    
    rm(Output_test, test_stat, Tracer_name, test_a, test_b, a1)
    
  } # fin boucle si N suffisant
}

Output_ALL_test[,6:10] <- round(Output_ALL_test[,6:10], 3)

rm(Tracer_data,i)
rm(Output_ALL_test)


# Statistical tests - AMONG SEASONS, FEMALES ONLY
Tracer_data <- as.data.frame(LUB_data %>% 
                               filter(sex %in% c("F"),
                                      !is.na(length),
                                      !is.na(d13C),
                                      !is.na(d15N)) %>% 
                               mutate(season = factor(season, levels = c("IMS","NWM","IMA"))) %>% 
                               select(season,d13C,d15N))

Output_ALL_test <- NULL

for (i in 2: length(colnames(Tracer_data))){ # metal loop
  #i = 3
  Tracer_name <- colnames(Tracer_data)[i]
  IMS <- subset(Tracer_data, season == "IMS")
  IMS <- IMS[,i]
  NWM <- subset(Tracer_data, season == "NWM")
  NWM <- NWM[,i]
  IMA <- subset(Tracer_data, season == "IMA")
  IMA <- IMA[,i]
  # SEM <- subset(Tracer_data, season == "SEM")
  # SEM <- SEM[,i]
  
  test_a <- fligner.test(Tracer_data[,i] ~ Tracer_data[,1]) # Homoscedasticity
  a1 <- aov(Tracer_data[,i] ~ Tracer_data[,1])
  test_b <- shapiro.test(resid(a1)) # Normality
  
  Output_test <- data.frame(Tracer = Tracer_name,
                            IMS = length(IMS[!is.na(IMS)]),
                            NWM = length(NWM[!is.na(NWM)]),
                            IMA = length(IMA[!is.na(IMA)]),
                            Test = NA,
                            P_val = NA,
                            stat_val = NA,
                            IMA_IMS = NA,
                            IMA_NWM = NA,
                            IMS_NWM = NA)
  
  rm(IMS,NWM,IMA)
  
  
  if (Output_test$IMS <= 2 | Output_test$NWM <= 2 | Output_test$IMA <= 2){ # pas de test pour 5 ech ou moins
    
    Output_ALL_test <- rbind(Output_ALL_test, Output_test)
    
  }else{
    
    if (Output_test$IMS <= 5 | Output_test$NWM <= 5 | Output_test$IMA <= 5 | test_a$p.value < 0.05 | test_b$p.value < 0.05) { # non paramétrique
      # Kruskal
      
      test_stat <- kruskal.test(Tracer_data[,i] ~ Tracer_data[,1])
      
      Output_test$Test <- "Kruskal"
      Output_test$P_val <- test_stat$p.value
      Output_test$stat_val <- test_stat$statistic
      
      if (Output_test$P_val < 0.05) { # post_hoc
        post_hoc <- dunnTest(Tracer_data[,i] ~ Tracer_data[,1], method = "bh")
        Output_test$IMA_IMS <- post_hoc[[2]][1,4]
        Output_test$IMA_NWM <- post_hoc[[2]][2,4]
        Output_test$IMS_NWM <- post_hoc[[2]][3,4]
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
        rm(post_hoc)
      }else{
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
      }
      
    }else{
      # test normalité
      
      lmFA <- lm(Tracer_data[,i] ~ Tracer_data$season)
      test_stat <- anova(lmFA)
      rm(lmFA)
      
      Output_test$Test <- "ANOVA"
      Output_test$P_val <- test_stat[1,5]
      Output_test$stat_val <- test_stat[1,4]
      
      if (Output_test$P_val < 0.05) {
        post_hoc <- TukeyHSD(a1, 'Tracer_data$season', conf.level = 0.95)
        Output_test$IMA_IMS <- post_hoc[[2]][1,4]
        Output_test$IMA_NWM <- post_hoc[[2]][2,4]
        Output_test$IMS_NWM <- post_hoc[[2]][3,4]
        Output_test$IMA_SEM <- post_hoc[[2]][4,4]
        Output_test$IMS_SEM <- post_hoc[[2]][5,4]
        Output_test$NWM_SEM <- post_hoc[[2]][6,4]
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
        rm(post_hoc)
      } else {
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
      }
      
    } # fin tests stat
    
    rm(Output_test, test_stat, Tracer_name, test_a, test_b, a1)
    
  } # fin boucle si N suffisant
}

Output_ALL_test[,6:10] <- round(Output_ALL_test[,6:10], 3)

rm(Tracer_data,i)
rm(Output_ALL_test)


# Statistical tests - AMONG SEASONS
Tracer_data <- as.data.frame(LUB_data %>% 
  filter(sex %in% c("F", "M"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(season = factor(season, levels = c("IMS","NWM","IMA"))) %>% 
  select(season,d13C,d15N))

Output_ALL_test <- NULL

for (i in 2: length(colnames(Tracer_data))){ # metal loop
  #i = 3
  Tracer_name <- colnames(Tracer_data)[i]
  IMS <- subset(Tracer_data, season == "IMS")
  IMS <- IMS[,i]
  NWM <- subset(Tracer_data, season == "NWM")
  NWM <- NWM[,i]
  IMA <- subset(Tracer_data, season == "IMA")
  IMA <- IMA[,i]
  # SEM <- subset(Tracer_data, season == "SEM")
  # SEM <- SEM[,i]
  
  test_a <- fligner.test(Tracer_data[,i] ~ Tracer_data[,1]) # Homoscedasticity
  a1 <- aov(Tracer_data[,i] ~ Tracer_data[,1])
  test_b <- shapiro.test(resid(a1)) # Normality
  
  Output_test <- data.frame(Tracer = Tracer_name,
                            IMS = length(IMS[!is.na(IMS)]),
                            NWM = length(NWM[!is.na(NWM)]),
                            IMA = length(IMA[!is.na(IMA)]),
                            Test = NA,
                            P_val = NA,
                            stat_val = NA,
                            IMA_IMS = NA,
                            IMA_NWM = NA,
                            IMS_NWM = NA)
  
  rm(IMS,NWM,IMA)
  
  
  if (Output_test$IMS <= 2 | Output_test$NWM <= 2 | Output_test$IMA <= 2){ # pas de test pour 5 ech ou moins
    
    Output_ALL_test <- rbind(Output_ALL_test, Output_test)
    
  }else{
    
    if (Output_test$IMS <= 5 | Output_test$NWM <= 5 | Output_test$IMA <= 5 | test_a$p.value < 0.05 | test_b$p.value < 0.05) { # non paramétrique
      # Kruskal
      
      test_stat <- kruskal.test(Tracer_data[,i] ~ Tracer_data[,1])
      
      Output_test$Test <- "Kruskal"
      Output_test$P_val <- test_stat$p.value
      Output_test$stat_val <- test_stat$statistic
      
      if (Output_test$P_val < 0.05) { # post_hoc
        post_hoc <- dunnTest(Tracer_data[,i] ~ Tracer_data[,1], method = "bh")
        Output_test$IMA_IMS <- post_hoc[[2]][1,4]
        Output_test$IMA_NWM <- post_hoc[[2]][2,4]
        Output_test$IMS_NWM <- post_hoc[[2]][3,4]
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
        rm(post_hoc)
      }else{
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
      }
      
    }else{
      # test normalité
      
      lmFA <- lm(Tracer_data[,i] ~ Tracer_data$season)
      test_stat <- anova(lmFA)
      rm(lmFA)
      
      Output_test$Test <- "ANOVA"
      Output_test$P_val <- test_stat[1,5]
      Output_test$stat_val <- test_stat[1,4]
      
      if (Output_test$P_val < 0.05) {
        post_hoc <- TukeyHSD(a1, 'Tracer_data$season', conf.level = 0.95)
        Output_test$IMA_IMS <- post_hoc[[2]][1,4]
        Output_test$IMA_NWM <- post_hoc[[2]][2,4]
        Output_test$IMS_NWM <- post_hoc[[2]][3,4]
        Output_test$IMA_SEM <- post_hoc[[2]][4,4]
        Output_test$IMS_SEM <- post_hoc[[2]][5,4]
        Output_test$NWM_SEM <- post_hoc[[2]][6,4]
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
        rm(post_hoc)
      } else {
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
      }
      
    } # fin tests stat
    
    rm(Output_test, test_stat, Tracer_name, test_a, test_b, a1)
    
  } # fin boucle si N suffisant
}

Output_ALL_test[,6:10] <- round(Output_ALL_test[,6:10], 3)

rm(Tracer_data,i)
rm(Output_ALL_test)



### III // Effect of sex and season on TE bioaccumulation ##############################################################################################

## 1 / TE profile ellipses ##############################################################################################

# Data selection
data_niche <- as.data.frame(LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(season = factor(season, levels = c("IMS","NWM","IMA")),
         sex = factor(sex, levels = c("M","F"))) %>% 
  select(organism_identifier,sex,season,As,Cd,Cr,Cu,Fe,Hg,Mn,Ni,Se,Zn))

rownames(data_niche) <- data_niche$organism_identifier

# Scaling data and data translation
scale_mtx <- as.matrix(data_niche[,4:13])
rownames(scale_mtx) <- data_niche$organism_identifier
scale_mtx <- as.data.frame(Linear_positiv_translation(scale(scale_mtx)))

data_niche <- data_niche %>% 
  select(organism_identifier,sex,season)
data_niche <- cbind(data_niche,scale_mtx)
rm(scale_mtx)

# nMDS computing
MDS_2017 <- metaMDS(data_niche[,4:13],distance = "bray",k = 2,try = 300)

# Get coordinates of the individuals on the nMDS
data_niche_TE <- as.data.frame(scores(MDS_2017))
data_niche_TE$sex <- data_niche$sex 
data_niche_TE$season <- data_niche$season
rm(data_niche)

# Get data for further plotting
FA.scores <- as.data.frame(scores(MDS_2017, "species"))  #Using the scores function from vegan to extract the sites scores and convert to a data.frame
FA.scores$FA <- rownames(FA.scores) # To add species names = FA

rm(MDS_2017)

# Plots
sex <- data_niche_TE %>% 
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_hline(yintercept = 0, linetype = "dashed", colour = 'grey')+
  geom_vline(xintercept = 0, linetype = "dashed", colour = 'grey')+
  stat_ellipse(aes(color = sex), level = 0.95)+
  geom_point(aes(fill = sex), size = 2.5, shape = 21, color = "grey50") + 
  geom_text(data = FA.scores, aes(x = NMDS1, y = NMDS2,label = FA), size = 3, fontface = "bold", color = "grey23") + 
  geom_segment(data = FA.scores, aes(x = 0, y = 0, 
                                     xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.22, "cm")), size = 0.7)+
  scale_color_manual(values = c("#257898","#E05F25"))+
  scale_fill_manual(values = c("#257898","#E05F25"))+
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", color = "Species")+
  theme(axis.title.x = element_text(face = "italic"), 
        axis.title.y = element_text(face = "italic"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.background = element_blank(),
        legend.title = element_blank())

season <- data_niche_TE %>% 
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_hline(yintercept = 0, linetype = "dashed", colour = 'grey')+
  geom_vline(xintercept = 0, linetype = "dashed", colour = 'grey')+
  stat_ellipse(aes(color = season), level = 0.95)+
  geom_point(aes(fill = season), size = 2.5, shape = 21, color = "grey50") + 
  geom_text(data = FA.scores, aes(x = NMDS1, y = NMDS2,label = FA), size = 3, fontface = "bold", color = "grey23") + 
  geom_segment(data = FA.scores, aes(x = 0, y = 0, 
                                     xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.22, "cm")), size = 0.7)+
  scale_color_manual(values = c("#3D4C53","#E64A45","#4DB3B3","#F2C249"))+
  scale_fill_manual(values = c("#3D4C53","#E64A45","#4DB3B3","#F2C249"))+
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", color = "Species")+
  theme(axis.title.x = element_text(face = "italic"), 
        axis.title.y = element_text(face = "italic"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.background = element_blank(),
        legend.title = element_blank())

ggarrange(sex, season, labels = c("A.","B."), ncol = 2, align = "hv")

rm(data_niche_TE, FA.scores)
rm(season,sex)


## 2 / Significant difference in TE mean concentrations between sexes and among seasons ##############################################################################################

# Statistical tests for factor sex
Tracer_data <- as.data.frame(LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(sex = factor(sex, levels = c("M","F"))) %>% 
  select(sex,As,Cd,Cr,Cu,Fe,Hg,Mn,Ni,Se,Zn))

Output_ALL_test <- NULL

for (i in 2: length(colnames(Tracer_data))){ # metal loop
  #i = 2
  
  Tracer_name <- colnames(Tracer_data)[i]
  Tracer_data_F1 <- subset(Tracer_data, sex == "M")
  Tracer_data_F1 <- Tracer_data_F1[,i]
  Tracer_data_F1 <- Tracer_data_F1[!is.na(Tracer_data_F1)]
  Tracer_data_F2 <- subset(Tracer_data, sex == "F")
  Tracer_data_F2 <- Tracer_data_F2[,i]
  Tracer_data_F2 <- Tracer_data_F2[!is.na(Tracer_data_F2)]
  
  test_a <- shapiro.test(Tracer_data_F1) # Normality
  test_b <- shapiro.test(Tracer_data_F2) # Normality
  test_c <- fligner.test(Tracer_data[,i] ~ Tracer_data[,1]) # Homoscedasticity
  
  Output_test <- data.frame(Trace_element = Tracer_name,
                            N1 = length(Tracer_data_F1),
                            N2 = length(Tracer_data_F2),
                            Test = NA,
                            P_val = NA,
                            stat_val = NA)
  
  
  if (Output_test$N1 <= 2 | Output_test$N2 <= 2){ # pas de test pour 5 ech ou moins
    
    Output_ALL_test <- rbind(Output_ALL_test, Output_test)
    
  }else{
    
    if (Output_test$N1 <= 10 | Output_test$N2 <= 10 | test_a$p.value < 0.05 | test_b$p.value < 0.05 | test_c$p.value < 0.05){ # non paramétrique
      # Wilcox
      
      test_stat <- wilcox.test(Tracer_data_F1, Tracer_data_F2, paired = F)
      
      Output_test$Test <- "Wilcoxon"
      Output_test$P_val <- test_stat$p.value
      Output_test$stat_val <- test_stat$statistic
      
    }else{
      # t.test
      test_stat <- t.test(Tracer_data_F1, Tracer_data_F2, paired = F)
      
      Output_test$Test <- "t test"
      Output_test$P_val <- test_stat$p.value
      Output_test$stat_val <- test_stat$statistic
      
    }
    
    Output_ALL_test <- rbind(Output_ALL_test, Output_test)
    
  } # fin test
  
  rm(Output_test, test_stat, Tracer_name, Tracer_data_F2,Tracer_data_F1,test_a, test_b, test_c)
  
} # fin boucle metal

Output_ALL_test[,5] <- round(Output_ALL_test[,5], 3)
rm(Tracer_data,Output_ALL_test,i)


# Statistical tests for factor season
Tracer_data <- as.data.frame(LUB_data %>% 
  filter(!is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(season = factor(season, levels = c("IMS","NWM","IMA"))) %>% 
  select(season,As,Cd,Cr,Cu,Fe,Hg,Mn,Ni,Se,Zn))

Output_ALL_test <- NULL

for (i in 2: length(colnames(Tracer_data))){ # metal loop
  #i = 3
  Tracer_name <- colnames(Tracer_data)[i]
  IMS <- subset(Tracer_data, season == "IMS")
  IMS <- IMS[,i]
  NWM <- subset(Tracer_data, season == "NWM")
  NWM <- NWM[,i]
  IMA <- subset(Tracer_data, season == "IMA")
  IMA <- IMA[,i]
  # SEM <- subset(Tracer_data, season == "SEM")
  # SEM <- SEM[,i]
  
  test_a <- fligner.test(Tracer_data[,i] ~ Tracer_data[,1]) # Homoscedasticity
  a1 <- aov(Tracer_data[,i] ~ Tracer_data[,1])
  test_b <- shapiro.test(resid(a1)) # Normality
  
  Output_test <- data.frame(Tracer = Tracer_name,
                            IMS = length(IMS[!is.na(IMS)]),
                            NWM = length(NWM[!is.na(NWM)]),
                            IMA = length(IMA[!is.na(IMA)]),
                            Test = NA,
                            P_val = NA,
                            stat_val = NA,
                            IMA_IMS = NA,
                            IMA_NWM = NA,
                            IMS_NWM = NA)
  
  rm(IMS,NWM,IMA)
  
  
  if (Output_test$IMS <= 2 | Output_test$NWM <= 2 | Output_test$IMA <= 2){ # pas de test pour 5 ech ou moins
    
    Output_ALL_test <- rbind(Output_ALL_test, Output_test)
    
  }else{
    
    if (Output_test$IMS <= 5 | Output_test$NWM <= 5 | Output_test$IMA <= 5 | test_a$p.value < 0.05 | test_b$p.value < 0.05) { # non paramétrique
      # Kruskal
      
      test_stat <- kruskal.test(Tracer_data[,i] ~ Tracer_data[,1])
      
      Output_test$Test <- "Kruskal"
      Output_test$P_val <- test_stat$p.value
      Output_test$stat_val <- test_stat$statistic
      
      if (Output_test$P_val < 0.05) { # post_hoc
        post_hoc <- dunnTest(Tracer_data[,i] ~ Tracer_data[,1], method = "bh")
        Output_test$IMA_IMS <- post_hoc[[2]][1,4]
        Output_test$IMA_NWM <- post_hoc[[2]][2,4]
        Output_test$IMS_NWM <- post_hoc[[2]][3,4]
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
        rm(post_hoc)
      }else{
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
      }
      
    }else{
      # test normalité
      
      lmFA <- lm(Tracer_data[,i] ~ Tracer_data$season)
      test_stat <- anova(lmFA)
      rm(lmFA)
      
      Output_test$Test <- "ANOVA"
      Output_test$P_val <- test_stat[1,5]
      Output_test$stat_val <- test_stat[1,4]
      
      if (Output_test$P_val < 0.05) {
        post_hoc <- TukeyHSD(a1, 'Tracer_data$season', conf.level = 0.95)
        Output_test$IMA_IMS <- post_hoc[[2]][1,4]
        Output_test$IMA_NWM <- post_hoc[[2]][2,4]
        Output_test$IMS_NWM <- post_hoc[[2]][3,4]
        Output_test$IMA_SEM <- post_hoc[[2]][4,4]
        Output_test$IMS_SEM <- post_hoc[[2]][5,4]
        Output_test$NWM_SEM <- post_hoc[[2]][6,4]
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
        rm(post_hoc)
      } else {
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
      }
      
    } # fin tests stat
    
    rm(Output_test, test_stat, Tracer_name, test_a, test_b, a1)
    
  } # fin boucle si N suffisant
}

Output_ALL_test[,6:10] <- round(Output_ALL_test[,6:10], 3)

rm(Tracer_data,i)
rm(Output_ALL_test)


## 3 / Plotting TE mean concentrations for each sex and each season ##############################################################################################

# Databases
TE_sex <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(sex = factor(sex, levels = c("M","F"))) %>% 
  select(sex,As,Cd,Cr,Cu,Fe,Hg,Mn,Ni,Se,Zn) %>% 
  rename(levels = sex) %>% 
  mutate(factor_type = "Sex") %>% 
  gather(metal, value, -factor_type, -levels)

TE_season <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(season = factor(season, levels = c("IMS","NWM","IMA"))) %>% 
  select(season,As,Cd,Cr,Cu,Fe,Hg,Mn,Ni,Se,Zn) %>% 
  rename(levels = season) %>% 
  mutate(factor_type = "Season") %>% 
  gather(metal, value, -factor_type, -levels)

TE_by_factor <- rbind(TE_sex, TE_season) %>% 
  mutate(factor_type = factor(factor_type, levels = c("Sex","Season")),
         levels = factor(levels, levels = c("M","F","IMS","NWM","IMA")))
rm(TE_sex, TE_season)

# Plot
TE_by_factor %>% 
  group_by(factor_type,levels,metal) %>% 
  summarise(mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE)) %>% 
  mutate(letter = NA,
         letter = ifelse(factor_type == "Sex" & metal == "Se" & levels == "M", "a", letter),
         letter = ifelse(factor_type == "Sex" & metal == "Se" & levels == "F", "b", letter),
         letter = ifelse(factor_type == "Season" & metal %in% c("Cd","Cr","Cu","Mn","Zn") & levels %in% c("NWM","IMA"), "a", letter),
         letter = ifelse(factor_type == "Season" & metal %in% c("Cd","Cr","Cu","Mn","Zn") & levels == "IMS", "b", letter),
         letter = ifelse(factor_type == "Season" & metal == "Hg" & levels == "IMS", "a", letter),
         letter = ifelse(factor_type == "Season" & metal == "Hg" & levels == "NWM", "ab", letter),
         letter = ifelse(factor_type == "Season" & metal == "Hg" & levels == "IMA", "b", letter),
         letter = ifelse(factor_type == "Season" & metal == "Se" & levels == "NWM", "a", letter),
         letter = ifelse(factor_type == "Season" & metal == "Se" & levels == "IMA", "ab", letter),
         letter = ifelse(factor_type == "Season" & metal == "Se" & levels == "IMS", "b", letter),
         letter = ifelse(factor_type == "Season" & metal == "Ni" & levels == "IMA", "a", letter),
         letter = ifelse(factor_type == "Season" & metal == "Ni" & levels == "NWM", "b", letter),
         letter = ifelse(factor_type == "Season" & metal == "Ni" & levels == "IMS", "c", letter)) %>%
  ggplot()+
  geom_errorbar(aes(x = levels, ymin = mean, ymax = mean+sd), color = "#404041", width = 0.3)+
  geom_bar(aes(x = levels, y = mean), color = "#404041", fill = "#BBBDC0", stat = "identity")+
  geom_text(aes(x = levels, y = mean+sd+sd, label = letter), color = "#58585B")+
  labs(y = "Concentrations (µg.g-1 ww)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(vjust = 1, hjust = 0.5),
        axis.title.y = element_text(face = "italic"),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.title=element_blank())+
  facet_grid(metal ~ factor_type, scales = "free")+
  theme(strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"))

rm(TE_by_factor)



### IV // Relationship between TE concentrations and size and d13C/d15N values ##############################################################################################

## 1 / Arsenic ##############################################################################################

database_As <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N,As)

# 1.1 / Step 1
mod_gam_As0 <- mgcv::gam(log10(As) ~ 1,
                         data = database_As, method = "REML")

mod_gam_As1 <- mgcv::gam(log10(As) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_As, method = "REML")
summary(mod_gam_As1) # Complete model

mod_gam_As2 <- mgcv::gam(log10(As) ~ s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_As, method = "REML")
summary(mod_gam_As2) # length

mod_gam_As3 <- mgcv::gam(log10(As) ~ s(length)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_As, method = "REML")
summary(mod_gam_As3) # d13C

mod_gam_As4 <- mgcv::gam(log10(As) ~ s(length)+
                           s(d13C)+
                           d13C:length+
                           d15N:length,
                         data = database_As, method = "REML")
summary(mod_gam_As4) # d15N

mod_gam_As5 <- mgcv::gam(log10(As) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d15N:length,
                         data = database_As, method = "REML")
summary(mod_gam_As5) # d13C:length

mod_gam_As6 <- mgcv::gam(log10(As) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length,
                         data = database_As, method = "REML")
summary(mod_gam_As6) # d15N:length

AIC(mod_gam_As0,mod_gam_As1,mod_gam_As2,mod_gam_As3,mod_gam_As4,
    mod_gam_As5,mod_gam_As6)


# 1.2 / Step 2

# Best = n°3
mod_gam_As7 <- mgcv::gam(log10(As) ~ s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_As, method = "REML")
summary(mod_gam_As7) # length

mod_gam_As8 <- mgcv::gam(log10(As) ~ s(length)+
                           d13C:length+
                           d15N:length,
                         data = database_As, method = "REML")
summary(mod_gam_As8) # d15N

mod_gam_As9 <- mgcv::gam(log10(As) ~ s(length)+
                           s(d15N)+
                           d15N:length,
                         data = database_As, method = "REML")
summary(mod_gam_As9) # d13C:length

mod_gam_As10 <- mgcv::gam(log10(As) ~ s(length)+
                            s(d15N)+
                            d13C:length,
                          data = database_As, method = "REML")
summary(mod_gam_As10) # d15N:length

AIC(mod_gam_As0,mod_gam_As1,mod_gam_As3,mod_gam_As7,mod_gam_As8,mod_gam_As9,
    mod_gam_As10)

# Best model is n°3

# 1.3 / Check GLM

par(mfrow=c(1,2))
plot(mod_gam_As3,se = TRUE, ylab = "", xlab = "",las = 1)

mod_glm_As3 <- glm(log10(As) ~ length+
                     d15N+
                     d13C:length+
                     d15N:length,
                   data = database_As)
summary(mod_glm_As3)
par(mfrow=c(2,2))
plot(mod_glm_As3)

rm(mod_gam_As0,mod_gam_As1,mod_gam_As2,mod_gam_As3,mod_gam_As4,mod_gam_As5,
   mod_gam_As6,mod_gam_As7,mod_gam_As8,mod_gam_As9,mod_gam_As10)
rm(database_As)
rm(mod_glm_As3)


## 2 / With residuals(logAs-size) ##############################################################################################

database_As <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N,As)

# Length-stand of [TE]
glm_logAs_length <- glm(log10(As) ~ length, data = database_As)
summary(glm_logAs_length)

# Diagnostic plots of linear model
par(mfrow=c(2,2))
plot(glm_logAs_length)

# Update database
database_As <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(residAs = residuals(glm_logAs_length)) %>% 
  select(length,d13C,d15N,As,residAs)

# GAM models
mod_gam_As0 <- mgcv::gam(residAs ~ 1,
                         data = database_As, method = "REML")

mod_gam_As1 <- mgcv::gam(residAs ~ s(d13C)+
                           s(d15N),
                         data = database_As, method = "REML")
summary(mod_gam_As1) # Complete model

AIC(mod_gam_As0,mod_gam_As1)

# Plot GAM results
par(mfrow=c(1,2))
plot(mod_gam_As1,se = TRUE, ylab = "", xlab = "", las = 1)

# Check GLM
mod_glm_As1 <- glm(residAs ~ d13C + d15N, data = database_As)
summary(mod_glm_As1)
par(mfrow=c(2,2))
plot(mod_glm_As1)

rm(database_As)
rm(glm_logAs_length)
rm(mod_gam_As0,mod_gam_As1)
rm(mod_glm_As1)


## 3 / Cadmium ##############################################################################################

database_Cd <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N,Cd)

# 3.1 / Step 1
mod_gam_Cd0 <- mgcv::gam(log10(Cd) ~ 1,
                         data = database_Cd, method = "REML")

mod_gam_Cd1 <- mgcv::gam(log10(Cd) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd1) # Complete model

mod_gam_Cd2 <- mgcv::gam(log10(Cd) ~ s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd2) # length

mod_gam_Cd3 <- mgcv::gam(log10(Cd) ~ s(length)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd3) # d13C

mod_gam_Cd4 <- mgcv::gam(log10(Cd) ~ s(length)+
                           s(d13C)+
                           d13C:length+
                           d15N:length,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd4) # d15N

mod_gam_Cd5 <- mgcv::gam(log10(Cd) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d15N:length,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd5) # d13C:length

mod_gam_Cd6 <- mgcv::gam(log10(Cd) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd6) # d15N:length

AIC(mod_gam_Cd0,mod_gam_Cd1,mod_gam_Cd2,mod_gam_Cd3,mod_gam_Cd4,
    mod_gam_Cd5,mod_gam_Cd6)

# 3.2 / Step 2

# Begin from n°6

mod_gam_Cd7 <- mgcv::gam(log10(Cd) ~ s(d13C)+
                           s(d15N)+
                           d13C:length,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd7) # length

mod_gam_Cd8 <- mgcv::gam(log10(Cd) ~ s(length)+
                           s(d15N)+
                           d13C:length,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd8) # d13C

mod_gam_Cd9 <- mgcv::gam(log10(Cd) ~ s(length)+
                           s(d13C)+
                           d13C:length,
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd9) #d15N

mod_gam_Cd10 <- mgcv::gam(log10(Cd) ~ s(length)+
                            s(d13C)+
                            s(d15N),
                          data = database_Cd, method = "REML")
summary(mod_gam_Cd10) #d13C:length

AIC(mod_gam_Cd0,mod_gam_Cd1,mod_gam_Cd6,mod_gam_Cd7,mod_gam_Cd8,mod_gam_Cd9,
    mod_gam_Cd10)

# Best model is n°6

# 3.3 / Check GLM

mod_glm_Cd6 <- glm(log10(Cd) ~ length + d13C + d15N + d13C:length,
                   data = database_Cd)
summary(mod_glm_Cd6)

par(mfrow=c(2,2))
plot(mod_glm_Cd6)

rm(database_Cd)
rm(mod_gam_Cd0,mod_gam_Cd1,mod_gam_Cd2,mod_gam_Cd3,mod_gam_Cd4,mod_gam_Cd5,
   mod_gam_Cd6,mod_gam_Cd7,mod_gam_Cd8,mod_gam_Cd9,mod_gam_Cd10)
rm(mod_glm_Cd6)


## 4 / With residuals(logCd-size) ##############################################################################################

database_Cd <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N,Cd)

# Length-stand of [TE]
glm_logCd_length <- glm(log10(Cd) ~ length, data = database_Cd)
summary(glm_logCd_length)

# Diagnostic plots of linear model
par(mfrow=c(2,2))
plot(glm_logCd_length)

# Database update
database_Cd <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(residCd = residuals(glm_logCd_length)) %>% 
  select(length,d13C,d15N,Cd,residCd)

# GAM models
mod_gam_Cd0 <- mgcv::gam(residCd ~ 1,
                         data = database_Cd, method = "REML")

mod_gam_Cd1 <- mgcv::gam(residCd ~ s(d13C)+
                           s(d15N),
                         data = database_Cd, method = "REML")
summary(mod_gam_Cd1) # Complete model

AIC(mod_gam_Cd0,mod_gam_Cd1)

# Plot GAM results
par(mfrow=c(1,2))
plot(mod_gam_Cd1,se = TRUE, ylab = "", xlab = "", las = 1)

# Check GLM
mod_glm_Cd1 <- glm(residCd ~ d13C + d15N, data = database_Cd)
summary(mod_glm_Cd1)
par(mfrow=c(2,2))
plot(mod_glm_Cd1)

rm(database_Cd)
rm(glm_logCd_length)
rm(mod_gam_Cd0,mod_gam_Cd1)
rm(mod_glm_Cd1)


## 5 / Chromium ##############################################################################################

database_Cr <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N,Cr)

# 5.1 / Step 1
mod_gam_Cr0 <- mgcv::gam(log10(Cr) ~ 1,
                         data = database_Cr, method = "REML")

mod_gam_Cr1 <- mgcv::gam(log10(Cr) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Cr, method = "REML")
summary(mod_gam_Cr1) # Complete model

mod_gam_Cr2 <- mgcv::gam(log10(Cr) ~ s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Cr, method = "REML")
summary(mod_gam_Cr2) # length

mod_gam_Cr3 <- mgcv::gam(log10(Cr) ~ s(length)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Cr, method = "REML")
summary(mod_gam_Cr3) # d13C

mod_gam_Cr4 <- mgcv::gam(log10(Cr) ~ s(length)+
                           s(d13C)+
                           d13C:length+
                           d15N:length,
                         data = database_Cr, method = "REML")
summary(mod_gam_Cr4) # d15N

mod_gam_Cr5 <- mgcv::gam(log10(Cr) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d15N:length,
                         data = database_Cr, method = "REML")
summary(mod_gam_Cr5) # d13C:length

mod_gam_Cr6 <- mgcv::gam(log10(Cr) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length,
                         data = database_Cr, method = "REML")
summary(mod_gam_Cr6) # d15N:length

AIC(mod_gam_Cr0,mod_gam_Cr1,mod_gam_Cr2,mod_gam_Cr3,mod_gam_Cr4,
    mod_gam_Cr5,mod_gam_Cr6)

# Best model is null model (n°0) => no relationship between Cr
# concentrations and length, d13C and d15N

rm(mod_gam_Cr0,mod_gam_Cr1,mod_gam_Cr1bis,mod_gam_Cr2,mod_gam_Cr3,
   mod_gam_Cr4,mod_gam_Cr5,mod_gam_Cr6)
rm(database_Cr)


## 6 / Copper ##############################################################################################

database_Cu <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N,Cu)

# 6.1 / Step 1
mod_gam_Cu0 <- mgcv::gam(log10(Cu) ~ 1,
                         data = database_Cu, method = "REML")

mod_gam_Cu1 <- mgcv::gam(log10(Cu) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu1) # Complete model

mod_gam_Cu2 <- mgcv::gam(log10(Cu) ~ s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu2) # length

mod_gam_Cu3 <- mgcv::gam(log10(Cu) ~ s(length)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu3) # d13C

mod_gam_Cu4 <- mgcv::gam(log10(Cu) ~ s(length)+
                           s(d13C)+
                           d13C:length+
                           d15N:length,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu4) # d15N

mod_gam_Cu5 <- mgcv::gam(log10(Cu) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d15N:length,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu5) # d13C:length

mod_gam_Cu6 <- mgcv::gam(log10(Cu) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length,
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu6) # d15N:length

AIC(mod_gam_Cu0,mod_gam_Cu1,mod_gam_Cu2,mod_gam_Cu3,mod_gam_Cu4,
    mod_gam_Cu5,mod_gam_Cu6)

# Best model still n°1 = complete model

# Check GLM
mod_glm_Cu1 <- glm(log10(Cu) ~ length + d13C + d15N + d13C:length+
                     d15N:length, data = database_Cu)
summary(mod_glm_Cu1)
par(mfrow=c(2,2))
plot(mod_glm_Cu1)

rm(mod_gam_Cu0,mod_gam_Cu1,mod_gam_Cu2,mod_gam_Cu3,mod_gam_Cu4,
   mod_gam_Cu5,mod_gam_Cu6)
rm(database_Cu)
rm(mod_glm_Cu1)


## 7 / With residuals(logCu-size) ##############################################################################################

database_Cu <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N,Cu)

# Length-stand of [TE]
glm_logCu_length <- glm(log10(Cu) ~ length, data = database_Cu)
summary(glm_logCu_length)

# Diagnostic plots of linear model
par(mfrow=c(2,2))
plot(glm_logCu_length)

# Database update
database_Cu <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(residCu = residuals(glm_logCu_length)) %>% 
  select(length,d13C,d15N,Cu,residCu)

# GAM models
mod_gam_Cu0 <- mgcv::gam(residCu ~ 1,
                         data = database_Cu, method = "REML")

mod_gam_Cu1 <- mgcv::gam(residCu ~ s(d13C, k = 3)+
                           s(d15N, k = 3),
                         data = database_Cu, method = "REML")
summary(mod_gam_Cu1) # Complete model

AIC(mod_gam_Cu0,mod_gam_Cu1)

# Plot GAM results
par(mfrow=c(1,2))
plot(mod_gam_Cu1,se = TRUE, ylab = "", xlab = "", las = 1)

# Check GLM
mod_glm_Cu1 <- glm(residCu ~ d13C + d15N, data = database_Cu)
summary(mod_glm_Cu1)
par(mfrow=c(2,2))
plot(mod_glm_Cu1)
# Non-linear, keep GAM

rm(database_Cu)
rm(glm_logCu_length)
rm(mod_gam_Cu0,mod_gam_Cu1)
rm(mod_glm_Cu1)


## 8 / Iron ##############################################################################################

database_Fe <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N,Fe)

# 8.1 / Step 1
mod_gam_Fe0 <- mgcv::gam(log10(Fe) ~ 1,
                         data = database_Fe, method = "REML")

mod_gam_Fe1 <- mgcv::gam(log10(Fe) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Fe, method = "REML")
summary(mod_gam_Fe1) # Complete model

mod_gam_Fe2 <- mgcv::gam(log10(Fe) ~ s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Fe, method = "REML")
summary(mod_gam_Fe2) # length

mod_gam_Fe3 <- mgcv::gam(log10(Fe) ~ s(length)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Fe, method = "REML")
summary(mod_gam_Fe3) # d13C

mod_gam_Fe4 <- mgcv::gam(log10(Fe) ~ s(length)+
                           s(d13C)+
                           d13C:length+
                           d15N:length,
                         data = database_Fe, method = "REML")
summary(mod_gam_Fe4) # d15N

mod_gam_Fe5 <- mgcv::gam(log10(Fe) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d15N:length,
                         data = database_Fe, method = "REML")
summary(mod_gam_Fe5) # d13C:length

mod_gam_Fe6 <- mgcv::gam(log10(Fe) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length,
                         data = database_Fe, method = "REML")
summary(mod_gam_Fe6) # d15N:length

AIC(mod_gam_Fe0,mod_gam_Fe1,mod_gam_Fe2,mod_gam_Fe3,mod_gam_Fe4,
    mod_gam_Fe5,mod_gam_Fe6)

# Best model is n°0, model NULL

rm(mod_gam_Fe0,mod_gam_Fe1,mod_gam_Fe2,mod_gam_Fe3,mod_gam_Fe4,
   mod_gam_Fe5,mod_gam_Fe6)
rm(database_Fe)


## 9 / Mercury ##############################################################################################

database_Hg <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N,Hg)

# 9.1 / Step 1
mod_gam_Hg0 <- mgcv::gam(log10(Hg) ~ 1,
                         data = database_Hg, method = "REML")

mod_gam_Hg1 <- mgcv::gam(log10(Hg) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Hg, method = "REML")
summary(mod_gam_Hg1) # Complete model

mod_gam_Hg2 <- mgcv::gam(log10(Hg) ~ s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Hg, method = "REML")
summary(mod_gam_Hg2) # length

mod_gam_Hg3 <- mgcv::gam(log10(Hg) ~ s(length)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Hg, method = "REML")
summary(mod_gam_Hg3) # d13C

mod_gam_Hg4 <- mgcv::gam(log10(Hg) ~ s(length)+
                           s(d13C)+
                           d13C:length+
                           d15N:length,
                         data = database_Hg, method = "REML")
summary(mod_gam_Hg4) # d15N

mod_gam_Hg5 <- mgcv::gam(log10(Hg) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d15N:length,
                         data = database_Hg, method = "REML")
summary(mod_gam_Hg5) # d13C:length

mod_gam_Hg6 <- mgcv::gam(log10(Hg) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length,
                         data = database_Hg, method = "REML")
summary(mod_gam_Hg6) # d15N:length

AIC(mod_gam_Hg0,mod_gam_Hg1,mod_gam_Hg2,mod_gam_Hg3,mod_gam_Hg4,
    mod_gam_Hg5,mod_gam_Hg6)

# 9.2 / Step 2

# Begin with model n°3
mod_gam_Hg7 <- mgcv::gam(log10(Hg) ~ s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Hg, method = "REML")
summary(mod_gam_Hg7) # length

mod_gam_Hg8 <- mgcv::gam(log10(Hg) ~ s(length)+
                           d13C:length+
                           d15N:length,
                         data = database_Hg, method = "REML")
summary(mod_gam_Hg8) # d15N

mod_gam_Hg9 <- mgcv::gam(log10(Hg) ~ s(length)+
                           s(d15N)+
                           d15N:length,
                         data = database_Hg, method = "REML")
summary(mod_gam_Hg9) #d13C:length

mod_gam_Hg10 <- mgcv::gam(log10(Hg) ~ s(length)+
                            s(d15N)+
                            d13C:length,
                          data = database_Hg, method = "REML")
summary(mod_gam_Hg10) #d15:length

AIC(mod_gam_Hg0,mod_gam_Hg3,mod_gam_Hg7,mod_gam_Hg8,mod_gam_Hg9,
    mod_gam_Hg10)

# Best model is still n°3

# 9.3 / Check GLM
mod_glm_Hg3 <- glm(log10(Hg) ~ length+
                     d15N+
                     d13C:length+
                     d15N:length,
                   data = database_Hg)
summary(mod_glm_Hg3)
par(mfrow=c(2,2))
plot(mod_glm_Hg3)

rm(mod_gam_Hg0,mod_gam_Hg1,mod_gam_Hg2,mod_gam_Hg3,mod_gam_Hg4,mod_gam_Hg5,
   mod_gam_Hg6,mod_gam_Hg7,mod_gam_Hg8,mod_gam_Hg9,mod_gam_Hg10)
rm(database_Hg)
rm(mod_glm_Hg3)


## 10 / With residuals(logHg-size) ##############################################################################################

database_Hg <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N,Hg)

# Length-stand of [TE]
glm_logHg_length <- glm(log10(Hg) ~ length, data = database_Hg)
summary(glm_logHg_length)

# Diagnostic plots of linear model
par(mfrow=c(2,2))
plot(glm_logHg_length)

# Database update
database_Hg <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(residHg = residuals(glm_logHg_length)) %>% 
  select(length,d13C,d15N,Hg,residHg)

# GAM models
mod_gam_Hg0 <- mgcv::gam(residHg ~ 1,
                         data = database_Hg, method = "REML")

mod_gam_Hg1 <- mgcv::gam(residHg ~ s(d13C)+
                           s(d15N),
                         data = database_Hg, method = "REML")
summary(mod_gam_Hg1) # Complete model

AIC(mod_gam_Hg0,mod_gam_Hg1)

# Check GLM
mod_glm_Hg1 <- glm(residHg ~ d13C + d15N, data = database_Hg)
summary(mod_glm_Hg1)
par(mfrow=c(2,2))
plot(mod_glm_Hg1)
# Non-linear = keep GAM

# Plot GAM results
par(mfrow=c(1,2))
plot(mod_gam_Hg1,se = TRUE, ylab = "", xlab = "", las = 1)

rm(database_Hg)
rm(glm_logHg_length)
rm(mod_gam_Hg0,mod_gam_Hg1)
rm(mod_glm_Hg1)


## 11 / Manganese ##############################################################################################

database_Mn <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N,Mn)

# 11.1 / Step 1
mod_gam_Mn0 <- mgcv::gam(log10(Mn) ~ 1,
                         data = database_Mn, method = "REML")

mod_gam_Mn1 <- mgcv::gam(log10(Mn) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn1) # Complete model

mod_gam_Mn2 <- mgcv::gam(log10(Mn) ~ s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn2) # length

mod_gam_Mn3 <- mgcv::gam(log10(Mn) ~ s(length)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn3) # d13C

mod_gam_Mn4 <- mgcv::gam(log10(Mn) ~ s(length)+
                           s(d13C)+
                           d13C:length+
                           d15N:length,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn4) # d15N

mod_gam_Mn5 <- mgcv::gam(log10(Mn) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d15N:length,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn5) # d13C:length

mod_gam_Mn6 <- mgcv::gam(log10(Mn) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn6) # d15N:length

AIC(mod_gam_Mn0,mod_gam_Mn1,mod_gam_Mn2,mod_gam_Mn3,mod_gam_Mn4,
    mod_gam_Mn5,mod_gam_Mn6)

# 11.2 / Step 2

# Begin from model n°3 and also remove d15N (n°4) and d13C:length (n°5)
mod_gam_Mn3bis <- mgcv::gam(log10(Mn) ~ s(length)+
                              d15N:length,
                            data = database_Mn, method = "REML")
summary(mod_gam_Mn3bis)

AIC(mod_gam_Mn0,mod_gam_Mn1,mod_gam_Mn3,mod_gam_Mn3bis)

# Begin from model n°3bis
mod_gam_Mn7 <- mgcv::gam(log10(Mn) ~ d15N:length,
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn7) # length

mod_gam_Mn8 <- mgcv::gam(log10(Mn) ~ s(length),
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn8) #d15N:length

AIC(mod_gam_Mn0,mod_gam_Mn1,mod_gam_Mn3,mod_gam_Mn3bis,mod_gam_Mn7,
    mod_gam_Mn8)

# Best model is still 3bis

# 11.3 / Check GLM
mod_glm_Mn3bis <- glm(log10(Mn) ~ length + d15N:length,
                      data = database_Mn)
summary(mod_glm_Mn3bis)

par(mfrow=c(2,2))
plot(mod_glm_Mn3bis)

rm(mod_gam_Mn0,mod_gam_Mn1,mod_gam_Mn2,mod_gam_Mn3,mod_gam_Mn4,mod_gam_Mn5,
   mod_gam_Mn6,mod_gam_Mn3bis,mod_gam_Mn7,mod_gam_Mn8)
rm(mod_glm_Mn3bis)
rm(database_Mn)


## 12 / With residuals(logMn-size) ##############################################################################################

database_Mn <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N,Mn)

# Length-stand of [TE]
glm_logMn_length <- glm(log10(Mn) ~ length, data = database_Mn)
summary(glm_logMn_length)

# Plot diagnostic plots of linear model
par(mfrow=c(2,2))
plot(glm_logMn_length)

# Database update
database_Mn <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(residMn = residuals(glm_logMn_length)) %>% 
  select(length,d13C,d15N,Mn,residMn)

# GAM models
mod_gam_Mn0 <- mgcv::gam(residMn ~ 1,
                         data = database_Mn, method = "REML")

mod_gam_Mn1 <- mgcv::gam(residMn ~ s(d13C)+
                           s(d15N),
                         data = database_Mn, method = "REML")
summary(mod_gam_Mn1) # Complete model

AIC(mod_gam_Mn0,mod_gam_Mn1)

# Check GLM
mod_glm_Mn1 <- glm(residMn ~ d13C + d15N, data = database_Mn)
summary(mod_glm_Mn1)
par(mfrow=c(2,2))
plot(mod_glm_Mn1)
# Linear, keep GLM

rm(database_Mn)
rm(glm_logMn_length)
rm(mod_gam_Mn0,mod_gam_Mn1)
rm(mod_glm_Mn1)


## 13 / Nickel ##############################################################################################

database_Ni <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N,Ni)

# 13.1 / Step 1
mod_gam_Ni0 <- mgcv::gam(log10(Ni) ~ 1,
                         data = database_Ni, method = "REML")

mod_gam_Ni1 <- mgcv::gam(log10(Ni) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Ni, method = "REML")
summary(mod_gam_Ni1) # Complete model

mod_gam_Ni2 <- mgcv::gam(log10(Ni) ~ s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Ni, method = "REML")
summary(mod_gam_Ni2) # length

mod_gam_Ni3 <- mgcv::gam(log10(Ni) ~ s(length)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Ni, method = "REML")
summary(mod_gam_Ni3) # d13C

mod_gam_Ni4 <- mgcv::gam(log10(Ni) ~ s(length)+
                           s(d13C)+
                           d13C:length+
                           d15N:length,
                         data = database_Ni, method = "REML")
summary(mod_gam_Ni4) # d15N

mod_gam_Ni5 <- mgcv::gam(log10(Ni) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d15N:length,
                         data = database_Ni, method = "REML")
summary(mod_gam_Ni5) # d13C:length

mod_gam_Ni6 <- mgcv::gam(log10(Ni) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length,
                         data = database_Ni, method = "REML")
summary(mod_gam_Ni6) # d15N:length

AIC(mod_gam_Ni0,mod_gam_Ni1,mod_gam_Ni2,mod_gam_Ni3,mod_gam_Ni4,
    mod_gam_Ni5,mod_gam_Ni6)

# 13.2 / Step 2

# Begin from n°6
mod_gam_Ni7 <- mgcv::gam(log10(Ni) ~ s(d13C)+
                           s(d15N)+
                           d13C:length,
                         data = database_Ni, method = "REML")
summary(mod_gam_Ni7) #length

mod_gam_Ni8 <- mgcv::gam(log10(Ni) ~ s(length)+
                           s(d15N)+
                           d13C:length,
                         data = database_Ni, method = "REML")
summary(mod_gam_Ni8) #d13C

mod_gam_Ni9 <- mgcv::gam(log10(Ni) ~ s(length)+
                           s(d13C)+
                           d13C:length,
                         data = database_Ni, method = "REML")
summary(mod_gam_Ni9) #d15?

mod_gam_Ni10 <- mgcv::gam(log10(Ni) ~ s(length)+
                            s(d13C)+
                            s(d15N),
                          data = database_Ni, method = "REML")
summary(mod_gam_Ni10) #d13C:length

AIC(mod_gam_Ni0,mod_gam_Ni1,mod_gam_Ni6,mod_gam_Ni7,mod_gam_Ni8,
    mod_gam_Ni9,mod_gam_Ni10)

# Best model is still n°6

# 13.3. / Check GLM
mod_glm_Ni6 <- glm(log10(Ni) ~ length + d13C + d15N + d13C:length,
                   data = database_Ni)
summary(mod_glm_Ni6)
par(mfrow=c(2,2))
plot(mod_glm_Ni6)

rm(mod_gam_Ni0,mod_gam_Ni1,mod_gam_Ni2,mod_gam_Ni3,mod_gam_Ni4,
   mod_gam_Ni5,mod_gam_Ni6,mod_gam_Ni7,mod_gam_Ni8,mod_gam_Ni9,
   mod_gam_Ni10)
rm(database_Ni)
rm(mod_glm_Ni6)


## 14 / With residuals(logNi-size) ##############################################################################################

database_Ni <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N,Ni)

# Length-stand of [TE]
glm_logNi_length <- glm(log10(Ni) ~ length, data = database_Ni)
summary(glm_logNi_length)

# Plot diagnostic plots of linear model
par(mfrow=c(2,2))
plot(glm_logNi_length)

# Database update
database_Ni <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(residNi = residuals(glm_logNi_length)) %>% 
  select(length,sex,season,d13C,d15N,Ni,residNi)

# GAM models
mod_gam_Ni0 <- mgcv::gam(residNi ~ 1,
                         data = database_Ni, method = "REML")

mod_gam_Ni1 <- mgcv::gam(residNi ~ s(d13C)+
                           s(d15N),
                         data = database_Ni, method = "REML")
summary(mod_gam_Ni1) # Complete model

AIC(mod_gam_Ni0,mod_gam_Ni1)

# Check GLM
mod_glm_Ni1 <- lm(residNi ~ d13C + d15N, data = database_Ni)
summary(mod_glm_Ni1)
par(mfrow=c(2,2))
plot(mod_glm_Ni1)
# Non-linear, keep GAM

# Plot GAM results
par(mfrow=c(1,2))
plot(mod_gam_Ni1,se = TRUE, ylab = "", xlab = "", las = 1)

rm(database_Ni)
rm(glm_logNi_length)
rm(mod_gam_Ni0,mod_gam_Ni1)
rm(mod_glm_Ni1)


## 15 / Selenium ##############################################################################################

database_Se <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N,Se)

# 15.1 / Step 1
mod_gam_Se0 <- mgcv::gam(log10(Se) ~ 1,
                         data = database_Se, method = "REML")

mod_gam_Se1 <- mgcv::gam(log10(Se) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Se, method = "REML")
summary(mod_gam_Se1) # Complete model

mod_gam_Se2 <- mgcv::gam(log10(Se) ~ s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Se, method = "REML")
summary(mod_gam_Se2) # length

mod_gam_Se3 <- mgcv::gam(log10(Se) ~ s(length)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Se, method = "REML")
summary(mod_gam_Se3) # d13C

mod_gam_Se4 <- mgcv::gam(log10(Se) ~ s(length)+
                           s(d13C)+
                           d13C:length+
                           d15N:length,
                         data = database_Se, method = "REML")
summary(mod_gam_Se4) # d15N

mod_gam_Se5 <- mgcv::gam(log10(Se) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d15N:length,
                         data = database_Se, method = "REML")
summary(mod_gam_Se5) # d13C:length

mod_gam_Se6 <- mgcv::gam(log10(Se) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length,
                         data = database_Se, method = "REML")
summary(mod_gam_Se6) # d15N:length

AIC(mod_gam_Se0,mod_gam_Se1,mod_gam_Se2,mod_gam_Se3,mod_gam_Se4,
    mod_gam_Se5,mod_gam_Se6)

# Best model is full model n°1

# Check GLM
mod_glm_Se1 <- glm(log10(Se) ~ length + d13C + d15N + d13C:length+
                     d15N:length,
                   data = database_Se)
summary(mod_glm_Se1)
par(mfrow=c(2,2))
plot(mod_glm_Se1)

rm(mod_gam_Se0,mod_gam_Se1,mod_gam_Se2,mod_gam_Se3,mod_gam_Se4,
   mod_gam_Se5,mod_gam_Se6)
rm(database_Se)
rm(mod_glm_Se1)


## 16 / With residuals(logSe-size) ##############################################################################################

database_Se <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N,Se)

# Length-stand of [TE]
glm_logSe_length <- glm(log10(Se) ~ length, data = database_Se)
summary(glm_logSe_length)

# Plot diagnostic plots of linear model
par(mfrow=c(2,2))
plot(glm_logSe_length)
# Non-linear, use GAM

# Length-stand of [TE] - GAM model
gam_logSe_length <- mgcv::gam(log10(Se) ~ s(length), data = database_Se)
summary(gam_logSe_length)
par(mfrow=c(1,1))
plot(gam_logSe_length)

# Update database
database_Se <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(residSe = residuals(gam_logSe_length)) %>% 
  select(length,d13C,d15N,Se,residSe)

# GAM models
mod_gam_Se0 <- mgcv::gam(residSe ~ 1,
                         data = database_Se, method = "REML")

mod_gam_Se1 <- mgcv::gam(residSe ~ s(d13C)+
                           s(d15N),
                         data = database_Se, method = "REML")
summary(mod_gam_Se1) # Complete model

AIC(mod_gam_Se0,mod_gam_Se1)

# Check GLM
mod_glm_Se1 <- glm(residSe ~ d13C + d15N, data = database_Se)
summary(mod_glm_Se1)
par(mfrow=c(2,2))
plot(mod_glm_Se1) # Non-linear, keep GAM

# Plot GAM results
par(mfrow=c(1,2))
plot(mod_gam_Se1,se = TRUE, ylab = "", xlab = "", las = 1)

rm(database_Se)
rm(glm_logSe_length)
rm(mod_gam_Se0,mod_gam_Se1)
rm(gam_logSe_length)
rm(mod_glm_Se1)


## 17 / Zinc ##############################################################################################

database_Zn <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N,Zn)

# 17.1 / Step 1
mod_gam_Zn0 <- mgcv::gam(log10(Zn) ~ 1,
                         data = database_Zn, method = "REML")

mod_gam_Zn1 <- mgcv::gam(log10(Zn) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn1) # Complete model

mod_gam_Zn2 <- mgcv::gam(log10(Zn) ~ s(d13C)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn2) # length

mod_gam_Zn3 <- mgcv::gam(log10(Zn) ~ s(length)+
                           s(d15N)+
                           d13C:length+
                           d15N:length,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn3) # d13C

mod_gam_Zn4 <- mgcv::gam(log10(Zn) ~ s(length)+
                           s(d13C)+
                           d13C:length+
                           d15N:length,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn4) # d15N

mod_gam_Zn5 <- mgcv::gam(log10(Zn) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d15N:length,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn5) # d13C:length

mod_gam_Zn6 <- mgcv::gam(log10(Zn) ~ s(length)+
                           s(d13C)+
                           s(d15N)+
                           d13C:length,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn6) # d15N:length

AIC(mod_gam_Zn0,mod_gam_Zn1,mod_gam_Zn2,mod_gam_Zn3,mod_gam_Zn4,
    mod_gam_Zn5,mod_gam_Zn6)

# 17.2 / Step 2

# Begin from n°6

mod_gam_Zn7 <- mgcv::gam(log10(Zn) ~ s(d13C)+
                           s(d15N)+
                           d13C:length,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn7) #length

mod_gam_Zn8 <- mgcv::gam(log10(Zn) ~ s(length)+
                           s(d15N)+
                           d13C:length,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn8) #d13C

mod_gam_Zn9 <- mgcv::gam(log10(Zn) ~ s(length)+
                           s(d13C)+
                           d13C:length,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn9) #d15N

mod_gam_Zn10 <- mgcv::gam(log10(Zn) ~ s(length)+
                            s(d13C)+
                            s(d15N),
                          data = database_Zn, method = "REML")
summary(mod_gam_Zn10) #d13C:length

AIC(mod_gam_Zn0,mod_gam_Zn1,mod_gam_Zn6,mod_gam_Zn7,mod_gam_Zn8,mod_gam_Zn9,
    mod_gam_Zn10)

# Best model is still n°6

# 17.3 / Check GLM
mod_glm_Zn6 <- glm(log10(Zn) ~ length + d13C + d15N + d13C:length,
                   data = database_Zn)
summary(mod_glm_Zn6)
par(mfrow=c(2,2))
plot(mod_glm_Zn6)

rm(mod_gam_Zn0,mod_gam_Zn1,mod_gam_Zn2,mod_gam_Zn3,mod_gam_Zn4,mod_gam_Zn5,
   mod_gam_Zn6,mod_gam_Zn7,mod_gam_Zn8,mod_gam_Zn9,mod_gam_Zn10)
rm(database_Zn)
rm(mod_glm_Zn6)


## 18 / With residuals(logZn-size) ##############################################################################################

database_Zn <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N,Zn)

# Length-stand of [TE]
glm_logZn_length <- glm(log10(Zn) ~ length, data = database_Zn)
summary(glm_logZn_length)

# Diagnostic plots of linear model
par(mfrow=c(2,2))
plot(glm_logZn_length)

# Update database
database_Zn <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  mutate(residZn = residuals(glm_logZn_length)) %>% 
  select(length,d13C,d15N,Zn,residZn)

# GAM models
mod_gam_Zn0 <- mgcv::gam(residZn ~ 1,
                         data = database_Zn, method = "REML")

mod_gam_Zn1 <- mgcv::gam(residZn ~ s(d13C, k = 3)+
                           s(d15N, k = 3),
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn1) # Complete model

mod_gam_Zn2 <- mgcv::gam(residZn ~ s(d13C)+
                           s(d15N)+
                           d13C:d15N,
                         data = database_Zn, method = "REML")
summary(mod_gam_Zn2) # Complete model

AIC(mod_gam_Zn0,mod_gam_Zn1,mod_gam_Zn2)

# Check GLM
mod_glm_Zn1 <- glm(residZn ~ d13C + d15N, data = database_Zn)
summary(mod_glm_Zn1)
par(mfrow=c(2,2))
plot(mod_glm_Zn1) # Non-linear, keep GAM

# Plot GAM results
par(mfrow=c(1,2))
plot(mod_gam_Zn1,se = TRUE, ylab = "", xlab = "", las = 1)

rm(database_Zn)
rm(glm_logZn_length)
rm(mod_gam_Zn0,mod_gam_Zn1,mod_gam_Zn2)
rm(mod_glm_Zn1)



### V // Appendix 5.1 - LUB fork length by sex and season ##############################################################################################

## Test for significant difference between sexes
data_sex <- LUB_data %>% 
  filter(sex %in% c("M","F")) %>% 
  mutate(sex = factor(sex, levels = c("M","F")))

fligner.test(data_sex$length ~ data_sex$sex) # Homoscedasticity
shapiro.test(data_sex$length)

wilcox.test(data_sex$length ~ data_sex$sex)
rm(data_sex)


## Calculation mean +/- SD fork length by sex
LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length)) %>% 
  mutate(sex = as.character(sex),
         sex = ifelse(sex == "M", "Males","Females"),
         sex = factor(sex, levels = c("Males","Females"))) %>% 
  group_by(sex) %>% 
  summarise(mean = mean(length),
            sd = sd(length))


## Test for significant difference among seasons
data_season <- LUB_data %>% 
  filter(sex %in% c("M","F")) %>% 
  mutate(season = factor(season, levels = c("IMS","NWM","IMA")))

fligner.test(data_season$length ~ data_season$season) # Homoscedasticity

a1 <- aov(data_season$length ~ data_season$season)
shapiro.test(resid(a1))

kruskal.test(data_season$length ~ data_season$season)
rm(data_season)


## Calculation mean +/- SD fork length by sex
LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length)) %>% 
  mutate(season = as.character(season),
         season = ifelse(season == "IMS", "Pre-NWM",season),
         season = ifelse(season == "IMA", "Pre-SEM",season),
         season = factor(season, levels = c("Pre-NWM","NWM","Pre-SEM"))) %>% 
  group_by(season) %>% 
  summarise(mean = mean(length),
            sd = sd(length))


## Plot fork length by sex and season
LUB_data %>% 
  filter(sex %in% c("M","F")) %>% 
  mutate(sex = ifelse(sex == "M", "Males", "Females"),
         sex = factor(sex, levels = c("Males","Females")),
         season = ifelse(season == "IMS", "Pre-NWM", season),
         season = ifelse(season == "IMA", "Pre-SEM", season),
         season = factor(season, levels = c("Pre-NWM","NWM","Pre-NWM"))) %>% 
  select(sex,season,)

sex <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length)) %>% 
  mutate(sex = as.character(sex),
         sex = ifelse(sex == "M", "Males","Females"),
         sex = factor(sex, levels = c("Males","Females")),
         title = "Sex") %>% 
  ggplot()+
  geom_boxplot(aes(x = sex, y = length), fill = "#BBBDC0" )+  
  theme_bw() + 
  labs(y = "Fork length (cm)")+
  ylim(30,95) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "italic"),
        legend.title = element_blank()) +
  facet_wrap(. ~title)+
  theme(strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"))

season <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length)) %>% 
  mutate(season = as.character(season),
         season = ifelse(season == "IMS", "Pre-NWM",season),
         season = ifelse(season == "IMA", "Pre-SEM",season),
         season = factor(season, levels = c("Pre-NWM","NWM","Pre-SEM")),
         title = "Season") %>% 
  ggplot()+
  geom_boxplot(aes(x = season, y = length), fill = "#BBBDC0")+
  theme_bw() + 
  ylim(30,95) +
  labs(y = "Fork length (cm)")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank())+
  facet_wrap(. ~title)+
  theme(strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"))

ggarrange(sex, season, ncol = 2, align = "hv", widths = c(0.8, 1))
rm(sex,season)



### VI // Appendix 5.2 - Estimation of length around which there is a shift in relationship between d13C and d15N ##############################################################################################

## 1. Plot the relationship between d13C and d15N
LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N) %>% 
  ggplot(aes(x = d13C, y = d15N))+
  geom_point(alpha = 0.3, size = 3)+
  geom_smooth(method = "gam", se = TRUE, color = "#DE8E06")+ # "#DE8E06","#48A4E3"
  geom_hline(yintercept = 14.5, color = "blue")+
  geom_hline(yintercept = 15, color = "blue")+
  geom_hline(yintercept = 14.4, color = "red", linetype = "dashed")+
  geom_vline(xintercept = -15.85, color = "red", linetype = "dashed")+
  geom_point(aes(x = -15.85, y = 14.4), size = 3, color = "red")+
  theme_bw()+
  theme(axis.title = element_text(face = "italic"),
        panel.grid = element_blank())


## 2. Calculation of d15N at shift point
database_test <- LUB_data %>% 
  filter(sex %in% c("M","F"),
         !is.na(length),
         !is.na(d13C),
         !is.na(d15N)) %>% 
  select(length,d13C,d15N)

# Relationship between d13C and d15N - GAM model
mod_gam_SI <- mgcv::gam(d15N ~ s(d13C), data = database_test, method = "REML")
par(mfrow=c(1,1))
plot(mod_gam_SI)

# Predict d15N value with shift point d13C = -15.8
newdata = data.frame(d13C = -15.8)
predict(mod_gam_SI, newdata, type = "response")


## 3. Calculation of size at shift point

# Check linear relationship with GLM
mod_glm_SI <- glm(length ~ d15N + d13C, data = database_test)
summary(mod_glm_SI)
par(mfrow=c(2,2))
plot(mod_glm_SI) # Non-linear relationqhip, use GAM

# Relationship between fork length and d13C + d15N - GAM model
mod_gam_SI <- mgcv::gam(length ~ s(d15N) + s(d13C),
                        data = database_test, method = "REML")
summary(mod_gam_SI)
par(mfrow=c(1,2))
plot(mod_gam_SI)
newdata = data.frame(d13C = -15.8, d15N = 14.4)
predict(mod_gam_SI, newdata, type = "response")
# Predicted size at shift point, for d13C = -15.8 and d15N = 14.4, fork length = 65.1 cm

rm(database_test)
rm(mod_gam_SI, mod_glm_SI, newdata)
