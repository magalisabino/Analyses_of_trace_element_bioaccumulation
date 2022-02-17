##--------------------------------------------------------------------------------------------------------
## SCRIPT : Allows to reproduce all analyses of Chapter 2 dedicated to the trace element bioaccumulation
##          in the three species of spiny lobsters (P. Penicillatus, P. longipes and P. Versicolor).
##          This script includes the following analyses :
##            - computing of correlation tests
##            - computing of trace element profile ellipses, and calculation of associated metrics
##            - computing of t-tests/Wilcoxon tests and of ANOVA/Kruskal-Wallis tests and associated post-hocs
##            - computing of isotopic and fatty acid niches for each sex, habitat reef type and time
##              period of habitat degradation within each species, whenever possible
##
## As part of :
##        Magali SABINO PhD - "Bioaccumulation of trace elements in Seychelles marine food webs"
##
## Author : Magali Sabino
## First created : 2022-01-13
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
lapply(c("tidyverse", "openxlsx", "vegan", "nicheROVER", "ggpubr",
         "FSA","corrplot","viridis"),
       library, character.only = TRUE)



### I // Creation of database with only spiny lobster samples ##############################################################################################

LOB_data <- data_TE_SI %>% 
  filter(c_sp_fao %in% c("NUP","NUV","LOJ")) %>% 
  select(organism_identifier,sex,d13C,d15N,Ag_stat,Cd_stat,Co_stat,
         Cr_stat,Cu_stat,Fe_stat,Mn_stat,Ni_stat,Pb_stat,Se_stat,TAs_stat,
         Zn_stat,THg_stat) %>% 
  rename(Ag = Ag_stat, Cd = Cd_stat, Co = Co_stat, Cr = Cr_stat,
         Cu = Cu_stat, Fe = Fe_stat, Mn = Mn_stat, Ni = Ni_stat,
         Pb = Pb_stat, Se = Se_stat, As = TAs_stat, Zn = Zn_stat,
         Hg = THg_stat) %>% 
  full_join(LOB_FA_percent, by = "organism_identifier")



### II // Size of sampled spiny lobsters ##############################################################################################

LOB_data %>% 
  select(c_sp_fao,sex,carapace_length) %>% 
  filter(!is.na(carapace_length)) %>% 
  gather(length,carapace_length,-c_sp_fao,-sex) %>% 
  select(-length) %>% 
  mutate(c_sp_fao = ifelse(c_sp_fao == "NUP", "Pronghorn spiny lobster", c_sp_fao),
         c_sp_fao = ifelse(c_sp_fao == "LOJ", "Longlegged spiny lobster", c_sp_fao),
         c_sp_fao = ifelse(c_sp_fao == "NUV", "Painted spiny lobster", c_sp_fao),
         c_sp_fao = factor(c_sp_fao, levels = c("Pronghorn spiny lobster","Longlegged spiny lobster","Painted spiny lobster")),
         sex = ifelse(sex == "M", "Male","Female"),
         sex = factor(sex, levels = c("Male","Female"))) %>% 
  group_by(c_sp_fao,sex) %>% 
  mutate(max_value = max(carapace_length,na.rm = TRUE),
         mean = mean(carapace_length, na.rm = TRUE),
         letter = NA,
         letter = ifelse(c_sp_fao %in% c("Pronghorn spiny lobster","Longlegged spiny lobster") & sex == "Male", "a", letter),
         letter = ifelse(c_sp_fao %in% c("Pronghorn spiny lobster","Longlegged spiny lobster") & sex == "Female", "b", letter)) %>% 
  ungroup() %>% 
  ggplot()+
  geom_boxplot(aes(x = sex, y = carapace_length))+
  geom_text(aes(x = sex, y = max_value+0.5, label = letter), color = "#58585D")+
  labs(x = "Sex", y = "Carapace length (cm)")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "italic"),
        panel.grid = element_blank())+
  facet_grid(.~ c_sp_fao)+
  theme(strip.text.x = element_text(face = "bold"))



### III // Interspecific variation in TE profile ##############################################################################################

## 1 / Computing TE profile ellipse for each spiny lobster species ##############################################################################################

# 1.1 / nMDS computing to obtain only 2 dimensions

data_niche <- LOB_data %>% 
  filter(!is.na(Cu),
         !is.na(Hg)) %>% 
  select(organism_identifier,c_sp_fao,biotope,As,Cd,Cr,Cu,Hg,Mn,Ni,Pb,Se,Zn)

data_niche$organism_identifier <- as.factor(data_niche$organism_identifier)
rownames(data_niche) <- data_niche$organism_identifier

# Scaling data and data translation
scale_mtx <- as.matrix(data_niche[,4:13])
rownames(scale_mtx) <- data_niche$organism_identifier
scale_mtx <- as.data.frame(Linear_positiv_translation(scale(scale_mtx)))

data_niche <- data_niche %>% 
  select(organism_identifier,c_sp_fao)
data_niche <- cbind(data_niche,scale_mtx)
rm(scale_mtx)

# nMDS computing
MDS_2017 <- metaMDS(data_niche[,3:12],distance = "bray",k = 2,try = 300)

data_niche_TE <- as.data.frame(scores(MDS_2017)) # Using the scores function from vegan to extract the site scores (here samples labels) and convert to a data.frame
data_niche_TE$c_sp_fao <- data_niche$c_sp_fao # To add sample numbs in data.frame for future merging
data_niche_TE$c_sp_fao <- factor(data_niche_TE$c_sp_fao, levels = c("NUP","LOJ","NUV"))
rm(data_niche)


# 1.2 / Calculation of posterior distribution of (mu, Sigma) for each species

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_TE), data_niche_TE$c_sp_fao,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_TE[ii,1:2]))


# 1.3 / Niche overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:3]
over.stat <- over.stat %>% 
  rename(NUP = `NUP.95%`, LOJ = `LOJ.95%`, NUV = `NUV.95%`) %>% 
  mutate(SpA = c("NUP","LOJ","NUV")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:9]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 3)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 5, 7))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("NUP","LOJ","NUV")),
         SpB = factor(SpB, levels = c("NUP","LOJ","NUV")))

rm(over.stat, over.mean, over.ICinf)


# 1.4 / Niche size calculation

# posterior distribution of niche size by species
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})

# Mise en forme data table
size.stat <- as.data.frame(t(rbind(est = colMeans(fish.size),
                                   quant = apply(fish.size, 2, quantile, prob = c(0.025,0.975)))))
size.stat$est <- round(size.stat$est, 3)
size.stat$`2.5%` <- round(size.stat$`2.5%`, 3)
size.stat$`97.5%` <- round(size.stat$`97.5%`, 3)
names(size.stat) <- c("Mean","Q1","Q2")
size.stat$Sp <- as.character(row.names(size.stat))

fish.probs <- fish.size
as.data.frame(t(rbind(est = colMeans(fish.probs),
                      quant = apply(fish.probs, 2, quantile, prob = c(0.025,0.975)))))
probs <- as.data.frame(fish.probs) %>% 
  mutate(NUV_NUP = ifelse(NUV > NUP,1,0),
         NUV_LOJ = ifelse(NUV > LOJ,1,0),
         LOJ_NUP = ifelse(LOJ > NUP,1,0))
mean(probs$LOJ_NUP)

rm(nsamples, fish.probs, probs)


# 1.5 / Extraction of ellipses for plotting - Mean and small/large

# Extraction of ellipses coordinates
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
data_niche <- data_niche_TE %>% 
  rename(SPP = c_sp_fao)
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

# Calculation of mean and CI95% ellipses for each species
ell.coord_MEAN <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("NUP","LOJ","NUV")))

ell.coord_CIinf <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("NUP","LOJ","NUV")))

ell.coord_CIsup <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("NUP","LOJ","NUV")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)


# 1.5 / Getting nMDS data for further plotting

FA.scores_2017 <- as.data.frame(scores(MDS_2017, "species"))  #Using the scores function from vegan to extract the sites scores and convert to a data.frame
FA.scores_2017$FA <- rownames(FA.scores_2017) # To add species names = FA

rm(MDS_2017)


# 1.2 / Plot ellipses and overlaps

plot1 <- data_niche_TE %>% 
  mutate(c_sp_fao = factor(c_sp_fao, levels = c("NUP","LOJ","NUV"))) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(fill = c_sp_fao), size = 2.5, shape = 21, color = "grey50") +
  geom_polygon(data = ell.coord_MEAN, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_text(data = FA.scores_2017, aes(x = NMDS1, y = NMDS2,label = FA), size = 3, color = "grey23") + # add the FA labels
  scale_color_manual(values = c("#6F94CC","#F2756D","#29B34A"))+
  scale_fill_manual(values = c("#6F94CC","#F2756D","#29B34A"))+
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", color = "Species")+
  theme(axis.title.x = element_text(size=10, face = "italic"),
        axis.title.y = element_text(size=10, face = "italic"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.title=element_blank())

plot2 <- output_overlap %>%
  filter(SpA == "NUP") %>% 
  mutate(title = "Pronghorn spiny lobster",
         SpB = ifelse(SpB == "LOJ","Longlegged","Painted"),
         SpB = factor(SpB, levels = c("Longlegged","Painted"))) %>% 
  ggplot(aes(x = SpB))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")+
  facet_grid(title ~.)+
  theme(strip.text.y = element_text(face = "bold"))

plot3 <- output_overlap %>%
  filter(SpA == "LOJ") %>% 
  mutate(title = "Longlegged spiny lobster",
         SpB = ifelse(SpB == "NUP","Pronghorn","Painted"),
         SpB = factor(SpB, levels = c("Pronghorn","Painted"))) %>% 
  ggplot(aes(x = SpB))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")+
  facet_grid(title ~.)+
  theme(strip.text.y = element_text(face = "bold"))

plot4 <- output_overlap %>%
  filter(SpA == "NUV") %>% 
  mutate(title = "Painted spiny lobster",
         SpB = ifelse(SpB == "NUP","Pronghorn","Longlegged"),
         SpB = factor(SpB, levels = c("Pronghorn","Longlegged"))) %>% 
  ggplot(aes(x = SpB))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")+
  facet_grid(title ~.)+
  theme(strip.text.y = element_text(face = "bold"))

ggarrange(plot1,ggarrange(plot2,plot3,plot4, nrow = 3, ncol = 1, align = "hv"),
          ncol = 2, nrow = 1,
          labels = c("A.","B."),
          align = "hv", widths = c(1,0.6))

rm(plot1,plot2,plot3,plot4,plot5)
rm(ell.coord_MEAN,ell.coord_CIinf, ell.coord_CIsup,FA.scores_2017)
rm(data_niche_TE)
rm(output_overlap)
rm(ii, size.stat)


## 2 / Computing ANOVA/Kruskal-Wallis to determine which species has highest/lowest mean TE concentration ##############################################################################################

# Statistical tests
Tracer_data <- LOB_data %>% 
  select(c_sp_fao,As,Cd,Cr,Cu,Hg,Mn,Ni,Pb,Se,Zn) %>% 
  mutate(c_sp_fao = factor(c_sp_fao, levels = c("NUP","LOJ","NUV")))
Tracer_data <- as.data.frame(Tracer_data)

Output_ALL_test <- NULL

for (i in 2: length(colnames(Tracer_data))){ # metal loop
  #i = 6
  Tracer_name <- colnames(Tracer_data)[i]
  NUP <- subset(Tracer_data, c_sp_fao == "NUP")
  NUP <- NUP[,i]
  LOJ <- subset(Tracer_data, c_sp_fao == "LOJ")
  LOJ <- LOJ[,i]
  NUV <- subset(Tracer_data, c_sp_fao == "NUV")
  NUV <- NUV[,i]
  
  test_a <- fligner.test(Tracer_data[,i] ~ Tracer_data[,1]) # Homoscedasticity
  a1 <- aov(Tracer_data[,i] ~ Tracer_data[,1])
  test_b <- shapiro.test(resid(a1)) # Normality
  
  Output_test <- data.frame(Tracer = Tracer_name,
                            NUP = length(NUP[!is.na(NUP)]),
                            LOJ = length(LOJ[!is.na(LOJ)]),
                            NUV = length(NUV[!is.na(NUV)]),
                            Test = NA,
                            P_val = NA,
                            stat_val = NA,
                            NUP_LOJ = NA,
                            NUP_NUV = NA,
                            LOJ_NUV = NA)
  
  rm(NUP, LOJ, NUV)
  
  
  if (Output_test$NUP <= 2 | Output_test$LOJ <= 2 | Output_test$NUV <= 2){ # pas de test pour 5 ech ou moins
    
    Output_ALL_test <- rbind(Output_ALL_test, Output_test)
    
  }else{
    
    if (Output_test$NUP <= 5 | Output_test$LOJ <= 5 | Output_test$NUV <= 5 | test_a$p.value < 0.05 | test_b$p.value < 0.05) { # non paramétrique
      # Kruskal
      
      test_stat <- kruskal.test(Tracer_data[,i] ~ Tracer_data[,1])
      
      Output_test$Test <- "Kruskal"
      Output_test$P_val <- test_stat$p.value
      Output_test$stat_val <- test_stat$statistic
      
      if (Output_test$P_val < 0.05) { # post_hoc
        post_hoc <- dunnTest(Tracer_data[,i] ~ Tracer_data[,1], method = "bh")
        Output_test$NUP_LOJ <- post_hoc[[2]][1,4]
        Output_test$NUP_NUV <- post_hoc[[2]][2,4]
        Output_test$LOJ_NUV <- post_hoc[[2]][3,4]
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
        rm(post_hoc)
      }else{
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
      }
      
    }else{
      # test normalité
      
      lmFA <- lm(Tracer_data[,i] ~ Tracer_data$c_sp_fao)
      test_stat <- anova(lmFA)
      rm(lmFA)
      
      Output_test$Test <- "ANOVA"
      Output_test$P_val <- test_stat[1,5]
      Output_test$stat_val <- test_stat[1,4]
      
      if (Output_test$P_val < 0.05) {
        a1 <- aov(Tracer_data[,i] ~ Tracer_data$c_sp_fao)
        post_hoc <- TukeyHSD(a1, 'Tracer_data$c_sp_fao', conf.level = 0.95)
        Output_test$NUP_LOJ <- post_hoc[[2]][1,4]
        Output_test$NUP_NUV <- post_hoc[[2]][2,4]
        Output_test$LOJ_NUV <- post_hoc[[2]][3,4]
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
        rm(post_hoc,a1)
      } else {
        Output_ALL_test <- rbind(Output_ALL_test, Output_test)
      }
      
    } # fin tests stat
    
    rm(Output_test, test_stat, Tracer_name, test_a, test_b, a1)
    
  } # fin boucle si N suffisant
} # fin boucle metal

Output_ALL_test[,8:10] <- round(Output_ALL_test[,8:10], 3)
rm(Tracer_data,i)

# Dataframe creation for further plotting
corr_val <- Output_ALL_test %>% 
  select(Tracer,Test,P_val) %>% 
  spread(Tracer,P_val)
rownames(corr_val) <- corr_val$Test
corr_val <- corr_val %>% 
  select(-Test)
corr_val <- as.matrix(corr_val)

sig_val <- Output_ALL_test %>% 
  select(Tracer,Test,P_val) %>% 
  spread(Tracer,P_val)
rownames(sig_val) <- sig_val$Test
sig_val <- sig_val %>% 
  select(-Test)
sig_val <- as.matrix(sig_val)

# Plotting ANOVA/Kruskal-Wallis test results
corrplot(corr_val,
         method = "color",
         addCoef.col = "black",
         tl.col="black", tl.srt = 0,
         p.mat = sig_val, sig.level = 0.05, insig = "blank",
         col = viridis(n = 8),
         diag = T)

rm(Output_ALL_test)
rm(corr_val,sig_val)


## 3 / Plotting TE concentrations for each species when there is a significant difference (p < 0.05) ##############################################################################################

LOB_data %>% 
  select(c_sp_fao,As,Hg,Mn,Ni,Pb,Zn) %>% 
  gather(metal,value,-c_sp_fao) %>% 
  group_by(c_sp_fao,metal) %>% 
  summarise(mean = mean(value),
            sd = sd(value)) %>% 
  ungroup() %>% 
  mutate(letter = NA,
         letter = ifelse(metal == "As" & c_sp_fao %in% c("NUP","NUV"), "a", letter),
         letter = ifelse(metal == "As" & c_sp_fao %in% c("LOJ"), "b", letter),
         letter = ifelse(metal == "Hg" & c_sp_fao %in% c("NUP","LOJ"), "a", letter),
         letter = ifelse(metal == "Hg" & c_sp_fao %in% c("NUV"), "b", letter),
         letter = ifelse(metal == "Mn" & c_sp_fao %in% c("NUP"), "b", letter),
         letter = ifelse(metal == "Mn" & c_sp_fao %in% c("NUV"), "a", letter),
         letter = ifelse(metal == "Pb" & c_sp_fao %in% c("NUP"), "a", letter),
         letter = ifelse(metal == "Pb" & c_sp_fao %in% c("NUV"), "b", letter),
         letter = ifelse(metal %in% c("Mn","Pb") & c_sp_fao == "LOJ", "ab", letter),
         letter = ifelse(metal == "Ni" & c_sp_fao == "NUP", "c", letter),
         letter = ifelse(metal == "Ni" & c_sp_fao == "NUV", "a", letter),
         letter = ifelse(metal == "Zn" & c_sp_fao == "NUP", "a", letter),
         letter = ifelse(metal == "Zn" & c_sp_fao == "NUV", "c", letter),
         letter = ifelse(metal %in% c("Ni","Zn") & c_sp_fao == "LOJ", "b", letter),
         c_sp_fao = ifelse(c_sp_fao == "NUP", "Pronghorn\nspiny lobster",c_sp_fao),
         c_sp_fao = ifelse(c_sp_fao == "LOJ", "Longlegged\nspiny lobster", c_sp_fao),
         c_sp_fao = ifelse(c_sp_fao == "NUV", "Painted\nspiny lobster", c_sp_fao),
         c_sp_fao = factor(c_sp_fao, levels = c("Pronghorn\nspiny lobster","Longlegged\nspiny lobster","Painted\nspiny lobster"))) %>% 
  ggplot()+
  geom_errorbar(aes(x = c_sp_fao, y = mean, ymin = mean-(sd/4), ymax = mean+sd), width = 0.5)+
  geom_bar(aes(x = c_sp_fao, y = mean, fill = c_sp_fao), color = "black", stat = "identity")+
  geom_text(aes(x = c_sp_fao, y = mean+sd+(sd/2), label = letter))+
  labs(y = "Concentration (µg.g-1 ww)")+
  scale_fill_manual(values = c("#6F94CC","#F2756D","#29B34A"))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(face = "italic"),
        axis.ticks.x = element_blank(),
        legend.title = element_blank())+
  facet_wrap(.~metal, ncol = 4, nrow = 3, scales = "free")+
  theme(strip.text.x = element_text(face = "bold"))


### IV // Intraspecific variation in TE profile according to physiology: effect of size ##############################################################################################

# Calculation of correlations
SPP <- sort(unique(LOB_data$c_sp_fao))
Output_ALL_cor <- NULL # création objet final

for (j in 1:length(SPP)) { # species loop
  
  species_name <- SPP[j]
  
  data_metaux <- as.data.frame(LOB_data %>% 
    filter(c_sp_fao == SPP[j]) %>% 
    select(As,Cd,Cr,Cu,Hg,Mn,Ni,Pb,Se,Zn))
  
  data_size <- as.data.frame(LOB_data %>% 
    filter(c_sp_fao == SPP[j]) %>% 
    select(carapace_length))
  
  
  for (i in 1: length(colnames(data_metaux))){ # metal loop
    
    Metal_name <- colnames(data_metaux)[i]
    Metal_data <- data_metaux[,i]
    
    test_a <- shapiro.test(data_size[,1])
    test_b <-shapiro.test(Metal_data)
    
    Output_cor <- data.frame(Species_name = species_name,
                             Metal = Metal_name,
                             N_val_used = sum((is.na(Metal_data)) == F),
                             Normality = NA,
                             Test_cor = NA,
                             P_val = NA,
                             Cor_val = NA)
    
    
    if (Output_cor$N_val_used <= 2){
      
      Output_ALL_cor <- rbind(Output_ALL_cor, Output_cor)
      
    }else{
      
      
      if(test_a$p.value >= 0.05 & test_b$p.value >= 0.05){# Pearson
        Output_cor$Normality <- T
        
        test_stat <- cor.test(data_size[,1], Metal_data, method = "pearson")
        
        Output_cor$Test_cor <- "pearson"
        Output_cor$P_val <- test_stat$p.value
        Output_cor$Cor_val <- test_stat$estimate
        
      }else{
        # Kendall
        Output_cor$Normality <- F
        
        test_stat <- cor.test(data_size[,1], Metal_data, method = "kendall")
        test_stat$p.value
        
        Output_cor$Test_cor <- "kendall"
        Output_cor$P_val <- test_stat$p.value
        Output_cor$Cor_val <- test_stat$estimate
        
      } # fin test correlation
      
      Output_ALL_cor <- rbind(Output_ALL_cor, Output_cor)
      rm(Output_cor, test_stat, test_a, test_b, Metal_data, Metal_name)
      
    } # fin boucle si N suffisant
  } # fin boucle metal
  rm(species_name,data_size,data_metaux)
  Output_ALL_cor[,6:7] <- round(Output_ALL_cor[,6:7], 3)
} # fin boucle 

rm(SPP,i,j)

# Dataframe creation for further plotting
corr_val <- Output_ALL_cor %>% 
  select(Metal,Species_name,Cor_val) %>% 
  mutate(Species_name = as.character(Species_name),
         Species_name = ifelse(Species_name == "LOJ", "Longlegged\nspiny lobster", Species_name),
         Species_name = ifelse(Species_name == "NUP", "Pronghorn\nspiny lobster", Species_name),
         Species_name = ifelse(Species_name == "NUV", "Painted\nspiny lobster", Species_name),
         Species_name = factor(Species_name, levels = c("Pronghorn\nspiny lobster","Longlegged\nspiny lobster","Painted\nspiny lobster"))) %>% 
  spread(Metal,Cor_val)
rownames(corr_val) <- corr_val$Species_name
corr_val <- corr_val %>% 
  select(-Species_name)
corr_val <- as.matrix(corr_val)

sig_val <- Output_ALL_cor %>% 
  select(Metal,Species_name,P_val) %>% 
  mutate(Species_name = as.character(Species_name),
         Species_name = ifelse(Species_name == "LOJ", "Longlegged\nspiny lobster", Species_name),
         Species_name = ifelse(Species_name == "NUP", "Pronghorn\nspiny lobster", Species_name),
         Species_name = ifelse(Species_name == "NUV", "Painted\nspiny lobster", Species_name),
         Species_name = factor(Species_name, levels = c("Pronghorn\nspiny lobster","Longlegged\nspiny lobster","Painted\nspiny lobster"))) %>% 
  spread(Metal,P_val)
rownames(sig_val) <- sig_val$Species_name
sig_val <- sig_val %>% 
  select(-Species_name)
sig_val <- as.matrix(sig_val)

# Correlation plot
corrplot(corr_val,
         method = "color",
         addCoef.col = "black",
         tl.col="black", tl.srt = 0,
         p.mat = sig_val, sig.level = 0.05, insig = "blank",
         diag = T)

rm(Output_ALL_cor)
rm(corr_val,sig_val)



### V // Intraspecific variation in TE profile according to physiology: effect of sex ##############################################################################################

## 1 / Computing TE profile ellipses for each sex and each species ##############################################################################################

# 1.1 / P. penicillatus - nMDS computing to obtain only 2 dimensions

data_niche <- LOB_data %>% 
  filter(c_sp_fao == "NUP") %>% 
  select(organism_identifier,sex,As,Cd,Cr,Cu,Hg,Mn,Ni,Pb,Se,Zn)

# Scaling data and data translation
scale_mtx <- as.matrix(data_niche[,3:12])
rownames(scale_mtx) <- data_niche$organism_identifier
scale_mtx <- as.data.frame(Linear_positiv_translation(scale(scale_mtx)))

data_niche <- data_niche %>% 
  select(organism_identifier,sex)
data_niche <- cbind(data_niche,scale_mtx)
rm(scale_mtx)

# nMDS computing
MDS_2017 <- metaMDS(data_niche[,3:12],distance = "bray",k = 2,try = 300)

# Getting coordinates of individuals on the nMDS
data_niche_NUP <- as.data.frame(scores(MDS_2017)) # Using the scores function from vegan to extract the site scores (here samples labels) and convert to a data.frame
data_niche_NUP$sex <- data_niche$sex # To add sample numbs in data.frame for future merging
data_niche_NUP$label <- data_niche$organism_identifier

rm(data_niche)


# 1.2 / P. penicillatus - Calculation of posterior distribution of (mu, Sigma) for each sex and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_NUP), data_niche_NUP$sex,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_NUP[ii,1:2]))


# Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


## 1.3 / P. penicillatus - Overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("F" = `F.95%`, "M" = `M.95%`) %>% 
  mutate(SpA = c("F","M")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 1)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 3, 3))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_NUP <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("M","F")),
         SpB = factor(SpB, levels = c("M","F")),
         Comp = factor(Comp, levels = c("M in F","F in M")))

rm(over.stat,over.mean,over.ICinf)


## 1.4 / P. penicillatus - Extraction of ellipses

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
data_niche <- data_niche_NUP %>% 
  rename(SPP = sex)
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

# Calculation of mean and CI95% ellipses for each sex
ell.coord_MEAN_NUP <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

ell.coord_CIinf_NUP <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

ell.coord_CIsup_NUP <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


## 1.5 / P. penicillatus - nMDS data for further nMDS plotting

FA.scores_NUP <- as.data.frame(scores(MDS_2017, "species"))  #Using the scores function from vegan to extract the sites scores and convert to a data.frame
FA.scores_NUP$TE <- rownames(FA.scores_NUP) # To add species names = FA

rm(MDS_2017)


## 1.6 / P. penicillatus - Plot ellipses and overlaps

NUP_ellipses <- data_niche_NUP %>% 
  mutate(sex = factor(sex, levels = c("M","F"))) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_polygon(data = ell.coord_MEAN_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_point(aes(fill = sex), size = 2.5, shape = 21, color = "grey50") + # add all individual points
  geom_text(data = FA.scores_NUP, aes(x = NMDS1, y = NMDS2,label = TE), size = 3, fontface = "bold", color = "grey23") + # add the FA labels
  scale_color_manual(values = c("#257898","#E05F25"))+
  scale_fill_manual(values = c("#257898","#E05F25"))+
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", title = "Pronghorn spiny lobsters")+
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

NUP_overlap <- output_overlap_NUP %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)", title = "Pronghorn spiny lobsters")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


## 1.7 / P. longipes - nMDS computing to obtain only 2 dimensions

data_niche <- LOB_data %>% 
  filter(c_sp_fao == "LOJ") %>% 
  select(organism_identifier,sex,As,Cd,Cr,Cu,Hg,Mn,Ni,Pb,Se,Zn)

# Scaling data and data translation
scale_mtx <- as.matrix(data_niche[,3:12])
rownames(scale_mtx) <- data_niche$organism_identifier
scale_mtx <- as.data.frame(Linear_positiv_translation(scale(scale_mtx)))

data_niche <- data_niche %>% 
  select(organism_identifier,sex)
data_niche <- cbind(data_niche,scale_mtx)
rm(scale_mtx)

# nMDS computing
MDS_2017 <- metaMDS(data_niche[,3:12],distance = "bray",k = 2,try = 300)

# Get coordinates of individuals on nMDS
data_niche_LOJ <- as.data.frame(scores(MDS_2017)) # Using the scores function from vegan to extract the site scores (here samples labels) and convert to a data.frame
data_niche_LOJ$sex <- data_niche$sex # To add sample numbs in data.frame for future merging
data_niche_LOJ$label <- data_niche$organism_identifier
rm(data_niche)


# 1.8 / P. longipes - Calculation of posterior distribution of (mu, Sigma) for each sex and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_LOJ), data_niche_LOJ$sex,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_LOJ[ii,1:2]))


# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# 1.9 / P. longipes - Niche overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("F" = `F.95%`, "M" = `M.95%`) %>% 
  mutate(SpA = c("F","M")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 1)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 3, 3))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_LOJ <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("M","F")),
         SpB = factor(SpB, levels = c("M","F")),
         Comp = factor(Comp, levels = c("M in F","F in M")))

rm(over.stat,over.mean,over.ICinf)


# 1.10 / P. longipes - Extraction of ellipses for plotting - Mean and small/large

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
data_niche <- data_niche_LOJ %>% 
  rename(SPP = sex)
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

# Calculation of mean and CI95% ellipses for each sex
ell.coord_MEAN_LOJ <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

ell.coord_CIinf_LOJ <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

ell.coord_CIsup_LOJ <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


# 1.11 / P. longipes - Getting nMDS data for further plotting

FA.scores_LOJ <- as.data.frame(scores(MDS_2017, "species"))  #Using the scores function from vegan to extract the sites scores and convert to a data.frame
FA.scores_LOJ$TE <- rownames(FA.scores_LOJ) # To add species names = FA

rm(MDS_2017)


# 1.12 / P. longipes - Plot ellipses and overlaps

LOJ_ellipses <- data_niche_LOJ %>% 
  mutate(sex = factor(sex, levels = c("M","F"))) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_polygon(data = ell.coord_MEAN_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_point(aes(fill = sex), size = 2.5, shape = 21, color = "grey50") +
  geom_text(data = FA.scores_LOJ, aes(x = NMDS1, y = NMDS2,label = TE), size = 3, fontface = "bold", color = "grey23") +
  scale_color_manual(values = c("#257898","#E05F25"))+
  scale_fill_manual(values = c("#257898","#E05F25"))+
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", color = "Species", title = "Longlegged spiny lobsters")+
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

LOJ_overlap <- output_overlap_LOJ %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)", title = "Longlegged spiny lobsters")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


# 1.13 / P. versicolor - nMDS computing to obtain only 2 dimensions

data_niche <- LOB_data %>% 
  filter(c_sp_fao == "NUV") %>% 
  select(organism_identifier,sex,As,Cd,Cr,Cu,Hg,Mn,Ni,Pb,Se,Zn)

# Scaling data and data translation
scale_mtx <- as.matrix(data_niche[,3:12])
rownames(scale_mtx) <- data_niche$organism_identifier
scale_mtx <- as.data.frame(Linear_positiv_translation(scale(scale_mtx)))

data_niche <- data_niche %>% 
  select(organism_identifier,sex)
data_niche <- cbind(data_niche,scale_mtx)
rm(scale_mtx)

# nMDS computing
MDS_2017 <- metaMDS(data_niche[,3:12],distance = "bray",k = 2,try = 300)

# Get coordinates of individuals on the nMDS
data_niche_NUV <- as.data.frame(scores(MDS_2017)) # Using the scores function from vegan to extract the site scores (here samples labels) and convert to a data.frame
data_niche_NUV$sex <- data_niche$sex # To add sample numbs in data.frame for future merging
data_niche_NUV$label <- data_niche$organism_identifier
rm(data_niche)


# 1.14 / P. versicolor - Calculation of posterior distribution of (mu, Sigma) for each sex and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_NUV), data_niche_NUV$sex,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_NUV[ii,1:2]))


# Niche sizes for each sex - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# 1.15 / P. versicolor - Niche overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("F" = `F.95%`, "M" = `M.95%`) %>% 
  mutate(SpA = c("F","M")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 1)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 3, 3))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_NUV <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("M","F")),
         SpB = factor(SpB, levels = c("M","F")),
         Comp = factor(Comp, levels = c("M in F","F in M")))

rm(over.stat,over.mean,over.ICinf)


# 1.16 / P. versicolor - Extraction of ellipses for plotting - Mean and small/large

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
data_niche <- data_niche_NUV %>% 
  rename(SPP = sex)
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
ell.coord_MEAN_NUV <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

ell.coord_CIinf_NUV <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

ell.coord_CIsup_NUV <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


# 1.17 / P. versicolor - Getting nMDS data for further plotting
FA.scores_NUV <- as.data.frame(scores(MDS_2017, "species"))  #Using the scores function from vegan to extract the sites scores and convert to a data.frame
FA.scores_NUV$TE <- rownames(FA.scores_NUV) # To add species names = FA

rm(MDS_2017)


# 1.18 / P. versicolor - Plot ellipses and overlaps

NUV_ellipses <- data_niche_NUV %>% 
  mutate(sex = factor(sex, levels = c("M","F"))) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_polygon(data = ell.coord_MEAN_NUV, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_NUV, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_NUV, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_point(aes(fill = sex), size = 2.5, shape = 21, color = "grey50") +
  geom_text(data = FA.scores_NUV, aes(x = NMDS1, y = NMDS2,label = TE), size = 3, fontface = "bold", color = "grey23") +
  scale_color_manual(values = c("#257898","#E05F25"))+
  scale_fill_manual(values = c("#257898","#E05F25"))+
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", title = "Painted spiny lobsters")+
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

NUV_overlap <- output_overlap_NUV %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)", title = "Painted spiny lobsters")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


# 1.19 / Plot ellipses and overlaps for the three species

# Ellipses
ggarrange(NUP_ellipses,LOJ_ellipses,NUV_ellipses,
          nrow = 3, align = "hv")

# Overlaps
ggarrange(NUP_overlap,LOJ_overlap,NUV_overlap,
          nrow = 3, align = "hv")

rm(NUP_ellipses,LOJ_ellipses,NUV_ellipses,NUP_overlap,LOJ_overlap,NUV_overlap)
rm(data_niche_NUP, FA.scores_NUP, ell.coord_MEAN_NUP, ell.coord_CIinf_NUP,
   ell.coord_CIsup_NUP, output_overlap_NUP)
rm(data_niche_LOJ, FA.scores_LOJ, ell.coord_MEAN_LOJ, ell.coord_CIinf_LOJ,
   ell.coord_CIsup_LOJ, output_overlap_LOJ)
rm(data_niche_NUV, FA.scores_NUV, ell.coord_MEAN_NUV, ell.coord_CIinf_NUV,
   ell.coord_CIsup_NUV, output_overlap_NUV)


## 2 / Statistical tests and plots to determine which sex has highest/lowest mean TE concentrations ##############################################################################################

# t-test/Wilcoxon test
SPP <- sort(unique(LOB_data$c_sp_fao))
Output_ALL_test <- NULL

for (j in 1:length(SPP)){ # Début boucle espèce
  #j = 3
  Metal_data <- as.data.frame(LOB_data %>% 
    filter(c_sp_fao == SPP[j]) %>% 
    select(c_sp_fao,sex,As,Cd,Cr,Cu,Hg,Mn,Ni,Pb,Se,Zn) %>% 
    mutate(sex = factor(sex, levels = c("M","F"))))
  Species_name = SPP[j]
  
  for (i in 3: length(colnames(Metal_data))){ # metal loop
    
    Metal_name <- colnames(Metal_data)[i]
    Metal_data_F1 <- subset(Metal_data, sex == "M")
    Metal_data_F1 <- Metal_data_F1[,i]
    Metal_data_F1 <- Metal_data_F1[!is.na(Metal_data_F1)]
    Metal_data_F2 <- subset(Metal_data, sex == "F")
    Metal_data_F2 <- Metal_data_F2[,i]
    Metal_data_F2 <- Metal_data_F2[!is.na(Metal_data_F2)]
    
    
    Output_test <- data.frame(Species = Species_name,
                              Metal = Metal_name,
                              N_Male = length(Metal_data_F1),
                              N_Female = length(Metal_data_F2),
                              Test = NA,
                              P_val = NA,
                              stat_val = NA)
    
    
    if (Output_test$N_Male <= 3 | Output_test$N_Female <= 3){ # pas de test pour 5 ech ou moins
      
      Output_ALL_test <- rbind(Output_ALL_test, Output_test)
      
    }else{
      
      if (Output_test$N_Male <= 10 | Output_test$N_Female <= 10){ # si 5 < n < 10, non paramétrique
        
        # Wilcox
        
        test_stat <- wilcox.test(Metal_data_F1, Metal_data_F2, paired = F)
        
        Output_test$Test <- "Wilcoxon"
        Output_test$P_val <- test_stat$p.value
        Output_test$stat_val <- test_stat$statistic
        
      }else{
        # Plus de 10 ech, test normalité & homogénéité
        
        test_a <- shapiro.test(Metal_data_F1) # Normality
        test_b <- shapiro.test(Metal_data_F2) # Normality
        test_c <- fligner.test(Metal_data[,i] ~ Metal_data[,2]) # Homoscedasticity
        
        if(test_a$p.value >= 0.05 & test_b$p.value >= 0.05 & test_c$p.value >= 0.05){# t.test
          
          test_stat <- t.test(Metal_data_F1, Metal_data_F2, paired = F)
          
          Output_test$Test <- "t test"
          Output_test$P_val <- test_stat$p.value
          Output_test$stat_val <- test_stat$statistic
          
        }else{
          # Wilcox
          
          test_stat <- wilcox.test(Metal_data_F1, Metal_data_F2, paired = F)
          
          Output_test$Test <- "Wilcoxon"
          Output_test$P_val <- test_stat$p.value
          Output_test$stat_val <- test_stat$statistic
          
        } # fin test stat
        
        rm(test_a,test_b,test_c)
        
      } # fin boucle si N suffisant
      
      Output_ALL_test <- rbind(Output_ALL_test, Output_test)
      rm(test_stat)
      
    } # fin test
    rm(Output_test, Metal_name, Metal_data_F2,Metal_data_F1)
  } # fin boucle metal
  Output_ALL_test[,6:7] <- round(Output_ALL_test[,6:7], 3)
  rm(Metal_data, Species_name)
} # Fin boucle espèce

rm(SPP,i,j)

# Dataframe creation for further plotting
corr_val <- Output_ALL_test %>% 
  select(Metal,Species,P_val) %>% 
  mutate(Species = as.character(Species),
         Species = ifelse(Species == "LOJ", "Longlegged\nspiny lobster", Species),
         Species = ifelse(Species == "NUP", "Pronghorn\nspiny lobster", Species),
         Species = ifelse(Species == "NUV", "Painted\nspiny lobster", Species),
         Species = factor(Species, levels = c("Pronghorn\nspiny lobster","Longlegged\nspiny lobster","Painted\nspiny lobster"))) %>% 
  spread(Metal,P_val)
rownames(corr_val) <- corr_val$Species
corr_val <- corr_val %>% 
  select(-Species)
corr_val <- as.matrix(corr_val)

sig_val <- Output_ALL_test %>% 
  select(Metal,Species,P_val) %>% 
  mutate(Species = as.character(Species),
         Species = ifelse(Species == "LOJ", "Longlegged\nspiny lobster", Species),
         Species = ifelse(Species == "NUP", "Pronghorn\nspiny lobster", Species),
         Species = ifelse(Species == "NUV", "Painted\nspiny lobster", Species),
         Species = factor(Species, levels = c("Pronghorn\nspiny lobster","Longlegged\nspiny lobster","Painted\nspiny lobster"))) %>% 
  spread(Metal,P_val)
rownames(sig_val) <- sig_val$Species
sig_val <- sig_val %>% 
  select(-Species)
sig_val <- as.matrix(sig_val)

# Plot
corrplot(corr_val,
         method = "color",
         addCoef.col = "black",
         tl.col="black", tl.srt = 0,
         p.mat = sig_val, sig.level = 0.05, insig = "blank",
         col = viridis(n = 8),
         diag = T)

rm(Output_ALL_test)
rm(corr_val,sig_val)


## 3 / Computing FA niches for each sex and each species ##############################################################################################

# 3.1 / P. penicillatus - nMDS computing to obtain only 2 dimensions

data_niche <- LOB_data %>% 
  filter(c_sp_fao == "NUP",
         !is.na(c22_6w3_p)) %>% 
  select(organism_identifier,sex,c14_p,c15_p,c16_1w7_p,c16_p,iso_c17_p,
         c17_p,c18_DMA_p,c18_1w7_p,c18_1w9_p,c18_2w6_p,c18_p,c20_1w11_p,
         c20_1w9_p,c20_5w3_p,c20_p,c22_2i_p,c22_2j_p,c22_4w6_p,c22_5w3_p,
         c22_6w3_p)

# nMDS computing
MDS_2017 <- metaMDS(data_niche[,3:22],distance = "bray",k = 2,try = 300)

# Get coordinates of individuals on the nMDS
data_niche_NUP <- as.data.frame(scores(MDS_2017))
data_niche_NUP$sex <- data_niche$sex
data_niche_NUP$label <- data_niche$organism_identifier
rm(data_niche)


# 3.2 / P. penicillatus - Calculation of posterior distribution of (mu, Sigma) for each sex and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_NUP), data_niche_NUP$sex,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_NUP[ii,1:2]))


# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# 3.3 / P. penicillatus - Niche overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("F" = `F.95%`, "M" = `M.95%`) %>% 
  mutate(SpA = c("F","M")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 1)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 3, 3))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_NUP <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("M","F")),
         SpB = factor(SpB, levels = c("M","F")),
         Comp = factor(Comp, levels = c("M in F","F in M")))

rm(over.stat,over.mean,over.ICinf)


# 3.4 / P. penicillatus - Extraction of ellipses for plotting - Mean and small/large

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
data_niche <- data_niche_NUP %>% 
  rename(SPP = sex)
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

# Calculation of mean and CI95% ellipses for each sex
ell.coord_MEAN_NUP <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

ell.coord_CIinf_NUP <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

ell.coord_CIsup_NUP <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


# 3.5 / P. penicillatus - Getting nMDS data for further plotting

FA.scores_NUP <- as.data.frame(scores(MDS_2017, "species"))  #Using the scores function from vegan to extract the sites scores and convert to a data.frame
FA.scores_NUP$TE <- rownames(FA.scores_NUP) # To add species names = FA

rm(MDS_2017)


# 3.6 / P. penicillatus - Plot ellipses and overlaps

plot1 <- data_niche_NUP %>% 
  mutate(sex = factor(sex, levels = c("M","F"))) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_polygon(data = ell.coord_MEAN_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_point(aes(fill = sex), size = 2.5, shape = 21, color = "grey50") +
  geom_text(data = FA.scores_NUP, aes(x = NMDS1, y = NMDS2,label = TE), size = 3, fontface = "bold", color = "grey23") + 
  scale_color_manual(values = c("#257898","#E05F25"))+
  scale_fill_manual(values = c("#257898","#E05F25"))+
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", color = "Species", title = "Pronghorn spiny lobster")+
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

plot2 <- output_overlap_NUP %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


# 3.7 / P. longipes - nMDS computing to obtain only 2 dimensions

data_niche <- LOB_data %>% 
  filter(c_sp_fao == "LOJ",
         !is.na(c22_6w3_p)) %>% 
  select(organism_identifier,sex,c14_p,c15_p,c16_1w7_p,c16_p,iso_c17_p,
         c17_p,c18_DMA_p,c18_1w7_p,c18_1w9_p,c18_2w6_p,c18_p,c20_1w11_p,
         c20_1w9_p,c20_5w3_p,c20_p,c22_2i_p,c22_2j_p,c22_4w6_p,c22_5w3_p,
         c22_6w3_p)

# nMDS computing
MDS_2017 <- metaMDS(data_niche[,3:22],distance = "bray",k = 2,try = 300)

# Get coordinates of the individuals on the nMDS
data_niche_LOJ <- as.data.frame(scores(MDS_2017)) # Using the scores function from vegan to extract the site scores (here samples labels) and convert to a data.frame
data_niche_LOJ$sex <- data_niche$sex # To add sample numbs in data.frame for future merging
data_niche_LOJ$label <- data_niche$organism_identifier
rm(data_niche)


# 3.8 / P. longipes - Calculation of posterior distribution of (mu, Sigma) for each sex and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_LOJ), data_niche_LOJ$sex,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_LOJ[ii,1:2]))


# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# 3.9 / P. longipes - Niche overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("F" = `F.95%`, "M" = `M.95%`) %>% 
  mutate(SpA = c("F","M")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 1)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 3, 3))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_LOJ <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("M","F")),
         SpB = factor(SpB, levels = c("M","F")),
         Comp = factor(Comp, levels = c("M in F","F in M")))

rm(over.stat,over.mean,over.ICinf)


# 3.10 / P. longipes - Extraction of ellipses for plotting - Mean and small/large

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
data_niche <- data_niche_LOJ %>% 
  rename(SPP = sex)
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

# Calculation of mean and CI95% ellipses for each sex
ell.coord_MEAN_LOJ <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

ell.coord_CIinf_LOJ <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

ell.coord_CIsup_LOJ <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


# 3.11 / P. longipes - Getting nMDS data for further plotting

FA.scores_LOJ <- as.data.frame(scores(MDS_2017, "species"))  #Using the scores function from vegan to extract the sites scores and convert to a data.frame
FA.scores_LOJ$TE <- rownames(FA.scores_LOJ) # To add species names = FA
head(FA.scores_LOJ)

rm(MDS_2017)


# 3.12 / P. longipes - Plot ellipses and overlaps

plot3 <- data_niche_LOJ %>% 
  mutate(sex = factor(sex, levels = c("M","F"))) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_polygon(data = ell.coord_MEAN_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_point(aes(fill = sex), size = 2.5, shape = 21, color = "grey50") + 
  geom_text(data = FA.scores_LOJ, aes(x = NMDS1, y = NMDS2,label = TE), size = 3, fontface = "bold", color = "grey23") + 
  scale_color_manual(values = c("#257898","#E05F25"))+
  scale_fill_manual(values = c("#257898","#E05F25"))+
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", color = "Species", title = "Longlegged spiny lobster")+
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

plot4 <- output_overlap_LOJ %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


# 3.13 / P. versicolor - nMDS computing to obtain only 2 dimensions

data_niche <- LOB_data %>% 
  filter(c_sp_fao == "NUV",
         !is.na(c22_6w3_p)) %>% 
  select(organism_identifier,sex,c14_p,c15_p,c16_1w7_p,c16_p,iso_c17_p,
         c17_p,c18_DMA_p,c18_1w7_p,c18_1w9_p,c18_2w6_p,c18_p,c20_1w11_p,
         c20_1w9_p,c20_5w3_p,c20_p,c22_2i_p,c22_2j_p,c22_4w6_p,c22_5w3_p,
         c22_6w3_p)

# nMDS computing
MDS_2017 <- metaMDS(data_niche[,3:22],distance = "bray",k = 2,try = 300)

# Get coordinates of individuals on the nMDS
data_niche_NUV <- as.data.frame(scores(MDS_2017)) # Using the scores function from vegan to extract the site scores (here samples labels) and convert to a data.frame
data_niche_NUV$sex <- data_niche$sex # To add sample numbs in data.frame for future merging
data_niche_NUV$label <- data_niche$organism_identifier
rm(data_niche)


# 3.14 / P. versicolor - Calculation of posterior distribution of (mu, Sigma) for each sex and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_NUV), data_niche_NUV$sex,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_NUV[ii,1:2]))

# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# 3.15 / P. versicolor - Niche overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("F" = `F.95%`, "M" = `M.95%`) %>% 
  mutate(SpA = c("F","M")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 1)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 3, 3))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_NUV <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("M","F")),
         SpB = factor(SpB, levels = c("M","F")),
         Comp = factor(Comp, levels = c("M in F","F in M")))

rm(over.stat,over.mean,over.ICinf)


# 3.16 / P. versicolor - Extraction of ellipses for plotting - Mean and small/large

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
data_niche <- data_niche_NUV %>% 
  rename(SPP = sex)
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

# Calculation of mean and CI95% ellipses for each sex
ell.coord_MEAN_NUV <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

ell.coord_CIinf_NUV <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

ell.coord_CIsup_NUV <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


# 3.17 / P. versicolor - Getting nMDS data for further plotting

FA.scores_NUV <- as.data.frame(scores(MDS_2017, "species"))  #Using the scores function from vegan to extract the sites scores and convert to a data.frame
FA.scores_NUV$TE <- rownames(FA.scores_NUV) # To add species names = FA
head(FA.scores_NUV)

rm(MDS_2017)


# 3.18 / P. versicolor - Plot ellipses and overlaps

plot5 <- data_niche_NUV %>% 
  mutate(sex = factor(sex, levels = c("M","F"))) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_polygon(data = ell.coord_MEAN_NUV, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_NUV, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_NUV, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_point(aes(fill = sex), size = 2.5, shape = 21, color = "grey50") +
  geom_text(data = FA.scores_NUV, aes(x = NMDS1, y = NMDS2,label = TE), size = 3, fontface = "bold", color = "grey23") + 
  scale_color_manual(values = c("#257898","#E05F25"))+
  scale_fill_manual(values = c("#257898","#E05F25"))+
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", color = "Species", title = "Painted spiny lobster")+
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

plot6 <- output_overlap_NUV %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


# 3.18 / Plot all on same figure

ggarrange(plot1, plot2, plot3, plot4, plot5, plot6,
          nrow = 3, ncol = 2, labels = c("A.","B.","C.","D.","E.","F."),
          widths = c(1,0.4,1,0.4,1,0.4), align = "hv")

rm(plot1, plot2, plot3, plot4, plot5, plot6)
rm(data_niche_NUP, FA.scores_NUP, ell.coord_MEAN_NUP, ell.coord_CIinf_NUP,
   ell.coord_CIsup_NUP, output_overlap_NUP)
rm(data_niche_LOJ, FA.scores_LOJ, ell.coord_MEAN_LOJ, ell.coord_CIinf_LOJ,
   ell.coord_CIsup_LOJ, output_overlap_LOJ)
rm(data_niche_NUV, FA.scores_NUV, ell.coord_MEAN_NUV, ell.coord_CIinf_NUV,
   ell.coord_CIsup_NUV, output_overlap_NUV)






## 4 / Computing isotopic niches for each sex and each species ##############################################################################################

# 4.1 / Difference in d13C and d15N

data_SI <- LOB_data %>% 
  filter(c_sp_fao == "NUV",
         !is.na(d13C)) %>% 
  select(organism_identifier,sex,d13C,d15N) %>% 
  mutate(sex = factor(sex, levels = c("M","F")))

# d13C
shapiro.test(data_SI$d13C) # Normality
fligner.test(data_SI$d13C ~ data_SI$sex) # Homoscedasticité

t.test(data_SI$d13C ~ data_SI$sex)
wilcox.test(data_SI$d13C ~ data_SI$sex)

# d15N
shapiro.test(data_SI$d15N) # Normality
fligner.test(data_SI$d15N ~ data_SI$sex) # Homoscedasticité

t.test(data_SI$d15N ~ data_SI$sex)
wilcox.test(data_SI$d15N ~ data_SI$sex)

rm(data_SI)

# 4.2 / P. penicillatus - Formatting data: data_niche

data_niche_NUP <- LOB_data %>% 
  filter(c_sp_fao == "NUP",
         !is.na(d13C)) %>% 
  select(organism_identifier,sex,d13C,d15N)


# 4.3 / P. penicillatus - Calculation of posterior distribution of (mu, Sigma) for each sex and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_NUP), data_niche_NUP$sex,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_NUP[ii,3:4]))

# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# 4.4 / P. penicillatus - Niche overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("F" = `F.95%`, "M" = `M.95%`) %>% 
  mutate(SpA = c("F","M")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 1)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 3, 3))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_NUP <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("M","F")),
         SpB = factor(SpB, levels = c("M","F")),
         Comp = factor(Comp, levels = c("M in F","F in M")))

rm(over.stat,over.mean,over.ICinf)


# 4.5 / P. penicillatus - Extraction of ellipses for plotting - Mean and small/large

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
data_niche <- data_niche_NUP %>% 
  rename(SPP = sex)
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
  MU1 <- mean(Ellipses_sel$d13C)
  MU2 <- mean(Ellipses_sel$d15N)
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

# Calculation of mean and CI95% ellipses for each sex
ell.coord_MEAN_NUP <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

ell.coord_CIinf_NUP <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

ell.coord_CIsup_NUP <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


# 4.6 / P. penicillatus - Plot ellipses and overlaps

plot1 <- data_niche_NUP %>% 
  group_by(sex) %>% 
  mutate(mean_d13C = mean(d13C, na.rm = TRUE),
         se_d13C = sd(d13C, na.rm = TRUE),
         mean_d15N = mean(d15N, na.rm = TRUE),
         se_d15N = sd(d15N, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(sex = factor(sex, levels = c("M","F"))) %>% 
  ggplot(aes(x = d13C, y = d15N))+
  geom_point(aes(fill = sex), size = 2.5, shape = 21, color = "grey50") + 
  geom_polygon(data = ell.coord_MEAN_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_errorbar(aes(x = mean_d13C, ymin = mean_d15N-se_d15N, ymax = mean_d15N+se_d15N), color = "grey50", width=.1)+
  geom_errorbarh(aes(y = mean_d15N, xmin = mean_d13C-se_d13C, xmax = mean_d13C+se_d13C), color = "grey50", height=.1) +
  geom_point(aes(x = mean_d13C, y = mean_d15N, fill = sex), shape = 23, color = 'grey50', size=5) +
  scale_color_manual(values = c("#257898","#E05F25"))+
  scale_fill_manual(values = c("#257898","#E05F25"))+
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", color = "Species", title = "Pronghorn spiny lobster")+
  theme(axis.title.x = element_text(face = "italic"), 
        axis.title.y = element_text(face = "italic"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

plot2 <- output_overlap_NUP %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


# 4.7 / P. longipes - Formatting data: data_niche

data_niche_LOJ <- LOB_data %>% 
  filter(c_sp_fao == "LOJ",
         !is.na(d13C)) %>% 
  select(organism_identifier,sex,d13C,d15N)


# 4.8 / P. longipes - Calculation of posterior distribution of (mu, Sigma) for each sex and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_LOJ), data_niche_LOJ$sex,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_LOJ[ii,3:4]))


# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# 4.9 / P. longipes - Niche overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("F" = `F.95%`, "M" = `M.95%`) %>% 
  mutate(SpA = c("F","M")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 1)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 3, 3))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_LOJ <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("M","F")),
         SpB = factor(SpB, levels = c("M","F")),
         Comp = factor(Comp, levels = c("M in F","F in M")))

rm(over.stat,over.mean,over.ICinf)


# 4.10 / P. longipes - Extraction of ellipses for plotting - Mean and small/large

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
data_niche <- data_niche_LOJ %>% 
  rename(SPP = sex)
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
  MU1 <- mean(Ellipses_sel$d13C)
  MU2 <- mean(Ellipses_sel$d15N)
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

# Calculation of mean and CI95% ellipses for each sex
ell.coord_MEAN_LOJ <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

ell.coord_CIinf_LOJ <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

ell.coord_CIsup_LOJ <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


# 4.11 / P. longipes - Plot ellipses and overlaps

plot3 <- data_niche_LOJ %>% 
  group_by(sex) %>% 
  mutate(mean_d13C = mean(d13C, na.rm = TRUE),
         se_d13C = sd(d13C, na.rm = TRUE),
         mean_d15N = mean(d15N, na.rm = TRUE),
         se_d15N = sd(d15N, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(sex = factor(sex, levels = c("M","F"))) %>% 
  ggplot(aes(x = d13C, y = d15N))+
  geom_point(aes(fill = sex), size = 2.5, shape = 21, color = "grey50") +
  geom_polygon(data = ell.coord_MEAN_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_errorbar(aes(x = mean_d13C, ymin = mean_d15N-se_d15N, ymax = mean_d15N+se_d15N), color = "grey50", width=.1)+
  geom_errorbarh(aes(y = mean_d15N, xmin = mean_d13C-se_d13C, xmax = mean_d13C+se_d13C), color = "grey50", height=.1) +
  geom_point(aes(x = mean_d13C, y = mean_d15N, fill = sex), shape = 23, color = 'grey50', size=5) +
  scale_color_manual(values = c("#257898","#E05F25"))+
  scale_fill_manual(values = c("#257898","#E05F25"))+
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", color = "Species", title = "Longlegged spiny lobster")+
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

plot4 <- output_overlap_LOJ %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


# 4.12 / P. versicolor - Formatting data: data_niche

data_niche_NUV <- LOB_data %>% 
  filter(c_sp_fao == "NUV",
         !is.na(d13C)) %>% 
  select(organism_identifier,sex,d13C,d15N)


# 4.13 / P. versicolor - Calculation of posterior distribution of (mu, Sigma) for each sex and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_NUV), data_niche_NUV$sex,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_NUV[ii,3:4]))


# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# 4.14 / P. versicolor - Niche overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("F" = `F.95%`, "M" = `M.95%`) %>% 
  mutate(SpA = c("F","M")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 1)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 3, 3))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_NUV <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("M","F")),
         SpB = factor(SpB, levels = c("M","F")),
         Comp = factor(Comp, levels = c("M in F","F in M")))

rm(over.stat,over.mean,over.ICinf)


# 4.15 / P. versicolor - Extraction of ellipses for plotting - Mean and small/large

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
data_niche <- data_niche_NUV %>% 
  rename(SPP = sex)
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
  MU1 <- mean(Ellipses_sel$d13C)
  MU2 <- mean(Ellipses_sel$d15N)
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
ell.coord_MEAN_NUV <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

ell.coord_CIinf_NUV <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

ell.coord_CIsup_NUV <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("M","F")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


# 4.16 / P. versicolor - Plot ellipses and overlaps

plot5 <- data_niche_NUV %>% 
  group_by(sex) %>% 
  mutate(mean_d13C = mean(d13C, na.rm = TRUE),
         se_d13C = sd(d13C, na.rm = TRUE),
         mean_d15N = mean(d15N, na.rm = TRUE),
         se_d15N = sd(d15N, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(sex = factor(sex, levels = c("M","F"))) %>% 
  ggplot(aes(x = d13C, y = d15N))+
  geom_point(aes(fill = sex), size = 2.5, shape = 21, color = "grey50") +
  geom_polygon(data = ell.coord_MEAN_NUV, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_NUV, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_NUV, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_errorbar(aes(x = mean_d13C, ymin = mean_d15N-se_d15N, ymax = mean_d15N+se_d15N), color = "grey50", width=.1)+
  geom_errorbarh(aes(y = mean_d15N, xmin = mean_d13C-se_d13C, xmax = mean_d13C+se_d13C), color = "grey50", height=.1) +
  geom_point(aes(x = mean_d13C, y = mean_d15N, fill = sex), shape = 23, color = 'grey50', size=5) +
  scale_color_manual(values = c("#257898","#E05F25"))+
  scale_fill_manual(values = c("#257898","#E05F25"))+
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", color = "Species", title = "Painted spiny lobster")+
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

plot6 <- output_overlap_NUV %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


# 4.17 / Plot all on same figure

ggarrange(plot1, plot2, plot3, plot4, plot5, plot6,
          nrow = 3, ncol = 2, labels = c("A.","B.","C.","D.","E.","F."),
          widths = c(1,0.4,1,0.4,1,0.4), align = "hv")

rm(plot1, plot2, plot3, plot4, plot5, plot6)
rm(data_niche_NUP, ell.coord_MEAN_NUP, ell.coord_CIinf_NUP,
   ell.coord_CIsup_NUP, output_overlap_NUP)
rm(data_niche_LOJ, ell.coord_MEAN_LOJ, ell.coord_CIinf_LOJ,
   ell.coord_CIsup_LOJ, output_overlap_LOJ)
rm(data_niche_NUV, ell.coord_MEAN_NUV, ell.coord_CIinf_NUV,
   ell.coord_CIsup_NUV, output_overlap_NUV)



### VI // Effect of reef habitat type on intraspecific variation in TE profile ##############################################################################################

## 1 / Computing TE profile ellipses for each reef type and each species ##############################################################################################

# 1.1 / P. penicillatus - nMDS computing to obtain only 2 dimensions

data_niche_NUP <- LOB_data %>% 
  filter(c_sp_fao == "NUP") %>% 
  select(organism_identifier,biotope,As,Cd,Cr,Cu,Hg,Mn,Ni,Pb,Se,Zn) %>% 
  mutate(biotope = ifelse(biotope == "Granite reef", "GRA","CAR"),
         biotope = factor(biotope, levels = c("CAR","GRA")))

rownames(data_niche_NUP) <- data_niche_NUP$organism_identifier

# Scaling data and data translation
scale_mtx <- as.matrix(data_niche_NUP[,3:12])
rownames(scale_mtx) <- data_niche_NUP$organism_identifier
scale_mtx <- as.data.frame(Linear_positiv_translation(scale(scale_mtx)))

data_niche_NUP <- data_niche_NUP %>% 
  select(organism_identifier,biotope)
data_niche_NUP <- cbind(data_niche_NUP,scale_mtx)
rm(scale_mtx)

# nMDS computing
MDS_2017 <- metaMDS(data_niche_NUP[,3:12],distance = "bray",k = 2,try = 300)

# Récupère les coordonnées des individus sur la nMDS
data_niche_TE_NUP <- as.data.frame(scores(MDS_2017)) # Using the scores function from vegan to extract the site scores (here samples labels) and convert to a data.frame
data_niche_TE_NUP$biotope <- data_niche_NUP$biotope # To add sample numbs in data.frame for future merging
data_niche_TE_NUP$label <- data_niche_NUP$organism_identifier
rm(data_niche_NUP)


# 1.2 / P. penicillatus - Calculation of posterior distribution of (mu, Sigma) for each habitat type and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_TE_NUP), data_niche_TE_NUP$biotope,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_TE_NUP[ii,1:2]))


# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# 1.3 / P. penicillatus - Niche size calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("CAR" = `CAR.95%`, "GRA" = `GRA.95%`) %>% 
  mutate(SpA = c("CAR","GRA")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 3)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 5, 7))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_TE_NUP <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("CAR","GRA")),
         SpB = factor(SpB, levels = c("CAR","GRA")),
         Comp = factor(Comp, levels = c("CAR in GRA","GRA in CAR")))

rm(over.stat,over.mean,over.ICinf)


# 1.4 / P. penicillatus - Extraction of ellipses for plotting - Mean and small/large

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
data_niche <- data_niche_TE_NUP %>% 
  rename(SPP = biotope)
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

rm(Sp_names,data_niche,ii)

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
ell.coord_MEAN_TE_NUP <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("CAR","GRA")))

ell.coord_CIinf_TE_NUP <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("CAR","GRA")))

ell.coord_CIsup_TE_NUP <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("CAR","GRA")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


# 1.5 / P. penicillatus - Getting nMDS data for further plotting

FA.scores_TE_NUP <- as.data.frame(scores(MDS_2017, "species"))  #Using the scores function from vegan to extract the sites scores and convert to a data.frame
FA.scores_TE_NUP$TE <- rownames(FA.scores_TE_NUP) # To add species names = FA

rm(MDS_2017)


# 1.6 / P. penicillatus - Plot ellipses and overlaps

plot1 <- data_niche_TE_NUP %>% 
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_polygon(data = ell.coord_MEAN_TE_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_TE_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_TE_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_point(aes(fill = biotope), size = 2.5, shape = 21, color = "grey50") + 
  geom_text(data = FA.scores_TE_NUP, aes(x = NMDS1, y = NMDS2,label = TE), size = 3, fontface = "bold", color = "grey23") +
  scale_color_manual(values = c("#8CC63E","#0A7B57"))+
  scale_fill_manual(values = c("#8CC63E","#0A7B57"))+
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", color = "Species", title = "Pronghorn spiny lobster")+
  theme(axis.title.x = element_text(face = "italic"), 
        axis.title.y = element_text(face = "italic"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

plot2 <- output_overlap_TE_NUP %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


# 1.7 / P. longipes - nMDS computing to obtain only 2 dimensions

data_niche_LOJ <- LOB_data %>% 
  filter(c_sp_fao == "LOJ") %>% 
  select(organism_identifier,biotope,As,Cd,Cr,Cu,Hg,Mn,Ni,Pb,Se,Zn) %>% 
  mutate(biotope = ifelse(biotope == "Granite reef", "GRA","CAR"),
         biotope = factor(biotope, levels = c("CAR","GRA")))

rownames(data_niche_LOJ) <- data_niche_LOJ$organism_identifier

# Scaling data and data translation
scale_mtx <- as.matrix(data_niche_LOJ[,3:12])
rownames(scale_mtx) <- data_niche_LOJ$organism_identifier
scale_mtx <- as.data.frame(Linear_positiv_translation(scale(scale_mtx)))

data_niche_LOJ <- data_niche_LOJ %>% 
  select(organism_identifier,biotope)
data_niche_LOJ <- cbind(data_niche_LOJ,scale_mtx)
rm(scale_mtx)

# nMDS computing
names(data_niche_LOJ)
MDS_2017 <- metaMDS(data_niche_LOJ[,3:12],distance = "bray",k = 2,try = 300)

# Get coordinates of individuals on the nMDS
data_niche_TE_LOJ <- as.data.frame(scores(MDS_2017)) # Using the scores function from vegan to extract the site scores (here samples labels) and convert to a data.frame
data_niche_TE_LOJ$biotope <- data_niche_LOJ$biotope # To add sample numbs in data.frame for future merging
data_niche_TE_LOJ$label <- data_niche_LOJ$organism_identifier
rm(data_niche_LOJ)


# 1.8 / P. longipes - Calculation of posterior distribution of (mu, Sigma) for each habitat type and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_TE_LOJ), data_niche_TE_LOJ$biotope,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_TE_LOJ[ii,1:2]))


# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# 1.9 / P. longipes - Niche overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("CAR" = `CAR.95%`, "GRA" = `GRA.95%`) %>% 
  mutate(SpA = c("CAR","GRA")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 3)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 5, 7))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_TE_LOJ <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("CAR","GRA")),
         SpB = factor(SpB, levels = c("CAR","GRA")),
         Comp = factor(Comp, levels = c("CAR in GRA","GRA in CAR")))

rm(over.stat,over.mean,over.ICinf)


# 1.10 / P. longipes - Extraction of ellipses for plotting - Mean and small/large

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
data_niche <- data_niche_TE_LOJ %>% 
  rename(SPP = biotope)
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

rm(Sp_names,data_niche,ii)

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
ell.coord_MEAN_TE_LOJ <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("CAR","GRA")))

ell.coord_CIinf_TE_LOJ <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("CAR","GRA")))

ell.coord_CIsup_TE_LOJ <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("CAR","GRA")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


# 1.11 / P. longipes - Getting nMDS data for further plotting
FA.scores_TE_LOJ <- as.data.frame(scores(MDS_2017, "species"))  #Using the scores function from vegan to extract the sites scores and convert to a data.frame
FA.scores_TE_LOJ$TE <- rownames(FA.scores_TE_LOJ) # To add species names = FA

rm(MDS_2017)


# 1.12 / P. longipes - Plot ellipses and overlaps

plot3 <- data_niche_TE_LOJ %>% 
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_polygon(data = ell.coord_MEAN_TE_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_TE_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_TE_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_point(aes(fill = biotope), size = 2.5, shape = 21, color = "grey50") + 
  geom_text(data = FA.scores_TE_LOJ, aes(x = NMDS1, y = NMDS2,label = TE), size = 3, fontface = "bold", color = "grey23") + 
  scale_color_manual(values = c("#8CC63E","#0A7B57"))+
  scale_fill_manual(values = c("#8CC63E","#0A7B57"))+
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", color = "Species", title = "Longlegged spiny lobster")+
  theme(axis.title.x = element_text(face = "italic"), 
        axis.title.y = element_text(face = "italic"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

plot4 <- output_overlap_TE_LOJ %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


# 1.13 / P. versicolor - nMDS computing to obtain only 2 dimensions

data_niche_NUV <- LOB_data %>% 
  filter(c_sp_fao == "NUV") %>% 
  select(organism_identifier,biotope,As,Cd,Cr,Cu,Hg,Mn,Ni,Pb,Se,Zn) %>% 
  mutate(biotope = ifelse(biotope == "Granite reef", "GRA","CAR"),
         biotope = factor(biotope, levels = c("CAR","GRA")))

rownames(data_niche_NUV) <- data_niche_NUV$organism_identifier

# Scaling data and data translation
scale_mtx <- as.matrix(data_niche_NUV[,3:12])
rownames(scale_mtx) <- data_niche_NUV$organism_identifier
scale_mtx <- as.data.frame(Linear_positiv_translation(scale(scale_mtx)))

data_niche_NUV <- data_niche_NUV %>% 
  select(organism_identifier,biotope)
data_niche_NUV <- cbind(data_niche_NUV,scale_mtx)
rm(scale_mtx)

# nMDS computing
names(data_niche_NUV)
MDS_2017 <- metaMDS(data_niche_NUV[,3:12],distance = "bray",k = 2,try = 300)

# Get coordinates of individuals on the nMDS
data_niche_TE_NUV <- as.data.frame(scores(MDS_2017)) # Using the scores function from vegan to extract the site scores (here samples labels) and convert to a data.frame
data_niche_TE_NUV$biotope <- data_niche_NUV$biotope # To add sample numbs in data.frame for future merging
data_niche_TE_NUV$label <- data_niche_NUV$organism_identifier
rm(data_niche_NUV)


# 1.14 / P. versicolor - Calculation of posterior distribution of (mu, Sigma) for each habitat type and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_TE_NUV), data_niche_TE_NUV$biotope,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_TE_NUV[ii,1:2]))


# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# 1.15 / P. versicolor - Niche overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("CAR" = `CAR.95%`, "GRA" = `GRA.95%`) %>% 
  mutate(SpA = c("CAR","GRA")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 3)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 5, 7))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_TE_NUV <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("CAR","GRA")),
         SpB = factor(SpB, levels = c("CAR","GRA")),
         Comp = factor(Comp, levels = c("CAR in GRA","GRA in CAR")))

rm(over.stat,over.mean,over.ICinf)


# 1.16 / P. versicolor - Extraction of ellipses for plotting - Mean and small/large

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
data_niche <- data_niche_TE_NUV %>% 
  rename(SPP = biotope)
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

rm(Sp_names,data_niche,ii)

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
ell.coord_MEAN_TE_NUV <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("CAR","GRA")))

ell.coord_CIinf_TE_NUV <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("CAR","GRA")))

ell.coord_CIsup_TE_NUV <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("CAR","GRA")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


# 1.17 / P. versicolor - Getting nMDS data for further plotting

FA.scores_TE_NUV <- as.data.frame(scores(MDS_2017, "species"))  #Using the scores function from vegan to extract the sites scores and convert to a data.frame
FA.scores_TE_NUV$TE <- rownames(FA.scores_TE_NUV) # To add species names = FA

rm(MDS_2017)

# 1.18 / P. versicolor - Plot ellipses and overlaps

plot5 <- data_niche_TE_NUV %>% 
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_polygon(data = ell.coord_MEAN_TE_NUV, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_TE_NUV, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_TE_NUV, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_point(aes(fill = biotope), size = 2.5, shape = 21, color = "grey50") + 
  geom_text(data = FA.scores_TE_NUV, aes(x = NMDS1, y = NMDS2,label = TE), size = 3, fontface = "bold", color = "grey23") + 
  scale_color_manual(values = c("#8CC63E","#0A7B57"))+
  scale_fill_manual(values = c("#8CC63E","#0A7B57"))+
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", color = "Species", title = "Painted spiny lobster")+
  theme(axis.title.x = element_text(face = "italic"), 
        axis.title.y = element_text(face = "italic"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

plot6 <- output_overlap_TE_NUV %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


# 1.19 / Plot all on same figure

ggarrange(plot1,plot2,plot3,plot4,plot5,plot6, ncol = 2, nrow = 3,
          align = "hv", widths = c(1,0.5,1,0.5,1,0.5),
          labels = c("A.","B.","C.","D.","E.","F."))

rm(plot1,plot2,plot3,plot4,plot5,plot6)
rm(output_overlap_TE_NUP,output_overlap_TE_LOJ,output_overlap_TE_NUV,
   FA.scores_TE_NUP,FA.scores_TE_LOJ,FA.scores_TE_NUV,
   ell.coord_MEAN_TE_NUP,ell.coord_MEAN_TE_LOJ,ell.coord_MEAN_TE_NUV,
   ell.coord_CIinf_TE_NUP,ell.coord_CIinf_TE_LOJ,ell.coord_CIinf_TE_NUV,
   ell.coord_CIsup_TE_NUP,ell.coord_CIsup_TE_LOJ,ell.coord_CIsup_TE_NUV)


## 2 / Statistical tests and plots to determine which reef type has highest/lowest mean TE concentrations ##############################################################################################

# t-test/Wilcoxon test
SPP <- sort(unique(LOB_data$c_sp_fao))
Output_ALL_test <- NULL

for (j in 1:length(SPP)){ # Début boucle espèce
  #j = 3
  Metal_data <- LOB_data %>% 
    filter(c_sp_fao == SPP[j]) %>% 
    select(c_sp_fao,biotope,As,Cd,Cr,Cu,Hg,Mn,Ni,Pb,Se,Zn) %>% 
    mutate(biotope = factor(biotope, levels = c("Carbonate reef","Granite reef")))
  Species_name = SPP[j]
  
  for (i in 3: length(colnames(Metal_data))){ # metal loop
    
    Metal_name <- colnames(Metal_data)[i]
    Metal_data_F1 <- subset(Metal_data, biotope == "Carbonate reef")
    Metal_data_F1 <- Metal_data_F1[,i]
    Metal_data_F1 <- Metal_data_F1[!is.na(Metal_data_F1)]
    Metal_data_F2 <- subset(Metal_data, biotope == "Granite reef")
    Metal_data_F2 <- Metal_data_F2[,i]
    Metal_data_F2 <- Metal_data_F2[!is.na(Metal_data_F2)]
    
    
    Output_test <- data.frame(Species = Species_name,
                              Metal = Metal_name,
                              N_CAR = length(Metal_data_F1),
                              N_GRA = length(Metal_data_F2),
                              Test = NA,
                              P_val = NA,
                              stat_val = NA)
    
    
    if (Output_test$N_CAR <= 2 | Output_test$N_GRA <= 2){ # pas de test pour 5 ech ou moins
      
      Output_ALL_test <- rbind(Output_ALL_test, Output_test)
      
    }else{
      
      if (Output_test$N_CAR <= 10 | Output_test$N_GRA <= 10){ # si 5 < n < 10, non paramétrique
        
        # Wilcox
        
        test_stat <- wilcox.test(Metal_data_F1, Metal_data_F2, paired = F)
        
        Output_test$Test <- "Wilcoxon"
        Output_test$P_val <- test_stat$p.value
        Output_test$stat_val <- test_stat$statistic
        
      }else{
        # Plus de 10 ech, test normalité & homogénéité
        
        test_a <- shapiro.test(Metal_data_F1) # Normality
        test_b <- shapiro.test(Metal_data_F2) # Normality
        test_c <- fligner.test(Metal_data[,i] ~ Metal_data[,2]) # Homoscedasticity
        
        if(test_a$p.value >= 0.05 & test_b$p.value >= 0.05 & test_c$p.value >= 0.05){# t.test
          
          test_stat <- t.test(Metal_data_F1, Metal_data_F2, paired = F)
          
          Output_test$Test <- "t test"
          Output_test$P_val <- test_stat$p.value
          Output_test$stat_val <- test_stat$statistic
          
        }else{
          # Wilcox
          
          test_stat <- wilcox.test(Metal_data_F1, Metal_data_F2, paired = F)
          
          Output_test$Test <- "Wilcoxon"
          Output_test$P_val <- test_stat$p.value
          Output_test$stat_val <- test_stat$statistic
          
        } # fin test stat
        
        rm(test_a,test_b,test_c)
        
      } # fin boucle si N suffisant
      
      Output_ALL_test <- rbind(Output_ALL_test, Output_test)
      rm(test_stat)
      
    } # fin test
    rm(Output_test, Metal_name, Metal_data_F2,Metal_data_F1)
  } # fin boucle metal
  Output_ALL_test[,6:7] <- round(Output_ALL_test[,6:7], 3)
  rm(Metal_data, Species_name)
} # Fin boucle espèce

rm(SPP,i,j)

# Plot significant differences
corr_val <- Output_ALL_test %>% 
  select(Metal,Species,P_val) %>% 
  mutate(Species = as.character(Species),
         Species = ifelse(Species == "LOJ", "Longlegged\nspiny lobster", Species),
         Species = ifelse(Species == "NUP", "Pronghorn\nspiny lobster", Species),
         Species = ifelse(Species == "NUV", "Painted\nspiny lobster", Species),
         Species = factor(Species, levels = c("Pronghorn\nspiny lobster","Longlegged\nspiny lobster","Painted\nspiny lobster"))) %>% 
  spread(Metal,P_val)
rownames(corr_val) <- corr_val$Species
corr_val <- corr_val %>% 
  select(-Species)
corr_val <- as.matrix(corr_val)

sig_val <- Output_ALL_test %>% 
  select(Metal,Species,P_val) %>% 
  mutate(Species = as.character(Species),
         Species = ifelse(Species == "LOJ", "Longlegged\nspiny lobster", Species),
         Species = ifelse(Species == "NUP", "Pronghorn\nspiny lobster", Species),
         Species = ifelse(Species == "NUV", "Painted\nspiny lobster", Species),
         Species = factor(Species, levels = c("Pronghorn\nspiny lobster","Longlegged\nspiny lobster","Painted\nspiny lobster"))) %>% 
  spread(Metal,P_val)
rownames(sig_val) <- sig_val$Species
sig_val <- sig_val %>% 
  select(-Species)
sig_val <- as.matrix(sig_val)

corrplot(corr_val,
         method = "color",
         addCoef.col = "black",
         tl.col="black", tl.srt = 0,
         p.mat = sig_val, sig.level = 0.05, insig = "blank",
         col = viridis(n = 8),
         diag = T)

rm(Output_ALL_test)
rm(corr_val,sig_val)


## 3 / Computing isotopic niches for each reef type and each species ##############################################################################################

# 3.1 / P. penicillatus - Formatting data: data_niche

data_niche_SI_NUP <- LOB_data %>% 
  filter(c_sp_fao == "NUP") %>% 
  select(biotope,d13C,d15N) %>% 
  mutate(biotope = ifelse(biotope == "Granite reef", "GRA","CAR"),
         biotope = factor(biotope, levels = c("CAR","GRA")))

# Calculate mean d13C and d15N values for each habitat type
data_niche_SI_NUP %>% 
  group_by(biotope) %>% 
  summarise(mean_d13C = mean(d13C),
            sd_d13C = sd(d13C),
            mean_d15N = mean(d15N),
            sd_d15N = sd(d15N))

# Test for significant difference in mean d13C or mean d15N between habitat types
fligner.test(data_niche_SI_NUP$d13C ~ data_niche_SI_NUP$biotope) # Homoscedasticity
shapiro.test(data_niche_SI_NUP$d13C) # Normality of residuals
t.test(data_niche_SI_NUP$d13C ~ data_niche_SI_NUP$biotope)

fligner.test(data_niche_SI_NUP$d15N ~ data_niche_SI_NUP$biotope) # Homoscedasticity
shapiro.test(data_niche_SI_NUP$d15N) # Normality of residuals
t.test(data_niche_SI_NUP$d15N ~ data_niche_SI_NUP$biotope)


# 3.2 / P. penicillatus - Calculation of posterior distribution of (mu, Sigma) for each habitat type and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_SI_NUP), data_niche_SI_NUP$biotope,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_SI_NUP[ii,2:3]))


# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# 3.3 / P. penicillatus - Niche overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("CAR" = `CAR.95%`, "GRA" = `GRA.95%`) %>% 
  mutate(SpA = c("CAR","GRA")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 3)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 5, 7))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_SI_NUP <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("CAR","GRA")),
         SpB = factor(SpB, levels = c("CAR","GRA")),
         Comp = factor(Comp, levels = c("CAR in GRA","GRA in CAR")))

rm(over.stat,over.mean,over.ICinf)


# 3.4 / P. penicillatus - Extraction of ellipses for plotting - Mean and small/large

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
data_niche <- data_niche_SI_NUP %>% 
  rename(SPP = biotope)
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
  MU1 <- mean(Ellipses_sel$d13C)
  MU2 <- mean(Ellipses_sel$d15N)
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

rm(Sp_names,data_niche,ii)

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
ell.coord_MEAN_SI_NUP <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("CAR","GRA")))

ell.coord_CIinf_SI_NUP <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("CAR","GRA")))

ell.coord_CIsup_SI_NUP <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("CAR","GRA")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


# 3.5 / P. penicillatus - Plot ellipses and overlaps

plot1 <- data_niche_SI_NUP %>% 
  group_by(biotope) %>% 
  mutate(mean_d13C = mean(d13C, na.rm = TRUE),
         se_d13C = sd(d13C, na.rm = TRUE),
         mean_d15N = mean(d15N, na.rm = TRUE),
         se_d15N = sd(d15N, na.rm = TRUE)) %>% 
  ungroup() %>% 
  ggplot(aes(x = d13C, y = d15N))+
  geom_point(aes(fill = biotope), size = 2.5, shape = 21, color = "grey50") + 
  geom_polygon(data = ell.coord_MEAN_SI_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_SI_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_SI_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_errorbar(aes(x = mean_d13C, ymin = mean_d15N-se_d15N, ymax = mean_d15N+se_d15N), color = "grey50", width=.1)+
  geom_errorbarh(aes(y = mean_d15N, xmin = mean_d13C-se_d13C, xmax = mean_d13C+se_d13C), color = "grey50", height=.1) +
  geom_point(aes(x = mean_d13C, y = mean_d15N, fill = biotope), shape = 23, color = 'grey50', size=5) +
  scale_color_manual(values = c("#8CC63E","#0A7B57"))+
  scale_fill_manual(values = c("#8CC63E","#0A7B57"))+
  theme_bw() + 
  labs(x = "d13C", y = "d15N", color = "Species", title = "Pronghorn spiny lobster")+
  theme(axis.title.x = element_text(face = "italic"), 
        axis.title.y = element_text(face = "italic"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

plot2 <- output_overlap_SI_NUP %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


# 3.6 / P. longipes - Formatting data: data_niche

data_niche_SI_LOJ <- LOB_data %>% 
  filter(c_sp_fao == "LOJ") %>% 
  select(biotope,d13C,d15N) %>% 
  mutate(biotope = ifelse(biotope == "Granite reef", "GRA","CAR"),
         biotope = factor(biotope, levels = c("CAR","GRA")))

# Calculate mean d13C and d15N values for each habitat type
data_niche_SI_LOJ %>% 
  group_by(biotope) %>% 
  summarise(mean_d13C = mean(d13C),
            sd_d13C = sd(d13C),
            mean_d15N = mean(d15N),
            sd_d15N = sd(d15N))

# Test for significant difference in mean d13C or mean d15N between habitat types
fligner.test(data_niche_SI_LOJ$d13C ~ data_niche_SI_LOJ$biotope) # Homoscedasticity
shapiro.test(data_niche_SI_LOJ$d13C) # Normality of residuals
t.test(data_niche_SI_LOJ$d13C ~ data_niche_SI_LOJ$biotope)

fligner.test(data_niche_SI_LOJ$d15N ~ data_niche_SI_LOJ$biotope) # Homoscedasticity
shapiro.test(data_niche_SI_LOJ$d15N) # Normality of residuals
t.test(data_niche_SI_LOJ$d15N ~ data_niche_SI_LOJ$biotope)


# 3.7 / P. longipes - Calculation of posterior distribution of (mu, Sigma) for each habitat type and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_SI_LOJ), data_niche_SI_LOJ$biotope,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_SI_LOJ[ii,2:3]))


# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# 3.8 / P. longipes - Niche overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("CAR" = `CAR.95%`, "GRA" = `GRA.95%`) %>% 
  mutate(SpA = c("CAR","GRA")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 3)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 5, 7))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_SI_LOJ <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("CAR","GRA")),
         SpB = factor(SpB, levels = c("CAR","GRA")),
         Comp = factor(Comp, levels = c("CAR in GRA","GRA in CAR")))

rm(over.stat,over.mean,over.ICinf)


# 3.9 / P. longipes - Extraction of ellipses for plotting - Mean and small/large

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
data_niche <- data_niche_SI_LOJ %>% 
  rename(SPP = biotope)
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
  MU1 <- mean(Ellipses_sel$d13C)
  MU2 <- mean(Ellipses_sel$d15N)
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

rm(Sp_names,data_niche,ii)

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

ell.coord_MEAN_SI_LOJ <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("CAR","GRA")))

ell.coord_CIinf_SI_LOJ <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("CAR","GRA")))

ell.coord_CIsup_SI_LOJ <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("CAR","GRA")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


# 3.10 / P. longipes - Plot ellipses and overlaps

plot3 <- data_niche_SI_LOJ %>% 
  group_by(biotope) %>% 
  mutate(mean_d13C = mean(d13C, na.rm = TRUE),
         se_d13C = sd(d13C, na.rm = TRUE),
         mean_d15N = mean(d15N, na.rm = TRUE),
         se_d15N = sd(d15N, na.rm = TRUE)) %>% 
  ungroup() %>% 
  ggplot(aes(x = d13C, y = d15N))+
  geom_point(aes(fill = biotope), size = 2.5, shape = 21, color = "grey50") + 
  geom_polygon(data = ell.coord_MEAN_SI_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_SI_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_SI_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_errorbar(aes(x = mean_d13C, ymin = mean_d15N-se_d15N, ymax = mean_d15N+se_d15N), color = "grey50", width=.1)+
  geom_errorbarh(aes(y = mean_d15N, xmin = mean_d13C-se_d13C, xmax = mean_d13C+se_d13C), color = "grey50", height=.1) +
  geom_point(aes(x = mean_d13C, y = mean_d15N, fill = biotope), shape = 23, color = 'grey50', size=5) +
  scale_color_manual(values = c("#8CC63E","#0A7B57"))+
  scale_fill_manual(values = c("#8CC63E","#0A7B57"))+
  theme_bw() + 
  labs(x = "d13C", y = "d15N", color = "Species", title = "Longlegged spiny lobster")+
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

plot4 <- output_overlap_SI_LOJ %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


# 3.11 / P. versicolor - Formatting data: data_niche

data_niche_SI_NUV <- LOB_data %>% 
  filter(c_sp_fao == "NUV") %>% 
  select(biotope,d13C,d15N) %>% 
  mutate(biotope = ifelse(biotope == "Granite reef", "GRA","CAR"),
         biotope = factor(biotope, levels = c("CAR","GRA")))

# Calculate mean d13C and d15N values for each habitat type
data_niche_SI_NUV %>% 
  group_by(biotope) %>% 
  summarise(mean_d13C = mean(d13C),
            sd_d13C = sd(d13C),
            mean_d15N = mean(d15N),
            sd_d15N = sd(d15N))

# Test for significant difference in mean d13C or mean d15N between habitat types
fligner.test(data_niche_SI_NUV$d13C ~ data_niche_SI_NUV$biotope) # Homoscedasticity
shapiro.test(data_niche_SI_NUV$d13C) # Normality of residuals
wilcox.test(data_niche_SI_NUV$d13C ~ data_niche_SI_NUV$biotope)

fligner.test(data_niche_SI_NUV$d15N ~ data_niche_SI_NUV$biotope) # Homoscedasticity
shapiro.test(data_niche_SI_NUV$d15N) # Normality of residuals
wilcox.test(data_niche_SI_NUV$d15N ~ data_niche_SI_NUV$biotope)


# 3.12 / P. versicolor - Calculation of posterior distribution of (mu, Sigma) for each habitat type and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_SI_NUV), data_niche_SI_NUV$biotope,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_SI_NUV[ii,2:3]))


# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# 3.13 / P. versicolor - Niche overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("CAR" = `CAR.95%`, "GRA" = `GRA.95%`) %>% 
  mutate(SpA = c("CAR","GRA")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 3)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 5, 7))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_SI_NUV <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("CAR","GRA")),
         SpB = factor(SpB, levels = c("CAR","GRA")),
         Comp = factor(Comp, levels = c("CAR in GRA","GRA in CAR")))

rm(over.stat,over.mean,over.ICinf)


# 3.14 / P. versicolor - Extraction of ellipses for plotting - Mean and small/large

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
data_niche <- data_niche_SI_NUV %>% 
  rename(SPP = biotope)
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
  MU1 <- mean(Ellipses_sel$d13C)
  MU2 <- mean(Ellipses_sel$d15N)
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

rm(Sp_names,data_niche,ii)

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
ell.coord_MEAN_SI_NUV <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("CAR","GRA")))

ell.coord_CIinf_SI_NUV <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("CAR","GRA")))

ell.coord_CIsup_SI_NUV <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("CAR","GRA")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


# 3.15 / P. versicolor - Plot ellipses and overlaps

plot5 <- data_niche_SI_NUV %>% 
  group_by(biotope) %>% 
  mutate(mean_d13C = mean(d13C, na.rm = TRUE),
         se_d13C = sd(d13C, na.rm = TRUE),
         mean_d15N = mean(d15N, na.rm = TRUE),
         se_d15N = sd(d15N, na.rm = TRUE)) %>% 
  ungroup() %>% 
  ggplot(aes(x = d13C, y = d15N))+
  geom_point(aes(fill = biotope), size = 2.5, shape = 21, color = "grey50") + 
  geom_polygon(data = ell.coord_MEAN_SI_NUV, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_SI_NUV, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_SI_NUV, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_errorbar(aes(x = mean_d13C, ymin = mean_d15N-se_d15N, ymax = mean_d15N+se_d15N), color = "grey50", width=.1)+
  geom_errorbarh(aes(y = mean_d15N, xmin = mean_d13C-se_d13C, xmax = mean_d13C+se_d13C), color = "grey50", height=.1) +
  geom_point(aes(x = mean_d13C, y = mean_d15N, fill = biotope), shape = 23, color = 'grey50', size=5) +
  scale_color_manual(values = c("#8CC63E","#0A7B57"))+
  scale_fill_manual(values = c("#8CC63E","#0A7B57"))+
  theme_bw() + 
  labs(x = "d13C", y = "d15N", color = "Species", title = "Painted spiny lobster")+
  theme(axis.title.x = element_text(face = "italic"), 
        axis.title.y = element_text(face = "italic"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

plot6 <- output_overlap_SI_NUV %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


# 3.16 / Plot all on same figure

ggarrange(plot1,plot2,plot3,plot4,plot5,plot6, ncol = 2, nrow = 3,
          align = "hv", widths = c(1,0.5,1,0.5,1,0.5),
          labels = c("A.","B.","C.","D.","E.","F."))

rm(plot1,plot2,plot3,plot4,plot5,plot6)
rm(data_niche_SI_NUP,data_niche_SI_LOJ,data_niche_SI_NUV)
rm(output_overlap_SI_NUP,output_overlap_SI_LOJ,output_overlap_SI_NUV,
   ell.coord_MEAN_SI_NUP,ell.coord_MEAN_SI_LOJ,ell.coord_MEAN_SI_NUV,
   ell.coord_CIinf_SI_NUP,ell.coord_CIinf_SI_LOJ,ell.coord_CIinf_SI_NUV,
   ell.coord_CIsup_SI_NUP,ell.coord_CIsup_SI_LOJ,ell.coord_CIsup_SI_NUV)



### VII // Effect of habitat degradation on intraspecific variation in TE profile ##############################################################################################

## 1 / Computing TE profile ellipses for period of habitat degradation and each species (only NUP & LOJ) ##############################################################################################

# 1.1 / P. penicillatus - nMDS computing to obtain only 2 dimensions

data_niche_NUP <- LOB_data %>% 
  filter(c_sp_fao == "NUP") %>% 
  select(organism_identifier,bleaching_event,As,Cd,Cr,Cu,Hg,Mn,Ni,Pb,Se,Zn) %>% 
  mutate(bleaching_event = ifelse(bleaching_event == "Before", "BEF","AFT"),
         bleaching_event = factor(bleaching_event, levels = c("BEF","AFT")))

rownames(data_niche_NUP) <- data_niche_NUP$organism_identifier

# Scaling data and data translation
scale_mtx <- as.matrix(data_niche_NUP[,3:12])
rownames(scale_mtx) <- data_niche_NUP$organism_identifier
scale_mtx <- as.data.frame(Linear_positiv_translation(scale(scale_mtx)))

data_niche_NUP <- data_niche_NUP %>% 
  select(organism_identifier,bleaching_event)
data_niche_NUP <- cbind(data_niche_NUP,scale_mtx)
rm(scale_mtx)

# nMDS computing
MDS_2017 <- metaMDS(data_niche_NUP[,3:12],distance = "bray",k = 2,try = 300)

# Récupère les coordonnées des individus sur la nMDS
data_niche_TE_NUP <- as.data.frame(scores(MDS_2017)) # Using the scores function from vegan to extract the site scores (here samples labels) and convert to a data.frame
data_niche_TE_NUP$bleaching_event <- data_niche_NUP$bleaching_event # To add sample numbs in data.frame for future merging
data_niche_TE_NUP$label <- data_niche_NUP$organism_identifier
rm(data_niche_NUP)


# 1.2 / P. penicillatus - Calculation of posterior distribution of (mu, Sigma) for each period of habitat degradation and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_TE_NUP), data_niche_TE_NUP$bleaching_event,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_TE_NUP[ii,1:2]))


# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# 1.3 / P. penicillatus - Niche overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("BEF" = `BEF.95%`, "AFT" = `AFT.95%`) %>% 
  mutate(SpA = c("BEF","AFT")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 3)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 5, 7))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_TE_NUP <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("BEF","AFT")),
         SpB = factor(SpB, levels = c("BEF","AFT")),
         Comp = factor(Comp, levels = c("BEF in AFT","AFT in BEF")))

rm(over.stat,over.mean,over.ICinf)


# 1.4 / P. penicillatus - Extraction of ellipses for plotting - Mean and small/large

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
data_niche <- data_niche_TE_NUP %>% 
  rename(SPP = bleaching_event)
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

rm(Sp_names,data_niche,ii)

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
ell.coord_MEAN_TE_NUP <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("BEF","AFT")))

ell.coord_CIinf_TE_NUP <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("BEF","AFT")))

ell.coord_CIsup_TE_NUP <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("BEF","AFT")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


# 1.5 / P. penicillatus - Getting nMDS data for further plotting

FA.scores_TE_NUP <- as.data.frame(scores(MDS_2017, "species"))  #Using the scores function from vegan to extract the sites scores and convert to a data.frame
FA.scores_TE_NUP$TE <- rownames(FA.scores_TE_NUP) # To add species names = FA

rm(MDS_2017)


# 1.6 / P. penicillatus - Plot ellipses and overlaps

plot1 <- data_niche_TE_NUP %>% 
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_polygon(data = ell.coord_MEAN_TE_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_TE_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_TE_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_point(aes(fill = bleaching_event), size = 2.5, shape = 21, color = "grey50") +
  geom_text(data = FA.scores_TE_NUP, aes(x = NMDS1, y = NMDS2,label = TE), size = 3, fontface = "bold", color = "grey23") +
  scale_color_manual(values = c("#157BBE","#113658"))+
  scale_fill_manual(values = c("#157BBE","#113658"))+
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", color = "Species", title = "Pronghorn spiny lobster")+
  theme(axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

plot2 <- output_overlap_TE_NUP %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


# 1.7 / P. longipes - nMDS computing to obtain only 2 dimensions

data_niche_LOJ <- LOB_data %>% 
  filter(c_sp_fao == "LOJ") %>% 
  select(organism_identifier,bleaching_event,As,Cd,Cr,Cu,Hg,Mn,Ni,Pb,Se,Zn) %>% 
  mutate(bleaching_event = ifelse(bleaching_event == "Before", "BEF","AFT"),
         bleaching_event = factor(bleaching_event, levels = c("BEF","AFT")))

rownames(data_niche_LOJ) <- data_niche_LOJ$organism_identifier

# Scaling data and data translation
scale_mtx <- as.matrix(data_niche_LOJ[,3:12])
rownames(scale_mtx) <- data_niche_LOJ$organism_identifier
scale_mtx <- as.data.frame(Linear_positiv_translation(scale(scale_mtx)))

data_niche_LOJ <- data_niche_LOJ %>% 
  select(organism_identifier,bleaching_event)
data_niche_LOJ <- cbind(data_niche_LOJ,scale_mtx)
rm(scale_mtx)

# nMDS computing
MDS_2017 <- metaMDS(data_niche_LOJ[,3:12],distance = "bray",k = 2,try = 300)

# Récupère les coordonnées des individus sur la nMDS
data_niche_TE_LOJ <- as.data.frame(scores(MDS_2017)) # Using the scores function from vegan to extract the site scores (here samples labels) and convert to a data.frame
data_niche_TE_LOJ$bleaching_event <- data_niche_LOJ$bleaching_event # To add sample numbs in data.frame for future merging
data_niche_TE_LOJ$label <- data_niche_LOJ$organism_identifier
rm(data_niche_LOJ)


# 1.8 / P. longipes - Calculation of posterior distribution of (mu, Sigma) for each period of habitat degradation and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_TE_LOJ), data_niche_TE_LOJ$bleaching_event,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_TE_LOJ[ii,1:2]))


# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# 1.9 / P. longipes - Niche overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("BEF" = `BEF.95%`, "AFT" = `AFT.95%`) %>% 
  mutate(SpA = c("BEF","AFT")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 3)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 5, 7))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_TE_LOJ <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("BEF","AFT")),
         SpB = factor(SpB, levels = c("BEF","AFT")),
         Comp = factor(Comp, levels = c("BEF in AFT","AFT in BEF")))

rm(over.stat,over.mean,over.ICinf)


# 1.10 / P. longipes - Extraction of ellipses for plotting - Mean and small/large

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
data_niche <- data_niche_TE_LOJ %>% 
  rename(SPP = bleaching_event)
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

rm(Sp_names,data_niche,ii)

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
ell.coord_MEAN_TE_LOJ <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("BEF","AFT")))

ell.coord_CIinf_TE_LOJ <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("BEF","AFT")))

ell.coord_CIsup_TE_LOJ <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("BEF","AFT")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


# 1.11 / P. longipes - Getting nMDS data for further plotting

FA.scores_TE_LOJ <- as.data.frame(scores(MDS_2017, "species"))  #Using the scores function from vegan to extract the sites scores and convert to a data.frame
FA.scores_TE_LOJ$TE <- rownames(FA.scores_TE_LOJ) # To add species names = FA

rm(MDS_2017)


# 1.12 / P. longipes - Plot ellipses and overlaps

plot3 <- data_niche_TE_LOJ %>% 
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_polygon(data = ell.coord_MEAN_TE_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_TE_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_TE_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_point(aes(fill = bleaching_event), size = 2.5, shape = 21, color = "grey50") + 
  geom_text(data = FA.scores_TE_LOJ, aes(x = NMDS1, y = NMDS2,label = TE), size = 3, fontface = "bold", color = "grey23") + 
  scale_color_manual(values = c("#157BBE","#113658"))+
  scale_fill_manual(values = c("#157BBE","#113658"))+
  theme_bw() + 
  labs(x = "nMDS1", y = "nMDS2", color = "Species", title = "Longlegged spiny lobster")+
  theme(axis.title.x = element_text(face = "italic"), 
        axis.title.y = element_text(face = "italic"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

plot4 <- output_overlap_TE_LOJ %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


# 1.13 / Plot all on same figure

ggarrange(plot1,plot2,plot3,plot4, ncol = 2, nrow = 2,
          align = "hv", widths = c(1,0.5,1,0.5),
          labels = c("A.","B.","C.","D."))

rm(plot1,plot2,plot3,plot4)
rm(output_overlap_TE_NUP,output_overlap_TE_LOJ,FA.scores_TE_NUP,
   FA.scores_TE_LOJ,ell.coord_MEAN_TE_NUP,ell.coord_MEAN_TE_LOJ,
   ell.coord_CIinf_TE_NUP,ell.coord_CIinf_TE_LOJ,
   ell.coord_CIsup_TE_NUP,ell.coord_CIsup_TE_LOJ)


## 2 / Statistical tests and plots to determine which period has highest/lowest mean TE concentrations ##############################################################################################

# Statistical tests
SPP <- c("NUP","LOJ")
Output_ALL_test <- NULL

for (j in 1:length(SPP)){ # Début boucle espèce
  #j = 3
  Metal_data <- LOB_data %>% 
    filter(c_sp_fao == SPP[j]) %>% 
    select(c_sp_fao,bleaching_event,As,Cd,Cr,Cu,Hg,Mn,Ni,Pb,Se,Zn) %>% 
    mutate(bleaching_event = factor(bleaching_event, levels = c("Before","After")))
  Species_name = SPP[j]
  
  for (i in 3: length(colnames(Metal_data))){ # metal loop
    
    Metal_name <- colnames(Metal_data)[i]
    Metal_data_F1 <- subset(Metal_data, bleaching_event == "Before")
    Metal_data_F1 <- Metal_data_F1[,i]
    Metal_data_F1 <- Metal_data_F1[!is.na(Metal_data_F1)]
    Metal_data_F2 <- subset(Metal_data, bleaching_event == "After")
    Metal_data_F2 <- Metal_data_F2[,i]
    Metal_data_F2 <- Metal_data_F2[!is.na(Metal_data_F2)]
    
    
    Output_test <- data.frame(Species = Species_name,
                              Metal = Metal_name,
                              N_BEF = length(Metal_data_F1),
                              N_AFT = length(Metal_data_F2),
                              Test = NA,
                              P_val = NA,
                              stat_val = NA)
    
    
    if (Output_test$N_BEF <= 3 | Output_test$N_AFT <= 3){ # pas de test pour 5 ech ou moins
      
      Output_ALL_test <- rbind(Output_ALL_test, Output_test)
      
    }else{
      
      if (Output_test$N_BEF <= 10 | Output_test$N_AFT <= 10){ # si 5 < n < 10, non paramétrique
        
        # Wilcox
        
        test_stat <- wilcox.test(Metal_data_F1, Metal_data_F2, paired = F)
        
        Output_test$Test <- "Wilcoxon"
        Output_test$P_val <- test_stat$p.value
        Output_test$stat_val <- test_stat$statistic
        
      }else{
        # Plus de 10 ech, test normalité & homogénéité
        
        test_a <- shapiro.test(Metal_data_F1) # Normality
        test_b <- shapiro.test(Metal_data_F2) # Normality
        test_c <- fligner.test(Metal_data[,i] ~ Metal_data[,2]) # Homoscedasticity
        
        if(test_a$p.value >= 0.05 & test_b$p.value >= 0.05 & test_c$p.value >= 0.05){# t.test
          
          test_stat <- t.test(Metal_data_F1, Metal_data_F2, paired = F)
          
          Output_test$Test <- "t test"
          Output_test$P_val <- test_stat$p.value
          Output_test$stat_val <- test_stat$statistic
          
        }else{
          # Wilcox
          
          test_stat <- wilcox.test(Metal_data_F1, Metal_data_F2, paired = F)
          
          Output_test$Test <- "Wilcoxon"
          Output_test$P_val <- test_stat$p.value
          Output_test$stat_val <- test_stat$statistic
          
        } # fin test stat
        
        rm(test_a,test_b,test_c)
        
      } # fin boucle si N suffisant
      
      Output_ALL_test <- rbind(Output_ALL_test, Output_test)
      rm(test_stat)
      
    } # fin test
    rm(Output_test, Metal_name, Metal_data_F2,Metal_data_F1)
  } # fin boucle metal
  Output_ALL_test[,6:7] <- round(Output_ALL_test[,6:7], 3)
  rm(Metal_data, Species_name)
} # Fin boucle espèce

rm(SPP,i,j)

# Plot significant differences
corr_val <- Output_ALL_test %>% 
  select(Metal,Species,P_val) %>% 
  mutate(Species = as.character(Species),
         Species = ifelse(Species == "LOJ", "Longlegged\nspiny lobster", Species),
         Species = ifelse(Species == "NUP", "Pronghorn\nspiny lobster", Species),
         Species = ifelse(Species == "NUV", "Painted\nspiny lobster", Species),
         Species = factor(Species, levels = c("Pronghorn\nspiny lobster","Longlegged\nspiny lobster","Painted\nspiny lobster"))) %>% 
  spread(Metal,P_val)
rownames(corr_val) <- corr_val$Species
corr_val <- corr_val %>% 
  select(-Species)
corr_val <- as.matrix(corr_val)

sig_val <- Output_ALL_test %>% 
  select(Metal,Species,P_val) %>% 
  mutate(Species = as.character(Species),
         Species = ifelse(Species == "LOJ", "Longlegged\nspiny lobster", Species),
         Species = ifelse(Species == "NUP", "Pronghorn\nspiny lobster", Species),
         Species = ifelse(Species == "NUV", "Painted\nspiny lobster", Species),
         Species = factor(Species, levels = c("Pronghorn\nspiny lobster","Longlegged\nspiny lobster","Painted\nspiny lobster"))) %>% 
  spread(Metal,P_val)
rownames(sig_val) <- sig_val$Species
sig_val <- sig_val %>% 
  select(-Species)
sig_val <- as.matrix(sig_val)

corrplot(corr_val,
         method = "color",
         addCoef.col = "black",
         tl.col="black", tl.srt = 0,
         p.mat = sig_val, sig.level = 0.05, insig = "blank",
         col = viridis(n = 8),
         diag = T)

rm(Output_ALL_test)
rm(corr_val,sig_val)


## 3 / Computing isotopic niches for period and each species (only NUP & LOJ) ##############################################################################################

# 3.1 / P. penicillatus - Formatting data: data_niche

data_niche_SI_NUP <- LOB_data %>% 
  filter(c_sp_fao == "NUP") %>% 
  select(bleaching_event,d13C,d15N) %>% 
  mutate(bleaching_event = ifelse(bleaching_event == "Before", "BEF","AFT"),
         bleaching_event = factor(bleaching_event, levels = c("BEF","AFT")))

# Calculate mean d13C and d15N values for each period of habitat degradation
data_niche_SI_NUP %>% 
  group_by(bleaching_event) %>% 
  summarise(mean_d13C = mean(d13C),
            sd_d13C = sd(d13C),
            mean_d15N = mean(d15N),
            sd_d15N = sd(d15N))

# Test for significant difference in mean d13C or mean d15N between periods of habitat degradation
fligner.test(data_niche_SI_NUP$d13C ~ data_niche_SI_NUP$bleaching_event) # Homoscedasticity
shapiro.test(data_niche_SI_NUP$d13C) # Normality of residuals
t.test(data_niche_SI_NUP$d13C ~ data_niche_SI_NUP$bleaching_event)

fligner.test(data_niche_SI_NUP$d15N ~ data_niche_SI_NUP$bleaching_event) # Homoscedasticity
shapiro.test(data_niche_SI_NUP$d15N) # Normality of residuals
t.test(data_niche_SI_NUP$d15N ~ data_niche_SI_NUP$bleaching_event)


# 3.2 / P. penicillatus - Calculation of posterior distribution of (mu, Sigma) for each period of habitat degradation and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_SI_NUP), data_niche_SI_NUP$bleaching_event,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_SI_NUP[ii,2:3]))


# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# 3.3 / P. penicillatus - Niche overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("BEF" = `BEF.95%`, "AFT" = `AFT.95%`) %>% 
  mutate(SpA = c("BEF","AFT")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 3)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 5, 7))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_SI_NUP <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("BEF","AFT")),
         SpB = factor(SpB, levels = c("BEF","AFT")),
         Comp = factor(Comp, levels = c("BEF in AFT","AFT in BEF")))

rm(over.stat,over.mean,over.ICinf)


# 3.4 / P. penicillatus - Extraction of ellipses for plotting - Mean and small/large

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
data_niche <- data_niche_SI_NUP %>% 
  rename(SPP = bleaching_event)
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
  MU1 <- mean(Ellipses_sel$d13C)
  MU2 <- mean(Ellipses_sel$d15N)
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

rm(Sp_names,data_niche,ii)

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
ell.coord_MEAN_SI_NUP <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("BEF","AFT")))

ell.coord_CIinf_SI_NUP <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("BEF","AFT")))

ell.coord_CIsup_SI_NUP <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("BEF","AFT")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


# 3.5 / P. penicillatus - Plot ellipses and overlaps

plot1 <- data_niche_SI_NUP %>% 
  group_by(bleaching_event) %>% 
  mutate(mean_d13C = mean(d13C, na.rm = TRUE),
         se_d13C = sd(d13C, na.rm = TRUE),
         mean_d15N = mean(d15N, na.rm = TRUE),
         se_d15N = sd(d15N, na.rm = TRUE)) %>% 
  ungroup() %>% 
  ggplot(aes(x = d13C, y = d15N))+
  geom_point(aes(fill = bleaching_event), size = 2.5, shape = 21, color = "grey50") + 
  geom_polygon(data = ell.coord_MEAN_SI_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_SI_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_SI_NUP, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_errorbar(aes(x = mean_d13C, ymin = mean_d15N-se_d15N, ymax = mean_d15N+se_d15N), color = "grey50", width=.1)+
  geom_errorbarh(aes(y = mean_d15N, xmin = mean_d13C-se_d13C, xmax = mean_d13C+se_d13C), color = "grey50", height=.1) +
  geom_point(aes(x = mean_d13C, y = mean_d15N, fill = bleaching_event), shape = 23, color = 'grey50', size=5) +
  scale_color_manual(values = c("#157BBE","#113658"))+
  scale_fill_manual(values = c("#157BBE","#113658"))+
  theme_bw() + 
  labs(x = "d13C", y = "d15N", color = "Species", title = "Pronghorn spiny lobster")+
  theme(axis.title.x = element_text(face = "italic"), 
        axis.title.y = element_text(face = "italic"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

plot2 <- output_overlap_SI_NUP %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


# 3.6 / P. longipes - Formatting data: data_niche

data_niche_SI_LOJ <- LOB_data %>% 
  filter(c_sp_fao == "LOJ") %>% 
  select(bleaching_event,d13C,d15N) %>% 
  mutate(bleaching_event = ifelse(bleaching_event == "Before", "BEF","AFT"),
         bleaching_event = factor(bleaching_event, levels = c("BEF","AFT")))

# Calculate mean d13C and d15N values for each period of habitat degradation
data_niche_SI_LOJ %>% 
  group_by(bleaching_event) %>% 
  summarise(mean_d13C = mean(d13C),
            sd_d13C = sd(d13C),
            mean_d15N = mean(d15N),
            sd_d15N = sd(d15N))

# Test for significant difference in mean d13C or mean d15N between periods of habitat degradation
fligner.test(data_niche_SI_LOJ$d13C ~ data_niche_SI_LOJ$bleaching_event) # Homoscedasticity
shapiro.test(data_niche_SI_LOJ$d13C) # Normality of residuals
t.test(data_niche_SI_LOJ$d13C ~ data_niche_SI_LOJ$bleaching_event)

fligner.test(data_niche_SI_LOJ$d15N ~ data_niche_SI_LOJ$bleaching_event) # Homoscedasticity
shapiro.test(data_niche_SI_LOJ$d15N) # Normality of residuals
t.test(data_niche_SI_LOJ$d15N ~ data_niche_SI_LOJ$bleaching_event)


# 3.7 / P. longipes - Calculation of posterior distribution of (mu, Sigma) for each period of habitat degradation and ellipse size calculation

nsamples <- 1000
fish.par <- tapply(1:nrow(data_niche_SI_LOJ), data_niche_SI_LOJ$bleaching_event,
                   function(ii) niw.post(nsamples = nsamples, X = data_niche_SI_LOJ[ii,2:3]))


# Niche sizes for each species - Calculation of the sizes of all ellipses
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


# 3.8 / P. longipes - Niche overlap calculation

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                         0.99))
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
over.ICinf <- apply(over.stat, c(1:2, 4), quantile, probs=c(0.025, 0.975), na.rm=T) * 100

over.stat <- as.data.frame(over.mean)
over.stat <- over.stat[,1:2]
over.stat <- over.stat %>% 
  rename("BEF" = `BEF.95%`, "AFT" = `AFT.95%`) %>% 
  mutate(SpA = c("BEF","AFT")) %>% 
  gather(SpB, overlap, -SpA) %>% 
  filter(!is.na(overlap)) %>% 
  mutate(Comp = paste0(SpA," in ",SpB)) %>% 
  mutate(overlap = round(overlap, 2))

over.ICinf <- as.data.frame(over.ICinf)
over.ICinf <- over.ICinf[,1:4]
over.ICinf <- as.data.frame(t(over.ICinf))
over.ICinf$SpA <- as.factor(substr(rownames(over.ICinf), 1, 3)) 
over.ICinf$SpB <- as.factor(substr(rownames(over.ICinf), 5, 7))
over.ICinf <- over.ICinf %>%
  rename(over_Q1 = `2.5%`, over_Q2 = `97.5%`) %>% 
  filter(!is.na(over_Q1)) %>%
  mutate(Comp = paste0(SpA," in ",SpB))

output_overlap_SI_LOJ <- over.stat %>%
  left_join(over.ICinf, by = c("SpA","SpB","Comp")) %>% 
  select(SpA,SpB,Comp,overlap,over_Q1,over_Q2) %>% 
  mutate(SpA = factor(SpA, levels = c("BEF","AFT")),
         SpB = factor(SpB, levels = c("BEF","AFT")),
         Comp = factor(Comp, levels = c("BEF in AFT","AFT in BEF")))

rm(over.stat,over.mean,over.ICinf)


# 3.9 / P. longipes - Extraction of ellipses for plotting - Mean and small/large

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
data_niche <- data_niche_SI_LOJ %>% 
  rename(SPP = bleaching_event)
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
  MU1 <- mean(Ellipses_sel$d13C)
  MU2 <- mean(Ellipses_sel$d15N)
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

rm(Sp_names,data_niche,ii)

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
ell.coord_MEAN_SI_LOJ <- as.data.frame(
  ell.coord_all_Ntrue %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("BEF","AFT")))

ell.coord_CIinf_SI_LOJ <- as.data.frame(
  ell.coord_inf_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("BEF","AFT")))

ell.coord_CIsup_SI_LOJ <- as.data.frame(
  ell.coord_sup_all %>%
    group_by(SPP,N_true) %>% 
    summarise(x_mean = mean(x), y_mean = mean(y)
    )
) %>% mutate(SPP = factor(SPP, levels = c("BEF","AFT")))

rm(ell.coord_all,ell.coord_inf_all,ell.coord_sup_all)
rm(ell.coord_mean95_all)
rm(ell.coord_all_Ntrue)
rm(nsamples)


# 3.10 / P. longipes - Plot ellipses and overlaps

plot3 <- data_niche_SI_LOJ %>% 
  group_by(bleaching_event) %>% 
  mutate(mean_d13C = mean(d13C, na.rm = TRUE),
         se_d13C = sd(d13C, na.rm = TRUE),
         mean_d15N = mean(d15N, na.rm = TRUE),
         se_d15N = sd(d15N, na.rm = TRUE)) %>% 
  ungroup() %>% 
  ggplot(aes(x = d13C, y = d15N))+
  geom_point(aes(fill = bleaching_event), size = 2.5, shape = 21, color = "grey50") + 
  geom_polygon(data = ell.coord_MEAN_SI_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, size = 1.2)+
  geom_polygon(data = ell.coord_CIinf_SI_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_polygon(data = ell.coord_CIsup_SI_LOJ, aes(x = x_mean, y = y_mean, color = SPP), fill = NA, linetype = "dashed", size = 0.8)+
  geom_errorbar(aes(x = mean_d13C, ymin = mean_d15N-se_d15N, ymax = mean_d15N+se_d15N), color = "grey50", width=.1)+
  geom_errorbarh(aes(y = mean_d15N, xmin = mean_d13C-se_d13C, xmax = mean_d13C+se_d13C), color = "grey50", height=.1) +
  geom_point(aes(x = mean_d13C, y = mean_d15N, fill = bleaching_event), shape = 23, color = 'grey50', size=5) +
  scale_color_manual(values = c("#157BBE","#113658"))+
  scale_fill_manual(values = c("#157BBE","#113658"))+
  theme_bw() + 
  labs(x = "d13C", y = "d15N", color = "Species", title = "Longlegged spiny lobster")+
  theme(axis.title.x = element_text(face = "italic"), 
        axis.title.y = element_text(face = "italic"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

plot4 <- output_overlap_SI_LOJ %>% 
  ggplot(aes(x = Comp))+
  geom_point(aes(y = overlap), size = 3, na.rm = TRUE)+
  geom_errorbar(aes(ymin = over_Q1, ymax = over_Q2), na.rm = TRUE, width = 0.2)+
  labs(y = "Probability of overlap (%)")+
  ylim(0,100)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")


# 3.11 / Plot all on same figure

ggarrange(plot1,plot2,plot3,plot4, ncol = 2, nrow = 2,
          align = "hv", widths = c(1,0.5,1,0.5),
          labels = c("A.","B.","C.","D."))

rm(plot1,plot2,plot3,plot4)
rm(data_niche_SI_NUP,data_niche_SI_LOJ)
rm(output_overlap_SI_NUP,output_overlap_SI_LOJ,
   ell.coord_MEAN_SI_NUP,ell.coord_MEAN_SI_LOJ,
   ell.coord_CIinf_SI_NUP,ell.coord_CIinf_SI_LOJ,
   ell.coord_CIsup_SI_NUP,ell.coord_CIsup_SI_LOJ)



### VIII // Effect of trophic ecology on TE bioaccumlation - Correlations between SI/FA and TE ##############################################################################################

# Selection of TE data
data_metaux <- LOB_data %>% 
  select(As,Cd,Cr,Cu,Hg,Mn,Ni,Pb,Se,Zn)

# Selection of tracers data
data_tracer <- LOB_data %>% 
  select(d13C,d15N,c18_1w9_p,c20_4w6_p,
         c20_5w3_p,c22_6w3_p,c20_1w7_p,c20_1w9_p,c22_2i_p,c22_2j_p,
         iso_c17_p,c18_DMA_p)

Output_ALL_cor <- NULL

# Loop for the calculation of correlations for all 3 species
for (j in 1 : length(colnames(data_tracer))){ # tracer loop
  
  # Sélection du traceur trophique
  Trophic_tracer <- data_tracer[,j]
  Tracer_name <- colnames(data_tracer)[j]
  
  
  for (i in 1: length(colnames(data_metaux))){ # metal loop
    
    Metal_name <- colnames(data_metaux)[i]
    Metal_data <- data_metaux[,i]
    
    test_a <- shapiro.test(Trophic_tracer)
    test_b <-shapiro.test(Metal_data)
    
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
}

rm(i,j)
rm(data_tracer,data_metaux)

# Creation dataframes with correlation values and significance values
cor_data <- Output_ALL_cor %>% 
  select(Metal,Trophic_tracer,Cor_val) %>% 
  spread(Metal,Cor_val)
rownames(cor_data) <- cor_data$Trophic_tracer
cor_data <- cor_data %>% 
  select(-Trophic_tracer)
cor_data <- as.matrix(cor_data)

sig_data <- Output_ALL_cor %>% 
  select(Metal,Trophic_tracer,P_val) %>% 
  spread(Metal,P_val)
rownames(sig_data) <- sig_data$Trophic_tracer
sig_data <- sig_data %>% 
  select(-Trophic_tracer)
sig_data <- as.matrix(sig_data)

# Correlation plot
corrplot(cor_data, method = "color",
         addCoef.col = "black", #add correlaction coeff
         tl.col="black", tl.srt = 0, #change text x and y
         p.mat = sig_data, sig.level = 0.05, insig = "blank") #add significance of correlations

rm(Output_ALL_cor)
rm(cor_data,sig_data)
