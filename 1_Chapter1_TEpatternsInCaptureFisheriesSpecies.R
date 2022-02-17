##--------------------------------------------------------------------------------------------------------
## SCRIPT : Allows to reproduce all analyses of Chapter 1 dedicated to the identification of trace
##          element concentrations patterns among Seychelles capture fisheries species.
##          This script includes the following analyses :
##            - clustering on trace element profiles
##            - computing of ANOVA/Kruskal-Wallis tests and associated post-hoc tests
##            - plotting the composition of each group of species inferred from trace element profiles
##            - plotting d13C and d15N values in each group of species
##            - plotting nMDS using trace element profiles, stable isotope values and organisms length to
##              investigate intragroup variability
##
## As part of :
##        Magali SABINO PhD - "Bioaccumulation of trace elements in Seychelles marine food webs"
##
## Note : Some of these results were used in Sabino et al. "The role of tropical capture fisheries in trace element delivery for a Small Island Developing State community, the Seychelles"
##
## Author : Magali Sabino
## First created : 2022-01-17
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
lapply(c("tidyverse", "openxlsx", "vegan", "ggdendro", "dendextend",
         "NbClust", "factoextra", "FSA", "ggpubr", "ggcorrplot",
         "corrplot","RColorBrewer","ggrepel"),
       library, character.only = TRUE)



### 1 // Identification of trace element concentration patterns among capture fisheries species ##########################################################################################

## 1 / Clustering using mean trace element concentrations

# Data preparation
data_cluster <- data_TE_SI %>% 
  filter(!is.na(Cu_stat),
         !is.na(THg_stat)) %>% 
  select(english_name,TAs_stat,Cd_stat,Cu_stat,Fe_stat,THg_stat,Mn_stat,
         Ni_stat,Se_stat,Zn_stat) %>% 
  rename(As = TAs_stat, Cd = Cd_stat, Cu = Cu_stat, Fe = Fe_stat,
         Hg = THg_stat, Mn = Mn_stat, Ni = Ni_stat,
         Se = Se_stat, Zn = Zn_stat) %>% 
  gather(metal, value, -english_name) %>% 
  group_by(english_name, metal) %>% 
  summarise(mean = mean(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  spread(metal, mean)

# Scaling data to reduce effect of value gaps between TE concentrations
clustering <- as.matrix(data_cluster[,2:10])
rownames(clustering) <- data_cluster$english_name
clustering <- scale(clustering)

# Need of data translation to avoid negative data
Linear_positiv_translation <- function(matrix){
  result <- matrix + abs(min(matrix, na.rm = TRUE)) + 0.1
}

clustering <- Linear_positiv_translation(clustering)
clustering_mtx <- clustering
clustering <- as.data.frame(clustering)

# Clustering
cluster_dist_mtx <- vegdist(clustering, method = "euclidean")
hc <- hclust(cluster_dist_mtx, method = "ward.D2")

# Formatting data for plotting
graph.cluster <- as.dendrogram(hc) # Build dendogram object from hclust results
graph.cluster <- dendro_data(graph.cluster, type = "rectangle") # Extract the data for rectangle lines
clusters <- data_cluster %>% 
  select(english_name) %>% 
  rename(label = english_name) # Extract infos on species
graph.cluster[["labels"]] <- merge(graph.cluster[["labels"]], clusters, by = "label") # Merge infos with labels of dendogram

# Determine the optimal number of clusters
res.nbclust <- NbClust(clustering_mtx, distance = "euclidean",
                       min.nc = 2, max.nc = 10,
                       method = "ward.D2", index = "all")

fviz_nbclust(res.nbclust, ggtheme = theme_minimal())

# Determination of clusters
clust_TE <- as.data.frame(cutree(hc, k = 4))
clust_TE$label <- row.names(clust_TE)
colnames(clust_TE)[1] <- c("cluster")
clust_TE$cluster <- as.character(clust_TE$cluster)
graph.cluster[["labels"]] <- merge(graph.cluster[["labels"]], clust_TE, by = "label") # Merge infos with labels of dendogram

# Plotting dendrogram
dendro.plot <- ggplot()+
  geom_segment(data = segment(graph.cluster), aes(x = y, y = x, xend = yend, yend = xend))+
  geom_text(data = label(graph.cluster), aes(x = y, y = x, label = label, hjust = 0,
                                             color = cluster,
                                             angle = 0), size = 3)+
  scale_color_manual(values = c("#ED7340","#FBAB19","#88BB9A","#D2062E","#048C7F"))+
  scale_x_reverse()+
  scale_y_reverse()+
  labs(x = "Weight", color = "Cluster")+
  theme(axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "italic"),
    legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank())

# Need to reorganise the clusters (cluster are not number according
# to their order of appearance on the dendrogram)
clust_TE <- clust_TE %>%
  mutate(new_cluster = cluster,
         new_cluster = ifelse(cluster == "3", "4", new_cluster),
         new_cluster = ifelse(cluster == "4", "3", new_cluster))

rm(data_cluster, clustering, hc, cluster_dist_mtx, clusters, 
   clustering_mtx, res.nbclust)


## 2 / Plotting trace element profiles and mean values for each cluster

# ANOVA/Kruskal-Wallis to determine which cluster has the highest or lowest concentration
data <- data_TE_SI %>% 
  filter(!is.na(Cu_stat),
         !is.na(THg_stat)) %>% 
  rename(label = english_name) %>% 
  left_join(clust_TE, by = "label") %>% 
  select(new_cluster,Cd_stat,Cu_stat,Fe_stat,Mn_stat,
         Ni_stat,Se_stat,TAs_stat,Zn_stat,THg_stat) %>% 
  rename(As = TAs_stat, Cd = Cd_stat, Cu = Cu_stat, Fe = Fe_stat,
         Hg = THg_stat, Mn = Mn_stat, Ni = Ni_stat,
         Se = Se_stat, Zn = Zn_stat) %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2","3","4")))

metals <- colnames(data)[2:10]
output_all_N <- as.data.frame(unique(data$new_cluster))
names(output_all_N) <- "new_cluster"
output_all_test <- NULL
output_all_posthoc <- NULL

for (ii in 1:length(metals)){
  #ii = 2
  # Création des variables
  metal_name <- metals[ii]
  data_test <- data %>% 
    select(new_cluster,metals[ii])
  data_test <- data_test[!is.na(data_test[2]),]
  data_test <- as.data.frame(data_test)
  
  output_test <- data.frame(Metal = metal_name,
                            Test = NA,
                            P_val = NA,
                            stat_val = NA)
  
  # Calcul du nombre de données par espèce
  output_N <- data_test %>% 
    group_by(new_cluster) %>% 
    summarise(N = n())
  colnames(output_N)[2] <- metal_name
  output_all_N <- output_all_N %>% 
    left_join(output_N, by = "new_cluster")
  
  # Tests stats en fonction du nbre de samples
  if (min(output_N[,2]) <= 5) {
    test_stat <- kruskal.test(data_test[,2] ~ data_test[,1])
    
    output_test$Test <- "Kruskal"
    output_test$P_val <- test_stat$p.value
    output_test$stat_val <- test_stat$statistic
    
    output_all_test <- rbind(output_all_test, output_test)
    
    # Post hoc
    if (output_test$P_val < 0.05) {
      test <- dunnTest(data_test[,2] ~ data_test[,1], method = "bh")
      
      test <- as.data.frame(test[[2]])
      #test <- test[!(test[,4]) >= 0.05,]
      test[,4] <- round(test[,4], 3)
      
      output_all_posthoc[[ii]] <- test
      names(output_all_posthoc)[ii] <- paste(metal_name,"Dunn", sep = ".")
      
      rm(test)
    } else {
      output_all_posthoc[[ii]] <- NULL
      names(output_all_posthoc)[ii] <- paste(metal_name,"Dunn", sep = ".")
    }
    
    rm(test_stat)
  } else {
    
    # Test homoscedasticité & normalité des résidus
    test_a <- fligner.test(data_test[,2] ~ data_test[,1]) # Homoscedasticity
    a1 <- aov(data_test[,2] ~ data_test[,1])
    test_b <- shapiro.test(resid(a1)) # Normality
    
    if (test_a$p.value < 0.05 | test_b$p.value < 0.05){
      # Kruskal
      test_stat <- kruskal.test(data_test[,2] ~ data_test[,1])
      
      output_test$Test <- "Kruskal"
      output_test$P_val <- test_stat$p.value
      output_test$stat_val <- test_stat$statistic
      
      output_all_test <- rbind(output_all_test, output_test)
      
      # Post hoc
      if (output_test$P_val < 0.05) {
        test <- dunnTest(data_test[,2] ~ data_test[,1], method = "bh")
        
        test <- as.data.frame(test[[2]])
        #test <- test[!(test[,4]) >= 0.05,]
        test[,4] <- round(test[,4], 3)
        
        output_all_posthoc[[ii]] <- test
        names(output_all_posthoc)[ii] <- paste(metal_name,"Dunn", sep = ".")
        
        rm(test)
      } else {
        output_all_posthoc[[ii]] <- NULL
        names(output_all_posthoc)[ii] <- paste(metal_name,"Dunn", sep = ".")
      }
      
      rm(test_stat)
    } else {
      # ANOVA
      lmFA <- lm(data_test[,2] ~ data_test[,1])
      test_stat <- anova(lmFA)
      rm(lmFA)
      
      output_test$Test <- "ANOVA"
      output_test$P_val <- test_stat[1,5]
      output_test$stat_val <- test_stat[1,4]
      
      output_all_test <- rbind(output_all_test, output_test)
      
      # Post-hoc
      if (output_test$P_val < 0.05) {
        test <- TukeyHSD(a1, 'data_test$new_cluster', conf.level = 0.95)
        output_all_posthoc[[ii]] <- as.data.frame(test[[2]])
        names(output_all_posthoc)[ii] <- paste(metal_name,"Tuke", sep = ".")
        
        rm(test)
      } else {
        output_all_posthoc[[ii]] <- NULL
        names(output_all_posthoc)[ii] <- paste(metal_name,"Tuke", sep = ".")
      }
      
      rm(test_stat)
    }
    
    rm(a1,test_a,test_b)
  }
  
  rm(metal_name,data_test,output_test,output_N)
}

output_all_test$P_val <- round(output_all_test$P_val, 5)

rm(ii)
rm(output_all_N,output_all_posthoc,output_all_test,mean)
rm(metals)


# Heatmap plot
heatmap.plot <- clust_TE %>% 
  # We define a value for each cluster and each trace element:
  # - if highest concentration among all clusters, then value = 1
  # - if lowest concentration among all clusters, then value = -1
  # - if intermediate concentration, then value = 0
  mutate(As = ifelse(new_cluster == "1", "1", "-1"),
         Cd = ifelse(new_cluster %in% c("1","3"), "1", "-1"),
         Cd = ifelse(new_cluster == "4", "0", Cd),
         Cu = ifelse(new_cluster == "1", "1", "-1"),
         Cu = ifelse(new_cluster == "4", "0", Cu),
         Fe = ifelse(new_cluster == "4", "1", "0"),
         Fe = ifelse(new_cluster == "1", "-1", Fe),
         Hg = ifelse(new_cluster == "3", "1", "0"),
         Hg = ifelse(new_cluster == "1", "-1", Hg),
         Mn = ifelse(new_cluster == "1", "1", "0"),
         Mn = ifelse(new_cluster == "3", "-1", Mn),
         Ni = ifelse(new_cluster %in% c("1","4"), "1", "0"),
         Ni = ifelse(new_cluster == "3", "-1", Ni),
         Se = ifelse(new_cluster == "3", "1", "0"),
         Se = ifelse(new_cluster == "1", "-1", Se),
         Zn = ifelse(new_cluster == "1", "1", "0"),
         Zn = ifelse(new_cluster == "2", "-1", Zn)) %>% 
  select(-new_cluster,-cluster) %>% 
  gather(metal, value, -label) %>% 
  # We define the levels of species and of values to plot
  mutate(label = factor(label, levels = c("Malabar trevally","Bluefin trevally",
                                          "Bludger","Blue-lined large-eye bream",
                                          "Yellowtail emperor","Two-spot red snapper",
                                          "Sky emperor","Slender emperor",
                                          "Golden trevally","Humpback red snapper",
                                          "Indian mackerel","Pickhandle barracuda",
                                          "Spangled emperor","Yellowspotted trevally",
                                          "Little tunny(=Atl.black skipj)","Bigeye trevally",
                                          "Pink ear emperor","Blackeye emperor",
                                          "Humphead snapper","Shoemaker spinefoot",
                                          "Grey reef shark","Dogtooth tuna",
                                          "Swordfish","Great hammerhead",
                                          "Tiger shark","Scalloped hammerhead",
                                          "Blacktip shark","Rosy goatfish",
                                          "Blue-barred parrotfish","Streamlined spinefoot",
                                          "Spot-tail shark","Spinner shark",
                                          "Big blue octopus","Brownspotted grouper",
                                          "Green jobfish","Yellow-edged lyretail",
                                          "Dash-and-dot goatfish","Honeycomb grouper",
                                          "Elongate surgeonfish","Peacock hind",
                                          "Blacktip grouper","Emperor red snapper",
                                          "Longspine grouper","White-blotched grouper",
                                          "Deepwater longtail red snapper","Eightbar grouper",
                                          "Brown-marbled grouper","Tomato hind",
                                          "Smalltooth emperor","Bigeye snapper",
                                          "Painted spiny lobster","Longlegged spiny lobster",
                                          "Pronghorn spiny lobster","Spanner crab")),
         value = factor(value, levels = c("1","0","-1")),
         metal = factor(metal, levels = c("Cd","Hg","As","Ni","Cu","Fe","Mn","Se","Zn"))) %>% 
  ggplot()+
  geom_tile(aes(x = metal, y = label, fill = value))+
  scale_fill_manual(values = c("#3B6078","#F1F1F2","#99BFC8"))+
  scale_x_discrete(position = "top")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

# Boxplots
mean_TE_boxplot <- data_TE_SI %>% 
  rename(label = english_name) %>% 
  left_join(clust_TE, by = "label") %>% 
  filter(!is.na(new_cluster)) %>% 
  select(label,new_cluster,TAs_stat,Cd_stat,Cu_stat,Fe_stat,THg_stat,Mn_stat,
         Ni_stat,Se_stat,Zn_stat) %>% 
  rename(As = TAs_stat, Cd = Cd_stat, Cu = Cu_stat, Fe = Fe_stat,
         Hg = THg_stat, Mn = Mn_stat, Ni = Ni_stat,
         Se = Se_stat, Zn = Zn_stat) %>% 
  gather(TE, value, -new_cluster, -label) %>% 
  mutate(TE = factor(TE, levels = c("As","Cd","Cu","Fe","Hg","Mn","Ni","Se","Zn"))) %>% 
  group_by(label,TE) %>% 
  mutate(mean_sp = mean(value, na.rm = TRUE),
         sd_sp = sd(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(new_cluster,TE) %>% 
  mutate(mean = mean(value, na.rm = TRUE),
         sd = sd(value, na.rm = TRUE),
         #max = max(value, na.rm = TRUE),
         max = max(mean_sp, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(mean2 = round(mean, 2),
         sd2 = round(sd, 2),
         mean_sd = paste0(mean2," ± ",sd2),
         letter = NA,
         letter = ifelse(TE == "Cd" & new_cluster %in% c("1","3"), "a", letter),
         letter = ifelse(TE == "Cd" & new_cluster == "4", "b", letter),
         letter = ifelse(TE == "Cd" & new_cluster == "2", "c", letter),
         letter = ifelse(TE == "Hg" & new_cluster == "3", "a", letter),
         letter = ifelse(TE == "Hg" & new_cluster == "4", "b", letter),
         letter = ifelse(TE == "Hg" & new_cluster == "2", "c", letter),
         letter = ifelse(TE == "Hg" & new_cluster == "1", "d", letter),
         letter = ifelse(TE == "As" & new_cluster == "1", "a", letter),
         letter = ifelse(TE == "As" & new_cluster %in% c("2","3","4"), "b", letter),
         letter = ifelse(TE == "Ni" & new_cluster %in% c("1","4"), "a", letter),
         letter = ifelse(TE == "Ni" & new_cluster == "2", "b", letter),
         letter = ifelse(TE == "Ni" & new_cluster == "3", "c", letter),
         letter = ifelse(TE == "Cu" & new_cluster == "1", "a", letter),
         letter = ifelse(TE == "Cu" & new_cluster == "4", "b", letter),
         letter = ifelse(TE == "Cu" & new_cluster %in% c("2","3"), "c", letter),
         letter = ifelse(TE == "Fe" & new_cluster == "4", "a", letter),
         letter = ifelse(TE == "Fe" & new_cluster == "3", "b", letter),
         letter = ifelse(TE == "Fe" & new_cluster == "2", "c", letter),
         letter = ifelse(TE == "Fe" & new_cluster == "1", "d", letter),
         letter = ifelse(TE == "Mn" & new_cluster == "1", "a", letter),
         letter = ifelse(TE == "Mn" & new_cluster %in% c("4","2"), "b", letter),
         letter = ifelse(TE == "Mn" & new_cluster == "3", "c", letter),
         letter = ifelse(TE == "Se" & new_cluster == "3", "a", letter),
         letter = ifelse(TE == "Se" & new_cluster == "4", "b", letter),
         letter = ifelse(TE == "Se" & new_cluster == "2", "c", letter),
         letter = ifelse(TE == "Se" & new_cluster == "1", "d", letter),
         letter = ifelse(TE == "Zn" & new_cluster == "1", "a", letter),
         letter = ifelse(TE == "Zn" & new_cluster == "3", "b", letter),
         letter = ifelse(TE == "Zn" & new_cluster == "4", "c", letter),
         letter = ifelse(TE == "Zn" & new_cluster == "2", "d", letter)) %>% 
  ggplot()+
  #geom_point(aes(x = new_cluster, y = value, fill = new_cluster), size = 2, shape = 21)+
  geom_boxplot(aes(x = new_cluster, y = value), fill = "#D1D2D4", outlier.shape = NA)+ #
  #geom_jitter(aes(x = new_cluster, y = value, fill = new_cluster), shape = 21, position=position_jitter(0.1))+
  geom_point(aes(x = new_cluster, y = mean_sp, fill = new_cluster), shape = 21, size = 2)+
  geom_text(aes(x = new_cluster, y = max+sd, label = mean_sd))+
  geom_text(aes(x = new_cluster, y = max+(2*sd), label = letter))+
  scale_fill_manual(values = c("#ED7340","#FBAB19","#D2062E","#88BB9A"))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")+
  facet_wrap(TE ~ ., ncol = 3, nrow = 3, scales = "free")+
  theme(strip.text.x = element_text(face = "bold"))

# Combine all plots in one
ggarrange(ggarrange(dendro.plot, heatmap.plot,align = "hv", widths = c(0.7, 1)),
          mean_TE_boxplot, align = "hv", nrow = 2)

rm(dendro.plot, heatmap.plot, mean_TE_boxplot, graph.cluster)



### 2 // Composition of each cluster in terms of functional group, habitat type and diet type ##########################################################################################

## 1 / Counting the total number of species in each cluster

ntot <- clust_TE %>% 
  group_by(new_cluster) %>% 
  summarise(ntot = n())


## 2 / Plot for functional group

data_functgroup <- data_TE_SI %>%
  mutate(funct_group = as.character(funct_group),
         funct_group = ifelse(funct_group == "benthic Crustacean", "Benthic crustacean", funct_group),
         funct_group = ifelse(funct_group == "benthic Cephalopod", "Benthic cephalopod", funct_group),
         funct_group = ifelse(funct_group == "benthic Teleost fish", "Benthic teleost fish", funct_group),
         funct_group = ifelse(funct_group == "demersal Teleost fish", "Demersal teleost fish", funct_group),
         #funct_group = ifelse(english_name == "Elongate surgeonfish", "Reef-associated teleost fish", funct_group),
         #funct_group = ifelse(english_name == "Slender emperor", "Demersal teleost fish", funct_group),
         funct_group = ifelse(funct_group == "pelagic-neritic Teleost fish", "Pelagic-neritic teleost fish", funct_group),
         funct_group = ifelse(funct_group == "pelagic-neritic Elasmobranchs", "Pelagic-neritic elasmobranch", funct_group),
         funct_group = ifelse(funct_group == "pelagic-oceanic Elasmobranchs", "Pelagic-oceanic elasmobranch", funct_group),
         funct_group = ifelse(funct_group == "pelagic-oceanic Teleost fish", "Pelagic-oceanic teleost fish", funct_group),
         funct_group = factor(funct_group, levels = c("Benthic crustacean","Benthic cephalopod",
                                                      "Benthic teleost fish","Demersal teleost fish",
                                                      "Pelagic-neritic teleost fish","Pelagic-neritic elasmobranch",
                                                      "Pelagic-oceanic elasmobranch","Pelagic-oceanic teleost fish"))) %>% 
  rename(label = english_name) %>% 
  left_join(clust_TE, by = "label") %>% 
  distinct(label, .keep_all = TRUE) %>% 
  group_by(new_cluster,funct_group) %>% 
  summarise(n = n()) %>% 
  left_join(ntot, by = "new_cluster") %>% 
  mutate(percent = n*100/ntot,
         percent = percent/100) %>% 
  select(new_cluster, funct_group, percent) %>% 
  spread(funct_group,percent)
data_functgroup <- data_functgroup[,-1]
data_functgroup[is.na(data_functgroup)] <- 0
data_functgroup <- as.matrix(data_functgroup)
rownames(data_functgroup) <- c("Clust 1","Clust 2","Clust 3","Clust 4")

corrplot(data_functgroup, method = "circle",
         cl.lim = c(0,1),
         col = brewer.pal(n = 8, name = "RdYlBu"),
         tl.col="black", tl.srt = 45)


## 3 / For habitat type

data_habitat <- data_TE_SI %>%
  mutate(Habitat.3 = as.character(Habitat.3),
         Habitat.3 = ifelse(english_name %in% c("Blacktip shark","Spangled emperor"), "Rocky and coral reefs",Habitat.3),
         Habitat.3 = ifelse(english_name %in% c("Brownspotted grouper","White-blotched grouper"), "Reefs and associated habitats",Habitat.3),
         Habitat.3 = ifelse(english_name == "Blue-barred parrotfish", "Coral reefs",Habitat.3),
         Habitat.3 = factor(Habitat.3, levels = c("Rocky reefs","Coral reefs","Rocky and coral reefs","Reefs and associated habitats",
                                                  "Sandy areas","Pelagic-neritic","Benthopelagic","Epipelagic","Bathypelagic"))) %>% 
  rename(label = english_name) %>% 
  left_join(clust_TE, by = "label") %>% 
  distinct(label, .keep_all = TRUE) %>% 
  group_by(new_cluster,Habitat.3) %>% 
  summarise(n = n()) %>% 
  left_join(ntot, by = "new_cluster") %>% 
  mutate(percent = n*100/ntot,
         percent = percent/100) %>% 
  select(new_cluster, Habitat.3, percent) %>% 
  spread(Habitat.3,percent)
rownames(data_habitat) <- data_habitat$new_cluster
data_habitat <- data_habitat[,-1]
data_habitat[is.na(data_habitat)] <- 0
data_habitat <- as.matrix(data_habitat)

corrplot(data_habitat, method = "circle",
         cl.lim = c(0,1),
         col = brewer.pal(n = 8, name = "RdYlBu"),
         tl.col="black", tl.srt = 45)


## 4 / For diet type

data_diet <- data_TE_SI %>%
  rename(label = english_name) %>% 
  left_join(clust_TE, by = "label") %>% 
  mutate(feeding2 = ifelse(feeding2 == "benthivore/scavanger", "Benthivore/scavenger",feeding2),
         feeding2 = ifelse(feeding2 == "benthivore", "Benthivore",feeding2),
         feeding2 = ifelse(feeding2 == "scraper", "Scraper",feeding2),
         feeding2 = ifelse(feeding2 == "grazer", "Grazer",feeding2),
         feeding2 = ifelse(feeding2 == "benthopelagivore", "Benthopelagivore",feeding2),
         feeding2 = ifelse(feeding2 == "pelagobenthivore", "Pelagobenthivore",feeding2),
         feeding2 = ifelse(feeding2 == "planktivore", "Planktivore",feeding2),
         feeding2 = ifelse(feeding2 == "pelagivore opportunist", "Pelagivore opportunist",feeding2),
         feeding2 = factor(feeding2, levels = c("Grazer","Scraper","Benthivore/scavenger",
                                                "Benthivore","Benthopelagivore","Pelagobenthivore",
                                                "Planktivore","Pelagivore opportunist"))) %>% 
  distinct(label, .keep_all = TRUE) %>% 
  group_by(new_cluster,feeding2) %>% 
  summarise(n = n()) %>% 
  left_join(ntot, by = "new_cluster") %>% 
  mutate(percent = n*100/ntot,
         percent = percent/100) %>% 
  select(new_cluster, feeding2, percent) %>% 
  spread(feeding2,percent)
rownames(data_diet) <- data_diet$new_cluster
data_diet <- data_diet[,-1]
data_diet[is.na(data_diet)] <- 0
data_diet <- as.matrix(data_diet)


corrplot(data_diet, method = "circle",
         cl.lim = c(0,1),
         col = brewer.pal(n = 8, name = "RdYlBu"),
         tl.col="black", tl.srt = 45)


rm(data_functgroup,data_habitat,data_diet,ntot)



### 3 // Trophic guild: SI values in each cluster ##########################################################################################

## 1 / Statistical tests to determine which cluster has highest/lowest mean d13C or d15N value

data <- data_TE_SI %>% 
  filter(!is.na(Cu_stat),
         !is.na(THg_stat)) %>% 
  rename(label = english_name) %>% 
  left_join(clust_TE, by = "label") %>% 
  select(new_cluster,d13C,d15N) %>% 
  mutate(new_cluster = factor(new_cluster, levels = c("1","2","3","4")))

metals <- colnames(data)[2:3]
output_all_N <- as.data.frame(unique(data$new_cluster))
names(output_all_N) <- "new_cluster"
output_all_test <- NULL
output_all_posthoc <- NULL

for (ii in 1:length(metals)){
  #ii = 2
  # Création des variables
  metal_name <- metals[ii]
  data_test <- data %>% 
    select(new_cluster,metals[ii])
  data_test <- data_test[!is.na(data_test[2]),]
  data_test <- as.data.frame(data_test)
  
  output_test <- data.frame(Metal = metal_name,
                            Test = NA,
                            P_val = NA,
                            stat_val = NA)
  
  # Calcul du nombre de données par espèce
  output_N <- data_test %>% 
    group_by(new_cluster) %>% 
    summarise(N = n())
  colnames(output_N)[2] <- metal_name
  output_all_N <- output_all_N %>% 
    left_join(output_N, by = "new_cluster")
  
  # Tests stats en fonction du nbre de samples
  if (min(output_N[,2]) <= 5) {
    test_stat <- kruskal.test(data_test[,2] ~ data_test[,1])
    
    output_test$Test <- "Kruskal"
    output_test$P_val <- test_stat$p.value
    output_test$stat_val <- test_stat$statistic
    
    output_all_test <- rbind(output_all_test, output_test)
    
    # Post hoc
    if (output_test$P_val < 0.05) {
      test <- dunnTest(data_test[,2] ~ data_test[,1], method = "bh")
      
      test <- as.data.frame(test[[2]])
      #test <- test[!(test[,4]) >= 0.05,]
      test[,4] <- round(test[,4], 3)
      
      output_all_posthoc[[ii]] <- test
      names(output_all_posthoc)[ii] <- paste(metal_name,"Dunn", sep = ".")
      
      rm(test)
    } else {
      output_all_posthoc[[ii]] <- NULL
      names(output_all_posthoc)[ii] <- paste(metal_name,"Dunn", sep = ".")
    }
    
    rm(test_stat)
  } else {
    
    # Test homoscedasticité & normalité des résidus
    test_a <- fligner.test(data_test[,2] ~ data_test[,1]) # Homoscedasticity
    a1 <- aov(data_test[,2] ~ data_test[,1])
    test_b <- shapiro.test(resid(a1)) # Normality
    
    if (test_a$p.value < 0.05 | test_b$p.value < 0.05){
      # Kruskal
      test_stat <- kruskal.test(data_test[,2] ~ data_test[,1])
      
      output_test$Test <- "Kruskal"
      output_test$P_val <- test_stat$p.value
      output_test$stat_val <- test_stat$statistic
      
      output_all_test <- rbind(output_all_test, output_test)
      
      # Post hoc
      if (output_test$P_val < 0.05) {
        test <- dunnTest(data_test[,2] ~ data_test[,1], method = "bh")
        
        test <- as.data.frame(test[[2]])
        #test <- test[!(test[,4]) >= 0.05,]
        test[,4] <- round(test[,4], 3)
        
        output_all_posthoc[[ii]] <- test
        names(output_all_posthoc)[ii] <- paste(metal_name,"Dunn", sep = ".")
        
        rm(test)
      } else {
        output_all_posthoc[[ii]] <- NULL
        names(output_all_posthoc)[ii] <- paste(metal_name,"Dunn", sep = ".")
      }
      
      rm(test_stat)
    } else {
      # ANOVA
      lmFA <- lm(data_test[,2] ~ data_test[,1])
      test_stat <- anova(lmFA)
      rm(lmFA)
      
      output_test$Test <- "ANOVA"
      output_test$P_val <- test_stat[1,5]
      output_test$stat_val <- test_stat[1,4]
      
      output_all_test <- rbind(output_all_test, output_test)
      
      # Post-hoc
      if (output_test$P_val < 0.05) {
        test <- TukeyHSD(a1, 'data_test$new_cluster', conf.level = 0.95)
        output_all_posthoc[[ii]] <- as.data.frame(test[[2]])
        names(output_all_posthoc)[ii] <- paste(metal_name,"Tuke", sep = ".")
        
        rm(test)
      } else {
        output_all_posthoc[[ii]] <- NULL
        names(output_all_posthoc)[ii] <- paste(metal_name,"Tuke", sep = ".")
      }
      
      rm(test_stat)
    }
    
    rm(a1,test_a,test_b)
  }
  
  rm(metal_name,data_test,output_test,output_N)
}

output_all_test$P_val <- round(output_all_test$P_val, 5)

rm(ii)
rm(data)
rm(output_all_N,output_all_posthoc,output_all_test)
rm(metals)


## 2 / Plotting SI values with statistical tests results

data_TE_SI %>% 
  filter(!is.na(Cu_stat),
         !is.na(THg_stat)) %>% 
  rename(label = english_name) %>% 
  left_join(clust_TE, by = "label") %>% 
  select(label,new_cluster,d13C, d15N) %>% 
  gather(SI, value, -new_cluster, -label) %>% 
  group_by(SI,label) %>% 
  mutate(mean = mean(value, na.rm = TRUE),
         se = sd(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(SI,new_cluster) %>% 
  mutate(mean_all = mean(value, na.rm = TRUE),
         se_all = sd(value, na.rm = TRUE),
         max_mean = max(mean, na.rm = TRUE),
         max_se = max(se, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(letter = NA,
         letter = ifelse(SI == "d13C" & new_cluster == "1", "a", letter),
         letter = ifelse(SI == "d13C" & new_cluster %in% c("2","3","4"), "b", letter),
         letter = ifelse(SI == "d15N" & new_cluster == "3", "a", letter),
         letter = ifelse(SI == "d15N" & new_cluster %in% c("2","4"), "b", letter),
         letter = ifelse(SI == "d15N" & new_cluster == "1", "c", letter)) %>% 
  ggplot()+
  geom_boxplot(aes(x = new_cluster, y = value), fill = "#D1D2D4", outlier.shape = NA)+
  geom_point(aes(x = new_cluster, y = mean, fill = new_cluster), shape = 21, color = "#404041", size = 3)+
  geom_text(aes(x = new_cluster, y = max_mean+max_se, label = letter), color = "#404041")+
  labs(y = "SI value (permil)")+
  scale_fill_manual(values = c("#ED7340","#FBAB19","#D2062E","#88BB9A"))+
  scale_color_manual(values = c("#ED7340","#FBAB19","#D2062E","#88BB9A"))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face = "italic"),
        legend.position = "none")+
  facet_wrap(SI ~ ., ncol = 4, nrow = 4, scales = "free")+
  theme(strip.text.x = element_text(face = "bold"))



### 4 // Intragroup variability in TE concentration, SI values and length ##########################################################################################

## 1 / nMDS for cluster 1

data_nMDS <- data_TE_SI %>% 
  filter(!is.na(Cu_stat),
         !is.na(THg_stat),
         !is.na(length)) %>% 
  select(english_name,length,d13C,d15N,TAs_stat,Cd_stat,Cu_stat,Fe_stat,
         THg_stat,Mn_stat,Ni_stat,Se_stat,Zn_stat) %>% 
  rename(As = TAs_stat, Cd = Cd_stat, Cu = Cu_stat, Fe = Fe_stat,
         Hg = THg_stat, Mn = Mn_stat, Ni = Ni_stat,
         Se = Se_stat, Zn = Zn_stat) %>% 
  rename(label = english_name) %>% 
  left_join(clust_TE, by = "label") %>% 
  rename(english_name = label) %>% 
  filter(new_cluster == "1")


# Scaling data to reduce effect of value gaps between variables
clustering <- as.matrix(data_nMDS[,2:13])
rownames(clustering) <- data_nMDS$english_name
clustering <- scale(clustering)

# Need of data translation to avoid negative data (i.e. Bray Curtis matrix)
clustering <- Linear_positiv_translation(clustering)
clustering <- as.data.frame(clustering)
data_nMDS <- data_nMDS %>% 
  select(english_name,new_cluster)
data_nMDS <- cbind(data_nMDS,clustering)
rm(clustering)

data_nMDS <- droplevels(data_nMDS)

# nMDS test
MDS_metals <- metaMDS(data_nMDS[,3:14],distance = "bray",k = 2,try = 300)

# Data manipulation for further plotting
data.scores_metals1 <- as.data.frame(scores(MDS_metals)) # Using the scores function from vegan to extract the site scores (here samples labels) and convert to a data.frame
data.scores_metals1$english_name <- data_nMDS$english_name # To add sample numbs in data.frame for future merging


# Get variables to plot them
metal.scores1 <- as.data.frame(scores(MDS_metals, "species"))  #Using the scores function from vegan to extract the sites scores and convert to a data.frame
metal.scores1$metals <- rownames(metal.scores1) # To add species names = FA

# Plot
Text_nudge <- 1.1

nMDS_clust1 <- data.scores_metals1 %>%
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_hline(yintercept = 0, linetype = "dashed", colour = 'grey')+
  geom_vline(xintercept = 0, linetype = "dashed", colour = 'grey')+
  geom_label_repel(data = metal.scores1,
                   aes(x = NMDS1*Text_nudge, y = NMDS2*Text_nudge, label = metals),
                   size = 4,
                   color = "black") +
  geom_segment(data = metal.scores1, aes(x = 0, y = 0, 
                                         xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.22, "cm")), size = 0.7)+
  labs(title = "Cluster 1", x = "nMDS1", y = "nMDS2")+
  theme_bw() +
  theme(axis.title  = element_text(face = "italic"),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())


## 2 / nMDS for Cluster 2

data_nMDS <- data_TE_SI %>% 
  filter(!is.na(Cu_stat),
         !is.na(THg_stat),
         !is.na(length),
         !is.na(d13C)) %>% 
  select(english_name,length,d13C,d15N,TAs_stat,Cd_stat,Cu_stat,Fe_stat,
         THg_stat,Mn_stat,Ni_stat,Se_stat,Zn_stat) %>% 
  rename(As = TAs_stat, Cd = Cd_stat, Cu = Cu_stat, Fe = Fe_stat,
         Hg = THg_stat, Mn = Mn_stat, Ni = Ni_stat,
         Se = Se_stat, Zn = Zn_stat) %>% 
  rename(label = english_name) %>% 
  left_join(clust_TE, by = "label") %>% 
  rename(english_name = label) %>% 
  select(-cluster) %>%
  filter(new_cluster == "2")


# Scaling data to reduce effect of value gaps between variables
clustering <- as.matrix(data_nMDS[,2:13])
rownames(clustering) <- data_nMDS$english_name
clustering <- scale(clustering)

# Need of data translation to avoid negative data (i.e. Bray Curtis matrix)
clustering <- Linear_positiv_translation(clustering)
clustering <- as.data.frame(clustering)
data_nMDS <- data_nMDS %>% 
  select(english_name,new_cluster)
data_nMDS <- cbind(data_nMDS,clustering)
rm(clustering)

data_nMDS <- droplevels(data_nMDS)

# nMDS test
MDS_metals <- metaMDS(data_nMDS[,3:14],distance = "bray",k = 2,try = 600)

# Data manipulation for further plotting
data.scores_metals2 <- as.data.frame(scores(MDS_metals)) # Using the scores function from vegan to extract the site scores (here samples labels) and convert to a data.frame
data.scores_metals2$sample_identifier <- data_nMDS$sample_identifier # To add sample numbs in data.frame for future merging
data.scores_metals2$english_name <- data_nMDS$english_name

# Get variables to plot them
metal.scores2 <- as.data.frame(scores(MDS_metals, "species"))  #Using the scores function from vegan to extract the sites scores and convert to a data.frame
metal.scores2$metals <- rownames(metal.scores2) # To add species names = FA
head(metal.scores2)

# Plot
Text_nudge <- 1.1

nMDS_clust2 <- data.scores_metals2 %>%
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_hline(yintercept = 0, linetype = "dashed", colour = 'grey')+
  geom_vline(xintercept = 0, linetype = "dashed", colour = 'grey')+
  geom_label_repel(data = metal.scores2,
                   aes(x = NMDS1*Text_nudge, y = NMDS2*Text_nudge, label = metals),
                   size = 4,
                   color = "black") +
  geom_segment(data = metal.scores2, aes(x = 0, y = 0, 
                                         xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.22, "cm")), size = 0.7)+
  labs(title = "Cluster 2", x = "nMDS1", y = "nMDS2")+
  theme_bw() +
  theme(axis.title  = element_text(face = "italic"),
    legend.title = element_blank(),
    legend.position = "none",
    panel.grid = element_blank())


## 3 / nMDS for Cluster 3

data_nMDS <- data_TE_SI %>% 
  filter(!is.na(Cu_stat),
         !is.na(THg_stat),
         !is.na(length),
         !is.na(d13C)) %>% 
  select(english_name,length,d13C,d15N,TAs_stat,Cd_stat,Cu_stat,Fe_stat,
         THg_stat,Mn_stat,Ni_stat,Se_stat,Zn_stat) %>% 
  rename(As = TAs_stat, Cd = Cd_stat, Cu = Cu_stat, Fe = Fe_stat,
         Hg = THg_stat, Mn = Mn_stat, Ni = Ni_stat,
         Se = Se_stat, Zn = Zn_stat) %>% 
  rename(label = english_name) %>% 
  left_join(clust_TE, by = "label") %>% 
  rename(english_name = label) %>% 
  select(-cluster) %>%
  filter(new_cluster == "3")


# Scaling data to reduce effect of value gaps between variables
clustering <- as.matrix(data_nMDS[,2:13])
rownames(clustering) <- data_nMDS$english_name
clustering <- scale(clustering)

# Need of data translation to avoid negative data (i.e. Bray Curtis matrix)
clustering <- Linear_positiv_translation(clustering)
clustering <- as.data.frame(clustering)
data_nMDS <- data_nMDS %>% 
  select(english_name,new_cluster)
data_nMDS <- cbind(data_nMDS,clustering)
rm(clustering)

data_nMDS <- droplevels(data_nMDS)

# nMDS test
MDS_metals <- metaMDS(data_nMDS[,3:14],distance = "bray",k = 2,try = 300)

# Data manipulation for further plotting
data.scores_metals3 <- as.data.frame(scores(MDS_metals)) # Using the scores function from vegan to extract the site scores (here samples labels) and convert to a data.frame
data.scores_metals3$sample_identifier <- data_nMDS$sample_identifier # To add sample numbs in data.frame for future merging
data.scores_metals3$english_name <- data_nMDS$english_name

# Get variable to plot them
metal.scores3 <- as.data.frame(scores(MDS_metals, "species"))  #Using the scores function from vegan to extract the sites scores and convert to a data.frame
metal.scores3$metals <- rownames(metal.scores3) # To add species names = FA

# Plot
Text_nudge <- 1.1

nMDS_clust3 <- data.scores_metals3 %>%
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_hline(yintercept = 0, linetype = "dashed", colour = 'grey')+
  geom_vline(xintercept = 0, linetype = "dashed", colour = 'grey')+
  geom_label_repel(data = metal.scores3,
                   aes(x = NMDS1*Text_nudge, y = NMDS2*Text_nudge, label = metals),
                   size = 4,
                   color = "black") +
  geom_segment(data = metal.scores3, aes(x = 0, y = 0, 
                                         xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.22, "cm")), size = 0.7)+
  labs(title = "Cluster 3", x = "nMDS1", y = "nMDS2")+
  theme_bw() +
  theme(axis.title  = element_text(face = "italic"),
    legend.title = element_blank(),
    legend.position = "none",
    panel.grid = element_blank())


## 4 / nMDS for Cluster 4

data_nMDS <- data_TE_SI %>% 
  filter(!is.na(Cu_stat),
         !is.na(THg_stat),
         !is.na(length),
         !is.na(d13C)) %>% 
  select(english_name,length,d13C,d15N,TAs_stat,Cd_stat,Cu_stat,Fe_stat,
         THg_stat,Mn_stat,Ni_stat,Se_stat,Zn_stat) %>% 
  rename(As = TAs_stat, Cd = Cd_stat, Cu = Cu_stat, Fe = Fe_stat,
         Hg = THg_stat, Mn = Mn_stat, Ni = Ni_stat,
         Se = Se_stat, Zn = Zn_stat) %>% 
  rename(label = english_name) %>% 
  left_join(clust_TE, by = "label") %>% 
  rename(english_name = label) %>% 
  select(-cluster) %>%
  filter(new_cluster == "4")


# Scaling data to reduce effect of value gaps between variables
clustering <- as.matrix(data_nMDS[,2:13])
rownames(clustering) <- data_nMDS$english_name
clustering <- scale(clustering)

# Need of data translation to avoid negative data (i.e. Bray Curtis matrix)
clustering <- Linear_positiv_translation(clustering)
clustering <- as.data.frame(clustering)
data_nMDS <- data_nMDS %>% 
  select(english_name,new_cluster)
data_nMDS <- cbind(data_nMDS,clustering)
rm(clustering)

data_nMDS <- droplevels(data_nMDS)

# nMDS test
MDS_metals <- metaMDS(data_nMDS[,3:14],distance = "bray",k = 2,try = 300)

# Data manipulation for further plotting
data.scores_metals4 <- as.data.frame(scores(MDS_metals)) # Using the scores function from vegan to extract the site scores (here samples labels) and convert to a data.frame
data.scores_metals4$sample_identifier <- data_nMDS$sample_identifier # To add sample numbs in data.frame for future merging
data.scores_metals4$english_name <- data_nMDS$english_name

# Get variables to plot them
metal.scores4 <- as.data.frame(scores(MDS_metals, "species"))  #Using the scores function from vegan to extract the sites scores and convert to a data.frame
metal.scores4$metals <- rownames(metal.scores4) # To add species names = FA

# Plot
Text_nudge <- 1.1

nMDS_clust4 <- data.scores_metals4 %>%
  ggplot(aes(x = NMDS1, y = NMDS2))+
  geom_hline(yintercept = 0, linetype = "dashed", colour = 'grey')+
  geom_vline(xintercept = 0, linetype = "dashed", colour = 'grey')+
  geom_label_repel(data = metal.scores4,
                   aes(x = NMDS1*Text_nudge, y = NMDS2*Text_nudge, label = metals),
                   size = 4, #fontface = "bold",
                   color = "black") +
  geom_segment(data = metal.scores4, aes(x = 0, y = 0, 
                                         xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.22, "cm")), size = 0.7)+
  labs(title = "Cluster 4", x = "nMDS1", y = "nMDS2")+
  theme_bw() +
  theme(axis.title  = element_text(face = "italic"),
    legend.title = element_blank(),
    legend.position = "none",
    panel.grid = element_blank())


## 5 / Combine all plots

ggarrange(nMDS_clust1,nMDS_clust2,nMDS_clust3,nMDS_clust4,
          align = "hv", ncol = 2, nrow = 2,
          labels = c("A.","B.","C.","D."))

rm(nMDS_clust1,nMDS_clust2,nMDS_clust3,nMDS_clust4)
rm(data.scores_metals1,data.scores_metals2,data.scores_metals3,data.scores_metals4)
rm(metal.scores1,metal.scores2,metal.scores3,metal.scores4)
rm(data_nMDS,MDS_metals,Text_nudge)

rm(clust_TE)