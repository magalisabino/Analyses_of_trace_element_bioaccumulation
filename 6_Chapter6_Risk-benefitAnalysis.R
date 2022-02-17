##--------------------------------------------------------------------------------------------------------
## SCRIPT : Allows to reproduce all analyses of Chapter 6 dedicated to the risk-benefit analysis.
##          This script includes the following analyses :
##            - Calculation of mean +/- SD TE concentrations for each species
##            - Calculation of MHg:MSe and HBVSe ratios to estimate interaction capacity between Hg and Se,
##              and estimation of theoretically bioavailable Hg and Se
##            - Estimation of risk-benefit or the consumption of Seychelles capture fisheries resources,
##              using %PTI and %RDI covered by one portion (calculated for each age class of interest)
##            - Estimation of iAs concentrations using TAs concentrations
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
lapply(c("tidyverse", "openxlsx", "sp", "sf", "rnaturalearth",
         "rnaturalearthdata", "ggspatial", "rgeos", "ggpubr",
         "FSA", "corrplot", "scatterpie"),
       library, character.only = TRUE)



### I // Database creation ##############################################################################################

data_risk_benefit <- data_TE_SI %>% 
  mutate(sp_group = NA,
         sp_group = ifelse(english_name == "Big blue octopus", "Octopuses",sp_group),
         sp_group = ifelse(english_name == "Spanner crab", "Crabs", sp_group),
         sp_group = ifelse(english_name %in% c("Pronghorn spiny lobster","Painted spiny lobster","Longlegged spiny lobster"),"Spiny lobsters",sp_group),
         sp_group = ifelse(english_name %in% c("Dash-and-dot goatfish","Rosy goatfish"), "Goatfishes",sp_group),
         sp_group = ifelse(english_name == "Blue-barred parrotfish", "Parrotfishes",sp_group),
         sp_group = ifelse(english_name %in% c("Shoemaker spinefoot","Streamlined spinefoot"), "Rabbitfishes",sp_group),
         sp_group = ifelse(english_name == "Elongate surgeonfish", "Surgeonfishes",sp_group),
         sp_group = ifelse(english_name %in% c("Blackeye emperor","Blue-lined large-eye bream",
                                                "Pink ear emperor","Slender emperor",
                                                "Sky emperor","Smalltooth emperor",
                                                "Spangled emperor","Yellowtail emperor"), "Emperors",sp_group),
         sp_group = ifelse(english_name %in% c("Peacock hind","Tomato hind",
                                               "Blacktip grouper","Honeycomb grouper",
                                               "Brownspotted grouper","Brown-marbled grouper",
                                               "Eightbar grouper","White-blotched grouper",
                                               "Longspine grouper","Yellow-edged lyretail"), "Groupers",sp_group),
         sp_group = ifelse(english_name %in% c("Bigeye trevally","Bluefin trevally",
                                               "Golden trevally","Malabar trevally",
                                               "Yellowspotted trevally","Bludger"), "Trevallies",sp_group),
         sp_group = ifelse(english_name %in% c("Green jobfish","Deepwater longtail red snapper",
                                               "Two-spot red snapper","Humpback red snapper",
                                               "Bigeye snapper","Emperor red snapper",
                                               "Humphead snapper"), "Snappers",sp_group),
         sp_group = ifelse(english_name == "Pickhandle barracuda","Barracudas",sp_group),
         sp_group = ifelse(english_name == "Indian mackerel","Mackerels",sp_group),
         sp_group = ifelse(english_name %in% c("Dogtooth tuna","Little tunny(=Atl.black skipj)"), "Tunas",sp_group),
         sp_group = ifelse(english_name %in% c("Grey reef shark","Spinner shark",
                                               "Blacktip shark","Spot-tail shark",
                                               "Tiger shark","Great hammerhead",
                                               "Scalloped hammerhead"), "Sharks",sp_group),
         sp_group = ifelse(english_name == "Swordfish","Billfishes",sp_group),
         sp_group = factor(sp_group, levels = c("Octopuses","Crabs","Spiny lobsters","Goatfishes",
                                                "Parrotfishes","Rabbitfishes","Surgeonfishes","Emperors",
                                                "Groupers","Trevallies","Snappers","Barracudas",
                                                "Mackerels","Tunas","Sharks","Billfishes")),
         english_name = factor(english_name, levels = c("Big blue octopus","Spanner crab",
                                                        "Pronghorn spiny lobster","Painted spiny lobster","Longlegged spiny lobster",
                                                        "Dash-and-dot goatfish","Rosy goatfish",
                                                        "Blue-barred parrotfish","Shoemaker spinefoot","Streamlined spinefoot",
                                                        "Elongate surgeonfish",
                                                        "Blackeye emperor","Blue-lined large-eye bream",
                                                        "Pink ear emperor","Slender emperor",
                                                        "Sky emperor","Smalltooth emperor",
                                                        "Spangled emperor","Yellowtail emperor",
                                                        "Peacock hind","Tomato hind",
                                                        "Blacktip grouper","Honeycomb grouper",
                                                        "Brownspotted grouper","Brown-marbled grouper",
                                                        "Eightbar grouper","White-blotched grouper",
                                                        "Longspine grouper","Yellow-edged lyretail",
                                                        "Bigeye trevally","Bluefin trevally",
                                                        "Golden trevally","Malabar trevally",
                                                        "Yellowspotted trevally","Bludger",
                                                        "Green jobfish","Deepwater longtail red snapper",
                                                        "Two-spot red snapper","Humpback red snapper",
                                                        "Bigeye snapper","Emperor red snapper",
                                                        "Humphead snapper",
                                                        "Pickhandle barracuda","Indian mackerel",
                                                        "Dogtooth tuna","Little tunny(=Atl.black skipj)",
                                                        "Grey reef shark","Spinner shark",
                                                        "Blacktip shark","Spot-tail shark",
                                                        "Tiger shark","Great hammerhead",
                                                        "Scalloped hammerhead","Swordfish")))


### II // Mean (+/- SD) TE concentrations by species (Appendix 8.1) ##############################################################################################

names(data_risk_benefit)

TE_concentrations <- data_risk_benefit %>% 
  select(sp_group,english_name,Ag_ww,Cd_ww,Co_ww,Cr_ww,
         Cu_ww,Fe_ww,Mn_ww,Ni_ww,Pb_ww,Se_ww,TAs_ww,Zn_ww,Hg_ww) %>% 
  rename(Ag = Ag_ww, Cd = Cd_ww, Co = Co_ww, Cr = Cr_ww, Cu = Cu_ww,
         Fe = Fe_ww, Mn = Mn_ww, Ni = Ni_ww, Pb = Pb_ww, Se = Se_ww,
         As = TAs_ww, Zn = Zn_ww, Hg = Hg_ww) %>% 
  group_by(sp_group,english_name) %>% 
  summarise(Ag_mean = mean(Ag, na.rm = T), Ag_sd = sd(Ag, na.rm = T),
            Cd_mean = mean(Cd, na.rm = T), Cd_sd = sd(Cd, na.rm = T),
            Co_mean = mean(Co, na.rm = T), Co_sd = sd(Co, na.rm = T),
            Cr_mean = mean(Cr, na.rm = T), Cr_sd = sd(Cr, na.rm = T),
            Cu_mean = mean(Cu, na.rm = T), Cu_sd = sd(Cu, na.rm = T),
            Fe_mean = mean(Fe, na.rm = T), Fe_sd = sd(Fe, na.rm = T),
            Mn_mean = mean(Mn, na.rm = T), Mn_sd = sd(Mn, na.rm = T),
            Ni_mean = mean(Ni, na.rm = T), Ni_sd = sd(Ni, na.rm = T),
            Pb_mean = mean(Pb, na.rm = T), Pb_sd = sd(Pb, na.rm = T),
            Se_mean = mean(Se, na.rm = T), Se_sd = sd(Se, na.rm = T),
            As_mean = mean(As, na.rm = T), As_sd = sd(As, na.rm = T),
            Zn_mean = mean(Zn, na.rm = T), Zn_sd = sd(Zn, na.rm = T),
            Hg_mean = mean(Hg, na.rm = T), Hg_sd = sd(Hg, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Ag_mean = round(Ag_mean, 3), Ag_sd = round(Ag_sd, 3),
         Cd_mean = round(Cd_mean, 3), Cd_sd = round(Cd_sd, 3),
         Hg_mean = round(Hg_mean, 3), Hg_sd = round(Hg_sd, 3),
         Pb_mean = round(Pb_mean, 3), Pb_sd = round(Pb_sd, 3),
         As_mean = round(As_mean, 1), As_sd = round(As_sd, 1),
         Cr_mean = round(Cr_mean, 2), Cr_sd = round(Cr_sd, 2),
         Ni_mean = round(Ni_mean, 2), Ni_sd = round(Ni_sd, 2),
         Co_mean = round(Co_mean, 3), Co_sd = round(Co_sd, 3),
         Cu_mean = round(Cu_mean, 2), Cu_sd = round(Cu_sd, 2),
         Fe_mean = round(Fe_mean, 1), Fe_sd = round(Fe_sd, 1),
         Mn_mean = round(Mn_mean, 2), Mn_sd = round(Mn_sd, 2),
         Se_mean = round(Se_mean, 2), Se_sd = round(Se_sd, 2),
         Zn_mean = round(Zn_mean, 1), Zn_sd = round(Zn_sd, 1),
         Ag = paste0(Ag_mean," ± ",Ag_sd),
         Cd = paste0(Cd_mean," ± ",Cd_sd),
         Hg = paste0(Hg_mean," ± ",Hg_sd),
         Pb = paste0(Pb_mean," ± ",Pb_sd),
         As = paste0(As_mean," ± ",As_sd),
         Cr = paste0(Cr_mean," ± ",Cr_sd),
         Ni = paste0(Ni_mean," ± ",Ni_sd),
         Co = paste0(Co_mean," ± ",Co_sd),
         Cu = paste0(Cu_mean," ± ",Cu_sd),
         Fe = paste0(Fe_mean," ± ",Fe_sd),
         Mn = paste0(Mn_mean," ± ",Mn_sd),
         Se = paste0(Se_mean," ± ",Se_sd),
         Zn = paste0(Zn_mean," ± ",Zn_sd),
         # Replace (NA) by mean value and NaN by <LOQ
         Ag = ifelse(is.na(Ag_sd) == TRUE, Ag_mean, Ag),
         Cd = ifelse(is.na(Cd_sd) == TRUE, Cd_mean, Cd),
         Hg = ifelse(is.na(Hg_sd) == TRUE, Hg_mean, Hg),
         Pb = ifelse(is.na(Pb_sd) == TRUE, Pb_mean, Pb),
         As = ifelse(is.na(As_sd) == TRUE, As_mean, As),
         Cr = ifelse(is.na(Cr_sd) == TRUE, Cr_mean, Cr),
         Ni = ifelse(is.na(Ni_sd) == TRUE, Ni_mean, Ni),
         Co = ifelse(is.na(Co_sd) == TRUE, Co_mean, Co),
         Cu = ifelse(is.na(Cu_sd) == TRUE, Cu_mean, Cu),
         Fe = ifelse(is.na(Fe_sd) == TRUE, Fe_mean, Fe),
         Mn = ifelse(is.na(Mn_sd) == TRUE, Mn_mean, Mn),
         Se = ifelse(is.na(Se_sd) == TRUE, Se_mean, Se),
         Zn = ifelse(is.na(Zn_sd) == TRUE, Zn_mean, Zn),
         Ag = ifelse(Ag == "NaN", "<LOQ", Ag),
         Cd = ifelse(Cd == "NaN", "<LOQ", Cd),
         Hg = ifelse(Hg == "NaN", "<LOQ", Hg),
         Pb = ifelse(Pb == "NaN", "<LOQ", Pb),
         As = ifelse(As == "NaN", "<LOQ", As),
         Cr = ifelse(Cr == "NaN", "<LOQ", Cr),
         Ni = ifelse(Ni == "NaN", "<LOQ", Ni),
         Co = ifelse(Co == "NaN", "<LOQ", Co),
         Cu = ifelse(Cu == "NaN", "<LOQ", Cu),
         Fe = ifelse(Fe == "NaN", "<LOQ", Fe),
         Mn = ifelse(Mn == "NaN", "<LOQ", Mn),
         Se = ifelse(Se == "NaN", "<LOQ", Se),
         Zn = ifelse(Zn == "NaN", "<LOQ", Zn),
         number = as.numeric(english_name)) %>% 
  arrange(number) %>% 
  select(sp_group,english_name,Co,Cu,Fe,Mn,Se,Zn,As,Cr,Ni,Ag,Cd,Hg,Pb)

write.table(TE_concentrations,
            paste0(direction,"1_MeanTEconcentrations_bySpecies.csv"),
            sep = ";", dec = ".", row.names = F)

rm(TE_concentrations)



### III // Interaction between Hg and Se ##############################################################################################

interHgSe <- data_risk_benefit %>% 
  select(organism_identifier,english_name,Hg_ww,Se_ww) %>% 
  mutate(MHg = Hg_ww/200590000,
         MSe = Se_ww/78960000,
         MHg_MSe = MHg/MSe,
         HBVSe = ((MSe-MHg)/MSe)*(MSe+MHg),
         Se_avail = (MSe-MHg)*78960000,
         Hg_avail = MHg-MSe,
         Hg_avail = ifelse(Hg_avail <= 0, 0, Hg_avail),
         Hg_avail = Hg_avail*200590000) %>% 
  group_by(english_name) %>% 
  mutate(mean_MHg_MSe = mean(MHg_MSe, na.rm = T),
         sd_MHg_MSe = sd(MHg_MSe, na.rm = T),
         se_MHg_MSe = std(MHg_MSe),
         CIinf_MHg_MSe = mean_MHg_MSe-(1.96*se_MHg_MSe),
         CIsup_MHg_MSe = mean_MHg_MSe+(1.96*se_MHg_MSe),
         mean_HBVSe = mean(HBVSe, na.rm = T),
         sd_HBVSe = sd(HBVSe, na.rm = T),
         se_HBVSe = std(HBVSe),
         CIinf_HBVSe = mean_HBVSe-(1.96*se_HBVSe),
         CIsup_HBVSe = mean_HBVSe+(1.96*se_HBVSe),
         mean_Se_avail = mean(Se_avail, na.rm = T),
         se_Se_avail = std(Se_avail),
         CIinf_Se_avail = mean_Se_avail-(1.96*se_Se_avail),
         CIsup_Se_avail = mean_Se_avail+(1.96*se_Se_avail),
         mean_Hg_avail = mean(Hg_avail, na.rm = T),
         se_Hg_avail = std(Hg_avail),
         CIinf_Hg_avail = mean_Hg_avail-(1.96*se_Hg_avail),
         CIsup_Hg_avail = mean_Hg_avail+(1.96*se_Hg_avail)) %>% 
         ungroup()

infos <- data_risk_benefit %>% select(sp_group,english_name) %>% distinct(english_name, .keep_all = TRUE)

names(interHgSe)

interHgSe_propre <- interHgSe %>% 
  left_join(infos, by = "english_name") %>% 
  mutate(MHg_MSe = paste0(round(mean_MHg_MSe,3)," ± ",round(sd_MHg_MSe,3)," [",round(CIinf_MHg_MSe,3),"-",round(CIsup_MHg_MSe,3),"]"),
         MHg_MSe = ifelse(is.na(sd_MHg_MSe) == TRUE, as.character(round(mean_MHg_MSe, 3)), MHg_MSe),
         HBVSe = paste0(round(mean_HBVSe,10)," ± ",round(sd_HBVSe,10)," [",round(CIinf_HBVSe,10),"-",round(CIsup_HBVSe,10),"]"),
         HBVSe = ifelse(is.na(sd_HBVSe) == TRUE, as.character(round(mean_HBVSe, 10)), HBVSe),
         Se_avail = paste0(round(mean_Se_avail,2)," [",round(CIinf_Se_avail,2),"-",round(CIsup_Se_avail,2),"]"),
         Se_avail = ifelse(is.na(CIinf_Se_avail) == TRUE, as.character(round(mean_Se_avail, 2)), Se_avail),
         Hg_avail = paste0(round(mean_Hg_avail,3)," [",round(CIinf_Hg_avail,3),"-",round(CIsup_Hg_avail,3),"]"),
         Hg_avail = ifelse(is.na(CIinf_Hg_avail) == TRUE, as.character(round(mean_Hg_avail, 3)), Hg_avail)) %>% 
  distinct(english_name, .keep_all = TRUE) %>% 
  mutate(number = as.numeric(english_name)) %>% 
  arrange(number) %>% 
  select(sp_group,english_name,MHg_MSe,HBVSe,Se_avail,Hg_avail)

write.table(interHgSe_propre,paste0(direction,"2_InteractionHgSe.csv"),
            sep = ";", dec = ".", row.names = F)

rm(infos,interHgSe_propre)



### IV // Risk ##############################################################################################

infos <- data_risk_benefit %>% select(sp_group,english_name) %>% distinct(english_name, .keep_all = TRUE)

## 1 / Calculation of %PTI for Se and Hg ##############################################################################################

portion_YoungChildren <- 30
portion_Children <- 72
portion_Adolescents <- 114
portion_YoungAdults_Adults <- 156

percentPTI_SeHg <- interHgSe %>% 
  ## Daily expo
  mutate(## For Se
    Se_day_YC = Se_avail*portion_YoungChildren,
    Se_day_C = Se_avail*portion_Children,
    Se_day_ADO = Se_avail*portion_Adolescents,
    Se_day_YA_A = Se_avail*portion_YoungAdults_Adults,
    ## For Hg
    Hg_day_YC = Hg_avail*portion_YoungChildren,
    Hg_day_C = Hg_avail*portion_Children,
    Hg_day_ADO = Hg_avail*portion_Adolescents,
    Hg_day_YA_A = Hg_avail*portion_YoungAdults_Adults) %>% 
  ## Daily %PTI
  mutate(## For Se
    percentPTI_Se_YC = Se_day_YC*100/90,
    percentPTI_Se_C = Se_day_C*100/150,
    percentPTI_Se_ADO = Se_day_ADO*100/280,
    percentPTI_Se_YA_A = Se_day_YA_A*100/400,
    ## For Hg
    percentPTI_Hg_YC = Hg_day_YC*100/7,
    percentPTI_Hg_C = Hg_day_C*100/12.6,
    percentPTI_Hg_ADO = Hg_day_ADO*100/20.6,
    percentPTI_Hg_YA = Hg_day_YA_A*100/29.7,
    percentPTI_Hg_A = Hg_day_YA_A*100/40) %>% 
  group_by(english_name) %>% 
  ## Mean [CI95%]
  summarise(## For Se
    meanPTI_Se_YC = mean(percentPTI_Se_YC, na.rm = TRUE),
    sePTI_Se_YC = std(percentPTI_Se_YC),
    CIinf_Se_YC = meanPTI_Se_YC-(1.96*sePTI_Se_YC),
    CIsup_Se_YC = meanPTI_Se_YC+(1.96*sePTI_Se_YC),
    CIinf_Se_YC = ifelse(CIinf_Se_YC < 0, 0, CIinf_Se_YC),
    meanPTI_Se_C = mean(percentPTI_Se_C, na.rm = TRUE),
    sePTI_Se_C = std(percentPTI_Se_C),
    CIinf_Se_C = meanPTI_Se_C-(1.96*sePTI_Se_C),
    CIsup_Se_C = meanPTI_Se_C+(1.96*sePTI_Se_C),
    CIinf_Se_C = ifelse(CIinf_Se_C < 0, 0, CIinf_Se_C),
    meanPTI_Se_ADO = mean(percentPTI_Se_ADO, na.rm = TRUE),
    sePTI_Se_ADO = std(percentPTI_Se_ADO),
    CIinf_Se_ADO = meanPTI_Se_ADO-(1.96*sePTI_Se_ADO),
    CIsup_Se_ADO = meanPTI_Se_ADO+(1.96*sePTI_Se_ADO),
    CIinf_Se_ADO = ifelse(CIinf_Se_ADO < 0, 0, CIinf_Se_ADO),
    meanPTI_Se_YA_A = mean(percentPTI_Se_YA_A, na.rm = TRUE),
    sePTI_Se_YA_A = std(percentPTI_Se_YA_A),
    CIinf_Se_YA_A = meanPTI_Se_YA_A-(1.96*sePTI_Se_YA_A),
    CIsup_Se_YA_A = meanPTI_Se_YA_A+(1.96*sePTI_Se_YA_A),
    CIinf_Se_YA_A = ifelse(CIinf_Se_YA_A < 0, 0, CIinf_Se_YA_A),
    ## For Hg
    meanPTI_Hg_YC = mean(percentPTI_Hg_YC, na.rm = TRUE),
    sePTI_Hg_YC = std(percentPTI_Hg_YC),
    CIinf_Hg_YC = meanPTI_Hg_YC-(1.96*sePTI_Hg_YC),
    CIsup_Hg_YC = meanPTI_Hg_YC+(1.96*sePTI_Hg_YC),
    CIinf_Hg_YC = ifelse(CIinf_Hg_YC < 0, 0, CIinf_Hg_YC),
    meanPTI_Hg_C = mean(percentPTI_Hg_C, na.rm = TRUE),
    sePTI_Hg_C = std(percentPTI_Hg_C),
    CIinf_Hg_C = meanPTI_Hg_C-(1.96*sePTI_Hg_C),
    CIsup_Hg_C = meanPTI_Hg_C+(1.96*sePTI_Hg_C),
    CIinf_Hg_C = ifelse(CIinf_Hg_C < 0, 0, CIinf_Hg_C),
    meanPTI_Hg_ADO = mean(percentPTI_Hg_ADO, na.rm = TRUE),
    sePTI_Hg_ADO = std(percentPTI_Hg_ADO),
    CIinf_Hg_ADO = meanPTI_Hg_ADO-(1.96*sePTI_Hg_ADO),
    CIsup_Hg_ADO = meanPTI_Hg_ADO+(1.96*sePTI_Hg_ADO),
    CIinf_Hg_ADO = ifelse(CIinf_Hg_ADO < 0, 0, CIinf_Hg_ADO),
    meanPTI_Hg_YA = mean(percentPTI_Hg_YA, na.rm = TRUE),
    sePTI_Hg_YA = std(percentPTI_Hg_YA),
    CIinf_Hg_YA = meanPTI_Hg_YA-(1.96*sePTI_Hg_YA),
    CIsup_Hg_YA = meanPTI_Hg_YA+(1.96*sePTI_Hg_YA),
    CIinf_Hg_YA = ifelse(CIinf_Hg_YA < 0, 0, CIinf_Hg_YA),
    meanPTI_Hg_A = mean(percentPTI_Hg_A, na.rm = TRUE),
    sePTI_Hg_A = std(percentPTI_Hg_A),
    CIinf_Hg_A = meanPTI_Hg_A-(1.96*sePTI_Hg_A),
    CIsup_Hg_A = meanPTI_Hg_A+(1.96*sePTI_Hg_A),
    CIinf_Hg_A = ifelse(CIinf_Hg_A < 0, 0, CIinf_Hg_A))

## 2 / Calculation of %PTI for Cd, Cu, Fe, Mn, Zn##############################################################################################

percentPTI <- data_risk_benefit %>% 
  select(english_name, Cu_ww, Fe_ww, Mn_ww, Zn_ww, Cd_ww) %>% 
  ## For Cd
  mutate(Cd_day_YC = Cd_ww*portion_YoungChildren, # According to EPA-USA, 1 serving = 1 ounce = 28.3g for children age 2
         Cd_day_C = Cd_ww*portion_Children, # According to EPA-USA, 1 serving progressively increases from 1 to 4 ounces from 2 years age to 11 years
         Cd_day_ADO = Cd_ww*portion_Adolescents,
         Cd_day_YA_A = Cd_ww*portion_YoungAdults_Adults,
         percentPTI_Cd_YC = Cd_day_YC*100/10,
         percentPTI_Cd_C = Cd_day_C*100/18.3,
         percentPTI_Cd_ADO = Cd_day_ADO*100/30,
         percentPTI_Cd_YA = Cd_day_YA_A*100/43.3,
         percentPTI_Cd_A = Cd_day_YA_A*100/58.3) %>% 
  ## For Cu
  mutate(Cu_day_YC = Cu_ww*portion_YoungChildren, # According to EPA-USA, 1 serving = 1 ounce = 28.3g for children age 2
         Cu_day_C = Cu_ww*portion_Children, # According to EPA-USA, 1 serving progressively increases from 1 to 4 ounces from 2 years age to 11 years
         Cu_day_ADO = Cu_ww*portion_Adolescents,
         Cu_day_YA_A = Cu_ww*portion_YoungAdults_Adults,
         percentPTI_Cu_YC = Cu_day_YC*100/1000,
         percentPTI_Cu_C = Cu_day_C*100/3000,
         percentPTI_Cu_ADO = Cu_day_ADO*100/5000,
         percentPTI_Cu_YA = Cu_day_YA_A*100/8000,
         percentPTI_Cu_A = Cu_day_YA_A*100/10000) %>% 
  ## For Fe
  mutate(Fe_day_YC = Fe_ww*portion_YoungChildren, # According to EPA-USA, 1 serving = 1 ounce = 28.3g for children age 2
         Fe_day_C = Fe_ww*portion_Children, # According to EPA-USA, 1 serving progressively increases from 1 to 4 ounces from 2 years age to 11 years
         Fe_day_ADO = Fe_ww*portion_Adolescents,
         Fe_day_YA_A = Fe_ww*portion_YoungAdults_Adults,
         percentPTI_Fe_YC = Fe_day_YC*100/40000,
         percentPTI_Fe_C = Fe_day_C*100/40000,
         percentPTI_Fe_ADO = Fe_day_ADO*100/40000,
         percentPTI_Fe_YA = Fe_day_YA_A*100/45000,
         percentPTI_Fe_A = Fe_day_YA_A*100/45000) %>% 
  ## For Mn
  mutate(Mn_day_YC = Mn_ww*portion_YoungChildren, # According to EPA-USA, 1 serving = 1 ounce = 28.3g for children age 2
         Mn_day_C = Mn_ww*portion_Children, # According to EPA-USA, 1 serving progressively increases from 1 to 4 ounces from 2 years age to 11 years
         Mn_day_ADO = Mn_ww*portion_Adolescents,
         Mn_day_YA_A = Mn_ww*portion_YoungAdults_Adults,
         percentPTI_Mn_YC = Mn_day_YC*100/2000,
         percentPTI_Mn_C = Mn_day_C*100/3000,
         percentPTI_Mn_ADO = Mn_day_ADO*100/6000,
         percentPTI_Mn_YA = Mn_day_YA_A*100/9000,
         percentPTI_Mn_A = Mn_day_YA_A*100/11000) %>% 
  ## For Zn
  mutate(Zn_day_YC = Zn_ww*portion_YoungChildren, # According to EPA-USA, 1 serving = 1 ounce = 28.3g for children age 2
         Zn_day_C = Zn_ww*portion_Children, # According to EPA-USA, 1 serving progressively increases from 1 to 4 ounces from 2 years age to 11 years
         Zn_day_ADO = Zn_ww*portion_Adolescents,
         Zn_day_YA_A = Zn_ww*portion_YoungAdults_Adults,
         percentPTI_Zn_YC = Zn_day_YC*100/7000,
         percentPTI_Zn_C = Zn_day_C*100/12000,
         percentPTI_Zn_ADO = Zn_day_ADO*100/23000,
         percentPTI_Zn_YA = Zn_day_YA_A*100/34000,
         percentPTI_Zn_A = Zn_day_YA_A*100/40000) %>% 
  group_by(english_name) %>% 
  ## For Cd
  summarise(meanPTI_Cd_YC = mean(percentPTI_Cd_YC, na.rm = TRUE),
            sePTI_Cd_YC = std(percentPTI_Cd_YC),
            CIinf_Cd_YC = meanPTI_Cd_YC-(1.96*sePTI_Cd_YC),
            CIsup_Cd_YC = meanPTI_Cd_YC+(1.96*sePTI_Cd_YC),
            meanPTI_Cd_C = mean(percentPTI_Cd_C, na.rm = TRUE),
            sePTI_Cd_C = std(percentPTI_Cd_C),
            CIinf_Cd_C = meanPTI_Cd_C-(1.96*sePTI_Cd_C),
            CIsup_Cd_C = meanPTI_Cd_C+(1.96*sePTI_Cd_C),
            meanPTI_Cd_ADO = mean(percentPTI_Cd_ADO, na.rm = TRUE),
            sePTI_Cd_ADO = std(percentPTI_Cd_ADO),
            CIinf_Cd_ADO = meanPTI_Cd_ADO-(1.96*sePTI_Cd_ADO),
            CIsup_Cd_ADO = meanPTI_Cd_ADO+(1.96*sePTI_Cd_ADO),
            meanPTI_Cd_YA = mean(percentPTI_Cd_YA, na.rm = TRUE),
            sePTI_Cd_YA = std(percentPTI_Cd_YA),
            CIinf_Cd_YA = meanPTI_Cd_YA-(1.96*sePTI_Cd_YA),
            CIsup_Cd_YA = meanPTI_Cd_YA+(1.96*sePTI_Cd_YA),
            meanPTI_Cd_A = mean(percentPTI_Cd_A, na.rm = TRUE),
            sePTI_Cd_A = std(percentPTI_Cd_A),
            CIinf_Cd_A = meanPTI_Cd_A-(1.96*sePTI_Cd_A),
            CIsup_Cd_A = meanPTI_Cd_A+(1.96*sePTI_Cd_A),
            ## For Cu
            meanPTI_Cu_YC = mean(percentPTI_Cu_YC, na.rm = TRUE),
            sePTI_Cu_YC = std(percentPTI_Cu_YC),
            CIinf_Cu_YC = meanPTI_Cu_YC-(1.96*sePTI_Cu_YC),
            CIsup_Cu_YC = meanPTI_Cu_YC+(1.96*sePTI_Cu_YC),
            meanPTI_Cu_C = mean(percentPTI_Cu_C, na.rm = TRUE),
            sePTI_Cu_C = std(percentPTI_Cu_C),
            CIinf_Cu_C = meanPTI_Cu_C-(1.96*sePTI_Cu_C),
            CIsup_Cu_C = meanPTI_Cu_C+(1.96*sePTI_Cu_C),
            meanPTI_Cu_ADO = mean(percentPTI_Cu_ADO, na.rm = TRUE),
            sePTI_Cu_ADO = std(percentPTI_Cu_ADO),
            CIinf_Cu_ADO = meanPTI_Cu_ADO-(1.96*sePTI_Cu_ADO),
            CIsup_Cu_ADO = meanPTI_Cu_ADO+(1.96*sePTI_Cu_ADO),
            meanPTI_Cu_YA = mean(percentPTI_Cu_YA, na.rm = TRUE),
            sePTI_Cu_YA = std(percentPTI_Cu_YA),
            CIinf_Cu_YA = meanPTI_Cu_YA-(1.96*sePTI_Cu_YA),
            CIsup_Cu_YA = meanPTI_Cu_YA+(1.96*sePTI_Cu_YA),
            meanPTI_Cu_A = mean(percentPTI_Cu_A, na.rm = TRUE),
            sePTI_Cu_A = std(percentPTI_Cu_A),
            CIinf_Cu_A = meanPTI_Cu_A-(1.96*sePTI_Cu_A),
            CIsup_Cu_A = meanPTI_Cu_A+(1.96*sePTI_Cu_A),
            ## For Fe
            meanPTI_Fe_YC = mean(percentPTI_Fe_YC, na.rm = TRUE),
            sePTI_Fe_YC = std(percentPTI_Fe_YC),
            CIinf_Fe_YC = meanPTI_Fe_YC-(1.96*sePTI_Fe_YC),
            CIsup_Fe_YC = meanPTI_Fe_YC+(1.96*sePTI_Fe_YC),
            meanPTI_Fe_C = mean(percentPTI_Fe_C, na.rm = TRUE),
            sePTI_Fe_C = std(percentPTI_Fe_C),
            CIinf_Fe_C = meanPTI_Fe_C-(1.96*sePTI_Fe_C),
            CIsup_Fe_C = meanPTI_Fe_C+(1.96*sePTI_Fe_C),
            meanPTI_Fe_ADO = mean(percentPTI_Fe_ADO, na.rm = TRUE),
            sePTI_Fe_ADO = std(percentPTI_Fe_ADO),
            CIinf_Fe_ADO = meanPTI_Fe_ADO-(1.96*sePTI_Fe_ADO),
            CIsup_Fe_ADO = meanPTI_Fe_ADO+(1.96*sePTI_Fe_ADO),
            meanPTI_Fe_YA = mean(percentPTI_Fe_YA, na.rm = TRUE),
            sePTI_Fe_YA = std(percentPTI_Fe_YA),
            CIinf_Fe_YA = meanPTI_Fe_YA-(1.96*sePTI_Fe_YA),
            CIsup_Fe_YA = meanPTI_Fe_YA+(1.96*sePTI_Fe_YA),
            meanPTI_Fe_A = mean(percentPTI_Fe_A, na.rm = TRUE),
            sePTI_Fe_A = std(percentPTI_Fe_A),
            CIinf_Fe_A = meanPTI_Fe_A-(1.96*sePTI_Fe_A),
            CIsup_Fe_A = meanPTI_Fe_A+(1.96*sePTI_Fe_A),
            ## For Mn
            meanPTI_Mn_YC = mean(percentPTI_Mn_YC, na.rm = TRUE),
            sePTI_Mn_YC = std(percentPTI_Mn_YC),
            CIinf_Mn_YC = meanPTI_Mn_YC-(1.96*sePTI_Mn_YC),
            CIsup_Mn_YC = meanPTI_Mn_YC+(1.96*sePTI_Mn_YC),
            meanPTI_Mn_C = mean(percentPTI_Mn_C, na.rm = TRUE),
            sePTI_Mn_C = std(percentPTI_Mn_C),
            CIinf_Mn_C = meanPTI_Mn_C-(1.96*sePTI_Mn_C),
            CIsup_Mn_C = meanPTI_Mn_C+(1.96*sePTI_Mn_C),
            meanPTI_Mn_ADO = mean(percentPTI_Mn_ADO, na.rm = TRUE),
            sePTI_Mn_ADO = std(percentPTI_Mn_ADO),
            CIinf_Mn_ADO = meanPTI_Mn_ADO-(1.96*sePTI_Mn_ADO),
            CIsup_Mn_ADO = meanPTI_Mn_ADO+(1.96*sePTI_Mn_ADO),
            meanPTI_Mn_YA = mean(percentPTI_Mn_YA, na.rm = TRUE),
            sePTI_Mn_YA = std(percentPTI_Mn_YA),
            CIinf_Mn_YA = meanPTI_Mn_YA-(1.96*sePTI_Mn_YA),
            CIsup_Mn_YA = meanPTI_Mn_YA+(1.96*sePTI_Mn_YA),
            meanPTI_Mn_A = mean(percentPTI_Mn_A, na.rm = TRUE),
            sePTI_Mn_A = std(percentPTI_Mn_A),
            CIinf_Mn_A = meanPTI_Mn_A-(1.96*sePTI_Mn_A),
            CIsup_Mn_A = meanPTI_Mn_A+(1.96*sePTI_Mn_A),
            ## For Zn
            meanPTI_Zn_YC = mean(percentPTI_Zn_YC, na.rm = TRUE),
            sePTI_Zn_YC = std(percentPTI_Zn_YC),
            CIinf_Zn_YC = meanPTI_Zn_YC-(1.96*sePTI_Zn_YC),
            CIsup_Zn_YC = meanPTI_Zn_YC+(1.96*sePTI_Zn_YC),
            meanPTI_Zn_C = mean(percentPTI_Zn_C, na.rm = TRUE),
            sePTI_Zn_C = std(percentPTI_Zn_C),
            CIinf_Zn_C = meanPTI_Zn_C-(1.96*sePTI_Zn_C),
            CIsup_Zn_C = meanPTI_Zn_C+(1.96*sePTI_Zn_C),
            meanPTI_Zn_ADO = mean(percentPTI_Zn_ADO, na.rm = TRUE),
            sePTI_Zn_ADO = std(percentPTI_Zn_ADO),
            CIinf_Zn_ADO = meanPTI_Zn_ADO-(1.96*sePTI_Zn_ADO),
            CIsup_Zn_ADO = meanPTI_Zn_ADO+(1.96*sePTI_Zn_ADO),
            meanPTI_Zn_YA = mean(percentPTI_Zn_YA, na.rm = TRUE),
            sePTI_Zn_YA = std(percentPTI_Zn_YA),
            CIinf_Zn_YA = meanPTI_Zn_YA-(1.96*sePTI_Zn_YA),
            CIsup_Zn_YA = meanPTI_Zn_YA+(1.96*sePTI_Zn_YA),
            meanPTI_Zn_A = mean(percentPTI_Zn_A, na.rm = TRUE),
            sePTI_Zn_A = std(percentPTI_Zn_A),
            CIinf_Zn_A = meanPTI_Zn_A-(1.96*sePTI_Zn_A),
            CIsup_Zn_A = meanPTI_Zn_A+(1.96*sePTI_Zn_A))

## 3 / Cleaning database ##############################################################################################

percentPTI_propre <- percentPTI %>% 
  left_join(percentPTI_SeHg, by = "english_name")
rm(percentPTI,percentPTI_SeHg)

names(percentPTI_propre)

ncol <- ncol(percentPTI_propre)
names(percentPTI_propre)

# Boucle pour Cd
percentPTI_propre$servCd_inf_Children <- as.numeric(NA)
percentPTI_propre$servCd_sup_Children <- as.numeric(NA)
for (i in 1:nrow(percentPTI_propre)) {
  
  ## 1. Children
  # 1.1. Min-max using min(CIinf) and max(CIsup)
  percentPTI_propre[i,ncol+1] <- min(percentPTI_propre[i,4],
                                     percentPTI_propre[i,8],
                                     percentPTI_propre[i,12])
  percentPTI_propre[i,ncol+2] <- max(percentPTI_propre[i,5],
                                     percentPTI_propre[i,9],
                                     percentPTI_propre[i,13])
  # 1.2. Dealing with unique measures and NAs
  percentPTI_propre[i,ncol+1] <- ifelse(is.na(percentPTI_propre[i,4]) == TRUE,
                                        min(percentPTI_propre[i,2], percentPTI_propre[i,6], percentPTI_propre[i,10]),
                                        percentPTI_propre[i,ncol+1])
  percentPTI_propre[i,ncol+2] <- ifelse(is.na(percentPTI_propre[i,5]) == TRUE,
                                        max(percentPTI_propre[i,2], percentPTI_propre[i,6], percentPTI_propre[i,10]),
                                        percentPTI_propre[i,ncol+2])
}

# Boucle pour Cu
percentPTI_propre$servCu_inf_Children <- as.numeric(NA)
percentPTI_propre$servCu_sup_Children <- as.numeric(NA)
for (i in 1:nrow(percentPTI_propre)) {
  
  ## 1. Children
  # 1.1. Min-max using min(CIinf) and max(CIsup)
  percentPTI_propre[i,ncol+3] <- min(percentPTI_propre[i,24],
                                     percentPTI_propre[i,28],
                                     percentPTI_propre[i,32])
  percentPTI_propre[i,ncol+4] <- max(percentPTI_propre[i,25],
                                     percentPTI_propre[i,29],
                                     percentPTI_propre[i,33])
  # 1.2. Dealing with unique measures and NAs
  percentPTI_propre[i,ncol+3] <- ifelse(is.na(percentPTI_propre[i,24]) == TRUE,
                                        min(percentPTI_propre[i,22], percentPTI_propre[i,26], percentPTI_propre[i,30]),
                                        percentPTI_propre[i,ncol+3])
  percentPTI_propre[i,ncol+4] <- ifelse(is.na(percentPTI_propre[i,25]) == TRUE,
                                        max(percentPTI_propre[i,22], percentPTI_propre[i,26], percentPTI_propre[i,30]),
                                        percentPTI_propre[i,ncol+4])
}

# Boucle pour Fe
percentPTI_propre$servFe_inf_Children <- as.numeric(NA)
percentPTI_propre$servFe_sup_Children <- as.numeric(NA)
for (i in 1:nrow(percentPTI_propre)) {
  
  ## 1. Children
  # 1.1. Min-max using min(CIinf) and max(CIsup)
  percentPTI_propre[i,ncol+5] <- min(percentPTI_propre[i,44],
                                     percentPTI_propre[i,48],
                                     percentPTI_propre[i,52])
  percentPTI_propre[i,ncol+6] <- max(percentPTI_propre[i,45],
                                     percentPTI_propre[i,49],
                                     percentPTI_propre[i,53])
  # 1.2. Dealing with unique measures and NAs
  percentPTI_propre[i,ncol+5] <- ifelse(is.na(percentPTI_propre[i,44]) == TRUE,
                                        min(percentPTI_propre[i,42], percentPTI_propre[i,46], percentPTI_propre[i,50]),
                                        percentPTI_propre[i,ncol+5])
  percentPTI_propre[i,ncol+6] <- ifelse(is.na(percentPTI_propre[i,45]) == TRUE,
                                        max(percentPTI_propre[i,42], percentPTI_propre[i,46], percentPTI_propre[i,50]),
                                        percentPTI_propre[i,ncol+6])
}

# Boucle pour Mn
percentPTI_propre$servMn_inf_Children <- as.numeric(NA)
percentPTI_propre$servMn_sup_Children <- as.numeric(NA)
for (i in 1:nrow(percentPTI_propre)) {
  
  ## 1. Children
  # 1.1. Min-max using min(CIinf) and max(CIsup)
  percentPTI_propre[i,ncol+7] <- min(percentPTI_propre[i,64],
                                     percentPTI_propre[i,68],
                                     percentPTI_propre[i,72])
  percentPTI_propre[i,ncol+8] <- max(percentPTI_propre[i,65],
                                     percentPTI_propre[i,69],
                                     percentPTI_propre[i,73])
  # 1.2. Dealing with unique measures and NAs
  percentPTI_propre[i,ncol+7] <- ifelse(is.na(percentPTI_propre[i,64]) == TRUE,
                                        min(percentPTI_propre[i,62], percentPTI_propre[i,66], percentPTI_propre[i,70]),
                                        percentPTI_propre[i,ncol+7])
  percentPTI_propre[i,ncol+8] <- ifelse(is.na(percentPTI_propre[i,64]) == TRUE,
                                        max(percentPTI_propre[i,62], percentPTI_propre[i,66], percentPTI_propre[i,70]),
                                        percentPTI_propre[i,ncol+8])
}

# Boucle pour Se
percentPTI_propre$servSe_inf_Children <- as.numeric(NA)
percentPTI_propre$servSe_sup_Children <- as.numeric(NA)
for (i in 1:nrow(percentPTI_propre)) {
  
  ## 1. Children
  # 1.1. Min-max using min(CIinf) and max(CIsup)
  percentPTI_propre[i,ncol+9] <- min(percentPTI_propre[i,104],
                                     percentPTI_propre[i,108],
                                     percentPTI_propre[i,112])
  percentPTI_propre[i,ncol+10] <- max(percentPTI_propre[i,105],
                                      percentPTI_propre[i,109],
                                      percentPTI_propre[i,113])
  # 1.2. Dealing with unique measures and NAs
  percentPTI_propre[i,ncol+9] <- ifelse(is.na(percentPTI_propre[i,104]) == TRUE,
                                        min(percentPTI_propre[i,102], percentPTI_propre[i,106], percentPTI_propre[i,110]),
                                        percentPTI_propre[i,ncol+9])
  percentPTI_propre[i,ncol+10] <- ifelse(is.na(percentPTI_propre[i,104]) == TRUE,
                                         max(percentPTI_propre[i,102], percentPTI_propre[i,106], percentPTI_propre[i,110]),
                                         percentPTI_propre[i,ncol+10])
}

# Boucle pour Zn
percentPTI_propre$servZn_inf_Children <- as.numeric(NA)
percentPTI_propre$servZn_sup_Children <- as.numeric(NA)
for (i in 1:nrow(percentPTI_propre)) {
  
  ## 1. Children
  # 1.1. Min-max using min(CIinf) and max(CIsup)
  percentPTI_propre[i,ncol+11] <- min(percentPTI_propre[i,84],
                                      percentPTI_propre[i,88],
                                      percentPTI_propre[i,92])
  percentPTI_propre[i,ncol+12] <- max(percentPTI_propre[i,85],
                                      percentPTI_propre[i,89],
                                      percentPTI_propre[i,93])
  # 1.2. Dealing with unique measures and NAs
  percentPTI_propre[i,ncol+11] <- ifelse(is.na(percentPTI_propre[i,84]) == TRUE,
                                         min(percentPTI_propre[i,82], percentPTI_propre[i,86], percentPTI_propre[i,90]),
                                         percentPTI_propre[i,ncol+11])
  percentPTI_propre[i,ncol+12] <- ifelse(is.na(percentPTI_propre[i,85]) == TRUE,
                                         max(percentPTI_propre[i,82], percentPTI_propre[i,86], percentPTI_propre[i,90]),
                                         percentPTI_propre[i,ncol+12])
}

# Boucle pour Hg
percentPTI_propre$servHg_inf_Children <- as.numeric(NA)
percentPTI_propre$servHg_sup_Children <- as.numeric(NA)
for (i in 1:nrow(percentPTI_propre)) {
  
  ## 1. Children
  # 1.1. Min-max using min(CIinf) and max(CIsup)
  percentPTI_propre[i,ncol+13] <- min(percentPTI_propre[i,120],
                                      percentPTI_propre[i,124],
                                      percentPTI_propre[i,128])
  percentPTI_propre[i,ncol+14] <- max(percentPTI_propre[i,121],
                                      percentPTI_propre[i,125],
                                      percentPTI_propre[i,129])
  # 1.2. Dealing with unique measures and NAs
  percentPTI_propre[i,ncol+13] <- ifelse(is.na(percentPTI_propre[i,120]) == TRUE,
                                         min(percentPTI_propre[i,118], percentPTI_propre[i,122], percentPTI_propre[i,126]),
                                         percentPTI_propre[i,ncol+13])
  percentPTI_propre[i,ncol+14] <- ifelse(is.na(percentPTI_propre[i,121]) == TRUE,
                                         max(percentPTI_propre[i,118], percentPTI_propre[i,122], percentPTI_propre[i,126]),
                                         percentPTI_propre[i,ncol+14])
}

rm(i,ncol)

names(percentPTI_propre)
percentPTI_clean <- percentPTI_propre %>% 
  mutate(meanPTI_Cd_Children = (meanPTI_Cd_YC+meanPTI_Cd_C+meanPTI_Cd_ADO)/3,
         meanPTI_Cu_Children = (meanPTI_Cu_YC+meanPTI_Cu_C+meanPTI_Cu_ADO)/3,
         meanPTI_Fe_Children = (meanPTI_Fe_YC+meanPTI_Fe_C+meanPTI_Fe_ADO)/3,
         meanPTI_Mn_Children = (meanPTI_Mn_YC+meanPTI_Mn_C+meanPTI_Mn_ADO)/3,
         meanPTI_Se_Children = (meanPTI_Se_YC+meanPTI_Se_C+meanPTI_Se_ADO)/3,
         meanPTI_Zn_Children = (meanPTI_Zn_YC+meanPTI_Zn_C+meanPTI_Zn_ADO)/3,
         meanPTI_Hg_Children = (meanPTI_Hg_YC+meanPTI_Hg_C+meanPTI_Hg_ADO)/3) %>% 
  mutate(## For Cu
    Cu_Children = paste0(round(meanPTI_Cu_Children, 1)," [",round(servCu_inf_Children, 1),"-",round(servCu_sup_Children, 1),"]"),
    Cu_Children = ifelse(is.na(servCu_inf_Children), as.character(round(meanPTI_Cu_Children,1)), Cu_Children),
    Cu_YoungAdults = paste0(round(meanPTI_Cu_YA,1)," [",round(CIinf_Cu_YA,1),"-",round(CIsup_Cu_YA,1),"]"),
    Cu_YoungAdults = ifelse(is.na(CIinf_Cu_YA), as.character(round(meanPTI_Cu_YA,1)), Cu_YoungAdults),
    Cu_Adults = paste0(round(meanPTI_Cu_A,1)," [",round(CIinf_Cu_A,1),"-",round(CIsup_Cu_A,1),"]"),
    Cu_Adults = ifelse(is.na(CIinf_Cu_A), as.character(round(meanPTI_Cu_A,1)), Cu_Adults),
    ## For Fe
    Fe_Children = paste0(round(meanPTI_Fe_Children,2)," [",round(servFe_inf_Children,2),"-",round(servFe_sup_Children,2),"]"),
    Fe_Children = ifelse(is.na(servFe_inf_Children), as.character(round(meanPTI_Fe_Children,2)), Fe_Children),
    Fe_YoungAdults = paste0(round(meanPTI_Fe_YA,2)," [",round(CIinf_Fe_YA,2),"-",round(CIsup_Fe_YA,2),"]"),
    Fe_YoungAdults = ifelse(is.na(CIinf_Fe_YA), as.character(round(meanPTI_Fe_YA,2)), Fe_YoungAdults),
    Fe_Adults = paste0(round(meanPTI_Fe_A,2)," [",round(CIinf_Fe_A,2),"-",round(CIsup_Fe_A,2),"]"),
    Fe_Adults = ifelse(is.na(CIinf_Fe_A), as.character(round(meanPTI_Fe_A,2)), Fe_Adults),
    ## For Mn
    Mn_Children = paste0(round(meanPTI_Mn_Children,2)," [",round(servMn_inf_Children,2),"-",round(servMn_sup_Children,2),"]"),
    Mn_Children = ifelse(is.na(servMn_inf_Children), as.character(round(meanPTI_Mn_Children,2)), Mn_Children),
    Mn_YoungAdults = paste0(round(meanPTI_Mn_YA,2)," [",round(CIinf_Mn_YA,2),"-",round(CIsup_Mn_YA,2),"]"),
    Mn_YoungAdults = ifelse(is.na(CIinf_Mn_YA), as.character(round(meanPTI_Mn_YA,2)), Mn_YoungAdults),
    Mn_Adults = paste0(round(meanPTI_Mn_A,2)," [",round(CIinf_Mn_A,2),"-",round(CIsup_Mn_A,2),"]"),
    Mn_Adults = ifelse(is.na(CIinf_Mn_A), as.character(round(meanPTI_Mn_A,2)), Mn_Adults),
    ## For Se
    Se_Children = paste0(round(meanPTI_Se_Children)," [",round(servSe_inf_Children),"-",round(servSe_sup_Children),"]"),
    Se_Children = ifelse(is.na(servSe_inf_Children), as.character(round(meanPTI_Se_Children)), Se_Children),
    Se_YoungAdults_Adults = paste0(round(meanPTI_Se_YA_A)," [",round(CIinf_Se_YA_A),"-",round(CIsup_Se_YA_A),"]"),
    Se_YoungAdults_Adults = ifelse(is.na(CIinf_Se_YA_A), as.character(round(meanPTI_Se_YA_A)), Se_YoungAdults_Adults),
    ## For Zn
    Zn_Children = paste0(round(meanPTI_Zn_Children)," [",round(servZn_inf_Children),"-",round(servZn_sup_Children),"]"),
    Zn_Children = ifelse(is.na(servZn_inf_Children), as.character(round(meanPTI_Zn_Children)), Zn_Children),
    Zn_YoungAdults = paste0(round(meanPTI_Zn_YA)," [",round(CIinf_Zn_YA),"-",round(CIsup_Zn_YA),"]"),
    Zn_YoungAdults = ifelse(is.na(CIinf_Zn_YA), as.character(round(meanPTI_Zn_YA)), Zn_YoungAdults),
    Zn_Adults = paste0(round(meanPTI_Zn_A)," [",round(CIinf_Zn_A),"-",round(CIsup_Zn_A),"]"),
    Zn_Adults = ifelse(is.na(CIinf_Zn_A), as.character(round(meanPTI_Zn_A)), Zn_Adults),
    ## For Cd
    Cd_Children = paste0(round(meanPTI_Cd_Children)," [",round(servCd_inf_Children),"-",round(servCd_sup_Children),"]"),
    Cd_Children = ifelse(is.na(servCd_inf_Children), as.character(round(meanPTI_Cd_Children)), Cd_Children),
    Cd_YoungAdults = paste0(round(meanPTI_Cd_YA)," [",round(CIinf_Cd_YA),"-",round(CIsup_Cd_YA),"]"),
    Cd_YoungAdults = ifelse(is.na(CIinf_Cd_YA), as.character(round(meanPTI_Cd_YA)), Cd_YoungAdults),
    Cd_Adults = paste0(round(meanPTI_Cd_A)," [",round(CIinf_Cd_A),"-",round(CIsup_Cd_A),"]"),
    Cd_Adults = ifelse(is.na(CIinf_Cd_A), as.character(round(meanPTI_Cd_A)), Cd_Adults),
    ## For Hg
    Hg_Children = paste0(round(meanPTI_Hg_Children,1)," [",round(servHg_inf_Children,1),"-",round(servHg_sup_Children,1),"]"),
    Hg_Children = ifelse(is.na(servHg_inf_Children), as.character(round(meanPTI_Hg_Children,1)), Hg_Children),
    Hg_YoungAdults = paste0(round(meanPTI_Hg_YA,1)," [",round(CIinf_Hg_YA,1),"-",round(CIsup_Hg_YA,1),"]"),
    Hg_YoungAdults = ifelse(is.na(CIinf_Hg_YA), as.character(round(meanPTI_Hg_YA,1)), Hg_YoungAdults),
    Hg_Adults = paste0(round(meanPTI_Hg_A,1)," [",round(CIinf_Hg_A,1),"-",round(CIsup_Hg_A,1),"]"),
    Hg_Adults = ifelse(is.na(CIinf_Hg_A), as.character(round(meanPTI_Hg_A,1)), Hg_Adults)) %>%
  left_join(infos, by = "english_name") %>%
  mutate(number = as.numeric(english_name)) %>%
  arrange(number) %>%
  select(sp_group,english_name,Cu_Children,Cu_YoungAdults,Cu_Adults,
         Fe_Children,Fe_YoungAdults,Fe_Adults,
         Mn_Children,Mn_YoungAdults,Mn_Adults,
         Se_Children,Se_YoungAdults_Adults,
         Zn_Children,Zn_YoungAdults,Zn_Adults,
         Cd_Children,Cd_YoungAdults,Cd_Adults,
         Hg_Children,Hg_YoungAdults,Hg_Adults)

write.table(percentPTI_clean,paste0(direction,"3_Risk_usingPercentPTI.csv"),
            sep = ";", dec = ".", row.names = F)

rm(percentPTI_clean)

## 4 / Grade each species according to %PTI (using CI95% sup) ####

Risk_step1 <- percentPTI_propre %>% 
  mutate(meanPTI_Cd_Children = (meanPTI_Cd_YC+meanPTI_Cd_C+meanPTI_Cd_ADO)/3,
         meanPTI_Cu_Children = (meanPTI_Cu_YC+meanPTI_Cu_C+meanPTI_Cu_ADO)/3,
         meanPTI_Fe_Children = (meanPTI_Fe_YC+meanPTI_Fe_C+meanPTI_Fe_ADO)/3,
         meanPTI_Mn_Children = (meanPTI_Mn_YC+meanPTI_Mn_C+meanPTI_Mn_ADO)/3,
         meanPTI_Se_Children = (meanPTI_Se_YC+meanPTI_Se_C+meanPTI_Se_ADO)/3,
         meanPTI_Zn_Children = (meanPTI_Zn_YC+meanPTI_Zn_C+meanPTI_Zn_ADO)/3,
         meanPTI_Hg_Children = (meanPTI_Hg_YC+meanPTI_Hg_C+meanPTI_Hg_ADO)/3) %>% 
  ## For Cu
  mutate(servCu_sup_Children = ifelse(is.na(servCu_sup_Children) == TRUE, meanPTI_Cu_Children, servCu_sup_Children),
         R_Cu_Children = "0",
         R_Cu_Children = ifelse(servCu_sup_Children > 75, "-1", R_Cu_Children),
         R_Cu_Children = ifelse(servCu_sup_Children > 90, "-2", R_Cu_Children),
         CIsup_Cu_YA = ifelse(is.na(CIsup_Cu_YA) == TRUE, meanPTI_Cu_YA, CIsup_Cu_YA),
         R_Cu_YoungAdults = "0",
         R_Cu_YoungAdults = ifelse(CIsup_Cu_YA > 75, "-1", R_Cu_YoungAdults),
         R_Cu_YoungAdults = ifelse(CIsup_Cu_YA > 90, "-2", R_Cu_YoungAdults),
         CIsup_Cu_A = ifelse(is.na(CIsup_Cu_A) == TRUE, meanPTI_Cu_A, CIsup_Cu_A),
         R_Cu_Adults = "0",
         R_Cu_Adults = ifelse(CIsup_Cu_A > 75, "-1", R_Cu_Adults),
         R_Cu_Adults = ifelse(CIsup_Cu_A > 90, "-2", R_Cu_Adults)) %>% 
  ## For Fe
  mutate(servFe_sup_Children = ifelse(is.na(servFe_sup_Children) == TRUE, meanPTI_Fe_Children, servFe_sup_Children),
         R_Fe_Children = "0",
         R_Fe_Children = ifelse(servFe_sup_Children > 75, "-1", R_Fe_Children),
         R_Fe_Children = ifelse(servFe_sup_Children > 90, "-2", R_Fe_Children),
         CIsup_Fe_YA = ifelse(is.na(CIsup_Fe_YA) == TRUE, meanPTI_Fe_YA, CIsup_Fe_YA),
         R_Fe_YoungAdults = "0",
         R_Fe_YoungAdults = ifelse(CIsup_Fe_YA > 75, "-1", R_Fe_YoungAdults),
         R_Fe_YoungAdults = ifelse(CIsup_Fe_YA > 90, "-2", R_Fe_YoungAdults),
         CIsup_Fe_A = ifelse(is.na(CIsup_Fe_A) == TRUE, meanPTI_Fe_A, CIsup_Fe_A),
         R_Fe_Adults = "0",
         R_Fe_Adults = ifelse(CIsup_Fe_A > 75, "-1", R_Fe_Adults),
         R_Fe_Adults = ifelse(CIsup_Fe_A > 90, "-2", R_Fe_Adults)) %>% 
  ## For Mn
  mutate(servMn_sup_Children = ifelse(is.na(servMn_sup_Children) == TRUE, meanPTI_Mn_Children, servMn_sup_Children),
         R_Mn_Children = "0",
         R_Mn_Children = ifelse(servMn_sup_Children > 75, "-1", R_Mn_Children),
         R_Mn_Children = ifelse(servMn_sup_Children > 90, "-2", R_Mn_Children),
         CIsup_Mn_YA = ifelse(is.na(CIsup_Mn_YA) == TRUE, meanPTI_Mn_YA, CIsup_Mn_YA),
         R_Mn_YoungAdults = "0",
         R_Mn_YoungAdults = ifelse(CIsup_Mn_YA > 75, "-1", R_Mn_YoungAdults),
         R_Mn_YoungAdults = ifelse(CIsup_Mn_YA > 90, "-2", R_Mn_YoungAdults),
         CIsup_Mn_A = ifelse(is.na(CIsup_Mn_A) == TRUE, meanPTI_Mn_A, CIsup_Mn_A),
         R_Mn_Adults = "0",
         R_Mn_Adults = ifelse(CIsup_Mn_A > 75, "-1", R_Mn_Adults),
         R_Mn_Adults = ifelse(CIsup_Mn_A > 90, "-2", R_Mn_Adults)) %>% 
  ## For Se
  mutate(servSe_sup_Children = ifelse(is.na(servSe_sup_Children) == TRUE, meanPTI_Se_Children, servSe_sup_Children),
         R_Se_Children = "0",
         R_Se_Children = ifelse(servSe_sup_Children > 75, "-1", R_Se_Children),
         R_Se_Children = ifelse(servSe_sup_Children > 90, "-2", R_Se_Children),
         CIsup_Se_YA_A = ifelse(is.na(CIsup_Se_YA_A) == TRUE, meanPTI_Se_YA_A, CIsup_Se_YA_A),
         R_Se_YoungAdults_Adults = "0",
         R_Se_YoungAdults_Adults = ifelse(CIsup_Se_YA_A > 75, "-1", R_Se_YoungAdults_Adults),
         R_Se_YoungAdults_Adults = ifelse(CIsup_Se_YA_A > 90, "-2", R_Se_YoungAdults_Adults)) %>% 
  ## For Zn
  mutate(servZn_sup_Children = ifelse(is.na(servZn_sup_Children) == TRUE, meanPTI_Zn_Children, servZn_sup_Children),
         R_Zn_Children = "0",
         R_Zn_Children = ifelse(servZn_sup_Children > 75, "-1", R_Zn_Children),
         R_Zn_Children = ifelse(servZn_sup_Children > 90, "-2", R_Zn_Children),
         CIsup_Zn_YA = ifelse(is.na(CIsup_Zn_YA) == TRUE, meanPTI_Zn_YA, CIsup_Zn_YA),
         R_Zn_YoungAdults = "0",
         R_Zn_YoungAdults = ifelse(CIsup_Zn_YA > 75, "-1", R_Zn_YoungAdults),
         R_Zn_YoungAdults = ifelse(CIsup_Zn_YA > 90, "-2", R_Zn_YoungAdults),
         CIsup_Zn_A = ifelse(is.na(CIsup_Zn_A) == TRUE, meanPTI_Zn_A, CIsup_Zn_A),
         R_Zn_Adults = "0",
         R_Zn_Adults = ifelse(CIsup_Zn_A > 75, "-1", R_Zn_Adults),
         R_Zn_Adults = ifelse(CIsup_Zn_A > 90, "-2", R_Zn_Adults)) %>% 
  ## For Cd
  mutate(servCd_sup_Children = ifelse(is.na(servCd_sup_Children) == TRUE, meanPTI_Cd_Children, servCd_sup_Children),
         R_Cd_Children = "0",
         R_Cd_Children = ifelse(servCd_sup_Children > 75, "-1", R_Cd_Children),
         R_Cd_Children = ifelse(servCd_sup_Children > 90, "-2", R_Cd_Children),
         CIsup_Cd_YA = ifelse(is.na(CIsup_Cd_YA) == TRUE, meanPTI_Cd_YA, CIsup_Cd_YA),
         R_Cd_YoungAdults = "0",
         R_Cd_YoungAdults = ifelse(CIsup_Cd_YA > 75, "-1", R_Cd_YoungAdults),
         R_Cd_YoungAdults = ifelse(CIsup_Cd_YA > 90, "-2", R_Cd_YoungAdults),
         CIsup_Cd_A = ifelse(is.na(CIsup_Cd_A) == TRUE, meanPTI_Cd_A, CIsup_Cd_A),
         R_Cd_Adults = "0",
         R_Cd_Adults = ifelse(CIsup_Cd_A > 75, "-1", R_Cd_Adults),
         R_Cd_Adults = ifelse(CIsup_Cd_A > 90, "-2", R_Cd_Adults)) %>%
  ## For Hg
  mutate(servHg_sup_Children = ifelse(is.na(servHg_sup_Children) == TRUE, meanPTI_Hg_Children, servHg_sup_Children),
         R_Hg_Children = "0",
         R_Hg_Children = ifelse(servHg_sup_Children > 75, "-1", R_Hg_Children),
         R_Hg_Children = ifelse(servHg_sup_Children > 90, "-2", R_Hg_Children),
         CIsup_Hg_YA = ifelse(is.na(CIsup_Hg_YA) == TRUE, meanPTI_Hg_YA, CIsup_Hg_YA),
         R_Hg_YoungAdults = "0",
         R_Hg_YoungAdults = ifelse(CIsup_Hg_YA > 75, "-1", R_Hg_YoungAdults),
         R_Hg_YoungAdults = ifelse(CIsup_Hg_YA > 90, "-2", R_Hg_YoungAdults),
         CIsup_Hg_A = ifelse(is.na(CIsup_Hg_A) == TRUE, meanPTI_Hg_A, CIsup_Hg_A),
         R_Hg_Adults = "0",
         R_Hg_Adults = ifelse(CIsup_Hg_A > 75, "-1", R_Hg_Adults),
         R_Hg_Adults = ifelse(CIsup_Hg_A > 90, "-2", R_Hg_Adults)) %>% 
  left_join(infos, by = "english_name") %>%
  mutate(number = as.numeric(english_name)) %>%
  arrange(number) %>%
  select(sp_group,english_name, R_Cu_Children, R_Cu_YoungAdults, R_Cu_Adults,
         R_Fe_Children, R_Fe_YoungAdults, R_Fe_Adults,
         R_Mn_Children, R_Mn_YoungAdults, R_Mn_Adults,
         R_Se_Children, R_Se_YoungAdults_Adults,
         R_Zn_Children, R_Zn_YoungAdults, R_Zn_Adults,
         R_Cd_Children, R_Cd_YoungAdults, R_Cd_Adults,
         R_Hg_Children, R_Hg_YoungAdults, R_Hg_Adults)
Risk_step1[is.na(Risk_step1)] = "0"

write.table(Risk_step1,paste0(direction,"4_percenPTI_ratingSp.csv"),
            sep = ";", dec = ".", row.names = F)

rm(Risk_step1)

rm(portion_Adolescents,portion_Children,portion_YoungAdults_Adults,
   portion_YoungChildren)
rm(percentPTI_propre,infos)



### V // Benefit ##############################################################################################

infos <- data_risk_benefit %>% select(sp_group,english_name) %>% distinct(english_name, .keep_all = TRUE)

## 1 / Calculation of %RDI for all essential and potentially essential TE ##############################################################################################

percentRDI <- data_risk_benefit %>% 
  left_join(interHgSe, by = "organism_identifier") %>% 
  rename(english_name = english_name.x) %>% 
  select(english_name, TAs_ww, Cu_ww, Fe_ww, Mn_ww, Se_avail, Zn_ww) %>% 
  ## For As
  mutate(As_day_YC = TAs_ww*30, # According to EPA-USA, 1 serving = 1 ounce = 28.3g for children age 2
         As_day_C = TAs_ww*70, # According to EPA-USA, 1 serving progressively increases from 1 to 4 ounces from 2 years age to 11 years
         As_day_ADO = TAs_ww*110,
         As_day_YA_A = TAs_ww*156,
         percentRDI_As_YC = As_day_YC*100/12,
         percentRDI_As_C = As_day_C*100/12,
         percentRDI_As_ADO = As_day_ADO*100/12,
         percentRDI_As_YA = As_day_YA_A*100/12,
         percentRDI_As_A = As_day_YA_A*100/12) %>% 
  ## For Cu
  mutate(Cu_day_YC = Cu_ww*30, # According to EPA-USA, 1 serving = 1 ounce = 28.3g for children age 2
         Cu_day_C = Cu_ww*70, # According to EPA-USA, 1 serving progressively increases from 1 to 4 ounces from 2 years age to 11 years
         Cu_day_ADO = Cu_ww*110,
         Cu_day_YA_A = Cu_ww*156,
         percentRDI_Cu_YC = Cu_day_YC*100/340,
         percentRDI_Cu_C = Cu_day_C*100/440,
         percentRDI_Cu_ADO = Cu_day_ADO*100/700,
         percentRDI_Cu_YA = Cu_day_YA_A*100/890,
         percentRDI_Cu_A = Cu_day_YA_A*100/900,
         percentRDI_Cu_Apre = Cu_day_YA_A*100/1000,
         percentRDI_Cu_Alact = Cu_day_YA_A*100/1300) %>% 
  ## For Fe
  mutate(Fe_day_YC = Fe_ww*30, # According to EPA-USA, 1 serving = 1 ounce = 28.3g for children age 2
         Fe_day_C = Fe_ww*70, # According to EPA-USA, 1 serving progressively increases from 1 to 4 ounces from 2 years age to 11 years
         Fe_day_ADO = Fe_ww*110,
         Fe_day_YA_A = Fe_ww*156,
         percentRDI_Fe_YC = Fe_day_YC*100/7000,
         percentRDI_Fe_C = Fe_day_C*100/10000,
         percentRDI_Fe_ADO = Fe_day_ADO*100/8000,
         percentRDI_Fe_YA = Fe_day_YA_A*100/15000,
         percentRDI_Fe_A = Fe_day_YA_A*100/18000,
         percentRDI_Fe_Apre = Fe_day_YA_A*100/27000,
         percentRDI_Fe_Alact = Fe_day_YA_A*100/9000) %>% 
  ## For Mn
  mutate(Mn_day_YC = Mn_ww*30, # According to EPA-USA, 1 serving = 1 ounce = 28.3g for children age 2
         Mn_day_C = Mn_ww*70, # According to EPA-USA, 1 serving progressively increases from 1 to 4 ounces from 2 years age to 11 years
         Mn_day_ADO = Mn_ww*110,
         Mn_day_YA_A = Mn_ww*156,
         percentRDI_Mn_YC = Mn_day_YC*100/1200,
         percentRDI_Mn_C = Mn_day_C*100/1500,
         percentRDI_Mn_ADO = Mn_day_ADO*100/1900,
         percentRDI_Mn_YA = Mn_day_YA_A*100/2200,
         percentRDI_Mn_A = Mn_day_YA_A*100/1800,
         percentRDI_Mn_Apre = Mn_day_YA_A*100/2000,
         percentRDI_Mn_Alact = Mn_day_YA_A*100/2600) %>% 
  ## For Se
  mutate(Se_day_YC = Se_avail*30, # According to EPA-USA, 1 serving = 1 ounce = 28.3g for children age 2
         Se_day_C = Se_avail*70, # According to EPA-USA, 1 serving progressively increases from 1 to 4 ounces from 2 years age to 11 years
         Se_day_ADO = Se_avail*110,
         Se_day_YA_A = Se_avail*156,
         percentRDI_Se_YC = Se_day_YC*100/20,
         percentRDI_Se_C = Se_day_C*100/30,
         percentRDI_Se_ADO = Se_day_ADO*100/40,
         percentRDI_Se_YA_A = Se_day_YA_A*100/55,
         percentRDI_Se_Apre = Se_day_YA_A*100/60,
         percentRDI_Se_Alact = Se_day_YA_A*100/70) %>% 
  ## For Zn
  mutate(Zn_day_YC = Zn_ww*30, # According to EPA-USA, 1 serving = 1 ounce = 28.3g for children age 2
         Zn_day_C = Zn_ww*70, # According to EPA-USA, 1 serving progressively increases from 1 to 4 ounces from 2 years age to 11 years
         Zn_day_ADO = Zn_ww*110,
         Zn_day_YA_A = Zn_ww*156,
         percentRDI_Zn_YC = Zn_day_YC*100/3000,
         percentRDI_Zn_C = Zn_day_C*100/5000,
         percentRDI_Zn_ADO_A = Zn_day_ADO*100/8000,
         percentRDI_Zn_YA_Apre = Zn_day_YA_A*100/11000,
         percentRDI_Zn_Alact = Zn_day_YA_A*100/12000) %>% 
  group_by(english_name) %>% 
  ## Calcul mean for each species
  summarise(## For As
    meanRDI_As_YC = mean(percentRDI_As_YC, na.rm = TRUE),
    seRDI_As_YC = std(percentRDI_As_YC),
    CIinf_As_YC = meanRDI_As_YC-(1.96*seRDI_As_YC),
    CIsup_As_YC = meanRDI_As_YC+(1.96*seRDI_As_YC),
    meanRDI_As_C = mean(percentRDI_As_C, na.rm = TRUE),
    seRDI_As_C = std(percentRDI_As_C),
    CIinf_As_C = meanRDI_As_C-(1.96*seRDI_As_C),
    CIsup_As_C = meanRDI_As_C+(1.96*seRDI_As_C),
    meanRDI_As_ADO = mean(percentRDI_As_ADO, na.rm = TRUE),
    seRDI_As_ADO = std(percentRDI_As_ADO),
    CIinf_As_ADO = meanRDI_As_ADO-(1.96*seRDI_As_ADO),
    CIsup_As_ADO = meanRDI_As_ADO+(1.96*seRDI_As_ADO),
    meanRDI_As_YA = mean(percentRDI_As_YA, na.rm = TRUE),
    seRDI_As_YA = std(percentRDI_As_YA),
    CIinf_As_YA = meanRDI_As_YA-(1.96*seRDI_As_YA),
    CIsup_As_YA = meanRDI_As_YA+(1.96*seRDI_As_YA),
    meanRDI_As_A = mean(percentRDI_As_A, na.rm = TRUE),
    seRDI_As_A = std(percentRDI_As_A),
    CIinf_As_A = meanRDI_As_A-(1.96*seRDI_As_A),
    CIsup_As_A = meanRDI_As_A+(1.96*seRDI_As_A),
    ## For Cu
    meanRDI_Cu_YC = mean(percentRDI_Cu_YC, na.rm = TRUE),
    seRDI_Cu_YC = std(percentRDI_Cu_YC),
    CIinf_Cu_YC = meanRDI_Cu_YC-(1.96*seRDI_Cu_YC),
    CIsup_Cu_YC = meanRDI_Cu_YC+(1.96*seRDI_Cu_YC),
    meanRDI_Cu_C = mean(percentRDI_Cu_C, na.rm = TRUE),
    seRDI_Cu_C = std(percentRDI_Cu_C),
    CIinf_Cu_C = meanRDI_Cu_C-(1.96*seRDI_Cu_C),
    CIsup_Cu_C = meanRDI_Cu_C+(1.96*seRDI_Cu_C),
    meanRDI_Cu_ADO = mean(percentRDI_Cu_ADO, na.rm = TRUE),
    seRDI_Cu_ADO = std(percentRDI_Cu_ADO),
    CIinf_Cu_ADO = meanRDI_Cu_ADO-(1.96*seRDI_Cu_ADO),
    CIsup_Cu_ADO = meanRDI_Cu_ADO+(1.96*seRDI_Cu_ADO),
    meanRDI_Cu_YA = mean(percentRDI_Cu_YA, na.rm = TRUE),
    seRDI_Cu_YA = std(percentRDI_Cu_YA),
    CIinf_Cu_YA = meanRDI_Cu_YA-(1.96*seRDI_Cu_YA),
    CIsup_Cu_YA = meanRDI_Cu_YA+(1.96*seRDI_Cu_YA),
    meanRDI_Cu_A = mean(percentRDI_Cu_A, na.rm = TRUE),
    seRDI_Cu_A = std(percentRDI_Cu_A),
    CIinf_Cu_A = meanRDI_Cu_A-(1.96*seRDI_Cu_A),
    CIsup_Cu_A = meanRDI_Cu_A+(1.96*seRDI_Cu_A),
    meanRDI_Cu_Apre = mean(percentRDI_Cu_Apre, na.rm = TRUE),
    seRDI_Cu_Apre = std(percentRDI_Cu_Apre),
    CIinf_Cu_Apre = meanRDI_Cu_Apre-(1.96*seRDI_Cu_Apre),
    CIsup_Cu_Apre = meanRDI_Cu_Apre+(1.96*seRDI_Cu_Apre),
    meanRDI_Cu_Alact = mean(percentRDI_Cu_Alact, na.rm = TRUE),
    seRDI_Cu_Alact = std(percentRDI_Cu_Alact),
    CIinf_Cu_Alact = meanRDI_Cu_A-(1.96*seRDI_Cu_Alact),
    CIsup_Cu_Alact = meanRDI_Cu_A+(1.96*seRDI_Cu_Alact),
    ## For Fe
    meanRDI_Fe_YC = mean(percentRDI_Fe_YC, na.rm = TRUE),
    seRDI_Fe_YC = std(percentRDI_Fe_YC),
    CIinf_Fe_YC = meanRDI_Fe_YC-(1.96*seRDI_Fe_YC),
    CIsup_Fe_YC = meanRDI_Fe_YC+(1.96*seRDI_Fe_YC),
    meanRDI_Fe_C = mean(percentRDI_Fe_C, na.rm = TRUE),
    seRDI_Fe_C = std(percentRDI_Fe_C),
    CIinf_Fe_C = meanRDI_Fe_C-(1.96*seRDI_Fe_C),
    CIsup_Fe_C = meanRDI_Fe_C+(1.96*seRDI_Fe_C),
    meanRDI_Fe_ADO = mean(percentRDI_Fe_ADO, na.rm = TRUE),
    seRDI_Fe_ADO = std(percentRDI_Fe_ADO),
    CIinf_Fe_ADO = meanRDI_Fe_ADO-(1.96*seRDI_Fe_ADO),
    CIsup_Fe_ADO = meanRDI_Fe_ADO+(1.96*seRDI_Fe_ADO),
    meanRDI_Fe_YA = mean(percentRDI_Fe_YA, na.rm = TRUE),
    seRDI_Fe_YA = std(percentRDI_Fe_YA),
    CIinf_Fe_YA = meanRDI_Fe_YA-(1.96*seRDI_Fe_YA),
    CIsup_Fe_YA = meanRDI_Fe_YA+(1.96*seRDI_Fe_YA),
    meanRDI_Fe_A = mean(percentRDI_Fe_A, na.rm = TRUE),
    seRDI_Fe_A = std(percentRDI_Fe_A),
    CIinf_Fe_A = meanRDI_Fe_A-(1.96*seRDI_Fe_A),
    CIsup_Fe_A = meanRDI_Fe_A+(1.96*seRDI_Fe_A),
    meanRDI_Fe_Apre = mean(percentRDI_Fe_Apre, na.rm = TRUE),
    seRDI_Fe_Apre = std(percentRDI_Fe_Apre),
    CIinf_Fe_Apre = meanRDI_Fe_Apre-(1.96*seRDI_Fe_Apre),
    CIsup_Fe_Apre = meanRDI_Fe_Apre+(1.96*seRDI_Fe_Apre),
    meanRDI_Fe_Alact = mean(percentRDI_Fe_Alact, na.rm = TRUE),
    seRDI_Fe_Alact = std(percentRDI_Fe_Alact),
    CIinf_Fe_Alact = meanRDI_Fe_A-(1.96*seRDI_Fe_Alact),
    CIsup_Fe_Alact = meanRDI_Fe_A+(1.96*seRDI_Fe_Alact),
    ## For Mn
    meanRDI_Mn_YC = mean(percentRDI_Mn_YC, na.rm = TRUE),
    seRDI_Mn_YC = std(percentRDI_Mn_YC),
    CIinf_Mn_YC = meanRDI_Mn_YC-(1.96*seRDI_Mn_YC),
    CIsup_Mn_YC = meanRDI_Mn_YC+(1.96*seRDI_Mn_YC),
    meanRDI_Mn_C = mean(percentRDI_Mn_C, na.rm = TRUE),
    seRDI_Mn_C = std(percentRDI_Mn_C),
    CIinf_Mn_C = meanRDI_Mn_C-(1.96*seRDI_Mn_C),
    CIsup_Mn_C = meanRDI_Mn_C+(1.96*seRDI_Mn_C),
    meanRDI_Mn_ADO = mean(percentRDI_Mn_ADO, na.rm = TRUE),
    seRDI_Mn_ADO = std(percentRDI_Mn_ADO),
    CIinf_Mn_ADO = meanRDI_Mn_ADO-(1.96*seRDI_Mn_ADO),
    CIsup_Mn_ADO = meanRDI_Mn_ADO+(1.96*seRDI_Mn_ADO),
    meanRDI_Mn_YA = mean(percentRDI_Mn_YA, na.rm = TRUE),
    seRDI_Mn_YA = std(percentRDI_Mn_YA),
    CIinf_Mn_YA = meanRDI_Mn_YA-(1.96*seRDI_Mn_YA),
    CIsup_Mn_YA = meanRDI_Mn_YA+(1.96*seRDI_Mn_YA),
    meanRDI_Mn_A = mean(percentRDI_Mn_A, na.rm = TRUE),
    seRDI_Mn_A = std(percentRDI_Mn_A),
    CIinf_Mn_A = meanRDI_Mn_A-(1.96*seRDI_Mn_A),
    CIsup_Mn_A = meanRDI_Mn_A+(1.96*seRDI_Mn_A),
    meanRDI_Mn_Apre = mean(percentRDI_Mn_Apre, na.rm = TRUE),
    seRDI_Mn_Apre = std(percentRDI_Mn_Apre),
    CIinf_Mn_Apre = meanRDI_Mn_Apre-(1.96*seRDI_Mn_Apre),
    CIsup_Mn_Apre = meanRDI_Mn_Apre+(1.96*seRDI_Mn_Apre),
    meanRDI_Mn_Alact = mean(percentRDI_Mn_Alact, na.rm = TRUE),
    seRDI_Mn_Alact = std(percentRDI_Mn_Alact),
    CIinf_Mn_Alact = meanRDI_Mn_A-(1.96*seRDI_Mn_Alact),
    CIsup_Mn_Alact = meanRDI_Mn_A+(1.96*seRDI_Mn_Alact),
    ## For Se
    meanRDI_Se_YC = mean(percentRDI_Se_YC, na.rm = TRUE),
    seRDI_Se_YC = std(percentRDI_Se_YC),
    CIinf_Se_YC = meanRDI_Se_YC-(1.96*seRDI_Se_YC),
    CIsup_Se_YC = meanRDI_Se_YC+(1.96*seRDI_Se_YC),
    meanRDI_Se_C = mean(percentRDI_Se_C, na.rm = TRUE),
    seRDI_Se_C = std(percentRDI_Se_C),
    CIinf_Se_C = meanRDI_Se_C-(1.96*seRDI_Se_C),
    CIsup_Se_C = meanRDI_Se_C+(1.96*seRDI_Se_C),
    meanRDI_Se_ADO = mean(percentRDI_Se_ADO, na.rm = TRUE),
    seRDI_Se_ADO = std(percentRDI_Se_ADO),
    CIinf_Se_ADO = meanRDI_Se_ADO-(1.96*seRDI_Se_ADO),
    CIsup_Se_ADO = meanRDI_Se_ADO+(1.96*seRDI_Se_ADO),
    meanRDI_Se_YA_A = mean(percentRDI_Se_YA_A, na.rm = TRUE),
    seRDI_Se_YA_A = std(percentRDI_Se_YA_A),
    CIinf_Se_YA_A = meanRDI_Se_YA_A-(1.96*seRDI_Se_YA_A),
    CIsup_Se_YA_A = meanRDI_Se_YA_A+(1.96*seRDI_Se_YA_A),
    meanRDI_Se_Apre = mean(percentRDI_Se_Apre, na.rm = TRUE),
    seRDI_Se_Apre = std(percentRDI_Se_Apre),
    CIinf_Se_Apre = meanRDI_Se_Apre-(1.96*seRDI_Se_Apre),
    CIsup_Se_Apre = meanRDI_Se_Apre+(1.96*seRDI_Se_Apre),
    meanRDI_Se_Alact = mean(percentRDI_Se_Alact, na.rm = TRUE),
    seRDI_Se_Alact = std(percentRDI_Se_Alact),
    CIinf_Se_Alact = meanRDI_Se_Alact-(1.96*seRDI_Se_Alact),
    CIsup_Se_Alact = meanRDI_Se_Alact+(1.96*seRDI_Se_Alact),
    ## For Zn
    meanRDI_Zn_YC = mean(percentRDI_Zn_YC, na.rm = TRUE),
    seRDI_Zn_YC = std(percentRDI_Zn_YC),
    CIinf_Zn_YC = meanRDI_Zn_YC-(1.96*seRDI_Zn_YC),
    CIsup_Zn_YC = meanRDI_Zn_YC+(1.96*seRDI_Zn_YC),
    meanRDI_Zn_C = mean(percentRDI_Zn_C, na.rm = TRUE),
    seRDI_Zn_C = std(percentRDI_Zn_C),
    CIinf_Zn_C = meanRDI_Zn_C-(1.96*seRDI_Zn_C),
    CIsup_Zn_C = meanRDI_Zn_C+(1.96*seRDI_Zn_C),
    meanRDI_Zn_ADO_A = mean(percentRDI_Zn_ADO_A, na.rm = TRUE),
    seRDI_Zn_ADO_A = std(percentRDI_Zn_ADO_A),
    CIinf_Zn_ADO_A = meanRDI_Zn_ADO_A-(1.96*seRDI_Zn_ADO_A),
    CIsup_Zn_ADO_A = meanRDI_Zn_ADO_A+(1.96*seRDI_Zn_ADO_A),
    meanRDI_Zn_YA_Apre = mean(percentRDI_Zn_YA_Apre, na.rm = TRUE),
    seRDI_Zn_YA_Apre = std(percentRDI_Zn_YA_Apre),
    CIinf_Zn_YA_Apre = meanRDI_Zn_YA_Apre-(1.96*seRDI_Zn_YA_Apre),
    CIsup_Zn_YA_Apre = meanRDI_Zn_YA_Apre+(1.96*seRDI_Zn_YA_Apre),
    meanRDI_Zn_Alact = mean(percentRDI_Zn_Alact, na.rm = TRUE),
    seRDI_Zn_Alact = std(percentRDI_Zn_Alact),
    CIinf_Zn_Alact = meanRDI_Zn_Alact-(1.96*seRDI_Zn_Alact),
    CIsup_Zn_Alact = meanRDI_Zn_Alact+(1.96*seRDI_Zn_Alact))

## 2 / Cleaning table ##############################################################################################

names(percentRDI)

percentRDI_propre <- percentRDI
ncol <- ncol(percentRDI_propre)
names(percentRDI_propre)

# Boucle pour As
percentRDI_propre$servAs_inf_Children <- as.numeric(NA)
percentRDI_propre$servAs_sup_Children <- as.numeric(NA)
for (i in 1:nrow(percentRDI_propre)) {
  
  ## 1. Children
  # 1.1. Min-max using min(CIinf) and max(CIsup)
  percentRDI_propre[i,ncol+1] <- min(percentRDI_propre[i,4],
                                     percentRDI_propre[i,8],
                                     percentRDI_propre[i,12])
  percentRDI_propre[i,ncol+2] <- max(percentRDI_propre[i,5],
                                     percentRDI_propre[i,9],
                                     percentRDI_propre[i,13])
  # 1.2. Dealing with unique measures and NAs
  percentRDI_propre[i,ncol+1] <- ifelse(is.na(percentRDI_propre[i,4]) == TRUE,
                                        min(percentRDI_propre[i,2], percentRDI_propre[i,6], percentRDI_propre[i,10]),
                                        percentRDI_propre[i,ncol+1])
  percentRDI_propre[i,ncol+2] <- ifelse(is.na(percentRDI_propre[i,5]) == TRUE,
                                        max(percentRDI_propre[i,2], percentRDI_propre[i,6], percentRDI_propre[i,10]),
                                        percentRDI_propre[i,ncol+2])
}

# Boucle pour Cu
percentRDI_propre$servCu_inf_Children <- as.numeric(NA)
percentRDI_propre$servCu_sup_Children <- as.numeric(NA)
percentRDI_propre$servCu_inf_Adults <- as.numeric(NA)
percentRDI_propre$servCu_sup_Adults <- as.numeric(NA)
for (i in 1:nrow(percentRDI_propre)) {
  
  ## 1. Children
  # 1.1. Min-max using min(CIinf) and max(CIsup)
  percentRDI_propre[i,ncol+3] <- min(percentRDI_propre[i,24],
                                     percentRDI_propre[i,28],
                                     percentRDI_propre[i,32])
  percentRDI_propre[i,ncol+4] <- max(percentRDI_propre[i,25],
                                     percentRDI_propre[i,29],
                                     percentRDI_propre[i,33])
  # 1.2. Dealing with unique measures and NAs
  percentRDI_propre[i,ncol+3] <- ifelse(is.na(percentRDI_propre[i,24]) == TRUE,
                                        min(percentRDI_propre[i,22], percentRDI_propre[i,26], percentRDI_propre[i,30]),
                                        percentRDI_propre[i,ncol+3])
  percentRDI_propre[i,ncol+4] <- ifelse(is.na(percentRDI_propre[i,25]) == TRUE,
                                        max(percentRDI_propre[i,22], percentRDI_propre[i,26], percentRDI_propre[i,30]),
                                        percentRDI_propre[i,ncol+4])
  
  ## 2. Adults
  # 2.1. Min-max using min(CIinf) and max(CIsup)
  percentRDI_propre[i,ncol+5] <- min(percentRDI_propre[i,40],
                                     percentRDI_propre[i,44],
                                     percentRDI_propre[i,48])
  percentRDI_propre[i,ncol+6] <- max(percentRDI_propre[i,41],
                                     percentRDI_propre[i,45],
                                     percentRDI_propre[i,49])
  # 2.2. Dealing with unique measures and NAs
  percentRDI_propre[i,ncol+5] <- ifelse(is.na(percentRDI_propre[i,40]) == TRUE,
                                        min(percentRDI_propre[i,38], percentRDI_propre[i,42], percentRDI_propre[i,46]),
                                        percentRDI_propre[i,ncol+5])
  percentRDI_propre[i,ncol+6] <- ifelse(is.na(percentRDI_propre[i,41]) == TRUE,
                                        max(percentRDI_propre[i,38], percentRDI_propre[i,42], percentRDI_propre[i,46]),
                                        percentRDI_propre[i,ncol+6])
}

# Boucle pour Fe
percentRDI_propre$servFe_inf_Children <- as.numeric(NA)
percentRDI_propre$servFe_sup_Children <- as.numeric(NA)
percentRDI_propre$servFe_inf_Adults <- as.numeric(NA)
percentRDI_propre$servFe_sup_Adults <- as.numeric(NA)
for (i in 1:nrow(percentRDI_propre)) {
  
  ## 1. Children
  # 1.1. Min-max using min(CIinf) and max(CIsup)
  percentRDI_propre[i,ncol+7] <- min(percentRDI_propre[i,52],
                                     percentRDI_propre[i,56],
                                     percentRDI_propre[i,60])
  percentRDI_propre[i,ncol+8] <- max(percentRDI_propre[i,53],
                                     percentRDI_propre[i,57],
                                     percentRDI_propre[i,61])
  # 1.2. Dealing with unique measures and NAs
  percentRDI_propre[i,ncol+7] <- ifelse(is.na(percentRDI_propre[i,52]) == TRUE,
                                        min(percentRDI_propre[i,50], percentRDI_propre[i,54], percentRDI_propre[i,58]),
                                        percentRDI_propre[i,ncol+7])
  percentRDI_propre[i,ncol+8] <- ifelse(is.na(percentRDI_propre[i,53]) == TRUE,
                                        max(percentRDI_propre[i,50], percentRDI_propre[i,54], percentRDI_propre[i,58]),
                                        percentRDI_propre[i,ncol+8])
  
  ## 2. Adults
  # 2.1. Min-max using min(CIinf) and max(CIsup)
  percentRDI_propre[i,ncol+9] <- min(percentRDI_propre[i,68],
                                     percentRDI_propre[i,72],
                                     percentRDI_propre[i,76])
  percentRDI_propre[i,ncol+10] <- max(percentRDI_propre[i,69],
                                      percentRDI_propre[i,73],
                                      percentRDI_propre[i,77])
  # 2.2. Dealing with unique measures and NAs
  percentRDI_propre[i,ncol+9] <- ifelse(is.na(percentRDI_propre[i,68]) == TRUE,
                                        min(percentRDI_propre[i,66], percentRDI_propre[i,70], percentRDI_propre[i,74]),
                                        percentRDI_propre[i,ncol+9])
  percentRDI_propre[i,ncol+10] <- ifelse(is.na(percentRDI_propre[i,69]) == TRUE,
                                         max(percentRDI_propre[i,66], percentRDI_propre[i,70], percentRDI_propre[i,74]),
                                         percentRDI_propre[i,ncol+10])
}

# Boucle pour Mn
percentRDI_propre$servMn_inf_Children <- as.numeric(NA)
percentRDI_propre$servMn_sup_Children <- as.numeric(NA)
percentRDI_propre$servMn_inf_Adults <- as.numeric(NA)
percentRDI_propre$servMn_sup_Adults <- as.numeric(NA)
for (i in 1:nrow(percentRDI_propre)) {
  
  ## 1. Children
  # 1.1. Min-max using min(CIinf) and max(CIsup)
  percentRDI_propre[i,ncol+11] <- min(percentRDI_propre[i,80],
                                      percentRDI_propre[i,84],
                                      percentRDI_propre[i,88])
  percentRDI_propre[i,ncol+12] <- max(percentRDI_propre[i,81],
                                      percentRDI_propre[i,85],
                                      percentRDI_propre[i,89])
  # 1.2. Dealing with unique measures and NAs
  percentRDI_propre[i,ncol+11] <- ifelse(is.na(percentRDI_propre[i,80]) == TRUE,
                                         min(percentRDI_propre[i,78], percentRDI_propre[i,82], percentRDI_propre[i,86]),
                                         percentRDI_propre[i,ncol+11])
  percentRDI_propre[i,ncol+12] <- ifelse(is.na(percentRDI_propre[i,80]) == TRUE,
                                         max(percentRDI_propre[i,78], percentRDI_propre[i,82], percentRDI_propre[i,86]),
                                         percentRDI_propre[i,ncol+12])
  
  ## 2. Adults
  # 2.1. Min-max using min(CIinf) and max(CIsup)
  percentRDI_propre[i,ncol+13] <- min(percentRDI_propre[i,96],
                                      percentRDI_propre[i,100],
                                      percentRDI_propre[i,104])
  percentRDI_propre[i,ncol+14] <- max(percentRDI_propre[i,97],
                                      percentRDI_propre[i,101],
                                      percentRDI_propre[i,105])
  # 2.2. Dealing with unique measures and NAs
  percentRDI_propre[i,ncol+13] <- ifelse(is.na(percentRDI_propre[i,96]) == TRUE,
                                         min(percentRDI_propre[i,94], percentRDI_propre[i,98], percentRDI_propre[i,102]),
                                         percentRDI_propre[i,ncol+13])
  percentRDI_propre[i,ncol+14] <- ifelse(is.na(percentRDI_propre[i,69]) == TRUE,
                                         max(percentRDI_propre[i,94], percentRDI_propre[i,98], percentRDI_propre[i,102]),
                                         percentRDI_propre[i,ncol+14])
}

# Boucle pour Se
percentRDI_propre$servSe_inf_Children <- as.numeric(NA)
percentRDI_propre$servSe_sup_Children <- as.numeric(NA)
percentRDI_propre$servSe_inf_YoungAdults_Adults <- as.numeric(NA)
percentRDI_propre$servSe_sup_YoungAdults_Adults <- as.numeric(NA)
for (i in 1:nrow(percentRDI_propre)) {
  
  ## 1. Children
  # 1.1. Min-max using min(CIinf) and max(CIsup)
  percentRDI_propre[i,ncol+15] <- min(percentRDI_propre[i,108],
                                      percentRDI_propre[i,112],
                                      percentRDI_propre[i,116])
  percentRDI_propre[i,ncol+16] <- max(percentRDI_propre[i,109],
                                      percentRDI_propre[i,113],
                                      percentRDI_propre[i,117])
  # 1.2. Dealing with unique measures and NAs
  percentRDI_propre[i,ncol+15] <- ifelse(is.na(percentRDI_propre[i,108]) == TRUE,
                                         min(percentRDI_propre[i,106], percentRDI_propre[i,110], percentRDI_propre[i,114]),
                                         percentRDI_propre[i,ncol+15])
  percentRDI_propre[i,ncol+16] <- ifelse(is.na(percentRDI_propre[i,109]) == TRUE,
                                         max(percentRDI_propre[i,106], percentRDI_propre[i,110], percentRDI_propre[i,114]),
                                         percentRDI_propre[i,ncol+16])
  
  ## 2. YoungAdults & Adults
  # 2.1. Min-max using min(CIinf) and max(CIsup)
  percentRDI_propre[i,ncol+17] <- min(percentRDI_propre[i,120],
                                      percentRDI_propre[i,124],
                                      percentRDI_propre[i,128])
  percentRDI_propre[i,ncol+18] <- max(percentRDI_propre[i,121],
                                      percentRDI_propre[i,125],
                                      percentRDI_propre[i,129])
  # 2.2. Dealing with unique measures and NAs
  percentRDI_propre[i,ncol+17] <- ifelse(is.na(percentRDI_propre[i,120]) == TRUE,
                                         min(percentRDI_propre[i,118], percentRDI_propre[i,122], percentRDI_propre[i,126]),
                                         percentRDI_propre[i,ncol+17])
  percentRDI_propre[i,ncol+18] <- ifelse(is.na(percentRDI_propre[i,121]) == TRUE,
                                         max(percentRDI_propre[i,118], percentRDI_propre[i,122], percentRDI_propre[i,126]),
                                         percentRDI_propre[i,ncol+18])
}

# Boucle pour Zn
percentRDI_propre$servZn_inf_Children <- as.numeric(NA)
percentRDI_propre$servZn_sup_Children <- as.numeric(NA)
percentRDI_propre$servZn_inf_Adults <- as.numeric(NA)
percentRDI_propre$servZn_sup_Adults <- as.numeric(NA)
for (i in 1:nrow(percentRDI_propre)) {
  
  ## 1. Children
  # 1.1. Min-max using min(CIinf) and max(CIsup)
  percentRDI_propre[i,ncol+19] <- min(percentRDI_propre[i,132],
                                      percentRDI_propre[i,136],
                                      percentRDI_propre[i,140])
  percentRDI_propre[i,ncol+20] <- max(percentRDI_propre[i,133],
                                      percentRDI_propre[i,137],
                                      percentRDI_propre[i,141])
  # 1.2. Dealing with unique measures and NAs
  percentRDI_propre[i,ncol+19] <- ifelse(is.na(percentRDI_propre[i,132]) == TRUE,
                                         min(percentRDI_propre[i,130], percentRDI_propre[i,134], percentRDI_propre[i,138]),
                                         percentRDI_propre[i,ncol+19])
  percentRDI_propre[i,ncol+20] <- ifelse(is.na(percentRDI_propre[i,133]) == TRUE,
                                         max(percentRDI_propre[i,130], percentRDI_propre[i,134], percentRDI_propre[i,138]),
                                         percentRDI_propre[i,ncol+20])
  
  ## 2. Adults
  # 2.1. Min-max using min(CIinf) and max(CIsup)
  percentRDI_propre[i,ncol+21] <- min(percentRDI_propre[i,140],
                                      percentRDI_propre[i,144],
                                      percentRDI_propre[i,148])
  percentRDI_propre[i,ncol+22] <- max(percentRDI_propre[i,141],
                                      percentRDI_propre[i,145],
                                      percentRDI_propre[i,149])
  # 2.2. Dealing with unique measures and NAs
  percentRDI_propre[i,ncol+21] <- ifelse(is.na(percentRDI_propre[i,140]) == TRUE,
                                         min(percentRDI_propre[i,138], percentRDI_propre[i,142], percentRDI_propre[i,146]),
                                         percentRDI_propre[i,ncol+21])
  percentRDI_propre[i,ncol+22] <- ifelse(is.na(percentRDI_propre[i,140]) == TRUE,
                                         max(percentRDI_propre[i,138], percentRDI_propre[i,142], percentRDI_propre[i,146]),
                                         percentRDI_propre[i,ncol+22])
}

rm(i,percentRDI, ncol)

names(percentRDI_propre)

percentRDI_propre <- percentRDI_propre %>% 
  mutate(meanRDI_Cu_Children = (meanRDI_Cu_YC+meanRDI_Cu_C+meanRDI_Cu_ADO)/3,
         meanRDI_Cu_Adults = (meanRDI_Cu_A+meanRDI_Cu_Apre+meanRDI_Cu_Alact)/3,
         meanRDI_Fe_Children = (meanRDI_Fe_YC+meanRDI_Fe_C+meanRDI_Fe_ADO)/3,
         meanRDI_Fe_Adults = (meanRDI_Fe_A+meanRDI_Fe_Apre+meanRDI_Fe_Alact)/3,
         meanRDI_Mn_Children = (meanRDI_Mn_YC+meanRDI_Mn_C+meanRDI_Mn_ADO)/3,
         meanRDI_Mn_Adults = (meanRDI_Mn_A+meanRDI_Mn_Apre+meanRDI_Mn_Alact)/3,
         meanRDI_Se_Children = (meanRDI_Se_YC+meanRDI_Se_C+meanRDI_Se_ADO)/3,
         meanRDI_Se_YoungAdults_Adults = (meanRDI_Se_YA_A+meanRDI_Se_Apre+meanRDI_Se_Alact)/3,
         meanRDI_Zn_Children = (meanRDI_Zn_YC+meanRDI_Zn_C+meanRDI_Zn_ADO_A)/3,
         meanRDI_Zn_Adults = (meanRDI_Zn_ADO_A+meanRDI_Zn_YA_Apre+meanRDI_Zn_Alact)/3) %>% 
  mutate(As_Adults = paste0(round(meanRDI_As_A,0)," [",round(CIinf_As_A,0),"-",round(CIsup_As_A,0),"]"),
         As_Adults = ifelse(is.na(CIinf_As_A), as.character(round(meanRDI_As_A,0)), As_Adults),
         Cu_Children = paste0(round(meanRDI_Cu_Children,0)," [",round(servCu_inf_Children,0),"-",round(servCu_sup_Children),"]"),
         Cu_Children = ifelse(is.na(servCu_inf_Children), as.character(round(meanRDI_Cu_Children,0)), Cu_Children),
         Cu_YoungAdults = paste0(round(meanRDI_Cu_YA,0)," [",round(CIinf_Cu_YA,0),"-",round(CIsup_Cu_YA,0),"]"),
         Cu_YoungAdults = ifelse(is.na(CIinf_Cu_YA), as.character(round(meanRDI_Cu_YA,0)), Cu_YoungAdults),
         Cu_Adults = paste0(round(meanRDI_Cu_Adults,0)," [",round(servCu_inf_Adults,0),"-",round(servCu_sup_Adults,0),"]"),
         Cu_Adults = ifelse(is.na(servCu_inf_Adults), as.character(round(meanRDI_Cu_Adults,0)), Cu_Adults),
         Fe_Children = paste0(round(meanRDI_Fe_Children,0)," [",round(servFe_inf_Children,0),"-",round(servFe_sup_Children),"]"),
         Fe_Children = ifelse(is.na(servFe_inf_Children), as.character(round(meanRDI_Fe_Children,0)), Fe_Children),
         Fe_YoungAdults = paste0(round(meanRDI_Fe_YA,0)," [",round(CIinf_Fe_YA,0),"-",round(CIsup_Fe_YA,0),"]"),
         Fe_YoungAdults = ifelse(is.na(CIinf_Fe_YA), as.character(round(meanRDI_Fe_YA,0)), Fe_YoungAdults),
         Fe_Adults = paste0(round(meanRDI_Fe_Adults,0)," [",round(servFe_inf_Adults,0),"-",round(servFe_sup_Adults,0),"]"),
         Fe_Adults = ifelse(is.na(servFe_inf_Adults), as.character(round(meanRDI_Fe_Adults,0)), Fe_Adults),
         Mn_Children = paste0(round(meanRDI_Mn_Children,1)," [",round(servMn_inf_Children,1),"-",round(servMn_sup_Children,1),"]"),
         Mn_Children = ifelse(is.na(servMn_inf_Children), as.character(round(meanRDI_Mn_Children,1)), Mn_Children),
         Mn_YoungAdults = paste0(round(meanRDI_Mn_YA,1)," [",round(CIinf_Mn_YA,1),"-",round(CIsup_Mn_YA,1),"]"),
         Mn_YoungAdults = ifelse(is.na(CIinf_Mn_YA), as.character(round(meanRDI_Mn_YA,1)), Mn_YoungAdults),
         Mn_Adults = paste0(round(meanRDI_Mn_Adults,1)," [",round(servMn_inf_Adults,1),"-",round(servMn_sup_Adults,1),"]"),
         Mn_Adults = ifelse(is.na(servMn_inf_Adults), as.character(round(meanRDI_Mn_Adults,1)), Mn_Adults),
         Se_Children = paste0(round(meanRDI_Se_Children,0)," [",round(servSe_inf_Children,0),"-",round(servSe_sup_Children),"]"),
         Se_Children = ifelse(is.na(servSe_inf_Children), as.character(round(meanRDI_Se_Children,0)), Se_Children),
         Se_YoungAdults_Adults = paste0(round(meanRDI_Se_YoungAdults_Adults,0)," [",round(servSe_inf_YoungAdults_Adults,0),"-",round(servSe_sup_YoungAdults_Adults,0),"]"),
         Se_YoungAdults_Adults = ifelse(is.na(servSe_inf_YoungAdults_Adults), as.character(round(meanRDI_Se_YoungAdults_Adults,0)), Se_YoungAdults_Adults),
         Zn_Children = paste0(round(meanRDI_Zn_Children,0)," [",round(servZn_inf_Children,0),"-",round(servZn_sup_Children),"]"),
         Zn_Children = ifelse(is.na(servZn_inf_Children), as.character(round(meanRDI_Zn_Children,0)), Zn_Children),
         Zn_YoungAdults = paste0(round(meanRDI_Zn_YA_Apre,0)," [",round(CIinf_Zn_YA_Apre,0),"-",round(CIsup_Zn_YA_Apre,0),"]"),
         Zn_YoungAdults = ifelse(is.na(CIinf_Zn_YA_Apre), as.character(round(meanRDI_Zn_YA_Apre,0)), Zn_YoungAdults),
         Zn_Adults = paste0(round(meanRDI_Zn_Adults,0)," [",round(servZn_inf_Adults,0),"-",round(servZn_sup_Adults,0),"]"),
         Zn_Adults = ifelse(is.na(servZn_inf_Adults), as.character(round(meanRDI_Zn_Adults,0)), Zn_Adults)) %>% 
  left_join(infos, by = "english_name") %>% 
  distinct(english_name, .keep_all = TRUE) %>% 
  mutate(number = as.numeric(english_name)) %>% 
  arrange(number) %>% 
  select(sp_group,english_name,Cu_Children,Cu_YoungAdults,Cu_Adults,
         Fe_Children,Fe_YoungAdults,Fe_Adults,
         Mn_Children,Mn_YoungAdults,Mn_Adults,
         Se_Children,Se_YoungAdults_Adults,
         Zn_Children,Zn_YoungAdults,Zn_Adults,As_Adults)

write.table(percentRDI_propre,paste0(direction,"5_Benefit_percenRDI.csv"),
            sep = ";", dec = ".", row.names = F)

rm(percentRDI_propre)
rm(infos)
rm(interHgSe)



### VI // iAs estimation ##############################################################################################

infos <- data_risk_benefit %>% select(sp_group,english_name) %>% distinct(english_name, .keep_all = TRUE)

iAs_theoretical <- data_risk_benefit %>% 
  mutate(iAs = ifelse(english_name %in% c("Spanner crab","Longlegged spiny lobster","Pronghorn spiny lobster","Painted spiny lobster"),
                      2.2 * TAs_ww /100, # If crustaceans
                      4.2 * TAs_ww /100)) %>%  # If no crustaceans
  group_by(english_name) %>% 
  summarise(mean_iAs = mean(iAs, na.rm = TRUE),
            sd_iAs = sd(iAs, na.rm = TRUE)) %>% 
  mutate(iAs = paste0(round(mean_iAs, 3)," ± ",round(sd_iAs, 3)),
         iAs = ifelse(is.na(sd_iAs), as.character(round(mean_iAs, 3)), iAs)) %>% 
  left_join(infos, by = "english_name") %>% 
  distinct(english_name, .keep_all = TRUE) %>% 
  mutate(number = as.numeric(english_name)) %>% 
  arrange(number) %>% 
  select(sp_group,english_name,iAs)

write.table(iAs_theoretical,paste0(direction,"6_Theoretical_iAS.csv"),
            sep = ";", dec = ".", row.names = F)

rm(iAs_theoretical,infos)

rm(data_risk_benefit)
