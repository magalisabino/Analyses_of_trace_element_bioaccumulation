##--------------------------------------------------------------------------------------------------------
## SCRIPT : Creation of the databases used in the PhD manuscript to compute all analyses. These
##          databases include all samples collected from all species of interest; when only
##          some (or just one) species were used for a chapter, a subset was created with only
##          the data of interest in the dedicated script.
##
## As part of :
##        Magali SABINO PhD - "Bioaccumulation of trace elements in Seychelles marine food webs"
##
## Author : Magali Sabino
## First created : 2022-01-11
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
lapply(c("tidyverse", "openxlsx", "lubridate"),
       library, character.only = TRUE)



### I // Creation of database with trace element concentrations and stable isotope values ##############################################################################################

## 1 / Dataset importation

data_TE_SI <- read.xlsx("C:/Users/msabin02/Nextcloud/PhD/R_projects/Final_scripts/seyfish_database_TE_SI.xlsx", detectDates = TRUE) # Modify with your own file path

# This step allows cleaning of the database and adding information about the species of interest
data_TE_SI <- data_TE_SI %>% 
  filter(!species_code_fao == "SKJ",
         !tm_old == "F") %>% 
  distinct(organism_identifier, .keep_all = TRUE) %>% 
  group_by(english_name) %>% 
  mutate(N = n()) %>% 
  ungroup() %>% 
  mutate(longitude_deg_dec = ifelse(is.na(longitude_deg_dec) == TRUE, (longitude_deg_dec_min+longitude_deg_dec_max)/2, longitude_deg_dec),
         latitude_deg_dec = ifelse(is.na(latitude_deg_dec) == TRUE, (latitude_deg_dec_min+latitude_deg_dec_max)/2, latitude_deg_dec)) %>% 
  # Add Habitat.3
  mutate(Habitat.3 = NA,
         Habitat.3 = ifelse(english_name %in% c("Big blue octopus","Streamlined spinefoot","Yellowtail emperor",
                                                "Sky emperor","Rainbow runner","Bigeye trevally","Green jobfish",
                                                "Two-spot red snapper","Humpback red snapper","Dogtooth tuna",
                                                "Pickhandle barracuda","Peacock hind","Bluefin trevally",
                                                "Blacktip grouper","Honeycomb grouper","Brown-marbled grouper",
                                                "Golden trevally","Shoemaker spinefoot","Smalltooth emperor",
                                                "Bigeye snapper","Blackeye emperor","Painted spiny lobster"),
                            "Coral reefs", Habitat.3),
         Habitat.3 = ifelse(english_name %in% c("Pronghorn spiny lobster","Longlegged spiny lobster",
                                                "Deepwater longtail red snapper","Eightbar grouper"),
                            "Rocky reefs",Habitat.3),
         Habitat.3 = ifelse(english_name %in% c("Blue-barred parrotfish","Spanner crab","Slender emperor",
                                                "Pink ear emperor","Rosy goatfish","Dash-and-dot goatfish"),
                            "Sandy areas",Habitat.3),
         Habitat.3 = ifelse(english_name %in% c("Bludger","Yellowspotted trevally","Emperor red snapper",
                                                "Tomato hind","Yellow-edged lyretail","Spinner shark",
                                                "Elongate surgeonfish","Longspine grouper",
                                                "Blue-lined large-eye bream","Humphead snapper","Malabar trevally",
                                                "Spot-tail shark"),
                            "Rocky and coral reefs",Habitat.3),
         Habitat.3 = ifelse(english_name %in% c("Brownspotted grouper","White-blotched grouper","Spangled emperor",
                                                "Blacktip shark"),"Diverse",Habitat.3),
         Habitat.3 = ifelse(english_name %in% c("Grey reef shark","Kawakawa","Little tunny(=Atl.black skipj)",
                                                "Scalloped hammerhead","Indian mackerel"),
                            "Pelagic-neritic", Habitat.3),
         Habitat.3 = ifelse(english_name %in% c("Common dolphinfish","Great hammerhead"),"Epipelagic",Habitat.3),
         #Habitat.3 = ifelse(english_name %in% c("Skipjack tuna"),"Mesopelagic",Habitat.3),
         Habitat.3 = ifelse(english_name %in% c("Tiger shark"),"Benthopelagic", Habitat.3),
         Habitat.3 = ifelse(english_name %in% c("Swordfish"),"Bathypelagic",Habitat.3)) %>% 
  # Add fishing_area
  mutate(fishing_area = "Plateau",
         fishing_area = ifelse(english_name %in% c("Pronghorn spiny lobster","Longlegged spiny lobster","Painted spiny lobster",
                                                   "Big blue octopus","Streamlined spinefoot","Blue-barred parrotfish",
                                                   "Shoemaker spinefoot","Elongate surgeonfish"),
                               "Reef", fishing_area),
         fishing_area = ifelse(habitat.1 %in% c("pelagic-oceanic"), "Offshore", fishing_area)) %>% 
  # Add vertical_habitat
  mutate(vertical_habitat = habitat.1,
         vertical_habitat = ifelse(species_code_fao %in% c("USY","IGA","IUU","LHV","QZH","RFP"), "benthic", vertical_habitat),
         vertical_habitat = ifelse(habitat.1 %in% c("pelagic-neritic","pelagic"), "pelagic (neritic)", vertical_habitat),
         vertical_habitat = ifelse(habitat.1 %in% c("pelagic-oceanic"), "pelagic (oceanic)", vertical_habitat)) %>% 
  # Add Classification
  mutate(Classification = "Teleost fish",
         Classification = ifelse(english_name %in% c("Longlegged spiny lobster","Pronghorn spiny lobster",       
                                                     "Painted spiny lobster","Spanner crab"), "Crustacean", Classification),
         Classification = ifelse(english_name %in% c("Big blue octopus"), "Cephalopod", Classification),
         Classification = ifelse(english_name %in% c("Grey reef shark","Spinner shark","Blacktip shark",
                                                     "Spot-tail shark","Scalloped hammerhead","Tiger shark",
                                                     "Great hammerhead"), "Elasmobranchs", Classification)) %>% 
  mutate(english_name = factor(english_name, levels = c("Longlegged spiny lobster","Pronghorn spiny lobster",
                                                        "Painted spiny lobster","Spanner crab","Big blue octopus",
                                                        "Blue-barred parrotfish","Streamlined spinefoot",
                                                        "Shoemaker spinefoot","Rosy goatfish",
                                                        "Dash-and-dot goatfish","Slender emperor",
                                                        "Elongate surgeonfish","Yellowtail emperor",
                                                        "Spangled emperor","Smalltooth emperor",
                                                        "Blackeye emperor","Sky emperor","Pink ear emperor",
                                                        "Blue-lined large-eye bream",
                                                        "Two-spot red snapper","Deepwater longtail red snapper",
                                                        "Humpback red snapper","Bigeye snapper","Emperor red snapper",
                                                        "Humphead snapper","Green jobfish","Blacktip grouper",
                                                        "Honeycomb grouper","Brownspotted grouper","Brown-marbled grouper",
                                                        "Eightbar grouper","White-blotched grouper","Longspine grouper",
                                                        "Tomato hind","Peacock hind","Yellow-edged lyretail","Bludger",
                                                        "Bluefin trevally","Bigeye trevally","Golden trevally",
                                                        "Malabar trevally","Yellowspotted trevally","Pickhandle barracuda",
                                                        "Kawakawa","Little tunny(=Atl.black skipj)","Indian mackerel","Dogtooth tuna",
                                                        "Grey reef shark","Spinner shark","Blacktip shark","Spot-tail shark",
                                                        "Scalloped hammerhead",
                                                        "Tiger shark","Great hammerhead","Swordfish")),
         Habitat.3 = factor(Habitat.3, levels = c("Rocky reefs","Coral reefs","Rocky and coral reefs",
                                                  "Sandy areas","Diverse","Pelagic-neritic","Benthopelagic","Epipelagic","Bathypelagic"))) %>%
  # Add location (nearshore vs offshore)
  mutate(loc = NA,
         loc = ifelse(habitat.1 %in% c("demersal","benthic","pelagic-neritic","pelagic"), "Nearshore", loc),
         loc = ifelse(habitat.1 %in% c("pelagic-oceanic"), "Offshore", loc),
         loc = factor(loc, levels = c("Nearshore","Offshore"))) %>% 
  # Add funct_group
  mutate(habitat = vertical_habitat,
         habitat = ifelse(habitat == "pelagic (neritic)", "pelagic-neritic",habitat),
         habitat = ifelse(habitat == "pelagic (oceanic)", "pelagic-oceanic",habitat),
         funct_group = paste(habitat,Classification, sep = " "),
         funct_group = factor(funct_group, levels = c("benthic Crustacean","benthic Cephalopod","benthic Teleost fish","demersal Teleost fish",
                                                      "pelagic-neritic Teleost fish","pelagic-neritic Elasmobranchs",
                                                      "pelagic-oceanic Elasmobranchs","pelagic-oceanic Teleost fish"))) %>% 
  # Add unique column length (all types of measures, i.e. lower jaw fork length, carapace length, fork length, etc)
  # are spread into several columns, but we need a unique column with these measures
  mutate(length = dorsal_mantle_length,
         length = ifelse(is.na(length) == TRUE, lowerjawfork_length, length),
         length = ifelse(is.na(length) == TRUE, carapace_length, length),
         length = ifelse(is.na(length) == TRUE, fork_length, length),
         length = ifelse(is.na(length) == TRUE, total_length, length),
         length = ifelse(is.na(length) == TRUE, standard_length, length)) %>% 
  select(organism_identifier,english_name, N, species_code_fao, scientific_name, french_name,
         seychelles_creole_name, family, order, Classification, latitude_deg_dec, longitude_deg_dec,loc, funct_group, fishery_type, fishing_area, gear_code,
         distribution, migration, habitat,habitat.1, habitat.2, vertical_habitat, Habitat.3, depth_range, 
         depth_preferred, feeding1, feeding2, fishery_type, date, season, 
         length,sex, water_p, d13C, d15N, CN,
         TAs, Cd, Co, Cr, Cu, Fe, Mn, Ni, Pb, Se, V, Zn, Ag, THg) %>% 
  rename(c_sp_fao = species_code_fao, longitude = longitude_deg_dec, latitude = latitude_deg_dec)

data_TE_SI <- droplevels(data_TE_SI)



## 2 / Correction of data < LOQ

# As data below the LOQ are concentrations that could not be measured (i.e. comprised between 0 and the LOQ),
# we need to replace these values either with "NA" for calculating mean TE concentrations in all species
# (nutritional purpose), or with a random value comprised between 0 and the LOQ for statistical analyses.

## 2.1. IF YOU DID NOT REPLACE NAs WITH RANDOM VALUES COMPRISED BETWEEN 0 AND LOQ YET

# Given LOQ have less decimals (usually no more than 3) than measured data -> we need to count decimals to
# determine which number is an LOQ and which is a data we need to keep

# (Step 0) Creation of a function to count decimals
count_decimals = function(x) {
  #length zero input
  if (length(x) == 0) return(numeric())

  #count decimals
  x_nchr = x %>% abs() %>% as.character() %>% nchar() %>% as.numeric()
  x_int = floor(x) %>% abs() %>% nchar()
  x_nchr = x_nchr - 1 - x_int
  x_nchr[x_nchr < 0] = 0

  x_nchr
}

# (Step 1) Creation of object with only LOQ
LOQ_data <- data_TE_SI %>%
  filter(!is.na(Cu), # Here, we need to suppress all lines that have not be filled, thus we need a TE column with a detection frequency = 100%
         !organism_identifier %in% c("VLO005","LZG012","NXM004","ETH008",
                                    "EFH008","LTQ011","KAR017","KAR023",
                                    "NXM002","KAK001","NXM001","IUU002",
                                    "RFP002","RFP005","zour02","zour11",
                                    "zour12","zour13","LTS003")) %>% # profiles that are not complete, probably contamination
  mutate(nb_TAs = count_decimals(TAs),
         nb_Cd = count_decimals(Cd),
         nb_Co = count_decimals(Co),
         nb_Cr = count_decimals(Cr),
         nb_Cu = count_decimals(Cu),
         nb_Fe = count_decimals(Fe),
         nb_Mn = count_decimals(Mn),
         nb_Ni = count_decimals(Ni),
         nb_Pb = count_decimals(Pb),
         nb_Se = count_decimals(Se),
         nb_V = count_decimals(V),
         nb_Zn = count_decimals(Zn),
         nb_Ag = count_decimals(Ag)) %>%
  mutate(TAs = ifelse(nb_TAs <= 3, TAs, NA),
         Cd = ifelse(nb_Cd <= 3, Cd, NA),
         Co = ifelse(nb_Co <= 3, Co, NA),
         Cr = ifelse(nb_Cr <= 3, Cr, NA),
         Cu = ifelse(nb_Cu <= 3, Cu, NA),
         Fe = ifelse(nb_Fe <= 3, Fe, NA),
         Mn = ifelse(nb_Mn <= 3, Mn, NA),
         Ni = ifelse(nb_Ni <= 3, Ni, NA),
         Pb = ifelse(nb_Pb <= 3, Pb, NA),
         Se = ifelse(nb_Se <= 3, Se, NA),
         V = ifelse(nb_V <= 3, V, NA),
         Zn = ifelse(nb_Zn <= 3, Zn, NA),
         Ag = ifelse(nb_Ag <= 3, Ag, NA)) %>%
  select(organism_identifier, TAs, Cd, Co, Cr, Cu, Fe, Mn, Ni, Pb, Se, V, Zn, Ag) %>%
  gather(TE,LOQ,-organism_identifier)

# (Step 2) Suppressing LOQ data in the database and replacing them by NAs
data_TE_SI_LOQ <- data_TE_SI %>%
  filter(!is.na(Cu),
         !organism_identifier %in% c("VLO005","LZG012","NXM004","ETH008",
                                    "EFH008","LTQ011","KAR017","KAR023",
                                    "NXM002","KAK001","NXM001","LTS003",
                                    "IUU002","RFP002","RFP005","zour02",
                                    "zour11","zour12","zour13")) %>%
  mutate(nb_TAs = count_decimals(TAs),
         nb_Cd = count_decimals(Cd),
         nb_Co = count_decimals(Co),
         nb_Cr = count_decimals(Cr),
         nb_Cu = count_decimals(Cu),
         nb_Fe = count_decimals(Fe),
         nb_Mn = count_decimals(Mn),
         nb_Ni = count_decimals(Ni),
         nb_Pb = count_decimals(Pb),
         nb_Se = count_decimals(Se),
         nb_V = count_decimals(V),
         nb_Zn = count_decimals(Zn),
         nb_Ag = count_decimals(Ag)) %>%
  mutate(TAs = ifelse(nb_TAs <= 3, NA, TAs),
         Cd = ifelse(nb_Cd <= 3, NA, Cd),
         Co = ifelse(nb_Co <= 3, NA, Co),
         Cr = ifelse(nb_Cr <= 3, NA, Cr),
         Cu = ifelse(nb_Cu <= 3, NA, Cu),
         Fe = ifelse(nb_Fe <= 3, NA, Fe),
         Mn = ifelse(nb_Mn <= 3, NA, Mn),
         Ni = ifelse(nb_Ni <= 3, NA, Ni),
         Pb = ifelse(nb_Pb <= 3, NA, Pb),
         Se = ifelse(nb_Se <= 3, NA, Se),
         V = ifelse(nb_V <= 3, NA, V),
         Zn = ifelse(nb_Zn <= 3, NA, Zn),
         Ag = ifelse(nb_Ag <= 3, NA, Ag)) %>%
  select(organism_identifier, TAs, Cd, Co, Cr, Cu, Fe, Mn, Ni, Pb, Se, V, Zn, Ag) %>%
  gather(TE,value,-organism_identifier) %>%
  left_join(LOQ_data, by = c("organism_identifier","TE"))

# (Step 3) Random generation of values in ]0;LOQ[
for (i in 1:nrow(data_TE_SI_LOQ)){
  if (is.na(data_TE_SI_LOQ[i,]$value)){
    data_TE_SI_LOQ[i,]$value <- runif(n = 1, min = 0, max = data_TE_SI_LOQ[i,]$LOQ)
  }
}

# (Step 4) Cleansing data
data_TE_SI_LOQ <- data_TE_SI_LOQ %>%
  select(-LOQ) %>%
  spread(TE, value) %>% 
  rename(Ag_stat = Ag, Cd_stat = Cd, Co_stat = Co, Cr_stat = Cr, Cu_stat = Cu, Fe_stat = Fe,
         Mn_stat = Mn, Ni_stat = Ni, Pb_stat = Pb, Se_stat = Se, TAs_stat = TAs, V_stat = V,
         Zn_stat = Zn)

# Calculating detection frequency for each trace element
data_TE_SI_LOQ %>%
  select(-organism_identifier) %>%
  gather(metal,value) %>%
  filter(!is.na(value)) %>%
  group_by(metal) %>%
  summarise(n = n()) %>%
  mutate(n_tot = 1039,
         percent_above_LOQ = n*100/n_tot) %>%
  select(metal,percent_above_LOQ) %>%
  spread(metal,percent_above_LOQ)
rm(data_TE_SI_LOQ)
  
# (Step 5) Saving generated TE data
write.csv(data_TE_SI_LOQ, file = "C:/Users/msabin02/Nextcloud/PhD/R_projects/Final_scripts/stat_TE_Data.csv")

rm(LOQ_data,data_TE_SI_LOQ,i)

# CAUTION! Now that TE data for statistical analyses have been completed and saved, you can
# skip steps 0 to 5 and just import your file whenever you want. This prevents the generation
# of unique data everytime you run this script.


## 2.2. IF YOU HAVE ALREADY GENERATED DATA FOR STATISTICAL ANALYSES

# Importation of TE data for statistical analyses
stat_TE_Data <- read.csv2("C:/Users/msabin02/Nextcloud/PhD/R_projects/Final_scripts/stat_TE_Data.csv", dec = ".")

# Replacing LOQ with NAs and adding TE data for statistical analyses
data_TE_SI <- data_TE_SI %>% 
  mutate(nb_TAs = count_decimals(TAs),
         nb_Cd = count_decimals(Cd),
         nb_Co = count_decimals(Co),
         nb_Cr = count_decimals(Cr),
         nb_Cu = count_decimals(Cu),
         nb_Fe = count_decimals(Fe),
         nb_Mn = count_decimals(Mn),
         nb_Ni = count_decimals(Ni),
         nb_Pb = count_decimals(Pb),
         nb_Se = count_decimals(Se),
         nb_V = count_decimals(V),
         nb_Zn = count_decimals(Zn),
         nb_Ag = count_decimals(Ag)) %>% 
  mutate(TAs = ifelse(nb_TAs <= 3 | is.na(nb_TAs) == TRUE, NA, TAs),
         Cd = ifelse(nb_Cd <= 3 | is.na(nb_Cd) == TRUE, NA, Cd),
         Co = ifelse(nb_Co <= 3 | is.na(nb_Co) == TRUE, NA, Co),
         Cr = ifelse(nb_Cr <= 3 | is.na(nb_Cr) == TRUE, NA, Cr),
         Cu = ifelse(nb_Cu <= 3 | is.na(nb_Cu) == TRUE, NA, Cu),
         Fe = ifelse(nb_Fe <= 3 | is.na(nb_Fe) == TRUE, NA, Fe),
         Mn = ifelse(nb_Mn <= 3 | is.na(nb_Mn) == TRUE, NA, Mn),
         Ni = ifelse(nb_Ni <= 3 | is.na(nb_Ni) == TRUE, NA, Ni),
         Pb = ifelse(nb_Pb <= 3 | is.na(nb_Pb) == TRUE, NA, Pb),
         Se = ifelse(nb_Se <= 3 | is.na(nb_Se) == TRUE, NA, Se),
         V = ifelse(nb_V <= 3 | is.na(nb_V) == TRUE, NA, V),
         Zn = ifelse(nb_Zn <= 3 | is.na(nb_Zn) == TRUE, NA, Zn),
         Ag = ifelse(nb_Ag <= 3 | is.na(nb_Ag) == TRUE, NA, Ag)) %>% 
  left_join(stat_TE_Data, by = "organism_identifier") %>% 
  select(-nb_Ag,-nb_TAs,-nb_Cd,-nb_Co,-nb_Cr,-nb_Cu,-nb_Fe,-nb_Mn,-nb_Ni,
         -nb_Pb,-nb_Se,-nb_V,-nb_Zn)

rm(stat_TE_Data)



## 3 / Conversion from dry weight (dw) to wet weight (ww)

# Calculation of conversion coefficient for each species
coeff <- data_TE_SI %>% 
  select(c_sp_fao, water_p) %>% 
  mutate(coeff = 100/(100 - water_p)) %>% 
  group_by(c_sp_fao) %>% 
  mutate(coeff_ww = mean(coeff, na.rm = TRUE)) %>% 
  ungroup() %>% 
  distinct(c_sp_fao, .keep_all = TRUE)
mean_coeff <- mean(coeff$coeff_ww, na.rm = TRUE)
coeff <- coeff %>%
  mutate(coeff_ww = replace_na(coeff_ww, 4.390106)) %>% 
  select(c_sp_fao, coeff_ww)

# Adding conversion coefficient in database then convert to ww
data_TE_SI <- data_TE_SI %>% 
  left_join(coeff, by = "c_sp_fao") %>% 
  mutate(Ag_ww = Ag/coeff_ww,
         Cd_ww = Cd/coeff_ww,
         Co_ww = Co/coeff_ww,
         Cr_ww = Cr/coeff_ww,
         Cu_ww = Cu/coeff_ww,
         Fe_ww = Fe/coeff_ww,
         Mn_ww = Mn/coeff_ww,
         Ni_ww = Ni/coeff_ww,
         Pb_ww = Pb/coeff_ww,
         Se_ww = Se/coeff_ww,
         TAs_ww = TAs/coeff_ww,
         V_ww = V/coeff_ww,
         Zn_ww = Zn/coeff_ww,
         Hg_ww = THg/coeff_ww) %>% 
  mutate(Ag_stat = Ag_stat/coeff_ww,
         Cd_stat = Cd_stat/coeff_ww,
         Co_stat = Co_stat/coeff_ww,
         Cr_stat = Cr_stat/coeff_ww,
         Cu_stat = Cu_stat/coeff_ww,
         Fe_stat = Fe_stat/coeff_ww,
         Mn_stat = Mn_stat/coeff_ww,
         Ni_stat = Ni_stat/coeff_ww,
         Pb_stat = Pb_stat/coeff_ww,
         Se_stat = Se_stat/coeff_ww,
         TAs_stat = TAs_stat/coeff_ww,
         V_stat = V_stat/coeff_ww,
         Zn_stat = Zn_stat/coeff_ww,
         THg_stat = Hg_ww)

rm(coeff, mean_coeff)



## 4 / Correction of d13C values for SWO samples for which lipids were not removed prior to SI analyses

# The equation used to correct d13C values and model constants were selected/calculated prior to
# this. Selection process and constant calculation are given in the related script entitled
# "model_choice_d13C_correction".

# Creation of vector with the names of non-lipid-free SWO samples (i.e. samples SWO001 to SWO099)
non_delip_samples <- "SWO001"
for (i in 2:99){
  if (i < 10){
    non_delip_samples <- c(non_delip_samples, paste0("SWO00",i))
  } else{
    non_delip_samples <- c(non_delip_samples, paste0("SWO0",i))
  }
}

# Correction of d13C values for non-lipid-free SWO samples
data_TE_SI <- data_TE_SI %>%
  mutate(d13C_corr = NA,
         d13C_corr = ifelse(organism_identifier %in% non_delip_samples, (d13C + ((7.05 * CN + -22.4)/(CN + -0.44))), d13C_corr)) %>% 
  mutate(d13C = ifelse(is.na(d13C_corr) == FALSE, d13C_corr, d13C)) %>% 
  select(-d13C_corr)

rm(non_delip_samples, i)


## 5 / Replacement of longitudes and latitudes for spiny lobsters

LOB_coordinates <- read.xlsx("C:/Users/msabin02/Nextcloud/PhD/R_projects/Final_scripts/LOB_coordinates.xlsx")

LOB_coordinates <- LOB_coordinates %>% 
  rename(organism_identifier = fish_identifier) %>% 
  dplyr::select(-site_name)

data_TE_SI <- data_TE_SI %>% 
  left_join(LOB_coordinates, by = "organism_identifier") %>% 
  mutate(longitude = ifelse(scientific_name %in% c("Panulirus penicillatus", "Panulirus longipes", "Panulirus versicolor"), long_centroid, longitude),
         latitude = ifelse(scientific_name %in% c("Panulirus penicillatus", "Panulirus longipes", "Panulirus versicolor"), lat_centroid, latitude)) %>% 
  dplyr::select(-long_centroid, -lat_centroid)

rm(LOB_coordinates)



### II // Creation of database with fatty acid profiles measured in spiny lobsters ##############################################################################################

## 1 / Dataset importation

LOB_FA_percent <- read.xlsx("C:/Users/msabin02/Nextcloud/PhD/R_projects/Final_scripts/LOB_FA_data.xlsx")

## 2 / Adding column bleaching event

LOB_FA_percent <- LOB_FA_percent %>% 
  mutate(sampling_year = as.character(sampling_year),
         bleaching_event = ifelse(sampling_year %in% c("2014","2015"), "Before","After"))

## 3 / Adding columns site name and substrate ####

LOB_FA_percent <- LOB_FA_percent %>% 
  # Add site names
  mutate(site_name = NA,
         site_name = ifelse(lat_centroid == -4.60 & long_centroid == 55.43, "Beau vallon",site_name),
         site_name = ifelse(lat_centroid == -4.62 & long_centroid == 55.39, "Belombre",site_name),
         site_name = ifelse(lat_centroid == -4.63 & long_centroid == 55.37, "Baie Ternay",site_name),
         site_name = ifelse(lat_centroid == -4.64 & long_centroid == 55.38, "Baie Ternay",site_name),
         site_name = ifelse(lat_centroid == -4.64 & long_centroid == 55.37, "Baie Ternay",site_name),
         site_name = ifelse(lat_centroid == -4.66 & long_centroid == 55.37, "Conception",site_name),
         site_name = ifelse(lat_centroid == -4.66 & long_centroid == 55.40, "Ile Let",site_name),
         site_name = ifelse(lat_centroid == -4.67 & long_centroid == 55.40, "Ile Therese East",site_name),
         site_name = ifelse(lat_centroid == -4.67 & long_centroid == 55.41, "Port Glaud",site_name),
         site_name = ifelse(lat_centroid == -4.68 & long_centroid == 55.40, "Ile Therese West",site_name),
         site_name = ifelse(lat_centroid == -4.68 & long_centroid == 55.43, "Ile Aux Vache",site_name),
         site_name = ifelse(lat_centroid == -4.71 & long_centroid == 55.47, "Anse Boileau - Anse Louis",site_name),
         site_name = ifelse(lat_centroid == -4.72 & long_centroid == 55.48, "Roche Canon",site_name),
         site_name = ifelse(lat_centroid == -4.74 & long_centroid == 55.47, "Ile Chauve Souris - Anse Soleil - Anse Poul Bleu",site_name),
         site_name = ifelse(lat_centroid == -4.75 & long_centroid == 55.46, "Four Season",site_name),
         site_name = ifelse(lat_centroid == -4.76 & long_centroid == 55.46, "Four Season",site_name),
         site_name = ifelse(lat_centroid == -4.77 & long_centroid == 55.49, "Maravi - Chez Baptista",site_name),
         site_name = ifelse(lat_centroid == -4.78 & long_centroid == 55.49, "Maravi - Chez Baptista",site_name)) %>% 
  # Add biotope
  mutate(biotope = "Granite reef",
         biotope = ifelse(site_name %in% c("Beau vallon","Ile Therese East","Port Glaud"),"Carbonate reef",biotope))



### III // Creation of database with fatty acid profiles measured in swordfish ##############################################################################################

## 1 / Dataset importation

SWO_FA <- read.xlsx("C:/Users/msabin02/Nextcloud/PhD/R_projects/Final_scripts/SWO_FA_data.xlsx")

SWO_FA <- SWO_FA %>% 
  select(-species_code_fao,-operator_name.y,-project,-organism_identifier_origin,
         -organism_length_unit,-dorsal_mantle_length,-gear_code,-vessel_name,
         -vessel_code,-landing_site,-remarks_capture,-scientific_name,-english_name,
         -french_name,-author,-family,-order,-analysis_lab,-fa_c_unit)


## 2 / Get sampling year in the date column

SWO_FA$date <- year(as.Date(SWO_FA$date,"%d/%m/%Y"))
class(SWO_FA$date)
SWO_FA$date <- as.character(SWO_FA$date)


## 3 / Conversion dw to ww

SWO_dwtoww <- SWO_FA %>% 
  filter(operator_name.x == "S Hollanda") %>% 
  mutate_at(vars(contains("_c")), funs(.*(100-67.1)/100))

SWO_ww <- SWO_FA %>% 
  filter(operator_name.x == "M Sabino")

SWO_FA_ww <- rbind(SWO_dwtoww,SWO_ww)

rm(SWO_dwtoww,SWO_ww)
rm(SWO_FA)


## 4 / Conversion content to percentage of total fatty acid (%TFA)

# Conversion
SWO_FA_percent <- SWO_FA_ww %>% 
  rename(TFA = tfa_c) %>% 
  mutate_at(vars(contains("_c")), funs(.*100/TFA)) %>% 
  select(-TFA)

# Nettoyage FA < 0.8%
FAmean <- colMeans(SWO_FA_percent[,9:81], na.rm = TRUE) # Calculating means of FA
others <- names(which(FAmean <0.8)) # Give names of FA which are less than 1%
SWO_FA_percent$others <- rowSums(SWO_FA_percent[,others], na.rm = TRUE)
SWO_FA_percent <- SWO_FA_percent[, !names(SWO_FA_percent)%in%(others)]
SWO_FA_percent <- droplevels(SWO_FA_percent)

rm(FAmean, others)

SWO_FA_percent <- SWO_FA_percent %>% 
  select(-c17_1w7_c,-c18_1w11_c,-sex,-season,-operator_name.x)

rm(SWO_FA_ww)
