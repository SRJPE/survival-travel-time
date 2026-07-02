# Data loading and data frame creation to run survival + travel time models 

#Clear the workspace
rm(list = ls())

# Prepare the environment ------------------------------------------------------------
library(tidyverse)
library(here) #better working directory management
library(SRJPEdata)

# If need to re-install SRJPEdata package
#remotes::install_github("SRJPE/SRJPEdata")

source('scripts/helper_fxt.R')

# Create Sacramento River's tagged fish detection history -------------------------------------------------------------
name <- "SacRiverAllYears"
path <- paste0("./data/", name, "/")
dir.create(path, showWarnings = FALSE) 

studyIDs <- c("ColemanFall_2013","ColemanFall_2016","ColemanFall_2017",
              "CNFH_FMR_2019","CNFH_FMR_2020","CNFH_FMR_2021",
              "ColemanAltRel_2021", "ColemanAltRel_2022",
              "RBDD_2017","RBDD_2018","Wild_stock_Chinook_Rbdd_2021",
              "Wild_stock_Chinook_RBDD_2022","Wild_stock_Chinook_Rbdd_2024",
              "SacRiverSpringJPE_2022","SacRiverSpringJPE_2023", "SacRiverSpringJPE_2024",
              "Seasonal_Survival_2024" ,"Spring_Pulse_2023","Spring_Pulse_2024",
              "DeerCk_Wild_CHK_2018","DeerCk_Wild_CHK_2020",
              "MillCk_Wild_CHK_2013","MillCk_Wild_CHK_2015","MillCk_Wild_CHK_2017",
              "MillCk_Wild_CHK_2022") 

## Retrieve ERDDAP data if first time
all_detections <- lapply(studyIDs, get_detections)

## Save detections to csv files
names <- lapply(studyIDs, function(x) paste0("./data/",name, "/", x, ".csv"))
pmap(list(all_detections, names), write_csv)

# ## Retrieve ERDDAP detection data saved as csv files
# names <- lapply(studyIDs, function(x) paste0("./data/",name, "/", x, ".csv"))
# all_detections <- lapply(names, vroom)

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections)

# Manually select receiver locations to use and combine for Sac River study
reach.meta <- reach.meta %>% 
  filter(GEN %in% c("BattleCk_CNFH_Rel","BattleCk_RST_Rel",
                    "RBDD_Rel","RBDD_Rel_Rec","RBDD1","RBDD2",
                    "Altube Island","Abv_Altube1","Abv_Altube2", 
                    "IrvineFinch_Rel","RB River Park_Rel",
                    "MillCk_RST_Rel", "MillCk2_Rel","DeerCk_RST_Rel",
                    "Abv_WoodsonBr","Blw_Woodson",
                    "ButteBr","BlwButteBr","AbvButteBr",
                     "I80-50_Br","TowerBridge", 
                    "ToeDrainBase","Hwy84Ferry",
                    "BeniciaE","BeniciaW",
                    "ChippsE","ChippsW","Chipps_1","Chipps_2")) %>% 
  mutate(Region = case_when(Region == 'Battle Ck' ~ 'Release',
                            Region == 'Deer Ck' ~ 'Release',
                            Region == 'Mill Ck' ~ 'Release',
                            GEN == 'RBDD_Rel'& Region == 'Upper Sac R' ~ 'Release',
                            GEN == 'RBDD_Rel_Rec'& Region == 'Upper Sac R' ~ 'Release',
                            GEN == "RBDD1" & Region == 'Upper Sac R' ~ 'Release',
                            GEN == "RBDD2" & Region == 'Upper Sac R' ~ 'Release',
                            GEN == "Altube Island" & Region == 'Upper Sac R' ~ 'Release',
                            GEN == "Abv_Altube1" & Region == 'Upper Sac R' ~ 'Release',
                            GEN == "Abv_Altube2" & Region == 'Upper Sac R' ~ 'Release',
                            GEN == "IrvineFinch_Rel" & Region == 'Upper Sac R' ~ 'Release',
                            GEN == "RB River Park_Rel" & Region == 'Upper Sac R' ~ 'Release',
                            Region == 'Yolo Bypass' ~ 'Lower Sac R',
                            Region == 'North Delta' ~ 'Lower Sac R',
                            Region == 'West Delta' ~ 'End',
                            Region == 'Carquinez Strait' ~ 'End',
                            TRUE ~ Region))

# Aggregate receiver locations and detections 
all_aggregated <- lapply(all_detections, aggregate_GEN_Sac)

# Create detection history
all_DH <- lapply(all_aggregated, make_DH)

# Avoid duplicate min_time for a given location and remove Nas
DH_df <- bind_rows(all_DH) %>% 
  finalize_DH() %>% 
  mutate(year=year(min_time),
         month=month(min_time))

# clean up DH by removing some fish
# Remove  2021  and 2022 non battle creek release location from Net pen studies
netpen_ids_to_remove <- c(
  sprintf("ColemanAltRel-2021-%03d", 1:300),
  sprintf("ColemanAltRel-2022-%03d", 301:600))

# Remove winter and late fall-run fish tagged
lfw_ids_to_remove <- c("RBDD_WS2022-028",sprintf("SS_2024-%03d", 1:766))

DH_df <- DH_df %>% 
  filter(FishID != "SP2023-1033") %>% # fish without size info
  filter(!FishID %in% lfw_ids_to_remove) %>% 
  filter(!FishID %in% netpen_ids_to_remove)

write.csv(DH_df,"./data/DetectionHistorySac.csv", row.names=FALSE)

# Create Sacramento River's tagged fish inp file  -------------------------------------------------------------
all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>%
  bind_rows()

# Add in fish information to inp file 
all.inp <- all.inp %>% 
  left_join(TaggedFish %>% select(fish_id, fish_length, fish_weight, fish_type,fish_release_date, 
                                  release_location),by = c("FishID" = "fish_id"))

# clean up inp by removing some fish 
all.inp <- all.inp %>% 
           filter(FishID != "SP2023-1033") %>% # fish without size info
           filter(!FishID %in% c("RBDD_WS2021-001","RBDD_WS2021-008")) %>%  # fish not found in the DH table
           filter(release_location != "SacRiver_ButteCity_Rel") %>% # non battle creek release location from Net pen studies
           filter(!fish_type %in% c("CNFH Late-fall Chinook","LSNFH Winter Chinook", "LAD_winter_Chinook")) # winter and late fall-run fish 

# calculate K factor and make sure fish cov are numeric factors
all.inp$fish_weight <- as.numeric(all.inp$fish_weight)
all.inp$fish_length <- as.numeric(all.inp$fish_length)
all.inp$fish_k <- (100 * all.inp$fish_weight) / (all.inp$fish_length^3)
all.inp$fish_k <- as.numeric(all.inp$fish_k)

# add year factor
all.inp$year <- as.factor(year(as.Date(all.inp$fish_release_date,format="%m/%d/%Y")))

write.csv(all.inp,"./data/SacInp.csv", row.names=FALSE)

# Create Feather River's tagged fish detection history -------------------------------------------------------------
name <- "FeatherRiverAllYears"
path <- paste0("./data/", name, "/")
dir.create(path, showWarnings = FALSE) 

studyIDs <- c("FR_Spring_2013","FR_Spring_2014","FR_Spring_2015","FR_Spring_2019","FR_Spring_2020",
              "FR_Spring_2021","FR_Spring_2023","FR_Spring_2024") 

## Retrieve ERDDAP data if first time
all_detections_F <- lapply(studyIDs, get_detections)
 
## Save detections to csv files
names <- lapply(studyIDs, function(x) paste0("./data/",name, "/", x, ".csv"))
pmap(list(all_detections_F, names), write_csv)

# ## Retrieve ERDDAP detection data saved as csv files
# names <- lapply(studyIDs, function(x) paste0("./data/",name, "/", x, ".csv"))
# all_detections_F <- lapply(names, vroom)

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections_F)

# Manually select receiver locations to use and combine for Feather River study
reach.meta <- reach.meta %>% 
  filter(GEN %in% c("FR_Gridley_Rel","FR_Boyds_Rel","FR_Boyds_Rel_Rec",
                    "I80-50_Br","TowerBridge", 
                    "ToeDrainBase","Hwy84Ferry",
                    "BeniciaE","BeniciaW","ChippsE","ChippsW"))%>% 
  mutate(Region = case_when(Region == 'Yolo Bypass' ~ 'Lower Sac R',
                            Region == 'North Delta' ~ 'Lower Sac R',
                            Region == 'West Delta' ~ 'End',
                            Region == 'Carquinez Strait' ~ 'End',
                            TRUE ~ Region))

# Aggregate receiver locations and detections
all_aggregated_F <- lapply(all_detections_F, aggregate_GEN_Feather)

# Create detection history 
all_DH_F <- lapply(all_aggregated_F, make_DH)

DH_F_df <- bind_rows(all_DH_F) %>% 
  finalize_DH() %>% 
  mutate(year=year(min_time))

write.csv(DH_F_df,"./data/DetectionHistoryFea.csv", row.names=FALSE)

# Create Feather River's tagged fish inp file  -------------------------------------------------------------
all_EH_F <- lapply(all_aggregated_F, make_EH)

all.inp_F <- pmap(list(all_detections_F,all_EH_F), create_inp) %>%
  bind_rows()

# Add in fish information to inp file 
# First add in fish info 
all.inp_F <- all.inp_F %>% 
  left_join(TaggedFish %>% select(fish_id, fish_length, fish_weight, fish_type,fish_release_date, 
                                  release_location),by = c("FishID" = "fish_id"))

# calculate K factor and make sure fish cov are numeric factors
all.inp_F$fish_weight <- as.numeric(all.inp_F$fish_weight)
all.inp_F$fish_length <- as.numeric(all.inp_F$fish_length)
all.inp_F$fish_k <- (100 * all.inp_F$fish_weight) / (all.inp_F$fish_length^3)
all.inp_F$fish_k <- as.numeric(all.inp_F$fish_k)

# add year and release location
all.inp_F$year <- as.factor(year(as.Date(all.inp_F$fish_release_date,format="%m/%d/%Y")))

write.csv(all.inp_F,"./data/FeatherInp.csv", row.names=FALSE)

# Create Butte Creek's tagged fish detection history -------------------------------------------------------------
name <- "ButteCreekAllYears"
path <- paste0("./data/", name, "/")
dir.create(path, showWarnings = FALSE) 

studyIDs <- c("SB_Spring_2015","SB_Spring_2016","SB_Spring_2017","SB_Spring_2018",
              "SB_Spring_2019","SB_Spring_2023", "Butte_Sink_2023","Butte_Sink_2024",
              "Butte_Sink_2021","Upper_Butte_2019","Upper_Butte_2020",'Upper_Butte_2021') 

## Retrieve ERDDAP data if first time
all_detections_B <- lapply(studyIDs, get_detections)

## Save detections to csv files
names <- lapply(studyIDs, function(x) paste0("./data/",name, "/", x, ".csv"))
pmap(list(all_detections_B, names), write_csv)

# ## Retrieve ERDDAP detection data saved as csv files
# names <- lapply(studyIDs, function(x) paste0("./data/",name, "/", x, ".csv"))
# all_detections_B <- lapply(names, vroom)

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections_B)

# Manually select receiver locations to use and combine for Butte Creek study
reach.meta <- reach.meta %>% 
  filter(GEN %in% c("UpperButte_RST_Rel","UpperButte_RST","UpperButte_SKWY",
                    "SutterBypass_Weir2_RST_Rel","SutterBypass Weir2 RST",
                    "Butte_Blw_Sanborn_Rel","Sanborn_Slough_Rel","North_Weir_Rel","Laux Rd",
                    "I80-50_Br","TowerBridge", 
                    "ToeDrainBase","Hwy84Ferry",
                    "BeniciaE","BeniciaW","ChippsE","ChippsW"))%>% 
  mutate(Region = case_when(Region == 'Sutter Bypass' ~ 'Butte Ck',
                            Region == 'Yolo Bypass' ~ 'Lower Sac R',
                            Region == 'North Delta' ~ 'Lower Sac R',
                            Region == 'West Delta' ~ 'End',
                            Region == 'Carquinez Strait' ~ 'End',
                            TRUE ~ Region))

# Aggregate receiver locations and detections
all_aggregated_B <- lapply(all_detections_B, aggregate_GEN_Butte)

# Create detection history 
all_DH_B <- lapply(all_aggregated_B, make_DH)

DH_B_df <- bind_rows(all_DH_B) %>% 
  finalize_DH() %>% 
  mutate(year=year(min_time))

write.csv(DH_B_df,"./data/DetectionHistoryBut.csv", row.names=FALSE)

# Create Butte Creek's tagged fish inp file  -------------------------------------------------------------
all_EH_B <- lapply(all_aggregated_B, make_EH)

all.inp_B <- pmap(list(all_detections_B,all_EH_B), create_inp) %>%
  bind_rows()

# Add in fish information to inp file 
# First add in fish info 
all.inp_B <- all.inp_B %>% 
  left_join(TaggedFish %>% select(fish_id, fish_length, fish_weight, fish_type,fish_release_date, 
                                  release_location),by = c("FishID" = "fish_id"))

# calculate K factor and make sure fish cov are numeric factors
all.inp_B$fish_weight <- as.numeric(all.inp_B$fish_weight)
all.inp_B$fish_length <- as.numeric(all.inp_B$fish_length)
all.inp_B$fish_k <- (100 * all.inp_B$fish_weight) / (all.inp_B$fish_length^3)
all.inp_B$fish_k <- as.numeric(all.inp_B$fish_k)

# add year and release location
all.inp_B$year <- as.factor(year(as.Date(all.inp_B$fish_release_date,format="%m/%d/%Y")))

write.csv(all.inp_B,"./data/ButteInp.csv", row.names=FALSE)

# Create Yuba River's tagged fish detection history  -----------------------------------------------------------------------------------------
name <- "YubaRiverAllYears"
path <- paste0("./data/", name, "/")
dir.create(path, showWarnings = FALSE) 

studyIDs <- c("Lower_Yuba_FRH_Chinook_2021", "Lower_Yuba_FRH_Chinook_2022",
              "Lower_Yuba_FRH_Chinook_2024") 

## Retrieve ERDDAP data if first time
all_detections_Y <- lapply(studyIDs, get_detections)

## Save detections to csv files
names <- lapply(studyIDs, function(x) paste0("./data/",name, "/", x, ".csv"))
pmap(list(all_detections_Y, names), write_csv)

# ## Retrieve ERDDAP detection data saved as csv files
# names <- lapply(studyIDs, function(x) paste0("./data/",name, "/", x, ".csv"))
# all_detections_Y <- lapply(names, vroom)

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections_Y)

# Manually select receiver locations to use and combine for Butte Creek study
reach.meta <- reach.meta %>% 
  filter(GEN %in% c("Daguerre_Dam_Rel" , "Englebright_Dam_Rel",
                    "I80-50_Br","TowerBridge", 
                    "ToeDrainBase","Hwy84Ferry",
                    "SacTrawl1","SacTrawl2",
                    "BeniciaE","BeniciaW","ChippsE","ChippsW"))%>% 
  mutate(Region = case_when(Region == 'Yolo Bypass' ~ 'Lower Sac R',
                            Region == 'North Delta' ~ 'Lower Sac R',
                            Region == 'West Delta' ~ 'End',
                            Region == 'Carquinez Strait' ~ 'End',
                            TRUE ~ Region))

# Aggregate receiver locations and detections
all_aggregated_Y <- lapply(all_detections_Y, aggregate_GEN_Yuba)

# Create detection history
all_DH_Y <- lapply(all_aggregated_Y, make_DH)

DH_Y_df <- bind_rows(all_DH_Y) %>% 
  finalize_DH() %>% 
  mutate(year=year(min_time))

write.csv(DH_Y_df,"./data/DetectionHistoryYuba.csv", row.names=FALSE)

# Create Yuba River's tagged fish inp file ------------------------------------------------------------------
all_EH_Y<- lapply(all_aggregated_Y, make_EH)

all.inp_Y <- pmap(list(all_detections_Y,all_EH_Y), create_inp) %>%
  bind_rows()

# Add in fish information to inp file 
# First add in fish info 
all.inp_Y <- all.inp_Y %>% 
  left_join(TaggedFish %>% select(fish_id, fish_length, fish_weight, fish_type,fish_release_date, 
                                  release_location),by = c("FishID" = "fish_id"))

# calculate K factor and make sure fish cov are numeric factors
all.inp_Y$fish_weight <- as.numeric(all.inp_Y$fish_weight)
all.inp_Y$fish_length <- as.numeric(all.inp_Y$fish_length)
all.inp_Y$fish_k <- (100 * all.inp_Y$fish_weight) / (all.inp_Y$fish_length^3)
all.inp_Y$fish_k <- as.numeric(all.inp_Y$fish_k)

# add year and release location
all.inp_Y$year <- as.factor(year(as.Date(all.inp_Y$fish_release_date,format="%m/%d/%Y")))

write.csv(all.inp_Y,"./data/YubaInp.csv", row.names=FALSE)

# Combine Feather and Butte (and eventually yuba) inp files ---------------------------------------------
all.inp_FB <- data.frame(rbind(all.inp_F,all.inp_B))

write.csv(all.inp_FB,"./data/FeatherButteInp.csv", row.names=FALSE)

all.inp_FBY <- data.frame(rbind(all.inp_F,all.inp_B, all.inp_Y))

write.csv(all.inp_FBY,"./data/FeatherButteYubaInp.csv", row.names=FALSE)

# Create full data set for CJS survival + travel time modelling -----------------------------------------------------------------------------
### Read required data 
## Inp files
d_Sac <- read.csv(here('data','SacInp.csv'),  stringsAsFactors = F)
d_FeaBut <-  read.csv(here('data','FeatherButteInp.csv'),  stringsAsFactors = F)

## Covariate data for from Flow West 
Envdata_flowwest <- SRJPEdata::forecast_covariates
unique(Envdata_flowwest$name)

maxflow_Sac <-Envdata_flowwest %>% 
  filter(name=="monthly_max_flow" & stream == "sacramento river" &
           site_group == "red bluff diversion dam") %>% 
  mutate(Maxflowsac = value,
         Maxflow.z = scale(Maxflowsac)) %>% 
  select(year,month,Maxflowsac,Maxflow.z)

maxflow_Mill <- Envdata_flowwest %>% 
  filter(name=="monthly_max_flow" & stream == "mill creek") %>% 
  mutate(Maxflowmil = value) %>% 
  select(year,month,Maxflowmil)

maxflow_Deer <- Envdata_flowwest %>% 
  filter(name=="monthly_max_flow" & stream == "deer creek") %>% 
  mutate(Maxflowdee = value) %>% 
  select(year,month,Maxflowdee)

maxflow_Battle <- Envdata_flowwest %>% 
  filter(name=="monthly_max_flow" & stream == "battle creek") %>% 
  mutate(Maxflowbat = value) %>% 
  select(year,month,Maxflowbat)

maxflow_uppersac <- left_join(maxflow_Sac,maxflow_Mill, by=c("year","month")) %>% 
  left_join(maxflow_Deer, by=c("year","month")) %>% 
  left_join(maxflow_Battle, by=c("year","month"))

maxflow_Butte <- Envdata_flowwest %>% 
  filter(name=="monthly_max_flow" & stream == "butte creek") %>% 
  mutate(Maxflowbut = value) %>% 
  select(year,month,Maxflowbut)

maxflow_Feather <- Envdata_flowwest %>% 
  filter(name=="monthly_max_flow" & stream == "feather river") %>% 
  mutate(Maxflowfea = value) %>% 
  select(year,month,Maxflowfea)

maxflow_ButteFeather <- left_join(maxflow_Butte,maxflow_Feather, by=c("year","month")) %>% 
  mutate(MeanflowB = mean(Maxflowbut,na.rm=TRUE),
         SdflowB = sd(Maxflowbut,na.rm=TRUE),
         MeanflowF = mean(Maxflowfea,na.rm=TRUE),
         SdflowF = sd(Maxflowfea,na.rm=TRUE))

flowexceed_Sac <-  Envdata_flowwest %>% 
  filter(name== "3_category_flow_exceedance_year_type" & 
           stream == "sacramento river") %>% 
  mutate(flowexceedtype_sac = text_value,
         year=water_year) %>% 
  select(year,flowexceedtype_sac)

flowexceed_Deer <-  Envdata_flowwest %>% 
  filter(name== "3_category_flow_exceedance_year_type" & stream == "deer creek") %>% 
  mutate(flowexceedtype_dee = text_value,
         year=water_year) %>% 
  select(year,flowexceedtype_dee)

flowexceed_Mill <-  Envdata_flowwest %>% 
  filter(name== "3_category_flow_exceedance_year_type" & stream == "mill creek") %>% 
  mutate(flowexceedtype_mil = text_value,
         year=water_year) %>% 
  select(year,flowexceedtype_mil)

flowexceed_Battle <-  Envdata_flowwest %>% 
  filter(name== "3_category_flow_exceedance_year_type" & stream == "battle creek") %>% 
  mutate(flowexceedtype_bat = text_value,
         year=water_year) %>% 
  select(year,flowexceedtype_bat)

flowexceed_uppersac <- left_join(flowexceed_Sac,flowexceed_Deer, by="year") %>%
  left_join(flowexceed_Mill, by="year") %>% 
  left_join(flowexceed_Battle, by="year") %>% 
  distinct(year, .keep_all = TRUE)

flowexceed_Butte <-  Envdata_flowwest %>% 
  filter(name== "3_category_flow_exceedance_year_type" & stream == "butte creek") %>% 
  mutate(flowexceedtype_but = text_value,
         year=water_year)%>% 
  select(year,flowexceedtype_but)

flowexceed_Feather <-  Envdata_flowwest %>% 
  filter(name== "3_category_flow_exceedance_year_type" & stream == "feather river") %>% 
  mutate(flowexceedtype_fea = text_value,
         year=water_year)%>% 
  select(year,flowexceedtype_fea) %>% 
  add_row(year=2024,flowexceedtype_fea="Average") # Added this row for now because no value in 2024

flowexceed_ButteFeather <- left_join(flowexceed_Butte,flowexceed_Feather, by="year") %>% 
  distinct(year, .keep_all = TRUE)

# flow data at Freeport from CDEC
FPTflow <- read.csv(here('data','FPTflow.csv'),  stringsAsFactors = F) %>%
  mutate(flow = as.numeric(gsub(",", "", VALUE)),
         date=as.Date(DATE.TIME,format = "%m/%d/%Y"),
         year=year(date),
         month=month(date))

MonthFPTflow <- FPTflow %>% 
  group_by(year,month) %>% 
  summarize(MaxFPTflow=max(flow,na.rm=TRUE)) %>% 
  mutate(MaxFPTflow.z=scale(MaxFPTflow)) %>% 
  ungroup()

## Final Sac data table
d_Sac_sort <- d_Sac %>% 
  filter(fish_type != "LAD_winter_Chinook") %>%  # remove winter-run fish tagged during seasonal tagging study
  # Add month number to match with flowwest enviro monthly data
  mutate(month = month(as.Date(fish_release_date,format="%m/%d/%Y"))) %>% 
  left_join(maxflow_Sac, by = c('year','month')) %>% 
  left_join(MonthFPTflow,by= c('year','month')) %>% 
  left_join(flowexceed_Sac, by = 'year') %>% 
  # Add dummy water year type variable,
  mutate(WY2 = case_when(year %in% c(2013,2015,2016,2018,2020, 2021,2022) ~ 0,
                         year %in% c(2017, 2019, 2023,2024) ~ 1), #0 or 1 for 2 water year type categories: dry (C,D,BN) and wet (AN,W) water year 
        WY3 = case_when(year %in% c(2015, 2021, 2022) ~ 0,
                         year %in% c(2013, 2016, 2018, 2020) ~ 1,
                         year %in% c(2017, 2019, 2023,2024) ~ 2), #0, 1 or 2 for 3 water year type categories: C, D-BN, AN-W water year
        Flowexceed_Sac = case_when(flowexceedtype_sac =="Dry" ~ 0,
                                 flowexceedtype_sac == "Average" ~ 1,
                                 flowexceedtype_sac == "Wet" ~ 2),
        rel_site = case_when(release_location == "BattleCk_CNFH_Rel" ~ "Battle",
                          release_location == "MillCk_RST_Rel" ~ "Mill",
                          release_location == "DeerCk_RST_Rel" ~ "Deer",
                          release_location %in% c("Altube Island","IrvineFinch_Rel","RBDD_Rel") ~ 
                            "Sac"),
         firstCap = 1,  # define first capture location, it is always the release location
         length.z = scale(fish_length), #standardized length
         weight.z = scale(fish_weight),#standardized weight
         k.z = scale(fish_k)) %>%  #standardized condition factor
  group_by(FishID) %>% 
  # find last capture location for each fish and each potential capture history ch
  mutate(lastCap = 
           case_when(ch == 10000 ~ 1,
                             ch == 11000 ~ 2,
                             ch %in% c(11100,10100) ~ 3,
                             ch %in% c(11110,10110,10010,11010) ~ 4 ,
                             ch %in% c(11111,10111,11011,11101,11001,10011,10101,10001) ~ 5),
         dist_rlwoodson = case_when(rel_site == "Battle" ~ 91.8, # dist from Battle Creek to Woodson Bridge
                                    rel_site  == "Sac" ~ 36.9, # dist from RBDD to Woodson Bridge
                                    rel_site  == "Mill" ~ 25.7, # dist from Mill Creek to Woodson Bridge
                                    rel_site  == "Deer" ~ 16.6), # dist from Deer Creek to Woodson Bridge
         dist_woodsonbutte = 88, # distance in km from Woodson Bridge to Butte Bridge
         dist_buttesac = 170, # distance in km from Butte Bridge to Sac
         dist_sacdelta = 110,
         dist_rlwoodson.z = dist_rlwoodson/100, # standardize distances per 100km
         dist_woodsonbutte.z = dist_woodsonbutte/100,# standardize distances per 100km
         dist_buttesac.z= dist_buttesac/100, # standardize distances per 100km
         dist_sacdelta.z =dist_sacdelta/100) %>% # standardize distances per 100km
  ungroup() %>% 
  arrange(year,release_location)

## Final Feather + Butte table
d_FeaBut_sort <- d_FeaBut %>% 
  # Add month number to match with flowwest enviro monthly data
  mutate(month = month(as.Date(fish_release_date,format="%m/%d/%Y"))) %>% 
  left_join(maxflow_ButteFeather, by = c('year','month')) %>% 
  left_join(MonthFPTflow,by= c('year','month')) %>% 
  left_join(flowexceed_ButteFeather, by = 'year') %>% 
  # Add dummy water year type variable
  mutate(WY2 = case_when(year %in% c(2013,2014,2015,2016,2018,2020, 2021) ~ 0,
                         year %in% c(2017, 2019, 2023,2024) ~ 1), #0 or 1 for 2 water year type categories: dry (C,D,BN) and wet (AN,W) water year 
         WY3 = case_when(year %in% c(2014,2015, 2021) ~ 0,
                         year %in% c(2013, 2016, 2018, 2020) ~ 1,
                         year %in% c(2017, 2019, 2023,2024) ~ 2), #0, 1 or 2 for 3 water year type categories: C, D-BN, AN-W water year
         rel_site = case_when(release_location == "UpperButte_RST_Rel" ~ "Butte_upper",
                              release_location == "Sanborn_Slough_Rel" ~ "Butte_sanborn",
                              release_location == "Butte_Blw_Sanborn_Rel" ~ "Butte_blwsanborn",
                              release_location == "SutterBypass_Weir2_RST_Rel" ~ "Butte_weir2",
                              release_location == "Laux Rd" ~ "Butte_laux",
                              release_location == "North_Weir_Rel" ~ "Butte_northweir",
                              release_location == 'FR_Boyds_Rel' ~ "Feather_boyds",
                              release_location == 'FR_Gridley_Rel' ~ "Feather_gridley"),
         release_trib = case_when(rel_site %in% c("Butte_upper","Butte_sanborn", "Butte_blwsanborn
                                                  Butte_weir2","Butte_laux","Butte_northweir") ~ 'Butte',
                                  TRUE ~ 'Feather'),
         firstCap = 1,  # define first capture location, it is always the release location
         length.z = scale(fish_length), #standardized length
         weight.z = scale(fish_weight),#standardized weight
         k.z = scale(fish_k)) %>%  #standardized condition factor)  
  group_by(FishID) %>% 
  # find last capture location for each fish and each potential capture history ch
  mutate(lastCap = case_when(ch == 100 ~ 1,
                             ch == 110 ~ 2,
                             ch == 111 ~ 3,
                             ch == 101 ~ 3),
         trib_ind = case_when(release_location %in% c('FR_Boyds_Rel','FR_Gridley_Rel') ~ 2,  # Feather = 2
                              TRUE ~ 1), # Butte = 1
         dist_rlsac = case_when(rel_site == "Butte_blwsanborn" ~ 120, # dist from Butte_Blw_Sanborn to Sac
                                rel_site == "Butte_laux" ~ 103.5, # dist from Laux Road to Sac
                                rel_site == "Butte_northweir" ~ 117, # dist from North Weir to Sac
                                rel_site == "Butte_sanborn" ~ 116.7, # dist from Sanborn Slough to Sac
                                rel_site == "Butte_weir2" ~ 78, # dist from Sutter Bypass Weir 2 to Sac
                                rel_site == "Butte_upper" ~ 168.5, # dist from Upper Butte to Sac
                                rel_site == "Feather_gridley" ~ 115, # dist from Gridley to Sac
                                rel_site == "Feather_boyds" ~ 69), # dist from Boyds to Sac
         dist_sacdelta = 110,
         dist_rlsac.z =dist_rlsac/100, # standardize distances per 100km
         dist_sacdelta.z = dist_sacdelta/100,
         Maxflow = case_when(release_location %in% c('FR_Boyds_Rel','FR_Gridley_Rel') ~ Maxflowfea,
                             TRUE ~ Maxflowbut),
         Maxflow.z = case_when(release_location %in% c('FR_Boyds_Rel','FR_Gridley_Rel') ~ 
                               (Maxflow-MeanflowF)/SdflowF,
                               TRUE ~  (Maxflow-MeanflowB)/SdflowB),
         flowexceedtype = case_when(release_location %in% c('FR_Boyds_Rel','FR_Gridley_Rel') ~ 
                                    flowexceedtype_fea,
                                    TRUE ~ flowexceedtype_but),
         FlowexceedT = case_when(flowexceedtype =="Dry" ~ 0,
                                flowexceedtype == "Average" ~ 1,
                                flowexceedtype == "Wet" ~ 2)) %>% 
  ungroup() %>% 
  arrange(year,release_location)

write.csv(d_Sac_sort,here("data", "Sac_data.csv"), row.names=FALSE)
write.csv(d_FeaBut_sort,here("data", "FeaBut_data.csv"), row.names=FALSE)

# Create data set for CWT survival + travel time model -----------------------------------

### Load CWT release data
drel0 <- hatchery_release %>% 
  mutate(mid_release_date = as.Date(mid_release_date),
         release_month = month(mid_release_date)) 

#only use release groups released in one day with an average forklength from RBDD or Battle
drel1 <- subset(drel0, date_span==1 & is.na(avg_length)==F &
               (release_location_name=="COLEMAN NFH"
                | release_location_name=="BATTLE CREEK BELOW CNFH"
                | release_location_name=="BATTLE CREEK NFK WILDAT"
                | release_location_name=="SAC R BEL RBDD"
                |release_location_name=="SAC R RED BLUFF DIV DAM"))

drel2 <- subset(drel1, run=="fall") #to exclude late fall and winter run tagged fish
drel2 <- drel2 %>% 
  filter(release_group_id != 873) %>% # Avg_length is inconsistent so remove for now
  distinct(release_group_id, .keep_all = TRUE)

unique(sort(drel2$year))

### Load CWT recapture data
drec0 <- rst_cwt_recaptures #RST recaptures at Knights Landing

drec1 <- subset(drec0,is.na(release_group_id)==F) %>% #drop any records without release id as can't link to release table
  mutate(Year= year(as.Date(date)))

drec2 <- drec1 %>% 
  select(-year) %>% 
  distinct(release_group_id, date, forklength, tag_code, .keep_all = TRUE)

unique(sort(drec2$Year))

### Combine release and recapture data by release group ID 
drerelrec <- left_join(drel2,drec2,by=c("release_group_id")) %>% 
  mutate(cwt_tt=as.numeric(date-mid_release_date), # estimate cwt travel time
         # define release location area to match later with forecast object
         relloc_area = case_when(release_location_name %in% c("SAC R BEL RBDD","SAC R RED BLUFF DIV DAM") ~
                                   "RBDD",
                                 release_location_name %in% c("COLEMAN NFH","BATTLE CREEK BELOW CNFH",
                                                              "BATTLE CREEK NFK WILDAT") ~
                                   "Battle")) %>% 
  filter(is.na(cwt_tt) | (cwt_tt > 0 & cwt_tt < 60)) %>% # filter out negative and too long travel time
  filter(!year %in% c(1981:1987)) %>% # remove early years with no flow covariate information
  select(release_group_id,release_location_name,run,avg_weight,avg_length,mid_release_date,group_total_marked_N,
         group_total_release_N,month,year,monthly_max_flow,date,forklength,tag_code,cwt_tt,relloc_area) %>% 
  distinct()

unique(sort(drerelrec$year))

write.csv(drerelrec,here("data", "hatchery_relrec.csv"), row.names=FALSE)

### Summary stats
summary_drerelrec_release <- drerelrec %>%
  distinct(release_group_id, .keep_all = TRUE) %>%
  group_by(year,month,relloc_area) %>%
  summarise(
    n_groups = n_distinct(release_group_id),
    total_released = sum(group_total_marked_N, na.rm = TRUE),
    avg_length = mean(avg_length, na.rm = TRUE),
    .groups = "drop"
  )

summary_drerelrec_tt <- drerelrec %>%
  filter(!is.na(cwt_tt)) %>%
  group_by(year,month,relloc_area) %>%
  summarise(
    mean_cwt_tt = mean(cwt_tt, na.rm = TRUE),
    .groups = "drop"
  )

summary_drerelrec <- summary_drerelrec_release %>%
  left_join(summary_drerelrec_tt, by = c("year", "month", "relloc_area"))
