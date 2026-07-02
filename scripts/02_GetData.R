# Dataframe loading and setting up model variables

#Clear the workspace
rm(list = ls())

# Prepare the environment ------------------------------------------------------------
library(tidyverse)
library(here) #better working directory management
library(ggpubr)
library(lubridate)
library(SRJPEdata)

# If SRJPEdata package needs to be updated
# remotes::install_github("SRJPE/SRJPEdata")

# Get water year type data from CDEC -------------------------------------
cdec_wyt <- readxl::read_excel(here("data", "cdec-water-year-type-jun-2025.xlsx")) %>%
  transmute(
    year = as.integer(WY),
    WYT = trimws(WYT)
  ) %>%
  distinct(year, .keep_all = TRUE) %>%
  mutate(
    dry_wet = case_when(
      WYT %in% c("C", "D", "BN") ~ "Dry",
      WYT %in% c("AN", "W") ~ "Wet",
      TRUE ~ NA_character_
    )
  ) %>%
  select(year, WYT, dry_wet)

# Load CJS data frame ---------------------------------------------------
d_Sac_sort <- read.csv(here("data", "Sac_data.csv"),  stringsAsFactors = F)
d_FeaBut_sort <- read.csv(here("data", "FeaBut_data.csv"),  stringsAsFactors = F)

# Set CJS variables for Sac mainstem survival model -------------------------------------------------------------------
Nind <- dim(d_Sac_sort)[1]
Nreaches <- 4 #1=release-woodson, 2=woodson-butte, 3=butte-sac, 4=sac-delta
Ndetlocs <- 5
Nsz <- 25 # Number of size classes to plot size effect over 
NsX <- 25  # Number of continuous variable values
Nrst <- 6
rch_covind=c(1,2,3,4)
(Year_sac <- sort(unique(d_Sac_sort$year)))
(Year_FeaBut <- sort(unique(d_FeaBut_sort$year)))
(Year_all <- sort(unique(c(Year_sac,Year_FeaBut))))
Nyrs <- length(Year_all)
site_sac <- c('Sac','Battle','Mill', 'Deer')
Nsites <- 6
RelGp <- unique(d_Sac_sort$StudyID)
Nrg <- length(RelGp)
firstCap <- d_Sac_sort$firstCap
lastCap <- d_Sac_sort$lastCap
WY2 <- d_Sac_sort$WY2
WY3 <- d_Sac_sort$WY3
MaxflowSac.z <- d_Sac_sort$Maxflow.z
MaxflowSac <- d_Sac_sort$Maxflowsac
MaxflowDelta.z <- d_Sac_sort$MaxFPTflow.z
MaxflowDelta <- d_Sac_sort$MaxFPTflow
FlowexceedSac <- d_Sac_sort$Flowexceed_Sac
FL <- d_Sac_sort$fish_length
WGT <- d_Sac_sort$fish_weight
CF <- d_Sac_sort$fish_k
CH <- data.frame(cbind(as.integer(substr(d_Sac_sort$ch, 1, 1)),as.integer(substr(d_Sac_sort$ch,2, 2)),
                       as.integer(substr(d_Sac_sort$ch, 3, 3)),as.integer(substr(d_Sac_sort$ch, 4, 4)),
                       as.integer(substr(d_Sac_sort$ch, 5, 5))))

Rmult <- data.frame(cbind(d_Sac_sort$dist_rlwoodson.z, d_Sac_sort$dist_woodsonbutte.z, # reach distances standardized per 100km
                          d_Sac_sort$dist_buttesac.z,d_Sac_sort$dist_sacdelta.z))
RmultSac <- 0.43 # distance for prediction model corresponds to the average distance from all release locations to Woodson Bridge per 100km
Rmult_Tis <- 0.78 # distance for prediction model from Butte Bridge to Tisdale per 100km
Rmult_Kni <- 1.27 # distance for prediction model from Butte Bridge to Knights Landing per 100km
dist_rstwoodson <- c(124.9512976,93.2391355490999,26.1273401272,18.3412666561999,123.0699,33.15917) # Distance from Upper Clear C, Battle C, Mill C, Deer C, Redding, Red Bluff to Woodson Bridge 
dist_rstwoodson.z <- dist_rstwoodson/100
dist_rsttisd <- c(283.37420088157,251.6620387843,184.5502433623,176.7641698913) # Distance from Upper Clear C, Battle C, Mill C, Deer C to Tisdale
dist_rstknights <- c(332.15522617927,300.443064082,233.33126866,225.545195189,330.2738,240.3631) # Distance from Upper Clear C, Battle C, Mill C, Deer C, Redding, Red Bluff, to Knights Landing
dist_rsttisd.z <- dist_rsttisd/100
dist_rstknights.z <- dist_rstknights/100
Rmultrst <- dist_rstwoodson.z # distance for prediction model corresponds to the distance from RST locations to Woodson Bridge

indsite <- vector(length=Nind)
for(i in 1:Nind){
  indsite[i] <- which(site_sac == d_Sac_sort$rel_site[i])
}

yrind <- vector(length=Nind)
rgind <- yrind
for(i in 1:Nind){
  yrind[i] <- which(Year_all == d_Sac_sort$year[i])
  rgind[i] <- which(RelGp == d_Sac_sort$StudyID[i])
}

### Add index for covariates
rgwy3.df <- d_Sac_sort %>% 
            group_by(StudyID) %>% 
            dplyr::reframe(ind = unique(WY3), 
                           year = unique(year)) %>% 
            arrange(year) %>% 
            ungroup()

rgwy3_ind <- rgwy3.df$ind

rgwy2.df <- d_Sac_sort %>% 
  group_by(StudyID) %>% 
  dplyr::reframe(ind = unique(WY2), 
                 year = unique(year)) %>% 
  arrange(year) %>% 
  ungroup()

rgwy2_ind <- rgwy2.df$ind

flowexceedSac.df <- d_Sac_sort %>% 
  group_by(StudyID) %>% 
  dplyr::reframe(ind = unique(Flowexceed_Sac), 
                 year = unique(year)) %>% 
  arrange(year) %>% 
  ungroup()

flowexceedSac_ind <- flowexceedSac.df$ind

maxflowSac.df <- d_Sac_sort %>% 
  group_by(FishID) %>% 
  dplyr::reframe(ind = unique(Maxflow.z), 
                 year = unique(year)) %>% 
  arrange(year) %>% 
  ungroup()

maxflowSac_ind <- maxflowSac.df$ind

maxflowDelta.df <- d_Sac_sort %>% 
  group_by(FishID) %>% 
  dplyr::reframe(ind = unique(MaxFPTflow.z), 
                 year = unique(year)) %>% 
  arrange(year) %>% 
  ungroup()

maxflowDelta_ind <- maxflowDelta.df$ind

# Set CJS variables for Sac mainstem travel time model -------------------------------------------------------------------
#relKM <- c(517.344,461.579,461.579,450.703,441.728) # 'BattleCk_CNFH_Rel','RBDD_Rel','Altube Island','Irvine_finch,'DeerCk_RST_Rel','MillCk_RST_Rel'
# Woodson bridge = 423.2790, Butte bridge = 343.19467, Sacramento = 150.0867
ReachKM <- c(45, 88, 170, 110)
ReachKM_Tis <- 78
ReachKM_Kni <- 127
ReachKM_ind <- data.frame(cbind(d_Sac_sort$dist_rlwoodson,
                                d_Sac_sort$dist_woodsonbutte,
                                d_Sac_sort$dist_buttesac,
                                d_Sac_sort$dist_sacdelta))
ReachKMrst <- dist_rstwoodson
  
#identify records with one or more detections after release and get their FishID. 
#Fish not seen after release provide no data for travel time
d_DHSac <- read.csv(here('data',"DetectionHistorySac.csv"),stringsAsFactors = F) %>%   #detection history file with times of date/time of detection
          arrange(FishID)

d_Sac_ord <- d_Sac_sort %>% 
            arrange(FishID)

# Check that Detection history and inp data have same fish IDs
missing_ids <- d_DHSac %>%
  anti_join(d_Sac_sort, by = "FishID") %>%
  pull(FishID)

release_times <- d_DHSac %>%
  filter(GEN == "Releasepoint") %>%
  select(FishID, release_time = min_time)

Woodson_times <- d_DHSac %>%
  filter(GEN == "WoodsonBridge") %>%
  select(FishID, woodson_time = min_time)

Butte_times <- d_DHSac %>%
  filter(GEN == "ButteBridge") %>%
  select(FishID, butte_time = min_time)

Sacramento_times <- d_DHSac %>%
  filter(GEN == "Sacramento") %>%
  select(FishID, sacramento_time = min_time)

Endpoint_times <- d_DHSac %>%
  filter(GEN == "Endpoint") %>%
  select(FishID, endpoint_time = min_time)

TTfR <- release_times %>%
  left_join(Butte_times, by = "FishID") %>%
  left_join(Woodson_times, by = "FishID") %>%
  left_join(Sacramento_times, by = "FishID") %>%
  left_join(Endpoint_times, by = "FishID") %>%
  mutate(TTfR1 = as.numeric(difftime(woodson_time, release_time, units = "days")),
         TTfR2 = as.numeric(difftime(butte_time, release_time, units = "days")),
         TTfR3 = as.numeric(difftime(sacramento_time, release_time, units = "days")),
         TTfR4 = as.numeric(difftime(endpoint_time, release_time, units = "days")),
         TTfR4 = if_else(TTfR4 < 0, NA_real_, TTfR4)) %>%# Replace negative TTfR4 with NA, example fish SP2023-810 
  arrange(FishID) %>% 
  select(TTfR1,TTfR2,TTfR3,TTfR4)

Nobs_df <- TTfR %>% 
  mutate(Obs = 4 - rowSums(is.na(across(everything()))))

Nobs <- sum(Nobs_df$Obs)

ObsTT <- TTfR %>% 
  rowwise() %>%
  mutate(row_vec = list(c_across(everything()))) %>%
  pull(row_vec) %>%
  unlist() %>% 
  na.omit()

vec <- as.vector(t(TTfR))                    # row-wise flattening
counter <- seq_len(sum(!is.na(vec)))       # counter for non-NA values
vec[!is.na(vec)] <- counter                # replace only non-NA

# Rebuild dataframe with same shape and column names
TTind <- matrix(vec, nrow = nrow(TTfR), byrow = TRUE) %>%
  as.data.frame() %>% 
  mutate(across(everything(), ~ replace_na(., 0)))

# Save travel time info in csv file
TTfinalSac_df <- data.frame(cbind(d_Sac_ord$FishID,d_Sac_ord$year,d_Sac_ord$fish_release_date, d_Sac_ord$StudyID,
                                  rgind,d_Sac_ord$fish_length, d_Sac_ord$fish_weight,d_Sac_ord$WY3,
                                  d_Sac_ord$ch,round(TTfR,digits=2))) %>%
  rename(c(FishID=d_Sac_ord.FishID,Year=d_Sac_ord.year,RelDate=d_Sac_ord.fish_release_date,
           StudyID=d_Sac_ord.StudyID,RelGp=rgind,FL=d_Sac_ord.fish_length,WGT=d_Sac_ord.fish_weight,
           WY = d_Sac_ord.WY3,CH=d_Sac_ord.ch))

write.csv(TTfinalSac_df ,here("data", "TravTime_Sac.csv"), row.names=FALSE)

summary_TT <-TTfinalSac_df %>%
  group_by(Year,StudyID) %>%
  summarise(
    TTfR1 = sum(!is.na(TTfR1)),
    TTfR2 = sum(!is.na(TTfR2)),
    TTfR3 = sum(!is.na(TTfR3)),
    TTfR4 = sum(!is.na(TTfR4))
  ) %>%
  ungroup()

summary_TT_SacDelta <-TTfinalSac_df %>%
                  group_by(FishID,StudyID) %>%
                  summarise(TTfR3_and_TTfR4_ind = sum(!is.na(TTfR3) & !is.na(TTfR4))) %>%
                  group_by(StudyID) %>%
                  summarise(TTfR3_and_TTfR4 = sum(TTfR3_and_TTfR4_ind))

summary_TT_final <- summary_TT %>%
                    left_join(summary_TT_SacDelta,by="StudyID")

# write.csv(summary_TT_final ,here("data", "TravTime_Sac_summary.csv"), row.names=FALSE)


# Set CJS variables for Feather/Butte survival model------------------------------------
trib_ind <- d_FeaBut_sort$trib_ind
site_T <- c("Butte","Feather")
NyrsT <- length(Year_FeaBut)
RelGpT <- unique(d_FeaBut_sort$StudyID)
NrgT <- length(RelGpT)
NindT <- dim(d_FeaBut_sort)[1]
Ntribs <- 2
NrstT <- 3
TribForRST <- c(1,2,2)
NreachesT <- 2
NdetlocsT <- 3  # of detection locations for tribs = release, sacramento and delta stations
firstCapT <- d_FeaBut_sort$firstCap
lastCapT <-d_FeaBut_sort$lastCap
WY2T <- d_FeaBut_sort$WY2
WY3T <- d_FeaBut_sort$WY3
MaxflowT.z <- d_FeaBut_sort$Maxflow.z
MaxflowT <- d_FeaBut_sort$Maxflow
MaxflowB <- d_FeaBut_sort$Maxflowbut
MaxflowF <- d_FeaBut_sort$Maxflowfea
MaxflowDeltaT.z <- d_FeaBut_sort$MaxFPTflow.z
MaxflowDeltaT <- d_FeaBut_sort$MaxFPTflow
FlowexceedT <- d_FeaBut_sort$FlowexceedT

# Upper Butte 2019 does not have weight information so use average values across all release group instead. 
d_FeaBut_sort <- d_FeaBut_sort %>% 
  mutate(across(fish_weight, ~ replace_na(., mean(., na.rm=TRUE))),
         across(fish_k, ~ replace_na(., mean(., na.rm=TRUE))))

FL_T <- d_FeaBut_sort$fish_length
WGT_T <- d_FeaBut_sort$fish_weight
CF_T <- d_FeaBut_sort$fish_k

CHT <- data.frame(cbind(as.integer(substr(d_FeaBut_sort$ch, 1, 1)),as.integer(substr(d_FeaBut_sort$ch,2, 2)),
                        as.integer(substr(d_FeaBut_sort$ch, 3, 3))))

RmultT <- d_FeaBut_sort$dist_rlsac.z # reach distances standardized per 100km
RmultTrib <- c(1.17,0.92) # distances for prediction model correspond to the average distance from all release locations for Butte and Feather respectively per 100km
dist_rstsac <- c(201.583633469999,89.5817927391999, 131.0731) # Distance from PPDD, Yuba Hallwood, and average distance of eye riffle, steep riffle, and gateway riffle to Delta entry
dist_rstsac.z =dist_rstsac/100 # standardize distances per 100km
RmultTrst <- dist_rstsac.z

yrindT <- vector(length=NindT)
rgindT <- yrindT
for(i in 1:NindT){
  yrindT[i] <- which(Year_all == d_FeaBut_sort$year[i])
  rgindT[i] <- which(RelGpT == d_FeaBut_sort$StudyID[i])
}

trib_rg <- vector(length=NrgT) #the tributary index for each release groups
for(irg in 1:NrgT){
  irecs <- which(rgindT == irg) #identify all records with the current release group index irg
  trib_rg[irg] <- d_FeaBut_sort$trib_ind[irecs[1]] #Get the tributary index. Only need first records as all individuals with same irg will be from same trib
}

indsitet <- vector(length=NindT)
for(i in 1:NindT){
  indsitet[i] <- which(site_T == d_FeaBut_sort$release_trib[i])
}

#### Add index for covariates
rgwy3_T.df <- d_FeaBut_sort %>% 
             group_by(StudyID) %>% 
              dplyr::reframe(ind = unique(WY3),
                      year = unique(year)) %>% 
              arrange(year) %>% 
              ungroup()

rgwy3_indT <- rgwy3_T.df$ind

rgwy2_T.df <- d_FeaBut_sort %>% 
  group_by(StudyID) %>% 
  dplyr::reframe(ind = unique(WY2),
                   year = unique(year)) %>% 
  arrange(year) %>% 
  ungroup()

rgwy2_indT <- rgwy2_T.df$ind

flowexceedT.df <- d_FeaBut_sort %>% 
  group_by(StudyID) %>% 
  dplyr::reframe(ind = unique(FlowexceedT), 
                 year = unique(year)) %>% 
  arrange(year) %>% 
  ungroup()

flowexceedT_ind <- flowexceedT.df$ind

maxflowT.df <- d_FeaBut_sort %>% 
  group_by(FishID) %>% 
  dplyr::reframe(ind = unique(Maxflow.z), 
                 year = unique(year),
                 month=unique(month)) %>% 
  arrange(year,month) %>% 
  ungroup()

maxflowT_ind <- maxflowT.df$ind

maxflowDeltaT.df <- d_FeaBut_sort %>% 
  group_by(FishID) %>% 
  dplyr::reframe(ind = unique(MaxFPTflow.z), 
                 year = unique(year)) %>% 
  arrange(year) %>% 
  ungroup()

maxflowDeltaT_ind <- maxflowDeltaT.df$ind

# Set CJS variables for Feather/Butte travel time model------------------------------------
#identify records with one or more detections after release and get their FishID. 
#Fish not seen after release provide no data for travel time
d_DHFea <- read.csv(here('data',"DetectionHistoryFea.csv"),stringsAsFactors = F) #detection history file with times of date/time of detection
d_DHBut <- read.csv(here('data',"DetectionHistoryBut.csv"),stringsAsFactors = F) #detection history file with times of date/time of detection
d_DHFeaBut <- data.frame(rbind(d_DHBut,d_DHFea)) %>% 
              arrange(FishID)

# relRKM <- 'UpperButte_RST_Rel' =340.854, 'Butte_Blw_Sanborn_Rel' =293.010, 'Laux Rd'=276.000 ,'North_Weir_Rel' =289.490
#'Sanborn_Slough_Rel' =288.650, 'SutterBypass_Weir2_RST_Rel' =249.541
#' 'FR_Gridley_Rel' = 287.387, 'FR_Boyds_Rel' = 240.755
#' # Sacramento = 150.0867

ReachKMT <- data.frame(rbind(cbind(117,110),cbind(92,110)))
ReachKMT_ind <- data.frame(cbind(d_FeaBut_sort$dist_rlsac,
                                 d_FeaBut_sort$dist_sacdelta))
ReachKMTrst <- dist_rstsac

release_timesT <- d_DHFeaBut %>%
  filter(GEN == "Releasepoint") %>%
  select(FishID, release_time = min_time)

Sacramento_timesT <- d_DHFeaBut %>%
  filter(GEN == "Sacramento") %>%
  select(FishID, sacramento_time = min_time)

Endpoint_timesT <- d_DHFeaBut %>%
  filter(GEN == "Endpoint") %>%
  select(FishID, endpoint_time = min_time)

TTfRT <- release_timesT %>%
  left_join(Sacramento_timesT, by = "FishID") %>%
  left_join(Endpoint_timesT, by = "FishID") %>%
  mutate(TTfR1 = as.numeric(difftime(sacramento_time, release_time, units = "days")),
         TTfR2 = as.numeric(difftime(endpoint_time, release_time, units = "days"))) %>% 
  arrange(FishID) %>% 
  select(TTfR1,TTfR2)

NobsT_df <- TTfRT %>% 
  mutate(Obs = 2 - rowSums(is.na(across(everything()))))
NobsT <- sum(NobsT_df$Obs)

ObsTTT <- TTfRT %>% 
  rowwise() %>%
  mutate(row_vec = list(c_across(everything()))) %>%
  pull(row_vec) %>%
  unlist() %>% 
  na.omit()

vecT <- as.vector(t(TTfRT))                    # row-wise flattening
counterT <- seq_len(sum(!is.na(vecT)))       # counter for non-NA values
vecT[!is.na(vecT)] <- counterT                # replace only non-NA

# Rebuild dataframe with same shape and column names
TTindT <- matrix(vecT, nrow = nrow(TTfRT), byrow = TRUE) %>%
  as.data.frame() %>% 
  mutate(across(everything(), ~ replace_na(., 0)))

# Save travel time info in csv file
d_FeaBut_ord <- d_FeaBut_sort %>%
  arrange(FishID)

TTfinalFeaBut_df <- data.frame(cbind(d_FeaBut_ord$FishID,d_FeaBut_ord$year,d_FeaBut_ord$fish_release_date,
                                     d_FeaBut_ord$StudyID,rgindT,d_FeaBut_ord$fish_length, d_FeaBut_ord$fish_weight,
                                     d_FeaBut_ord$WY3,d_FeaBut_ord$ch,round(TTfRT,digits=2))) %>%
  rename(c(FishID=d_FeaBut_ord.FishID,Year=d_FeaBut_ord.year,RelDate=d_FeaBut_ord.fish_release_date,
           StudyID=d_FeaBut_ord.StudyID,RelGp=rgindT,FL=d_FeaBut_ord.fish_length,WGT=d_FeaBut_ord.fish_weight,
           WY = d_FeaBut_ord.WY3,CH=d_FeaBut_ord.ch))

write.csv(TTfinalFeaBut_df ,here("data", "TravTime_FeaBut.csv"), row.names=FALSE)

summary_TTT <-TTfinalFeaBut_df %>%
  group_by(Year,StudyID) %>%
  summarise(
    TTfR1 = sum(!is.na(TTfR1)),
    TTfR2 = sum(!is.na(TTfR2))
  ) %>%
  ungroup()

summary_TTT_SacDelta <-TTfinalFeaBut_df %>%
  group_by(FishID,StudyID) %>%
  summarise(TTfR1_and_TTfR2_ind = sum(!is.na(TTfR1) & !is.na(TTfR2))) %>%
  group_by(StudyID) %>%
  summarise(TTfR1_and_TTfR2 = sum(TTfR1_and_TTfR2_ind))

summary_TTT_final <- summary_TTT %>%
  left_join(summary_TTT_SacDelta,by="StudyID")


#write.csv(summary_TTT_final ,here("data", "TravTime_FeaBut_summary.csv"), row.names=FALSE)

# Summary stats -----------------------------------------------------------------------------
d_Sac_summary <- d_Sac_sort %>% 
  group_by(StudyID) %>% 
  rename(studyid = StudyID) %>% 
  dplyr::reframe(Mean_FL = mean(fish_length, na.rm=TRUE),
                   Mean_Wt = mean(fish_weight, na.rm=TRUE),
                   Mean_k = mean(fish_k, na.rm=TRUE),
                   Min_FL = min(fish_length, na.rm=TRUE),
                   Min_Wt = min(fish_weight, na.rm=TRUE),
                   Max_FL = max(fish_length, na.rm=TRUE),
                   Max_Wt = max(fish_weight, na.rm=TRUE),
                   CV_FL = sd(fish_length, na.rm=TRUE) / mean(fish_length, na.rm=TRUE),
                   CV_Wt = sd(fish_weight, na.rm=TRUE) / mean(fish_weight, na.rm=TRUE),
                   CV_k = sd(fish_k, na.rm=TRUE) / mean(fish_k, na.rm=TRUE),
                   ntot = n(),
                   ndetect_Woodson = sum(as.numeric(substr(as.character(ch), 2, 2))),
                   ndetect_Butte = sum(as.numeric(substr(as.character(ch), 3, 3))),
                   ndetect_Sac = sum(as.numeric(substr(as.character(ch), 4, 4))),
                   ndetect_Delta = sum(as.numeric(substr(as.character(ch), 5, 5)))) %>% 
  mutate(watershed ="Sac")

d_FeaBut_summary <- d_FeaBut_sort %>% 
  group_by(StudyID) %>% 
  rename(studyid = StudyID) %>% 
  dplyr::reframe(Mean_FL = mean(fish_length, na.rm=TRUE),
                   Mean_Wt = mean(fish_weight, na.rm=TRUE),
                   Mean_k = mean(fish_k, na.rm=TRUE),
                   Min_FL = min(fish_length, na.rm=TRUE),
                   Min_Wt = min(fish_weight, na.rm=TRUE),
                   Max_FL = max(fish_length, na.rm=TRUE),
                   Max_Wt = max(fish_weight, na.rm=TRUE),
                   CV_FL = sd(fish_length, na.rm=TRUE) / mean(fish_length, na.rm=TRUE),
                   CV_Wt = sd(fish_weight, na.rm=TRUE) / mean(fish_weight, na.rm=TRUE),
                   CV_k = sd(fish_k, na.rm=TRUE) / mean(fish_k, na.rm=TRUE),
                   ntot = n(),
                   ndetect_Sac = sum(as.numeric(substr(as.character(ch), 2, 2))),
                   ndetect_Delta = sum(as.numeric(substr(as.character(ch), 3, 3)))) %>% 
  mutate(watershed ="FeaBut")

write.csv(d_Sac_summary,here("data", "Sac_summary.csv"), row.names=FALSE)
write.csv(d_FeaBut_summary,here("data", "FeaBut_summary.csv"), row.names=FALSE)

# Summary figures -------------------------------------------------------------
### mainstem vs trib fish size 
d_SacFeaBut_summary <- data.frame(rbind(d_Sac_summary[,c(1:12,17)],d_FeaBut_summary[,c(1:12,15)]))

(fishsizedens <- ggplot()+
    geom_density(data=d_SacFeaBut_summary, aes(x=Mean_FL, fill = watershed,color=watershed,linetype = watershed),
                 alpha=0.4,adjust=1.3, size=1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    scale_color_manual(values = c("#000000", "#E69F00"))+
    scale_fill_manual(values = c("#000000", "#E69F00"))+
    scale_linetype_manual(values = c( "dashed","solid")) +
    scale_x_continuous(limits=c(70,140), breaks=seq(70,140,20))
)

 #ggsave('results/figures/fishsizedens.jpg',plot=fishsizedens, dpi = 350, height = 4, width = 5)
 
 # Load CWT data frame ------------------------------------------------
drerelrec <- read.csv(here("data", "hatchery_relrec.csv"),  stringsAsFactors = F)

