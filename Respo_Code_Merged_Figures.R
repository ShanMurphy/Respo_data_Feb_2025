# Load libraries
library(effects)
library(lme4)
library(emmeans)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(reshape2)
library(data.table)
library(rstatix)
library(car)
library(cowplot)
library(ggsignif)

# Read in the data called respo_data
respo_data <- read.csv("/Users/shannonmurphy/Desktop/Respirometry_modified_data_Jan_15_2025.csv", stringsAsFactors=TRUE)

# Triple-checking there are no NAs in the data 
respo_data <- respo_data[!is.na(respo_data$NPP_SA_Corrected_cm2_hr) & 
                           !is.na(respo_data$GPP_SA_Corrected_cm2_hr) &
                           !is.na(respo_data$R_SA_Corrected_cm2_hr) &
                           !is.na(respo_data$Nutrient_Treatment) & 
                           !is.na(respo_data$Time_Point) & 
                           !is.na(respo_data$Temp_Treatment) & 
                           !is.na(respo_data$Colony_ID), ]

# turning respo_data into a data frame
respo_data <- data.frame(respo_data)
respo_data$R_SA_Corrected_cm2_hr <- abs(respo_data$R_SA_Corrected_cm2_hr)

# Making sure there are only four nutrient treatments 
respo_data$Nutrient_Treatment <- trimws(respo_data$Nutrient_Treatment)

# visually inspecting the dataframe 
head(respo_data)

# Changing the corresponding columns to factors
respo_data$Nutrient_Treatment <- as.factor(respo_data$Nutrient_Treatment)
respo_data$Colony_ID <- as.factor(respo_data$Colony_ID)
respo_data$Run_Number <- as.factor(respo_data$Run_Number)
respo_data$Channel <- as.factor(respo_data$Channel)
respo_data$Sample_Day <- as.factor(respo_data$Sample_Day)

# Double-checking the levels of each column
levels(respo_data$Nutrient_Treatment)
levels(respo_data$Temp_Treatment)
levels(respo_data$Time_Point)
levels(respo_data$Colony_ID)
levels(respo_data$Sample_Day)

#################################################################################
# Calculating and removing outliers
# MCAP - T0 - NPP
filtered_data_MCAP_T0_NPP  <- respo_data %>%
  filter(Species == "MCAP", Time_Point == "T0") %>%
  mutate(
    Q1 = quantile(NPP_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(NPP_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = NPP_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      NPP_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

# MCAP - T0 - GPP
filtered_data_MCAP_T0_GPP  <- respo_data %>%
  filter(Species == "MCAP", Time_Point == "T0") %>%
  mutate(
    Q1 = quantile(GPP_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(GPP_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = GPP_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      GPP_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

# MCAP - T0 - R
filtered_data_MCAP_T0_R  <- respo_data %>%
  filter(Species == "MCAP", Time_Point == "T0") %>%
  mutate(
    Q1 = quantile(R_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(R_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = R_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      R_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 
#--------------------------------------------------------------------------------
# MCAP - T1 - NPP 
filtered_data_MCAP_T1_NPP  <- respo_data %>%
  filter(Species == "MCAP", Time_Point == "T1") %>%
  mutate(
    Q1 = quantile(NPP_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(NPP_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = NPP_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      NPP_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

# MCAP - T1 - GPP
filtered_data_MCAP_T1_GPP  <- respo_data %>%
  filter(Species == "MCAP", Time_Point == "T1") %>%
  mutate(
    Q1 = quantile(GPP_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(GPP_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = GPP_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      GPP_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

# MCAP - T1 - R
filtered_data_MCAP_T1_R  <- respo_data %>%
  filter(Species == "MCAP", Time_Point == "T1") %>%
  mutate(
    Q1 = quantile(R_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(R_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = R_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      R_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

#--------------------------------------------------------------------------------
# MCAP - T2 - NPP 
filtered_data_MCAP_T2_NPP  <- respo_data %>%
  filter(Species == "MCAP", Time_Point == "T2") %>%
  mutate(
    Q1 = quantile(NPP_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(NPP_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = NPP_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      NPP_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

# MCAP - T2 - GPP
filtered_data_MCAP_T2_GPP  <- respo_data %>%
  filter(Species == "MCAP", Time_Point == "T2") %>%
  mutate(
    Q1 = quantile(GPP_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(GPP_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = GPP_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      GPP_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

# MCAP - T2 - R
filtered_data_MCAP_T2_R <- respo_data %>%
  filter(Species == "MCAP", Time_Point == "T2") %>%
  mutate(
    Q1 = quantile(R_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(R_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = R_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      R_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

table(filtered_data_MCAP_T2_R$is_outlier)
#--------------------------------------------------------------------------------
# MCAP - T3 - NPP 
filtered_data_MCAP_T3_NPP <- respo_data%>%
  filter(Species == "MCAP", Time_Point == "T3") %>%
  mutate(
    Q1 = quantile(NPP_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(NPP_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = NPP_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) | 
      NPP_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier)  

# MCAP - T3 - GPP
filtered_data_MCAP_T3_GPP  <- respo_data %>%
  filter(Species == "MCAP", Time_Point == "T3") %>%
  mutate(
    Q1 = quantile(GPP_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(GPP_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = GPP_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      GPP_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

# MCAP - T3 - R
filtered_data_MCAP_T3_R  <- respo_data %>%
  filter(Species == "MCAP", Time_Point == "T3") %>%
  mutate(
    Q1 = quantile(R_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(R_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = R_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      R_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

#################################################################################
# Calculating and removing outliers
# PCOMP - T0 - NPP
filtered_data_PCOMP_T0_NPP  <- respo_data %>%
  filter(Species == "PCOMP", Time_Point == "T0") %>%
  mutate(
    Q1 = quantile(NPP_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(NPP_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = NPP_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      NPP_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

# PCOMP - T0 - GPP
filtered_data_PCOMP_T0_GPP  <- respo_data %>%
  filter(Species == "PCOMP", Time_Point == "T0") %>%
  mutate(
    Q1 = quantile(GPP_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(GPP_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = GPP_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      GPP_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

# PCOMP - T0 - R
filtered_data_PCOMP_T0_R  <- respo_data %>%
  filter(Species == "PCOMP", Time_Point == "T0") %>%
  mutate(
    Q1 = quantile(R_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(R_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = R_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      R_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 
#--------------------------------------------------------------------------------
# PCOMP - T1 - NPP 
filtered_data_PCOMP_T1_NPP  <- respo_data %>%
  filter(Species == "PCOMP", Time_Point == "T1") %>%
  mutate(
    Q1 = quantile(NPP_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(NPP_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = NPP_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      NPP_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

# PCOMP - T1 - GPP
filtered_data_PCOMP_T1_GPP  <- respo_data %>%
  filter(Species == "PCOMP", Time_Point == "T1") %>%
  mutate(
    Q1 = quantile(GPP_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(GPP_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = GPP_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      GPP_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

# PCOMP - T1 - R
filtered_data_PCOMP_T1_R  <- respo_data %>%
  filter(Species == "PCOMP", Time_Point == "T1") %>%
  mutate(
    Q1 = quantile(R_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(R_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = R_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      R_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

#--------------------------------------------------------------------------------
# PCOMP - T2 - NPP 
filtered_data_PCOMP_T2_NPP  <- respo_data %>%
  filter(Species == "PCOMP", Time_Point == "T2") %>%
  mutate(
    Q1 = quantile(NPP_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(NPP_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = NPP_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      NPP_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

# PCOMP - T2 - GPP
filtered_data_PCOMP_T2_GPP  <- respo_data %>%
  filter(Species == "PCOMP", Time_Point == "T2") %>%
  mutate(
    Q1 = quantile(GPP_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(GPP_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = GPP_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      GPP_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

# PCOMP - T2 - R
filtered_data_PCOMP_T2_R  <- respo_data %>%
  filter(Species == "PCOMP", Time_Point == "T2") %>%
  mutate(
    Q1 = quantile(R_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(R_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = R_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      R_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 
#--------------------------------------------------------------------------------
# PCOMP - T3 - NPP 
filtered_data_PCOMP_T3_NPP  <- respo_data %>%
  filter(Species == "PCOMP", Time_Point == "T3") %>%
  mutate(
    Q1 = quantile(NPP_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(NPP_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = NPP_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      NPP_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

# PCOMP - T3 - GPP
filtered_data_PCOMP_T3_GPP  <- respo_data %>%
  filter(Species == "PCOMP", Time_Point == "T3") %>%
  mutate(
    Q1 = quantile(GPP_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(GPP_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = GPP_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      GPP_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

# PCOMP - T3 - R
filtered_data_PCOMP_T3_R  <- respo_data %>%
  filter(Species == "PCOMP", Time_Point == "T3") %>%
  mutate(
    Q1 = quantile(R_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(R_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = R_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      R_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

##################################################################################
# count how many data points are in each nutrient treatment
# MCAP - T0 
counts_MCAP_T0_NPP <- filtered_data_MCAP_T0_NPP %>% 
  group_by(Nutrient_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_MCAP_T0_NPP

counts_MCAP_T0_GPP <- filtered_data_MCAP_T0_GPP %>% 
  group_by(Nutrient_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_MCAP_T0_GPP

counts_MCAP_T0_R <- filtered_data_MCAP_T0_R %>% 
  group_by(Nutrient_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_MCAP_T0_R

# MCAP - T1
counts_MCAP_T1_NPP <- filtered_data_MCAP_T1_NPP %>% 
  group_by(Nutrient_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_MCAP_T1_NPP

counts_MCAP_T1_GPP <- filtered_data_MCAP_T1_GPP %>% 
  group_by(Nutrient_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_MCAP_T1_GPP

counts_MCAP_T1_R <- filtered_data_MCAP_T1_R %>% 
  group_by(Nutrient_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_MCAP_T1_R

# MCAP - T2
counts_MCAP_T2_NPP <- filtered_data_MCAP_T2_NPP %>% 
  group_by(Nutrient_Treatment, Temp_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_MCAP_T2_NPP

counts_MCAP_T2_GPP <- filtered_data_MCAP_T2_GPP %>% 
  group_by(Nutrient_Treatment, Temp_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_MCAP_T2_GPP

counts_MCAP_T2_R <- filtered_data_MCAP_T2_R %>% 
  group_by(Nutrient_Treatment, Temp_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_MCAP_T2_R

# MCAP - T3
counts_MCAP_T3_NPP <- filtered_data_MCAP_T3_NPP %>%
  group_by(Nutrient_Treatment, Temp_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_MCAP_T3_NPP

counts_MCAP_T3_GPP <- filtered_data_MCAP_T3_GPP %>% 
  group_by(Nutrient_Treatment, Temp_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_MCAP_T3_GPP

counts_MCAP_T3_R <- filtered_data_MCAP_T3_R %>% 
  group_by(Nutrient_Treatment, Temp_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_MCAP_T3_R
#--------------------------------------------------------------------------------
# count how many data points are in each nutrient treatment
# PCOMP - T0 
counts_PCOMP_T0_NPP <- filtered_data_PCOMP_T0_NPP %>% 
  group_by(Nutrient_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_PCOMP_T0_NPP

counts_PCOMP_T0_GPP <- filtered_data_PCOMP_T0_GPP %>% 
  group_by(Nutrient_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_PCOMP_T0_GPP

counts_PCOMP_T0_R <- filtered_data_PCOMP_T0_R %>% 
  group_by(Nutrient_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_PCOMP_T0_R

# PCOMP - T1
counts_PCOMP_T1_NPP <- filtered_data_PCOMP_T1_NPP %>% 
  group_by(Nutrient_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_PCOMP_T1_NPP

counts_PCOMP_T1_GPP <- filtered_data_PCOMP_T1_GPP %>% 
  group_by(Nutrient_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_PCOMP_T1_GPP

counts_PCOMP_T1_R <- filtered_data_PCOMP_T1_R %>% 
  group_by(Nutrient_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_PCOMP_T1_R

# PCOMP - T2
counts_PCOMP_T2_NPP <- filtered_data_PCOMP_T2_NPP %>% 
  group_by(Nutrient_Treatment, Temp_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_PCOMP_T2_NPP

counts_PCOMP_T2_GPP <- filtered_data_PCOMP_T2_GPP %>% 
  group_by(Nutrient_Treatment, Temp_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_PCOMP_T2_GPP

counts_PCOMP_T2_R <- filtered_data_PCOMP_T2_R %>% 
  group_by(Nutrient_Treatment, Temp_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_PCOMP_T2_R

# PCOMP - T3
counts_PCOMP_T3_NPP <- filtered_data_PCOMP_T3_NPP %>% 
  group_by(Nutrient_Treatment, Temp_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_PCOMP_T3_NPP

counts_PCOMP_T3_GPP <- filtered_data_PCOMP_T3_GPP %>% 
  group_by(Nutrient_Treatment, Temp_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_PCOMP_T3_GPP

counts_PCOMP_T3_R <- filtered_data_PCOMP_T3_R %>% 
  group_by(Nutrient_Treatment, Temp_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts_PCOMP_T3_R

##################################################################################
# T0 Plots 
T0_M_plot_NPP <- ggplot(data = filtered_data_MCAP_T0_NPP, 
                    aes(x = Nutrient_Treatment, 
                        y = NPP_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_MCAP_T0_NPP, aes(x = Nutrient_Treatment, y = 0.15, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) + 
  labs(title = expression(italic("Montipora capitata") ~ "Net Primary Production at T0"),
       x = "Nutrient Treatment",
       y = expression("NPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T0_M_plot_NPP

# Adding significant lines over boxplots 
T0_M_plot_NPP <- T0_M_plot_NPP +
  geom_signif(
    comparisons = list(c("Control", "Effluent"), 
                       c("Effluent", "Guano"), 
                       c("Effluent", "Inorganic")),
    annotations = c("*", "*", "*"),  
    y_position = c(1.6, 1.8, 2.0),  
    tip_length = 0.02,
    size = 0.5,
    textsize = 7,
    color = "black"
  )
T0_M_plot_NPP

T0_M_plot_GPP <- ggplot(data = filtered_data_MCAP_T0_GPP, 
                    aes(x = Nutrient_Treatment, 
                        y = GPP_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_MCAP_T0_GPP, aes(x = Nutrient_Treatment, y = 0.15, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) + 
  labs(title = expression(italic("Montipora capitata") ~ "Gross Primary Production at T0"),
       x = "Nutrient Treatment",
       y = expression("GPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T0_M_plot_GPP

# Adding significant lines over boxplots 
T0_M_plot_GPP <- T0_M_plot_GPP +
  geom_signif(
    comparisons = list(c("Control", "Effluent"), 
                       c("Effluent", "Guano"), 
                       c("Effluent", "Inorganic")),
    annotations = c("*", "*", "*"),  
    y_position = c(2.0, 2.2, 2.4),  
    tip_length = 0.02,
    size = 0.5,
    textsize = 7,
    color = "black")

T0_M_plot_GPP

T0_M_plot_R <- ggplot(data = filtered_data_MCAP_T0_R, 
                    aes(x = Nutrient_Treatment, 
                        y = R_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_MCAP_T0_R, aes(x = Nutrient_Treatment, y = 0.15, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) + 
  labs(title = expression(italic("Montipora capitata") ~ "Net Primary Production at T0"),
       x = "Nutrient Treatment",
       y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T0_M_plot_R

# Adding significant lines over boxplots 
T0_M_plot_R <- T0_M_plot_R +
  geom_signif(
    comparisons = list(c("Control", "Effluent"), 
                       c("Effluent", "Guano"), 
                       c("Effluent", "Inorganic")),
    annotations = c("*", "*", "**"),  
    y_position = c(0.8, 1.0, 1.2),  
    tip_length = 0.02,
    size = 0.5,
    textsize = 7,
    color = "black")

T0_M_plot_R

##################################################################################
# T1 plots
T1_M_plot_NPP <- ggplot(data= filtered_data_MCAP_T1_NPP, 
                    aes(x = Nutrient_Treatment, 
                        y = NPP_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_MCAP_T1_NPP, aes(x = Nutrient_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  labs(title = expression(italic("Montipora capitata") ~ "Net Primary Production at T1"),
       x = "Nutrient Treatment",
       y = expression("NPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T1_M_plot_NPP

# Adding significant lines over boxplots 
T1_M_plot_NPP <- T1_M_plot_NPP +
  geom_signif(
    comparisons = list(c("Control", "Guano"), 
                       c("Guano", "Inorganic")),
    annotations = c("**", "*"),  
    y_position = c(1.5, 1.7),  # Adjust as needed for better placement
    tip_length = 0.02,
    size = 0.5,
    textsize = 7,
    color = "black"
  )
T1_M_plot_NPP

T1_M_plot_GPP <- ggplot(data= filtered_data_MCAP_T1_GPP, 
                    aes(x = Nutrient_Treatment, 
                        y = GPP_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_MCAP_T1_GPP, aes(x = Nutrient_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  labs(title = expression(italic("Montipora capitata") ~ "Gross Primary Production at T1"),
       x = "Nutrient Treatment",
       y = expression("GPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T1_M_plot_GPP

# Adding significant lines over boxplots 
T1_M_plot_GPP <- T1_M_plot_GPP +
  geom_signif(
    comparisons = list(c("Control", "Guano")),
    annotations = c("*"),  
    y_position = c(2.1),  # Adjust as needed for better placement
    tip_length = 0.02,
    size = 0.5,
    textsize = 7,
    color = "black")
T1_M_plot_GPP

T1_M_plot_R <- ggplot(data= filtered_data_MCAP_T1_R, 
                    aes(x = Nutrient_Treatment, 
                        y = R_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_MCAP_T1_R, aes(x = Nutrient_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  labs(title = expression(italic("Montipora capitata") ~ "Net Primary Production at T1"),
       x = "Nutrient Treatment",
       y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T1_M_plot_R

##################################################################################
# T2 plots 
T2_M_plot_NPP <- ggplot(filtered_data_MCAP_T2_NPP, 
                    aes(x = Temp_Treatment, 
                        y = NPP_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  facet_wrap(~ Nutrient_Treatment, ncol = 4) +
  geom_text(data = counts_MCAP_T2_NPP, aes(x = Temp_Treatment, y = -0.15, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  labs(title = expression(italic("Montipora capitata") ~ "Net Primary Production at T2"),
       x = "Temperature Treatment",
       y = expression("NPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T2_M_plot_NPP

# Add notation and adjust slightly to avoid overlap
annotations_M_T2_NPP <- data.frame(
  Nutrient_Treatment = c("Control", "Inorganic"),
  Temp_Treatment = "Ambient",
  y_position = c(1.4, 1.4),  # Slightly adjust y_position to prevent overlap
  label = c("***", "**")  # For simplicity, we just use the * label here
)

# Adding significant lines over boxplots 
T2_M_plot_NPP <- T2_M_plot_NPP +
  geom_text(
    data = annotations_M_T2_NPP,  
    aes(x = 1.5, y = y_position + 0.1, label = label), 
    size = 7, color = "black") +
  geom_segment(
    data = annotations_M_T2_NPP,  
    aes(x = 1.0, xend = 2.0, y = y_position, yend = y_position),
    color = "black", size = 0.5, lineend = "round") +
  geom_segment(
    data = annotations_M_T2_NPP,  
    aes(x = 1.0, xend = 1.0, y = y_position - 0.02, yend = y_position + 0),
    color = "black", size = 0.5) +
  geom_segment(
    data = annotations_M_T2_NPP,  
    aes(x = 2.0, xend = 2.0, y = y_position - 0.02, yend = y_position + 0),
    color = "black", size = 0.5)

T2_M_plot_NPP

T2_M_plot_GPP <- ggplot(filtered_data_MCAP_T2_GPP, 
                    aes(x = Temp_Treatment, 
                        y = GPP_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  facet_wrap(~ Nutrient_Treatment, ncol = 4) +
  geom_text(data = counts_MCAP_T2_GPP, aes(x = Temp_Treatment, y = -0.15, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  labs(title = expression(italic("Montipora capitata") ~ "Gross Primary Production at T2"),
       x = "Temperature Treatment",
       y = expression("GPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T2_M_plot_GPP

T2_M_plot_R <- ggplot(filtered_data_MCAP_T2_R, 
                    aes(x = Temp_Treatment, 
                        y = R_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  facet_wrap(~ Nutrient_Treatment, ncol = 4) +
  geom_text(data = counts_MCAP_T2_R, aes(x = Temp_Treatment, y = -0.15, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  labs(title = expression(italic("Montipora capitata") ~ "Oxygen Consumption at T2"),
       x = "Temperature Treatment",
       y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T2_M_plot_R

##################################################################################
# T3 plots 
T3_M_plot_NPP <- ggplot(data=filtered_data_MCAP_T3_NPP, 
                    aes(x =Temp_Treatment, 
                        y = NPP_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_MCAP_T3_NPP, aes(x = Temp_Treatment, y = -0.15, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  facet_wrap(~ Nutrient_Treatment, ncol=4) +
  labs(title = expression(italic("Montipora capitata") ~ "Net Primary Production at T3"),
       x = "Temperature Treatment",
       y = expression("NPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14)) 

T3_M_plot_NPP

T3_M_plot_GPP <- ggplot(data=filtered_data_MCAP_T3_GPP, 
                    aes(x =Temp_Treatment, 
                        y = GPP_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_MCAP_T3_GPP, aes(x = Temp_Treatment, y = -0.15, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  facet_wrap(~ Nutrient_Treatment, ncol=4) +
  labs(title = expression(italic("Montipora capitata") ~ "Gross Primary Production at T3"),
       x = "Temperature Treatment",
       y = expression("GPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T3_M_plot_GPP

T3_M_plot_R <- ggplot(data=filtered_data_MCAP_T3_R, 
                    aes(x =Temp_Treatment, 
                        y = R_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_MCAP_T3_R, aes(x = Temp_Treatment, y = -0.15, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  facet_wrap(~ Nutrient_Treatment, ncol=4) +
  labs(title = expression(italic("Montipora capitata") ~ "Oxygen Consumption at T3"),
       x = "Temperature Treatment",
       y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T3_M_plot_R

#################################################################################
#################################################################################

T0_P_plot_NPP <- ggplot(data= filtered_data_PCOMP_T0_NPP, 
                    aes(x = Nutrient_Treatment, 
                        y = NPP_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_PCOMP_T0_NPP, aes(x = Nutrient_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) + 
  labs(title = expression(italic("Porites compressa") ~ "Net Primary Production at T0"),
       x = "Nutrient Treatment",
       y = expression("NPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T0_P_plot_NPP

T0_P_plot_GPP <- ggplot(data= filtered_data_PCOMP_T0_GPP, 
                    aes(x = Nutrient_Treatment, 
                        y = GPP_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_PCOMP_T0_GPP, aes(x = Nutrient_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) + 
  labs(title = expression(italic("Porites compressa") ~ "Gross Primary Production at T0"),
       x = "Nutrient Treatment",
       y = expression("GPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T0_P_plot_GPP

T0_P_plot_R <- ggplot(data= filtered_data_PCOMP_T0_R, 
                    aes(x = Nutrient_Treatment, 
                        y = R_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_PCOMP_T0_R, aes(x = Nutrient_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) + 
  labs(title = expression(italic("Porites compressa") ~ "Oxygen Consumption at T0"),
       x = "Nutrient Treatment",
       y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T0_P_plot_R

T1_P_plot_NPP <- ggplot(data= filtered_data_PCOMP_T1_NPP, 
                    aes(x = Nutrient_Treatment, 
                        y = NPP_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_PCOMP_T1_NPP, aes(x = Nutrient_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  labs(title = expression(italic("Porites compressa") ~ "Net Primary Production at T1"),
       x = "Nutrient Treatment",
       y = expression("NPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T1_P_plot_NPP

T1_P_plot_GPP <- ggplot(data= filtered_data_PCOMP_T1_GPP, 
                    aes(x = Nutrient_Treatment, 
                        y = GPP_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_PCOMP_T1_GPP, aes(x = Nutrient_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  labs(title = expression(italic("Porites compressa") ~ "Gross Primary Production at T1"),
       x = "Nutrient Treatment",
       y = expression("GPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T1_P_plot_GPP

T1_P_plot_R <- ggplot(data= filtered_data_PCOMP_T1_R, 
                    aes(x = Nutrient_Treatment, 
                        y = R_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_PCOMP_T1_R, aes(x = Nutrient_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  labs(title = expression(italic("Porites compressa") ~ "Oxygen Consumption at T1"),
       x = "Nutrient Treatment",
       y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T1_P_plot_R

T2_P_plot_NPP <- ggplot(filtered_data_PCOMP_T2_NPP, 
                    aes(x = Temp_Treatment, 
                        y = NPP_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_PCOMP_T2_NPP, aes(x = Temp_Treatment, y = -0.25, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  facet_wrap(~ Nutrient_Treatment, ncol=4) +
  labs(title = expression(italic("Porites compressa") ~ "Net Primary Production at T2"),
       x = "Temperature Treatment",
       y = expression("NPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T2_P_plot_NPP

# Adding significant lines over boxplots 
annotations_P_T2_NPP <- data.frame(
  Nutrient_Treatment = c("Control", "Effluent", "Guano", "Inorganic"),
  Temp_Treatment = "Ambient",
  y_position = c(3.2, 3.2, 2.4, 2.4),  
  label = c("***", "***", "***", "***"))

T2_P_plot_NPP <- T2_P_plot_NPP +
  geom_text(data = annotations_P_T2_NPP, 
            aes(x = 1.5, y = y_position + 0.1, label = label), 
            size = 7, color = "black") +
  
  geom_segment(data = annotations_P_T2_NPP,
               aes(x = 1.0, xend = 2.0, y = y_position, yend = y_position),
               color = "black", size = 0.5, lineend = "round") +
  geom_segment(data = annotations_P_T2_NPP,
               aes(x = 1.0, xend = 1.0, y = y_position - 0.02, yend = y_position + 0),
               color = "black", size = 0.5) +
  geom_segment(data = annotations_P_T2_NPP,
               aes(x = 2.0, xend = 2.0, y = y_position - 0.02, yend = y_position + 0),
               color = "black", size = 0.5)

T2_P_plot_NPP

T2_P_plot_GPP <- ggplot(filtered_data_PCOMP_T2_GPP, 
                    aes(x = Temp_Treatment, 
                        y = GPP_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_PCOMP_T2_GPP, aes(x = Temp_Treatment, y = -0.25, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  facet_wrap(~ Nutrient_Treatment, ncol=4) +
  labs(title = expression(italic("Porites compressa") ~ "Gross Primary Production at T2"),
       x = "Temperature Treatment",
       y = expression("GPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T2_P_plot_GPP

# Adding significant lines over boxplots 
annotations_P_T2_GPP <- data.frame(
  Nutrient_Treatment = c("Control", "Effluent", "Guano", "Inorganic"),
  Temp_Treatment = "Ambient",
  y_position = c(4.4, 4.4, 3.4, 3.4), 
  label = c("***", "***", "***", "***"))

T2_P_plot_GPP <- T2_P_plot_GPP +
  geom_text(data = annotations_P_T2_GPP, 
            aes(x = 1.5, y = y_position + 0.1, label = label), 
            size = 7, color = "black") +
  
  geom_segment(data = annotations_P_T2_GPP,
               aes(x = 1.0, xend = 2.0, y = y_position, yend = y_position),
               color = "black", size = 0.5, lineend = "round") +
  geom_segment(data = annotations_P_T2_GPP,
               aes(x = 1.0, xend = 1.0, y = y_position - 0.02, yend = y_position + 0),
               color = "black", size = 0.5) +
  geom_segment(data = annotations_P_T2_GPP,
               aes(x = 2.0, xend = 2.0, y = y_position - 0.02, yend = y_position + 0),
               color = "black", size = 0.5)

T2_P_plot_GPP

T2_P_plot_R <- ggplot(filtered_data_PCOMP_T2_R, 
                    aes(x = Temp_Treatment, 
                        y = R_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_PCOMP_T2_R, aes(x = Temp_Treatment, y = -0.25, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  facet_wrap(~ Nutrient_Treatment, ncol=4) +
  labs(title = expression(italic("Porites compressa") ~ "Oxygen Consumption at T2"),
       x = "Temperature Treatment",
       y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T2_P_plot_R

T3_P_plot_NPP <- ggplot(data=filtered_data_PCOMP_T3_NPP, 
                    aes(x = Temp_Treatment, 
                        y = NPP_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_PCOMP_T3_NPP, aes(x = Temp_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  facet_wrap(~ Nutrient_Treatment, ncol=4) +
  labs(title = expression(italic("Porites compressa") ~ "Net Primary Production at T3"),
       x = "Temperature Treatment",
       y = expression("NPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T3_P_plot_NPP

T3_P_plot_GPP <- ggplot(data=filtered_data_PCOMP_T3_GPP, 
                    aes(x = Temp_Treatment, 
                        y = GPP_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_PCOMP_T3_GPP, aes(x = Temp_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  facet_wrap(~ Nutrient_Treatment, ncol=4) +
  labs(title = expression(italic("Porites compressa") ~ "Gross Primary Production at T3"),
       x = "Temperature Treatment",
       y = expression("GPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T3_P_plot_GPP

T3_P_plot_R <- ggplot(data=filtered_data_PCOMP_T3_R, 
                    aes(x = Temp_Treatment, 
                        y = R_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts_PCOMP_T3_R, aes(x = Temp_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  facet_wrap(~ Nutrient_Treatment, ncol=4) +
  labs(title = expression(italic("Porites compressa") ~ "Oxygen Consumption at T3"),
       x = "Temperature Treatment",
       y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T3_P_plot_R

#################################################################################
T0_M_plot_modified_NPP <- T0_M_plot_NPP + 
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(-0, 2.5, by = 0.25)) + 
  labs(subtitle = expression(bold("A) Net Primary Production")), y = expression("NPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  ggtitle(expression(italic("Montipora capitata") ~ "T0 Metabolic Rates")) +  
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", size = 16))

T0_M_plot_modified_GPP <- T0_M_plot_GPP + 
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, by = 0.25)) + 
  labs(subtitle = expression(bold("B) Gross Primary Production")), y = expression("GPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  theme(legend.position = "none", plot.title = element_blank())

T0_M_plot_modified_R <- T0_M_plot_R + 
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, by = 0.25)) + 
  labs(subtitle = expression(bold("C) [Respiration]")), y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_blank())

# Combine all four modified plots
combined_plot <- plot_grid(
  T0_M_plot_modified_NPP, T0_M_plot_modified_GPP, T0_M_plot_modified_R,
  nrow = 1,  
  align = "hv",  
  rel_widths = c(1, 1),  
  rel_heights = c(1, 1))

combined_plot

#################################################################################
T1_M_plot_modified_NPP <- T1_M_plot_NPP + 
  scale_y_continuous(limits = c(0, 2.35), breaks = seq(-0, 2.35, by = 0.25)) + 
  labs(subtitle = expression(bold("A) Net Primary Production")), y = expression("NPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  ggtitle(expression(italic("Montipora capitata") ~ "T1 Metabolic Rates")) +  
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16))

T1_M_plot_modified_GPP <- T1_M_plot_GPP + 
  scale_y_continuous(limits = c(0, 2.35), breaks = seq(0, 2.35, by = 0.25)) + 
  labs(subtitle = expression(bold("B) Gross Primary Production")), y = expression("GPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  theme(legend.position = "none", plot.title = element_blank())

T1_M_plot_modified_R <- T1_M_plot_R + 
  scale_y_continuous(limits = c(0, 2.35), breaks = seq(0, 2.35, by = 0.25)) + 
  labs(subtitle = expression(bold("C) [Respiration]")), y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_blank())

combined_plot_M_T1 <- plot_grid(
  T1_M_plot_modified_NPP, T1_M_plot_modified_GPP, T1_M_plot_modified_R,
  nrow = 1,  
  align = "hv",  
  rel_widths = c(1, 1),  
  rel_heights = c(1, 1))

combined_plot_M_T1
#################################################################################
T2_M_plot_modified_NPP <- T2_M_plot_NPP + 
  scale_y_continuous(limits = c(-0.25, 1.75), breaks = seq(-0.25, 1.75, by = 0.25)) + 
  labs(subtitle = expression(bold("A) Net Primary Production")), y = expression("NPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_x_discrete(labels = c("Ambient" = "Amb", "Heated" = "Heat")) +
  ggtitle(expression(italic("Montipora capitata") ~ "T2 Metabolic Rates")) +  
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16))

T2_M_plot_modified_GPP <- T2_M_plot_GPP + 
  scale_y_continuous(limits = c(-0.25, 1.75), breaks = seq(-0.25, 1.75, by = 0.25)) + 
  labs(subtitle = expression(bold("B) Gross Primary Production")), y = expression("GPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_x_discrete(labels = c("Ambient" = "Amb", "Heated" = "Heat")) +
  theme(legend.position = "none", plot.title = element_blank())

T2_M_plot_modified_R <- T2_M_plot_R + 
  scale_y_continuous(limits = c(-0.25, 1.75), breaks = seq(-0.25, 1.75, by = 0.25)) + 
  labs(subtitle = expression(bold("C) [Respiration]")), y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_x_discrete(labels = c("Ambient" = "Amb", "Heated" = "Heat")) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_blank())

combined_plot_M_T2 <- plot_grid(
  T2_M_plot_modified_NPP, T2_M_plot_modified_GPP, T2_M_plot_modified_R,
  nrow = 1,  
  align = "hv",  
  rel_widths = c(1, 1),  
  rel_heights = c(1, 1))

combined_plot_M_T2

#################################################################################
T3_M_plot_modified_NPP <- T3_M_plot_NPP + 
  scale_y_continuous(limits = c(-0.25, 2.0), breaks = seq(-0.25, 2.0, by = 0.25)) + 
  labs(subtitle = expression(bold("A) Net Primary Production")), y = expression("NPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_x_discrete(labels = c("Ambient" = "Amb", "Heated" = "Heat")) +
  ggtitle(expression(italic("Montipora capitata") ~ "T3 Metabolic Rates")) +  
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16))

T3_M_plot_modified_GPP <- T3_M_plot_GPP + 
  scale_y_continuous(limits = c(-0.25, 2.0), breaks = seq(-0.25, 2.0, by = 0.25)) + 
  labs(subtitle = expression(bold("B) Gross Primary Production")), y = expression("GPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_x_discrete(labels = c("Ambient" = "Amb", "Heated" = "Heat")) +
  theme(legend.position = "none", plot.title = element_blank())

T3_M_plot_modified_R <- T3_M_plot_R + 
  scale_y_continuous(limits = c(-0.25, 2.0), breaks = seq(-0.25, 2.0, by = 0.25)) + 
  labs(subtitle = expression(bold("C) [Respiration]")), y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_x_discrete(labels = c("Ambient" = "Amb", "Heated" = "Heat")) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_blank())

combined_plot_M_T3 <- plot_grid(
  T3_M_plot_modified_NPP, T3_M_plot_modified_GPP, T3_M_plot_modified_R,
  nrow = 1,  
  align = "hv",  
  rel_widths = c(1, 1),  
  rel_heights = c(1, 1))

combined_plot_M_T3

#################################################################################
#################################################################################
T0_P_plot_modified_NPP <- T0_P_plot_NPP + 
  scale_y_continuous(limits = c(0, 4.0), breaks = seq(-0, 4.0, by = 0.25)) + 
  labs(subtitle = expression(bold("A) Net Primary Production")), y = expression("NPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  ggtitle(expression(italic("Porites compressa") ~ "T0 Metabolic Rates")) +  
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16))

T0_P_plot_modified_GPP <- T0_P_plot_GPP + 
  scale_y_continuous(limits = c(0, 4.0), breaks = seq(0, 4.0, by = 0.25)) + 
  labs(subtitle = expression(bold("B) Gross Primary Production")), y = expression("GPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  theme(legend.position = "none", plot.title = element_blank())

T0_P_plot_modified_R <- T0_P_plot_R + 
  scale_y_continuous(limits = c(0, 4.0), breaks = seq(0, 4.0, by = 0.25)) + 
  labs(subtitle = expression(bold("C) [Respiration]")), y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_blank())

# Combine all four modified plots
combined_plot_P_T0 <- plot_grid(
  T0_P_plot_modified_NPP, T0_P_plot_modified_GPP, T0_P_plot_modified_R,
  nrow = 1,  
  align = "hv",  
  rel_widths = c(1, 1),  
  rel_heights = c(1, 1))

combined_plot_P_T0

#################################################################################
T1_P_plot_modified_NPP <- T1_P_plot_NPP + 
  scale_y_continuous(limits = c(0, 4.5), breaks = seq(0, 4.5, by = 0.25)) + 
  labs(subtitle = expression(bold("A) Net Primary Production")), y = expression("NPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  ggtitle(expression(italic("Porites compressa") ~ "T1 Metabolic Rates")) +  
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16))

T1_P_plot_modified_GPP <- T1_P_plot_GPP + 
  scale_y_continuous(limits = c(0, 4.5), breaks = seq(0, 4.5, by = 0.25)) + 
  labs(subtitle = expression(bold("B) Gross Primary Production")), y = expression("GPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  theme(legend.position = "none", plot.title = element_blank())

T1_P_plot_modified_R <- T1_P_plot_R + 
  scale_y_continuous(limits = c(0, 4.5), breaks = seq(0, 4.5, by = 0.25)) + 
  labs(subtitle = expression(bold("C) [Respiration]")), y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_blank())

combined_plot_P_T1 <- plot_grid(
  T1_P_plot_modified_NPP, T1_P_plot_modified_GPP, T1_P_plot_modified_R,
  nrow = 1,  
  align = "hv",  
  rel_widths = c(1, 1),  
  rel_heights = c(1, 1))

combined_plot_P_T1
#################################################################################
T2_P_plot_modified_NPP <- T2_P_plot_NPP + 
  scale_y_continuous(limits = c(-0.25, 4.5), breaks = seq(-0.25, 4.5, by = 0.25)) + 
  labs(subtitle = expression(bold("A) Net Primary Production")), y = expression("NPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_x_discrete(labels = c("Ambient" = "Amb", "Heated" = "Heat")) +
  ggtitle(expression(italic("Porites compressa") ~ "T2 Metabolic Rates")) +  
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16))

T2_P_plot_modified_GPP <- T2_P_plot_GPP + 
  scale_y_continuous(limits = c(-0.25, 4.5), breaks = seq(-0.25, 4.5, by = 0.25)) + 
  labs(subtitle = expression(bold("B) Gross Primary Production")), y = expression("GPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_x_discrete(labels = c("Ambient" = "Amb", "Heated" = "Heat")) +
  theme(legend.position = "none", plot.title = element_blank())

T2_P_plot_modified_R <- T2_P_plot_R + 
  scale_y_continuous(limits = c(-0.25, 4.5), breaks = seq(-0.25, 4.5, by = 0.25)) + 
  labs(subtitle = expression(bold("C) [Respiration]")), y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_x_discrete(labels = c("Ambient" = "Amb", "Heated" = "Heat")) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_blank())

combined_plot_P_T2 <- plot_grid(
  T2_P_plot_modified_NPP, T2_P_plot_modified_GPP, T2_P_plot_modified_R,
  nrow = 1,  
  align = "hv",  
  rel_widths = c(1, 1),  
  rel_heights = c(1, 1))

combined_plot_P_T2

#################################################################################
T3_P_plot_modified_NPP <- T3_P_plot_NPP + 
  scale_y_continuous(limits = c(0.0, 4.0), breaks = seq(0.0, 4.0, by = 0.25)) + 
  labs(subtitle = expression(bold("A) Net Primary Production")), y = expression("NPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_x_discrete(labels = c("Ambient" = "Amb", "Heated" = "Heat")) +
  ggtitle(expression(italic("Porites compressa") ~ "T3 Metabolic Rates")) +  
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16))

T3_P_plot_modified_GPP <- T3_P_plot_GPP + 
  scale_y_continuous(limits = c(0.0, 4.0), breaks = seq(0.0, 4.0, by = 0.25)) + 
  labs(subtitle = expression(bold("B) Gross Primary Production")), y = expression("GPP [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_x_discrete(labels = c("Ambient" = "Amb", "Heated" = "Heat")) +
  theme(legend.position = "none", plot.title = element_blank())

T3_P_plot_modified_R <- T3_P_plot_R + 
  scale_y_continuous(limits = c(0.0, 4.0), breaks = seq(0.0, 4.0, by = 0.25)) + 
  labs(subtitle = expression(bold("C) [Respiration]")), y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_x_discrete(labels = c("Ambient" = "Amb", "Heated" = "Heat")) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.title = element_blank())

combined_plot_P_T3 <- plot_grid(
  T3_P_plot_modified_NPP, T3_P_plot_modified_GPP, T3_P_plot_modified_R,
  nrow = 1,  
  align = "hv",  
  rel_widths = c(1, 1),  
  rel_heights = c(1, 1))

combined_plot_P_T3