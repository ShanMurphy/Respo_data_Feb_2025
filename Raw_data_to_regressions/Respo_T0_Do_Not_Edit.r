# Import new packages 
install.packages("magrittr")
install.packages("knitr")

# Read in the necessary packages 

library(readxl)
library(dplyr)
library(ggplot2)
library(magrittr)
library(knitr)

# Find the data file on my computer and rename it 

file_path <- '/Users/shannonmurphy/Desktop/Respirometry/T0/P_Compressa/T0_PCOM_run_5_light_7_3_2023.xlsx'

# Create a new vector for the document
channel_list <- list()

# Combine all of the sheet names together from the individual excel document 

sheet_names <- c('SABD0000000003, Ch 2', 'SABD0000000003, Ch 3', 'SABD0000000003, Ch 4',
                 'SABD0000000003, Ch 5', 'SABD0000000003, Ch 6', 'SABD0000000003, Ch 7',
                 'SABD0000000003, Ch 8', 'SABD0000000003, Ch 9', 'SABD0000000003, Ch 10')

# For loop for reading the excel file and combining the sheet names and renaming it "channel_list"

for (sheet_name in sheet_names) {
  df <- read_excel(file_path, sheet = sheet_name)
  channel_list[[sheet_name]] <- df
    }

# Printing the channel list to double check all of the sheets appear 

#channel_list

# For loop assigns each channel name from the original data file 

for (i in 1:length(sheet_names)) {
  assign(paste0("channel_", i), channel_list[[sheet_names[i]]], envir = .GlobalEnv)
}

# Create a new column containing the list of the channel numbers

channel_list <- list(channel_2, channel_3, channel_4, channel_5, channel_6, channel_7, channel_8, channel_9, channel_10)

# For loop changes the column "Delta T" to "Delta_Time"

for (i in 1:length(channel_list)) {
  names(channel_list[[i]])[3] <- "Delta_Time" 
}

#channel_list

# For loop that creates a new column in the data frame for the different channels

for (i in 1:length(channel_list)) {
  channel_name <- paste0("channel_", i + 1)
  channel_list[[i]]$Channel <- as.factor(channel_name)
}

# For loop that merges all the data sheets together into one data frame

combined_df <- channel_list[[1]]
for (i in 2:length(channel_list)) {
  combined_df <- merge(combined_df, channel_list[[i]], all = TRUE)
}

# Selecting only the columns in the data frame that we want to work with
# The computer program collects more data than we need for this analysis

selected_columns <- combined_df %>%
  select(Delta_Time, Oxygen, Temperature, Channel)

head(selected_columns)

# Plotting the data full oxygen versus time data set categorized by channel 

ggplot(selected_columns, aes(x = Delta_Time, y = Oxygen, color = Channel, group = Channel)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Set color to black
  labs(title = "Oxygen vs. Time 0-60 Minutes",
       x = "Time (Minutes)",
       y = "Oxygen") +
  theme_minimal()

# Now we want to narrow down the selected data for 5-15 minutes
# We also want to remove any NA values from the data set

start_time <- 15.0
end_time <- 25.0
selected_data <- selected_columns[selected_columns$Delta_Time >= start_time & selected_columns$Delta_Time <= end_time, ]
selected_data <- na.omit(selected_data)

# Double checking the first and last data sections, which should contain 
# 5 to 15 minutes
head(selected_data)
tail(selected_data)

# We want to plot the designated 5 to 15 minutes and eventually take the slope value 
ggplot(selected_data, 
       aes(x = Delta_Time, y = Oxygen, color = Channel, group = Channel)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Set color to black
  labs(title = "Oxygen vs. Time 5-15 Minutes",
       x = "Time (Minutes)",
       y = "Oxygen") +
  theme_minimal()

# Applying the linear model to find the slope of each channel 

lm_models <- lapply(split(selected_data, selected_data$Channel), function(sub_data) {
  model <- lm(Oxygen ~ Delta_Time, data = sub_data)
  coefficients <- coef(model)
  data.frame(
    Channel = unique(sub_data$Channel),
    Slope = coefficients[2],
    Intercept = coefficients[1]
  )
})

# Creating a vector containing the slopes from the light runs 

slope_table_light <- do.call(rbind, lm_models)

# Need to define the "blank" channel that does not contain any coral 
slope_channel_8 <- slope_table_light$Slope[slope_table_light$Channel == "channel_8"]

intercept_channel_8 <- slope_table_light$Intercept[slope_table_light$Channel == "channel_8"]

# Finding the difference between the channels with coral and the blank channel 
# This removes any background photosynthesis by phytoplankton

slope_table_light$Light_Slope_Difference <- slope_table_light$Slope - slope_channel_8

# Adding the average temperature of each channel to the data frame 

Light_avg_temp <- selected_data %>% 
    group_by(Channel) %>%
    summarize(Light_avg_temp = mean(Temperature, na.rm = TRUE))

head(Light_avg_temp)

slope_table_light <- merge(slope_table_light, Light_avg_temp, by = "Channel", all.x = TRUE)

# Printing the slope table 

print(slope_table_light)

# Now repeating the entire process but importing the dark run data

file_path <- '/Users/shannonmurphy/Desktop/Respirometry/T0/P_Compressa/T0_PCOM_run_5_dark_7_3_2023.xlsx'

channel_list <- list()

sheet_names <- c('SABD0000000003, Ch 2', 'SABD0000000003, Ch 3', 'SABD0000000003, Ch 4',
                 'SABD0000000003, Ch 5', 'SABD0000000003, Ch 6', 'SABD0000000003, Ch 7',
                 'SABD0000000003, Ch 8', 'SABD0000000003, Ch 9', 'SABD0000000003, Ch 10')

for (sheet_name in sheet_names) {
  df <- read_excel(file_path, sheet = sheet_name)
  channel_list[[sheet_name]] <- df
    }

for (i in 1:length(sheet_names)) {
  assign(paste0("channel_", i), channel_list[[sheet_names[i]]], envir = .GlobalEnv)
}

channel_list <- list(channel_2, channel_3, channel_4, channel_5, channel_6, channel_7, channel_8, channel_9, channel_10)

for (i in 1:length(channel_list)) {
  names(channel_list[[i]])[3] <- "Delta_Time" 
}

for (i in 1:length(channel_list)) {
  channel_name <- paste0("channel_", i + 1)
  channel_list[[i]]$Channel <- as.factor(channel_name)
}

combined_df <- channel_list[[1]]
for (i in 2:length(channel_list)) {
  combined_df <- merge(combined_df, channel_list[[i]], all = TRUE)
}

#channel_list

selected_columns <- combined_df %>%
  select(Delta_Time, Oxygen, Temperature, Channel)

head(selected_data)

library(ggplot2)
ggplot(selected_columns, aes(x = Delta_Time, y = Oxygen, color = Channel, group = Channel)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Set color to black
  labs(title = "Oxygen vs. Full Time Period",
       x = "Time (Minutes)",
       y = "Oxygen") +
  theme_minimal()

start_time <- 5.0
end_time <- 15.0
selected_data <- selected_columns[selected_columns$Delta_Time >= start_time & selected_columns$Delta_Time <= end_time, ]
selected_data <- na.omit(selected_data)

head(selected_data)
tail(selected_data)

library(ggplot2)
ggplot(selected_data, aes(x = Delta_Time, y = Oxygen, color = Channel, group = Channel)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Set color to black
  labs(title = "Oxygen vs. Time 5-15 Minutes",
       x = "Time (Minutes)",
       y = "Oxygen") +
  theme_minimal()

lm_models <- lapply(split(selected_data, selected_data$Channel), function(sub_data) {
  model <- lm(Oxygen ~ Delta_Time, data = sub_data)
  coefficients <- coef(model)
  data.frame(
    Channel = unique(sub_data$Channel),
    Slope = coefficients[2],
    Intercept = coefficients[1]
  )
})

# Extract slopes

slope_table_dark <- do.call(rbind, lm_models)

slope_channel_8 <- slope_table_dark$Slope[slope_table_dark$Channel == "channel_8"]

intercept_channel_8 <- slope_table_dark$Intercept[slope_table_dark$Channel == "channel_8"]

slope_table_dark$Dark_Slope_Difference <- slope_table_dark$Slope - slope_channel_8

Dark_avg_temp <- selected_data %>%
    group_by(Channel) %>%
    summarize(Dark_avg_temp = mean(Temperature, na.rm = TRUE))

head(Dark_avg_temp)

slope_table_dark <- merge(slope_table_dark, Dark_avg_temp, by = "Channel", all.x = TRUE)

print(slope_table_dark)

# Merging the slope tables for the light and dark run 

merged_slope_table <- merge(slope_table_light, slope_table_dark, by = "Channel", suffixes = c("_Light", "_Dark" ))

# Creating a new column in the newly merged table called NEP (Photosynthesis - Respiration) 

merged_slope_table$NEP <- merged_slope_table$Light_Slope_Difference - merged_slope_table$Dark_Slope_Difference

# Arranging the merged slope table by channel to keep channel 10 at the bottom 

merged_slope_table <- merged_slope_table %>%
  arrange(as.numeric(gsub("channel_", "", Channel)))

# Printing the merged table to have all the data put together 

merged_slope_table

file_path <- '/Users/Shannonmurphy/Desktop/Respirometry/T0/T0_PCOMP_Frag_ID_List.xlsx'

sheet_name <- 'Run5'

frag_ids <- read_excel(file_path, sheet = sheet_name)

frag_ids <- frag_ids[, c("Channel", "Frag_ID")]

merged_slope_table <- merge(merged_slope_table, frag_ids, by = "Channel", all.x = TRUE)

merged_slope_table <- merged_slope_table %>%
  arrange(as.numeric(gsub("channel_", "", Channel)))

merged_slope_table

write.csv(merged_slope_table, file = '/Users/Shannonmurphy/Desktop/Respirometry/R_Analyzed/T0/P_Comp/T0_PCOMP_run_5_analyzed.csv', row.names = FALSE)




