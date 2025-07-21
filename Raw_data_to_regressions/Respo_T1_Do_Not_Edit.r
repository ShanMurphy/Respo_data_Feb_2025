# Read in the necessary packages 

library(readxl)
library(dplyr)
library(ggplot2)
library(magrittr)
library(knitr)

file <- read.csv("/Users/shannonmurphy/Desktop/Respirometry/T1/MCAP/T1_MCAP_run_1_light_8_26_2023_O2.csv", header = FALSE)

head(file)

# Assuming your data frame is named 'data'
# Extract the values from the second row
col_names <- as.character(file[2, , drop = TRUE])

# Set the column names manually
colnames(file) <- col_names

# Remove the second row
data <- file[-c(1, 2), , drop = FALSE]
data <- data[-nrow(data), , drop = FALSE]
data

data$Channel <- as.factor(data$Channel)
data$delta_t <- as.numeric(data$delta_t)
data$Value <- as.numeric(data$Value)
data$Temp <- as.numeric(data$Temp)

selected_columns <- data %>%
  select(Channel, delta_t, Value, Temp)

selected_columns <- selected_columns %>%
  filter(Channel != 4)

head(selected_columns)

# Assuming your data frame is named 'data'
# Arrange the data by 'Channel'
ordered_data <- selected_columns %>% 
  arrange(desc(Channel))

head(ordered_data)

ggplot(ordered_data, aes(x = delta_t, y = Value, color = Channel, group = Channel)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Set color to black
  labs(title = "Oxygen vs. Time 0-60 Minutes",
       x = "Time (Minutes)",
       y = "Oxygen") +
  theme_minimal()

start_time <- 5.0
end_time <- 15.0
selected_data <- ordered_data[ordered_data$delta_t >= start_time & ordered_data$delta_t <= end_time, ]
selected_data <- na.omit(selected_data)

head(selected_data)
tail(selected_data)

ggplot(selected_data, 
       aes(x = delta_t, y = Value, color = Channel, group = Channel)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Set color to black
  labs(title = "Oxygen vs. Time 5-15 Minutes",
       x = "Time (Minutes)",
       y = "Oxygen") +
  theme_minimal()

# Define the channels to skip
channels_to_skip <- c(4)

# Filter out data for channels to skip
filtered_data <- selected_data[!(selected_data$Channel %in% channels_to_skip), ]

# Split and apply linear regression models
lm_models <- lapply(split(filtered_data, filtered_data$Channel), function(sub_data) {
  if (sum(!is.na(sub_data$Value)) > 1 && sum(!is.na(sub_data$delta_t)) > 1) {
    # Check for at least 2 non-NA values in both Value and delta_t
    model <- lm(Value ~ delta_t, data = sub_data)
    coefficients <- coef(model)
    data.frame(
      Channel = unique(sub_data$Channel),
      Slope = coefficients[2],
      Intercept = coefficients[1]
    )
  } else {
    # Return a placeholder data frame when there are not enough non-NA values
    data.frame(Channel = NA, Slope = NA, Intercept = NA)
  }
})

slope_table_light <- do.call(rbind, lm_models)

slope_table_light

slope_channel_e_blank <- na.omit(slope_table_light$Slope[slope_table_light$Channel == 10])
slope_channel_c_blank <- na.omit(slope_table_light$Slope[slope_table_light$Channel == 7])

slope_channel_e_blank
slope_channel_c_blank

# Calculate differences for channels 2, 5, 8, and 10
difference_channel_2 <- na.omit(slope_table_light$Slope[slope_table_light$Channel == 2]) - slope_channel_e_blank
difference_channel_5 <- na.omit(slope_table_light$Slope[slope_table_light$Channel == 5]) - slope_channel_e_blank
difference_channel_8 <- na.omit(slope_table_light$Slope[slope_table_light$Channel == 8]) - slope_channel_e_blank
difference_channel_10 <- na.omit(slope_table_light$Slope[slope_table_light$Channel == 10]) - slope_channel_e_blank

# Add differences to the existing data frame
slope_table_light$Light_Slope_Difference[slope_table_light$Channel == 2] <- difference_channel_2
slope_table_light$Light_Slope_Difference[slope_table_light$Channel == 5] <- difference_channel_5
slope_table_light$Light_Slope_Difference[slope_table_light$Channel == 8] <- difference_channel_8
slope_table_light$Light_Slope_Difference[slope_table_light$Channel == 10] <- difference_channel_10

# Calculate differences for channels 3, 6, 9, and 7
difference_channel_3 <- na.omit(slope_table_light$Slope[slope_table_light$Channel == 3]) - slope_channel_c_blank
difference_channel_6 <- na.omit(slope_table_light$Slope[slope_table_light$Channel == 6]) - slope_channel_c_blank
difference_channel_9 <- na.omit(slope_table_light$Slope[slope_table_light$Channel == 9]) - slope_channel_c_blank
difference_channel_7 <- na.omit(slope_table_light$Slope[slope_table_light$Channel == 7]) - slope_channel_c_blank

# Add differences to the existing data frame
slope_table_light$Light_Slope_Difference[slope_table_light$Channel == 3] <- difference_channel_3
slope_table_light$Light_Slope_Difference[slope_table_light$Channel == 6] <- difference_channel_6
slope_table_light$Light_Slope_Difference[slope_table_light$Channel == 9] <- difference_channel_9
slope_table_light$Light_Slope_Difference[slope_table_light$Channel == 7] <- difference_channel_7

slope_table_light

# Adding the average temperature of each channel to the data frame 

Light_avg_temp <- selected_data %>% 
    group_by(Channel) %>%
    summarize(Light_avg_temp = mean(Temp, na.rm = TRUE))

Light_avg_temp

slope_table_light <- merge(slope_table_light, Light_avg_temp, by = "Channel", all.x = TRUE)

slope_table_light

file <- read.csv("/Users/shannonmurphy/Desktop/Respirometry/T1/MCAP/T1_MCAP_run_1_dark_8_26_2023_O2.csv", header = FALSE)

head(file)

# Assuming your data frame is named 'data'
# Extract the values from the second row
col_names <- as.character(file[2, , drop = TRUE])

# Set the column names manually
colnames(file) <- col_names

# Remove the second row
data <- file[-c(1, 2), , drop = FALSE]
data <- data[-nrow(data), , drop = FALSE]
data

data$Channel <- as.factor(data$Channel)
data$delta_t <- as.numeric(data$delta_t)
data$Value <- as.numeric(data$Value)
data$Temp <- as.numeric(data$Temp)

selected_columns <- data %>%
  select(Channel, delta_t, Value, Temp)

selected_columns <- selected_columns %>%
  filter(Channel != 4)

head(selected_columns)

ordered_data <- selected_columns %>% 
  arrange(desc(Channel))

head(ordered_data)

ggplot(ordered_data, aes(x = delta_t, y = Value, color = Channel, group = Channel)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Set color to black
  labs(title = "Oxygen vs. Time 0-60 Minutes",
       x = "Time (Minutes)",
       y = "Oxygen") +
  theme_minimal()

start_time <- 10.0
end_time <- 20.0
selected_data <- ordered_data[ordered_data$delta_t >= start_time & ordered_data$delta_t <= end_time, ]
selected_data <- na.omit(selected_data)

head(selected_data)
tail(selected_data)

ggplot(selected_data, 
       aes(x = delta_t, y = Value, color = Channel, group = Channel)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Set color to black
  labs(title = "Oxygen vs. Time 5-15 Minutes",
       x = "Time (Minutes)",
       y = "Oxygen") +
  theme_minimal()

# Define the channels to skip
channels_to_skip <- c(4)

# Filter out data for channels to skip
filtered_data <- selected_data[!(selected_data$Channel %in% channels_to_skip), ]

# Split and apply linear regression models
lm_models <- lapply(split(filtered_data, filtered_data$Channel), function(sub_data) {
  if (sum(!is.na(sub_data$Value)) > 1 && sum(!is.na(sub_data$delta_t)) > 1) {
    # Check for at least 2 non-NA values in both Value and delta_t
    model <- lm(Value ~ delta_t, data = sub_data)
    coefficients <- coef(model)
    data.frame(
      Channel = unique(sub_data$Channel),
      Slope = coefficients[2],
      Intercept = coefficients[1]
    )
  } else {
    # Return a placeholder data frame when there are not enough non-NA values
    data.frame(Channel = NA, Slope = NA, Intercept = NA)
  }
})

slope_table_dark <- do.call(rbind, lm_models)

slope_table_dark

slope_channel_e_blank <- na.omit(slope_table_dark$Slope[slope_table_dark$Channel == 10])
slope_channel_c_blank <- na.omit(slope_table_dark$Slope[slope_table_dark$Channel == 7])

slope_channel_e_blank
slope_channel_c_blank

# Calculate differences for channels 2, 5, 8, and 10
difference_channel_2 <- na.omit(slope_table_dark$Slope[slope_table_dark$Channel == 2]) - slope_channel_e_blank
difference_channel_5 <- na.omit(slope_table_dark$Slope[slope_table_dark$Channel == 5]) - slope_channel_e_blank
difference_channel_8 <- na.omit(slope_table_dark$Slope[slope_table_dark$Channel == 8]) - slope_channel_e_blank
difference_channel_10 <- na.omit(slope_table_dark$Slope[slope_table_dark$Channel == 10]) - slope_channel_e_blank

# Add differences to the existing data frame
slope_table_dark$Dark_Slope_Difference[slope_table_dark$Channel == 2] <- difference_channel_2
slope_table_dark$Dark_Slope_Difference[slope_table_dark$Channel == 5] <- difference_channel_5
slope_table_dark$Dark_Slope_Difference[slope_table_dark$Channel == 8] <- difference_channel_8
slope_table_dark$Dark_Slope_Difference[slope_table_dark$Channel == 10] <- difference_channel_10

# Calculate differences for channels 3, 6, 9, and 7
difference_channel_3 <- na.omit(slope_table_dark$Slope[slope_table_dark$Channel == 3]) - slope_channel_c_blank
difference_channel_6 <- na.omit(slope_table_dark$Slope[slope_table_dark$Channel == 6]) - slope_channel_c_blank
difference_channel_9 <- na.omit(slope_table_dark$Slope[slope_table_dark$Channel == 9]) - slope_channel_c_blank
difference_channel_7 <- na.omit(slope_table_dark$Slope[slope_table_dark$Channel == 7]) - slope_channel_c_blank

# Add differences to the existing data frame
slope_table_dark$Dark_Slope_Difference[slope_table_dark$Channel == 3] <- difference_channel_3
slope_table_dark$Dark_Slope_Difference[slope_table_dark$Channel == 6] <- difference_channel_6
slope_table_dark$Dark_Slope_Difference[slope_table_dark$Channel == 9] <- difference_channel_9
slope_table_dark$Dark_Slope_Difference[slope_table_dark$Channel == 7] <- difference_channel_7

slope_table_dark

Dark_avg_temp <- selected_data %>%
    group_by(Channel) %>%
    summarize(Dark_avg_temp = mean(Temp, na.rm = TRUE))

Dark_avg_temp

slope_table_dark <- merge(slope_table_dark, Dark_avg_temp, by = "Channel", all.x = TRUE)

slope_table_dark

merged_slope_table <- merge(slope_table_light, slope_table_dark, by = "Channel", suffixes = c("_Light", "_Dark" ))

merged_slope_table

# Calculate NEP by subtracting Dark_Slope_Difference from Light_Slope_Difference
merged_slope_table$NEP <- merged_slope_table$Light_Slope_Difference - merged_slope_table$Dark_Slope_Difference

# If you want to fill missing values with NA, you can use:
# merged_slope_table$NEP <- merged_slope_table$Light_Slope_Difference - merged_slope_table$Dark_Slope_Difference

#merged_slope_table$NEP <- na.omit(merged_slope_table$Light_Slope_Difference - merged_slope_table$Dark_Slope_Difference)

merged_slope_table <- merged_slope_table %>%
  arrange(as.numeric(gsub("channel_", "", Channel)))

file_path <- '/Users/Shannonmurphy/Desktop/Respirometry/T1/T1_MCAP_Frag_ID_List.xlsx'

sheet_name <- 'Run1'

frag_ids <- read_excel(file_path, sheet = sheet_name)

frag_ids <- frag_ids[, c("Channel", "Frag_ID")]

merged_slope_table <- merge(merged_slope_table, frag_ids, by = "Channel", all.x = TRUE)

# Printing the merged table to have all the data put together 

merged_slope_table

write.csv(merged_slope_table, file = '/Users/Shannonmurphy/Desktop/Respirometry/R_Analyzed/T1/M_Cap/T1_MCAP_run_1_analyzed.csv', row.names = FALSE)


















