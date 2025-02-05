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
respo_data <- respo_data[!is.na(respo_data$R_SA_Corrected_cm2_hr) & 
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
filtered_data_MCAP_T0 <- respo_data %>%
  filter(Species == "MCAP", Time_Point == "T0") %>%
  mutate(
    Q1 = quantile(R_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(R_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = R_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) |
      R_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier) 

# Fit the model without the outliers
T0_M_lmer <- lmer(R_SA_Corrected_cm2_hr ~ Nutrient_Treatment + (1|Colony_ID), 
                  data = filtered_data_MCAP_T0)
summary(T0_M_lmer)

# Perform ANOVA and Tukey posthoc test
T0_M_anova <- anova(T0_M_lmer)
T0_M_anova

T0_M_emmeans <- emmeans(T0_M_lmer, ~ Nutrient_Treatment)
T0_M_emmeans

T0_M_tukey <- contrast(T0_M_emmeans, method = "pairwise", adjust = "tukey")
T0_M_tukey

# Print significant comparisons
significant_comparisons <- summary(T0_M_tukey) %>%
  dplyr::filter(p.value < 0.05)
significant_comparisons

# Checking residuals and qq-plot 
# Do not want any real patterns 
residuals <- resid(T0_M_lmer)
fitted_values <- fitted(T0_M_lmer)

# Residuals vs Fitted Plot
ggplot(data = data.frame(fitted = fitted_values, residuals = residuals), 
       aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "T0 MCAP Residuals vs Fitted Values",
       x = "Fitted Values",
       y = "Residuals") +
  theme_bw()

# QQ-Plot
qqnorm(residuals, main = "T0 MCAP QQ Plot of Residuals")
qqline(residuals, col = "red")

# Testing normality and homogeneity
shapiro.test(residuals)
leveneTest(R_SA_Corrected_cm2_hr ~ Nutrient_Treatment, 
           data = filtered_data_MCAP_T0)

# count how many data points are in each nutrient treatment
counts <- filtered_data_MCAP_T0 %>% 
  group_by(Nutrient_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts

# Plotting the data 
T0_M_plot <- ggplot(data = filtered_data_MCAP_T0, 
                    aes(x = Nutrient_Treatment, 
                        y = R_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts, aes(x = Nutrient_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) + 
  labs(title = expression(italic("Montipora capitata") ~ "Net Primary Production at T0"),
       x = "Nutrient Treatment",
       y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T0_M_plot

# Adding significant lines over boxplots 
T0_M_plot <- T0_M_plot +
  geom_signif(
    comparisons = list(c("Control", "Effluent"), 
                       c("Effluent", "Guano"), 
                       c("Effluent", "Inorganic")),
    annotations = c("*", "*", "**"),  
    y_position = c(0.7, 0.8, 1.0),  
    tip_length = 0.02,
    size = 0.5,
    textsize = 7,
    color = "black")

T0_M_plot
################################################################################
# Calculating and removing outliers
filtered_data_PCOMP_T0 <- respo_data %>%
  filter(Species == "PCOMP", Time_Point == "T0") %>%
  mutate(
    Q1 = quantile(R_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(R_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = R_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) | 
      R_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>% 
  select(-Q1, -Q3, -IQR, -is_outlier) 

# Fit the model without the outliers
T0_P_lmer <- lmer(R_SA_Corrected_cm2_hr ~ Nutrient_Treatment + (1|Colony_ID), 
                  data = filtered_data_PCOMP_T0)
T0_P_lmer

# Perform ANOVA. No post-hoc comparison since ANOVA is insignificant
T0_P_anova <- anova(T0_P_lmer)
T0_P_anova

# Checking residuals and qq-plot 
residuals <- resid(T0_P_lmer)
fitted_values <- fitted(T0_P_lmer)

# Residuals vs Fitted Plot
ggplot(data = data.frame(fitted = fitted_values, residuals = residuals), 
       aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "T0 PCOMP Residuals vs Fitted Values",
       x = "Fitted Values",
       y = "Residuals") +
  theme_bw()

# QQ-Plot
qqnorm(residuals, main = "T0 PCOMP QQ Plot of Residuals")
qqline(residuals, col = "red")

# Testing normality and homogeneity
shapiro.test(residuals)
leveneTest(R_SA_Corrected_cm2_hr ~ Nutrient_Treatment, 
           data = filtered_data_PCOMP_T0)

# count how many data points are in each nutrient treatment
counts <- filtered_data_PCOMP_T0 %>% 
  group_by(Nutrient_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts

# Plotting the data 
T0_P_plot <- ggplot(data= filtered_data_PCOMP_T0, 
                    aes(x = Nutrient_Treatment, 
                        y = R_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts, aes(x = Nutrient_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) + 
  labs(title = expression(italic("Porites compressa") ~ "Net Primary Production at T0"),
       x = "Nutrient Treatment",
       y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T0_P_plot

#################################################################################
# Calculating and removing outliers
filtered_data_MCAP_T1 <- respo_data %>%
  filter(Species == "MCAP", Time_Point == "T1") %>%
  mutate(
    Q1 = quantile(R_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(R_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = R_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) | 
      R_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>% 
  select(-Q1, -Q3, -IQR, -is_outlier) 

# Fit the model without outliers
T1_M_lmer <- lmer(R_SA_Corrected_cm2_hr ~ Nutrient_Treatment + (1|Colony_ID), 
                  data = filtered_data_MCAP_T1)
summary(T1_M_lmer)

# Perform ANOVA. No post-hoc test with insignificant ANOVA. 
T1_M_anova <- anova(T1_M_lmer)
T1_M_anova

#T1_M_emm <- emmeans(T1_M_lmer, pairwise ~ Nutrient_Treatment)
#T1_M_emm

#T1_M_tukey <- contrast(T1_M_emm, method = "pairwise", adjust = "tukey")
#T1_M_tukey

# Print significant comparisons
#sig_com_M_T1 <- summary(T1_M_tukey) %>%
 # dplyr::filter(p.value < 0.05)
#sig_com_M_T1

# Checking residuals and qq-plot 
residuals <- resid(T1_M_lmer)
fitted_values <- fitted(T1_M_lmer)

# Residuals vs Fitted Plot
ggplot(data = data.frame(fitted = fitted_values, residuals = residuals), 
       aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "T1 MCAP Residuals vs Fitted Values",
       x = "Fitted Values",
       y = "Residuals") +
  theme_bw()

# QQ-Plot
qqnorm(residuals, main = "T1 MCAP QQ Plot of Residuals")
qqline(residuals, col = "red")

# Testing normality and homogeneity
shapiro.test(residuals)
leveneTest(R_SA_Corrected_cm2_hr ~ Nutrient_Treatment, 
           data = filtered_data_MCAP_T1)

# count how many data points are in each nutrient treatment
counts <- filtered_data_MCAP_T1 %>% 
  group_by(Nutrient_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts

# Plotting the data 
T1_M_plot <- ggplot(data= filtered_data_MCAP_T1, 
                    aes(x = Nutrient_Treatment, 
                        y = R_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts, aes(x = Nutrient_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  labs(title = expression(italic("Montipora capitata") ~ "Net Primary Production at T1"),
       x = "Nutrient Treatment",
       y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T1_M_plot

#################################################################################
# Calculating and removing outliers
filtered_data_PCOMP_T1 <- respo_data %>%
  filter(Species == "PCOMP", Time_Point == "T1") %>%
  mutate(
    Q1 = quantile(R_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(R_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = R_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) | 
      R_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier)  

# Fit the model without the outliers
T1_P_lmer <- lmer(R_SA_Corrected_cm2_hr ~ Nutrient_Treatment + (1|Colony_ID), 
                  data = filtered_data_PCOMP_T1)
summary(T1_P_lmer)

# Perform ANOVA. Since the ANOVA is not significant, no posthoc tests were conducted
T1_P_anova <- anova(T1_P_lmer)
T1_P_anova

# Checking residuals and qq-plot 
residuals <- resid(T1_P_lmer)
fitted_values <- fitted(T1_P_lmer)

# Residuals vs Fitted Plot
ggplot(data = data.frame(fitted = fitted_values, residuals = residuals), 
       aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "T1 PCOMP Residuals vs Fitted Values",
       x = "Fitted Values",
       y = "Residuals") +
  theme_bw() 

# QQ-Plot
qqnorm(residuals, main = "T1 PCOMP QQ Plot of Residuals")
qqline(residuals, col = "red")

# Testing normality and homogeneity
shapiro.test(residuals)
leveneTest(R_SA_Corrected_cm2_hr ~ Nutrient_Treatment, 
           data = filtered_data_PCOMP_T1)

# count how many data points are in each nutrient treatment
counts <- filtered_data_PCOMP_T1 %>% 
  group_by(Nutrient_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts

# Plotting the data 
T1_P_plot <- ggplot(data= filtered_data_PCOMP_T1, 
                    aes(x = Nutrient_Treatment, 
                        y = R_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts, aes(x = Nutrient_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  labs(title = expression(italic("Porites compressa") ~ "Net Primary Production at T1"),
       x = "Nutrient Treatment",
       y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T1_P_plot

#################################################################################
# Calculating and removing outliers
filtered_data_MCAP_T2 <- respo_data %>%
  filter(Species == "MCAP", Time_Point == "T2") %>%
  mutate(
    Q1 = quantile(R_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(R_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = R_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) | 
      R_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier)  

# Fit the model without the outliers
T2_M_lmer <- lmer(R_SA_Corrected_cm2_hr ~ Nutrient_Treatment * Temp_Treatment + (1|Colony_ID), 
                  data = filtered_data_MCAP_T2)
summary(T2_M_lmer)

# Perform ANOVA. No post hoc test with insignificant ANOVA 
T2_M_anova <- anova(T2_M_lmer)
T2_M_anova

#T2_M_emm <- emmeans(T2_M_lmer, ~ Nutrient_Treatment * Temp_Treatment)
#T2_M_emm 

# Define contrasts
#contrasts <- list(
#  "Control Ambient - Control Heated" = c(-1, 0, 0, 0, 1, 0, 0, 0),
#  "Effluent Ambient - Effluent Heated" = c(0, -1, 0, 0, 0, 1, 0, 0),
#  "Guano Ambient - Guano Heated" = c(0, 0, 1, 0, 0, 0, -1, 0),
#  "Inorganic Ambient - Inorganic Heated" = c(0, 0, 0, 1, 0, 0, 0, -1))

# Perform contrasts with FDR adjustment
#contrast_M_T2 <- contrast(T2_M_emm, contrasts, adjust = "fdr")
#contrast_M_T2

# Print significant comparisons
#significant_comparisons_M_T2 <- summary(contrast_M_T2) %>%
#  dplyr::filter(p.value < 0.05)
#significant_comparisons_M_T2

#--------------------------------------------------------------------------------
#tukey_M_T2 <- contrast(T2_M_emm, method = "pairwise", adjust = "tukey")
#tukey_M_T2
#significant_comp_M_T2 <- summary(tukey_M_T2) %>%
#  dplyr::filter(p.value < 0.05)
#significant_comp_M_T2
#---------------------------------------------------------------------------------

# Checking residuals and qq-plot 
# Do not want any real patterns 
residuals <- resid(T2_M_lmer)
fitted_values <- fitted(T2_M_lmer)

# Residuals vs Fitted Plot
ggplot(data = data.frame(fitted = fitted_values, residuals = residuals), 
       aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "T2 MCAP Residuals vs Fitted Values",
       x = "Fitted Values",
       y = "Residuals") +
  theme_bw()

# QQ-Plot
qqnorm(residuals, main = "T2 MCAP QQ Plot of Residuals")
qqline(residuals, col = "red")

# Testing normality and homogeneity
shapiro.test(residuals)
leveneTest(R_SA_Corrected_cm2_hr ~ interaction(Nutrient_Treatment, Temp_Treatment), 
           data = filtered_data_MCAP_T2)

# count how many data points are in each nutrient treatment
counts <- filtered_data_MCAP_T2 %>% 
  group_by(Nutrient_Treatment, Temp_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts

# Plotting the data 
T2_M_plot <- ggplot(filtered_data_MCAP_T2, 
                    aes(x = Temp_Treatment, 
                        y = R_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  facet_wrap(~ Nutrient_Treatment, ncol = 4) +
  geom_text(data = counts, aes(x = Temp_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  labs(title = expression(italic("Montipora capitata") ~ "Net Primary Production at T2"),
       x = "Temperature Treatment",
       y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T2_M_plot

#################################################################################
# Calculating and removing outliers
filtered_data_PCOMP_T2 <- respo_data %>%
  filter(Species == "PCOMP", Time_Point == "T2") %>%
  mutate(
    Q1 = quantile(R_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(R_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = R_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) | 
      R_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>% 
  select(-Q1, -Q3, -IQR, -is_outlier) 

#Fit the model without outliers
T2_P_lmer <- lmer(R_SA_Corrected_cm2_hr ~ Nutrient_Treatment * Temp_Treatment + (1|Colony_ID), 
                  data = filtered_data_PCOMP_T2)
summary(T2_P_lmer)

# Step 3: Perform ANOVA. Since insignificant, no post hoc tests done. 
T2_P_anova <- anova(T2_P_lmer)
T2_P_anova

#T2_P_emm <- emmeans(T2_P_lmer, ~ Nutrient_Treatment * Temp_Treatment)

#Define contrasts
#contrasts <- list(
#  "Control Ambient - Control Heated" = c(-1, 0, 0, 0, 1, 0, 0, 0),
#  "Effluent Ambient - Effluent Heated" = c(0, -1, 0, 0, 0, 1, 0, 0),
#  "Guano Ambient - Guano Heated" = c(0, 0, 1, 0, 0, 0, -1, 0),
# "Inorganic Ambient - Inorganic Heated" = c(0, 0, 0, 1, 0, 0, 0, -1))

# Perform contrasts with FDR adjustment
#contrast_P_T2 <- contrast(T2_P_emm, contrasts, adjust = "fdr")
#contrast_P_T2

#significant_comparisons_P_T2 <- summary(contrast_P_T2) %>%
#  dplyr::filter(p.value < 0.05)
#significant_comparisons_P_T2

#-------------------------------------------------------------------------------
#tukey_P_T2 <- contrast(T2_P_emm, method = "pairwise", adjust = "tukey")
#tukey_P_T2
#significant_comp_P_T2 <- summary(tukey_P_T2) %>%
#  dplyr::filter(p.value < 0.05)
#significant_comp_P_T2
#-------------------------------------------------------------------------------

# Checking residuals and qq-plot 
# Do not want any real patterns 
residuals <- resid(T2_P_lmer)
fitted_values <- fitted(T2_P_lmer)

# Residuals vs Fitted Plot
ggplot(data = data.frame(fitted = fitted_values, residuals = residuals), 
       aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "T2 PCOMP Residuals vs Fitted Values",
       x = "Fitted Values",
       y = "Residuals") +
  theme_bw()

# QQ-Plot
qqnorm(residuals, main = "T2 PCOMP QQ Plot of Residuals")
qqline(residuals, col = "red")

# Testing normality and homogeneity
shapiro.test(residuals)
leveneTest(R_SA_Corrected_cm2_hr ~ interaction(Nutrient_Treatment, Temp_Treatment), 
           data = filtered_data_PCOMP_T2)

# count how many data points are in each nutrient treatment
counts <- filtered_data_PCOMP_T2 %>% 
  group_by(Nutrient_Treatment, Temp_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts

# Plotting the data 
T2_P_plot <- ggplot(filtered_data_PCOMP_T2, 
                    aes(x = Temp_Treatment, 
                        y = R_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts, aes(x = Temp_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  facet_wrap(~ Nutrient_Treatment, ncol=4) +
  labs(title = expression(italic("Porites compressa") ~ "Net Primary Production at T2"),
       x = "Temperature Treatment",
       y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T2_P_plot

#################################################################################
# Calculating and removing outliers
filtered_data_MCAP_T3 <- respo_data%>%
  filter(Species == "MCAP", Time_Point == "T3") %>%
  mutate(
    Q1 = quantile(R_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(R_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = R_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) | 
      R_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier)  

# Fit the model without the outliers
T3_M_lmer <- lmer(R_SA_Corrected_cm2_hr ~ Nutrient_Treatment * Temp_Treatment + (1|Colony_ID), 
                  data = filtered_data_MCAP_T3)
summary(T3_M_lmer)

# Perform ANOVA. Since ANOVA is insignificant, no planned contrasts take place. 
T3_M_anova <- anova(T3_M_lmer)
T3_M_anova

# Checking residuals and qq-plot 
# Do not want any real patterns 
residuals <- resid(T3_M_lmer)
fitted_values <- fitted(T3_M_lmer)

# Residuals vs Fitted Plot
ggplot(data = data.frame(fitted = fitted_values, residuals = residuals), 
       aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "T3 MCAP Residuals vs Fitted Values",
       x = "Fitted Values",
       y = "Residuals") +
  theme_bw()

# QQ-Plot
qqnorm(residuals, main = "T3 MCAP QQ Plot of Residuals")
qqline(residuals, col = "red")

# Testing normality and homogeneity
shapiro.test(residuals)
leveneTest(R_SA_Corrected_cm2_hr ~ interaction(Nutrient_Treatment, Temp_Treatment), 
           data = filtered_data_MCAP_T3)

# count how many data points are in each nutrient treatment
counts <- filtered_data_MCAP_T3 %>% 
  group_by(Nutrient_Treatment, Temp_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts

# Plotting the data 
T3_M_plot <- ggplot(data=filtered_data_MCAP_T3, 
                    aes(x =Temp_Treatment, 
                        y = R_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts, aes(x = Temp_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  facet_wrap(~ Nutrient_Treatment, ncol=4) +
  labs(title = expression(italic("Montipora capitata") ~ "Net Primary Production at T3"),
       x = "Temperature Treatment",
       y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~"]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T3_M_plot

####################################################################################
# Calculating and removing outliers
filtered_data_PCOMP_T3 <- respo_data %>%
  filter(Species == "PCOMP", Time_Point == "T3") %>%
  mutate(
    Q1 = quantile(R_SA_Corrected_cm2_hr, 0.25, na.rm = TRUE),
    Q3 = quantile(R_SA_Corrected_cm2_hr, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = R_SA_Corrected_cm2_hr < (Q1 - 1.5 * IQR) | 
      R_SA_Corrected_cm2_hr > (Q3 + 1.5 * IQR)) %>%
  filter(!is_outlier) %>%  
  select(-Q1, -Q3, -IQR, -is_outlier)  

# Fit the model without the outliers
T3_P_lmer <- lmer(R_SA_Corrected_cm2_hr ~ Nutrient_Treatment * Temp_Treatment + (1|Colony_ID), 
                  data = filtered_data_PCOMP_T3)
summary(T3_P_lmer)
T3_P_lmer_gen <- lmer(R_SA_Corrected_cm2_hr ~ Temp_Treatment + (1|Colony_ID), 
                      data = filtered_data_PCOMP_T3)

# Perform ANOVA. Since insignificant, no post hoc test done.
T3_P_anova <- anova(T3_P_lmer)
T3_P_anova

#T3_P_emm <- emmeans(T3_P_lmer, ~ Nutrient_Treatment * Temp_Treatment)

#Define contrasts
#contrasts <- list(
#  "Control Ambient - Control Heated" = c(-1, 0, 0, 0, 1, 0, 0, 0),
#  "Effluent Ambient - Effluent Heated" = c(0, -1, 0, 0, 0, 1, 0, 0),
#  "Guano Ambient - Guano Heated" = c(0, 0, 1, 0, 0, 0, -1, 0),
#  "Inorganic Ambient - Inorganic Heated" = c(0, 0, 0, 1, 0, 0, 0, -1))

#Perform contrasts with FDR adjustment
#contrast_P_T3 <- contrast(T3_P_emm, contrasts, adjust = "fdr")
#contrast_P_T3

# Checking residuals and qq-plot 
# Do not want any real patterns 
residuals <- resid(T3_P_lmer)
fitted_values <- fitted(T3_P_lmer)

# Residuals vs Fitted Plot
ggplot(data = data.frame(fitted = fitted_values, residuals = residuals), 
       aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "T3 PCOMP Residuals vs Fitted Values",
       x = "Fitted Values",
       y = "Residuals") +
  theme_bw() 

# QQ-Plot
qqnorm(residuals, main = "T3 PCOMP QQ Plot of Residuals")
qqline(residuals, col = "red")

# Testing normality and homogeneity
shapiro.test(residuals)
leveneTest(R_SA_Corrected_cm2_hr ~ interaction(Nutrient_Treatment, Temp_Treatment), 
           data = filtered_data_PCOMP_T3)

# count how many data points are in each nutrient treatment
counts <- filtered_data_PCOMP_T3 %>% 
  group_by(Nutrient_Treatment, Temp_Treatment) %>%
  summarise(n = n(), .groups = "drop")
counts

# Plotting the data 
T3_P_plot <- ggplot(data=filtered_data_PCOMP_T3, 
                    aes(x = Temp_Treatment, 
                        y = R_SA_Corrected_cm2_hr, 
                        color = Nutrient_Treatment, 
                        fill = Nutrient_Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Nutrient_Treatment), position = position_dodge(width = 0.8), alpha = 1) + 
  geom_text(data = counts, aes(x = Temp_Treatment, y = 0.0, label = paste0("n=", n)), 
            inherit.aes = FALSE, vjust = 0) +
  facet_wrap(~ Nutrient_Treatment, ncol=4) +
  labs(title = expression(italic("Porites compressa") ~ "Net Primary Production at T3"),
       x = "Temperature Treatment",
       y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_color_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

T3_P_plot

############################################################################################

# Modify the individual plots with new subtitles and y-axis title on the left only
T0_M_plot_modified <- T0_M_plot + 
  scale_y_continuous(limits = c(0.25, 2.0), breaks = seq(0.25, 2.0, by = 0.25)) + 
  labs(subtitle = "T0", y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +  
  theme(legend.position = "none", plot.title = element_blank())  
T0_M_plot_modified

T1_M_plot_modified <- T1_M_plot + 
  scale_y_continuous(limits = c(0.25, 2.0), breaks = seq(0.25, 2.0, by = 0.25)) + 
  labs(subtitle = "T1", y = NULL) + 
  theme(legend.position = "none", plot.title = element_blank())  

# Combine the modified plots side by side
combined_plot <- plot_grid(
  T0_M_plot_modified,
  T1_M_plot_modified,
  ncol = 2,
  align = "v",
  rel_widths = c(1, 1))

# Extract legend from one of the plots
legend <- get_legend(T0_M_plot + theme(legend.position = "bottom"))

# Create a joint title as a separate plot
joint_title <- ggdraw() + 
  draw_label(expression(italic("Montipora capitata") ~ "Net Primary Production", fontface = "bold", size = 14))

# Combine everything: title, combined plots, and legend
final_plot_M_T0 <- plot_grid(
  joint_title,  
  combined_plot, 
  legend, 
  ncol = 1,
  rel_heights = c(0.1, 1, 0.1))

final_plot_M_T0

########################################################################################

# Modify the individual plots with new subtitles and y-axis title on the left only
#T0_P_plot_modified <- T0_P_plot + 
#  scale_y_continuous(limits = c(0.25, 3.5), breaks = seq(0.25, 3.5, by = 0.25)) + 
#  labs(subtitle = "T0", y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) + 
#  theme(legend.position = "none", plot.title = element_blank())

T1_P_plot_modified <- T1_P_plot + 
  scale_y_continuous(limits = c(0.25, 3.5), breaks = seq(0.25, 3.5, by = 0.25)) + 
  labs(subtitle = "T1", y = NULL) + 
  theme(legend.position = "none", plot.title = element_blank())  

# Combine the modified plots side by side
combined_plot_P <- plot_grid(
  T0_P_plot_modified,
  T1_P_plot_modified,
  ncol = 2,
  align = "v",
  rel_widths = c(1, 1)
)

# Extract legend from one of the plots
legend <- get_legend(T0_P_plot + theme(legend.position = "bottom"))

# Create a joint title as a separate plot
joint_title_P <- ggdraw() + 
  draw_label(expression(italic("Porites compressa") ~ "Net Primary Production", fontface = "bold", size = 14))

# Combine everything: title, combined plots, and legend
final_plot_P_T0 <- plot_grid(
  joint_title_P,  
  combined_plot_P,  
  legend, 
  ncol = 1,
  rel_heights = c(0.1, 1, 0.1)  
)

final_plot_P_T0

########################################################################################
T2_M_plot_modified <- T2_M_plot + 
  scale_y_continuous(limits = c(-0.25, 1.75), breaks = seq(-0.25, 1.75, by = 0.25)) + 
  labs(subtitle = "T2", y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) + 
  theme(legend.position = "none", plot.title = element_blank()) 

T3_M_plot_modified <- T3_M_plot + 
  scale_y_continuous(limits = c(-0.25, 1.75), breaks = seq(-0.25, 1.75, by = 0.25)) + 
  labs(subtitle = "T3", y = NULL) +  
  theme(legend.position = "none", plot.title = element_blank()) 

# Combine the modified plots side by side
combined_plot_T2 <- plot_grid(
  T2_M_plot_modified,
  T3_M_plot_modified,
  ncol = 2,
  align = "v",
  rel_widths = c(1, 1))

# Extract legend from one of the plots
legend <- get_legend(T2_M_plot + theme(legend.position = "bottom"))

# Combine everything: title, combined plots, and legend
final_plot_M_T2 <- plot_grid(
  joint_title,  
  combined_plot_T2, 
  legend,  
  ncol = 1,
  rel_heights = c(0.1, 1, 0.1))

final_plot_M_T2
#########################################################################################

T2_P_plot_modified <- T2_P_plot + 
  scale_y_continuous(limits = c(-0.25, 3.5), breaks = seq(-0.25, 3.5, by = 0.25)) + 
  labs(subtitle = "T2", y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +  
  theme(legend.position = "none", plot.title = element_blank()) 

T3_P_plot_modified <- T3_P_plot + 
  scale_y_continuous(limits = c(-0.25, 3.5), breaks = seq(-0.25, 3.5, by = 0.25)) + 
  labs(subtitle = "T3", y = NULL) +  
  theme(legend.position = "none", plot.title = element_blank())  

# Combine the modified plots side by side
combined_plot_T3_P <- plot_grid(
  T2_P_plot_modified,
  T3_P_plot_modified,
  ncol = 2,
  align = "v",
  rel_widths = c(1, 1))

# Extract legend from one of the plots
legend <- get_legend(T2_P_plot + theme(legend.position = "bottom"))

# Combine everything: title, combined plots, and legend
final_plot_P_T2 <- plot_grid(
  joint_title_P,  
  combined_plot_T3_P, 
  legend,  
  ncol = 1,
  rel_heights = c(0.1, 1, 0.1))

final_plot_P_T2

####################################################################################
T0_M_plot_modified <- T0_M_plot + 
  scale_y_continuous(limits = c(-0.25, 1.25), breaks = seq(-0.25, 1.25, by = 0.25)) + 
  labs(subtitle = expression(bold("A) T0")), y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  theme(legend.position = "none", plot.title = element_blank())

T1_M_plot_modified <- T1_M_plot + 
  scale_y_continuous(limits = c(-0.25, 1.25), breaks = seq(-0.25, 1.25, by = 0.25)) + 
  labs(subtitle = expression(bold("B) T1")), y = NULL) +
  theme(legend.position = "none", plot.title = element_blank())

T2_M_plot_modified <- T2_M_plot + 
  scale_y_continuous(limits = c(-0.25, 1.25), breaks = seq(-0.25, 1.25, by = 0.25)) + 
  labs(subtitle = expression(bold("C) T2")), y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  theme(legend.position = "none", plot.title = element_blank())

T3_M_plot_modified <- T3_M_plot + 
  scale_y_continuous(limits = c(-0.25, 1.25), breaks = seq(-0.25, 1.25, by = 0.25)) + 
  labs(subtitle = expression(bold("D) T3")), y = NULL) +
  theme(legend.position = "none", plot.title = element_blank())

# Combine all four modified plots into a 2x2 grid
combined_plot_2x2 <- plot_grid(
  T0_M_plot_modified, T1_M_plot_modified,  
  T2_M_plot_modified, T3_M_plot_modified, 
  ncol = 2,  
  align = "hv",  
  rel_widths = c(1, 1),  
  rel_heights = c(1, 1))

# Extract legend from one of the plots
legend <- get_legend(T0_M_plot + theme(legend.position = "bottom"))

# Combine everything: title, combined plots, and legend
final_plot_ <- plot_grid(
  joint_title,       
  combined_plot_2x2,  
  legend,             
  ncol = 1,          
  rel_heights = c(0.1, 1, 0.1)  )

# Display the final plot
final_plot_
#####################################################################################

# Modify the individual plots with subtitles and y-axis titles
T0_P_plot_modified <- T0_P_plot + 
  scale_y_continuous(limits = c(-0.25, 2.0), breaks = seq(-0.25, 2.0, by = 0.25)) + 
  labs(subtitle = expression(bold("A) T0")), y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  theme(legend.position = "none", plot.title = element_blank())

T1_P_plot_modified <- T1_P_plot + 
  scale_y_continuous(limits = c(-0.25, 2.0), breaks = seq(-0.25, 2.0, by = 0.25)) + 
  labs(subtitle = expression(bold("B) T1")), y = NULL) +
  theme(legend.position = "none", plot.title = element_blank())

T2_P_plot_modified <- T2_P_plot + 
  scale_y_continuous(limits = c(-0.25, 2.0), breaks = seq(-0.25, 2.0, by = 0.25)) + 
  labs(subtitle = expression(bold("C) T2")), y = expression("Oxygen Consumption [µmol O"[2] ~ " cm"^"-2" ~ " hr"^"-1" ~ "]")) +
  theme(legend.position = "none", plot.title = element_blank())

T3_P_plot_modified <- T3_P_plot + 
  scale_y_continuous(limits = c(-0.25, 3.5), breaks = seq(-0.25, 3.5, by = 0.25)) + 
  labs(subtitle = expression(bold("D) T3")), y = NULL) +
  theme(legend.position = "none", plot.title = element_blank())

# Combine all four modified plots into a 2x2 grid
combined_plot_2x2 <- plot_grid(
  T0_P_plot_modified, T1_P_plot_modified,  # Top row
  T2_P_plot_modified, T3_P_plot_modified,  # Bottom row
  ncol = 2,  
  align = "hv",  
  rel_widths = c(1, 1), 
  rel_heights = c(1, 1))

# Extract legend from one of the plots
legend <- get_legend(T0_P_plot + theme(legend.position = "bottom"))

# Combine everything: title, combined plots, and legend
final_plot <- plot_grid(
  joint_title_P,        
  combined_plot_2x2,  
  legend,             
  ncol = 1,          
  rel_heights = c(0.1, 1, 0.1))

# Display the final plot
final_plot

