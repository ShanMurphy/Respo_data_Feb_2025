library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

nut_data <- read.csv("/Users/shannonmurphy/Desktop/Respirometry_Extras/Respo_Nutrient_Conc_Nov_7_2024.csv")
nut_data <- as.data.frame(nut_data)

nut_data$Nutrient_Treatment <- as.factor(nut_data$Nutrient_Treatment)
nut_data$Species <- as.factor(nut_data$Species)
nut_data$Time_Point <- as.factor(nut_data$Time_Point)

# Trim whitespace from the levels of Nutrient_Treatment
levels(nut_data$Nutrient_Treatment) <- trimws(levels(nut_data$Nutrient_Treatment))
levels(nut_data$Nutrient_Treatment)


nut_data <- nut_data %>%
  pivot_longer(cols = c(Dosing_Tot_N, 
                        Dosing_Tot_P, 
                        Dosing_N_P, 
                        Dosing_TIN,
                        Dosing_PO4, 
                        Dosing_NH4,
                        Dosing_N_N, 
                        Dosing_SiO4),
               names_to = "Dosing_Type",
               values_to = "Dosing_Value")

summary_data <- nut_data %>%
  group_by(Dosing_Type, Nutrient_Treatment, Species) %>%  
  summarise(
    mean_value = mean(Dosing_Value, na.rm = TRUE),  
    sd_value = sd(Dosing_Value, na.rm = TRUE) 
  ) %>%
  ungroup()

nut_plot <- ggplot(data = summary_data %>% filter(Species == "MCAP"),
                   aes(x =fct_reorder(Dosing_Type, mean_value, .desc = TRUE),
                       y = mean_value, 
                       fill = Nutrient_Treatment)) +
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
                position = position_dodge(width = 0.9), width = 0.2) +  
  theme_bw() +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") 
nut_plot

nut_plot_A <- ggplot(data = summary_data, 
       aes(x = Dosing_Type, 
            y = mean_value, 
           fill = Nutrient_Treatment)) +
geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
                position = position_dodge(width = 0.9), width = 0.2) +  # Error bars
  facet_wrap(~ Species, labeller = labeller(Species = c("MCAP" = "Montipora capitata",
                                                      "PCOMP" = "Porites compressa"))) +
  theme_bw() +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_x_discrete(limits = c("Dosing_TIN", "Dosing_N_P", 
                              "Dosing_Tot_N", "Dosing_Tot_P"),
                   labels = c("Dosing_Tot_N" = "Tot N", 
                              "Dosing_Tot_P" = "Tot P", 
                              "Dosing_N_P" = "N:P", 
                              "Dosing_TIN" = "Tot Inorg N")) +
  labs(title = "Nutrient Concentrations from Dosing Containers",
       x = "Nutrient Type", 
       y = expression("Nutrient Concentration [µmol L" ^"-1]")) +
  theme(plot.title = element_text(face = "bold", size = 14), strip.text = element_text(face = "italic", size = 12))
nut_plot_A

nut_plot_B <- ggplot(data = summary_data, 
       aes(x = Dosing_Type, 
           y = mean_value, 
           fill = Nutrient_Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
                position = position_dodge(width = 0.9), width = 0.2) +  
  facet_wrap(~ Species, labeller = labeller(Species = c("MCAP" = "Montipora capitata",
                                                        "PCOMP" = "Porites compressa"))) +
  theme_bw() +
  scale_fill_manual(values = c("Effluent" = "red1", "Control" = "cyan3", "Guano" = "palegreen3", "Inorganic" = "#FFCC00"), name = "Nutrient Treatment") +
  scale_x_discrete(limits = c("Dosing_NH4", "Dosing_N_N",
                              "Dosing_PO4", "Dosing_SiO4"), 
                   labels = c("Dosing_PO4" = expression(PO[4]), 
                              "Dosing_NH4" = expression(NH[4]), 
                              "Dosing_SiO4" = expression(SiO[4]), 
                              "Dosing_N_N" = bquote(NO[2] + NO[3]))) +
  labs(title = "Nutrient Concentrations from Dosing Containers",
       x = "Nutrient Type", 
       y = expression("Nutrient Concentration [µmol L" ^"-1]")) +
  theme(plot.title = element_text(face = "bold", size = 14), strip.text = element_text(face = "italic", size = 12))
nut_plot_B

