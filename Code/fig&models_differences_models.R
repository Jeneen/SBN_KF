# Load required libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(ggpubr)
library(purrr)
library(patchwork)
library(cowplot)
library(brms)
library(tidybayes)
library(ggplot2)

# Helper function for standard error
se <- function(x) sd(x) / sqrt(length(x))


#make ggplot theme
theme <-  theme(panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.text = element_text(size =20),
                axis.title = element_text(size =20),
                axis.line = element_line(size = 0.5, linetype = "solid",
                                         colour = "black"),
                legend.title=element_text(size=20), 
                legend.text=element_text(size=20),
                plot.title = element_text(size=30))

# Set color palette for bleaching stages
cols_bleaching <- brewer.pal(3, "Dark2")

# Load data
coral <- read.csv("Data/SBN_KF_surveys_averages.csv") %>% na.omit()

#clean data
coral_long <- gather(coral, label, value, Acanthastrea.sp:Turf.algae)
coral_long$label <-gsub("\\..*","",coral_long$label)
coral_long$label <- as.factor(coral_long$label)
coral_long$Site <- as.factor(coral_long$Site)
coral_long$Transect <- as.factor(coral_long$Transect)
coral_long$Date <- as.factor(coral_long$Date)
coral_long$Month <- as.factor(coral_long$Month)
coral_long$label <- gsub("Transect", "Tape", coral_long$label)

#other = other invertebrate, all = other (fix names and merge columns)
Other <- c("Tape", "Unidentified", "Shadow", "All")
Macroalgae <- c("Algae", "Macroalgae") #add "Turf"
Dead <- c("Broken", "Dead" )
Urchin <- c("Pencil", "Diadema")
OtherInvertebrate <- c("Other", "Anemone", "Ascidian", "Bivalve", "Bryozoan", "Holothuria",
                       "Sponge")
coral_long <- coral_long %>% mutate(label2 = case_when(
  label %in% Other ~ "Other",
  label %in% Macroalgae ~ "Macroalgae",
  label %in% Dead ~ "Dead",
  label %in% Urchin ~ "Urchin",
  label %in% OtherInvertebrate ~ "Other Invertebrate"))
coral_long$label <- ifelse(is.na(coral_long$label2), coral_long$label, coral_long$label2)
coral_long <- coral_long %>% group_by(Site, Transect, Date, Month, Year, label) %>%
  summarise(value = sum(value))

#make matrix again
coral <- coral_long %>% pivot_wider(names_from= label, values_from = value, 
                                    names_repair = "unique")
colnames(coral)

#add region
levels(as.factor(coral$Site))
East <- c("MartiniBay", "MartiniWall",  "SharkIsland")
SBN <- c("SBN-N","SBN-S", "SBN-W" )

coral <- coral %>% mutate(Region = case_when(
  Site %in% East ~ "East coast",
  Site %in% SBN ~ "SBN"))

#add bleaching period
subset(coral, select = c(Site, Date)) %>% unique()
pre2021 <- c("11/05/2021", "18/05/2021")
post2021 <- c("04/09/2021", "19/09/2021")
recovery2022 <- c("04/05/2022", "09/06/2022")

coral <- coral %>% mutate(time_period = case_when(
  Date %in% pre2021 ~ "Pre bleaching",
  Date %in% post2021 ~ "During bleaching",
  Date %in% recovery2022 ~ "Post bleaching"))

#filter out non-coral benthos
coral <- coral %>%
  select(-any_of(c("Rock_Pavement", "Sand", "Shell", "Silt", 
                   "Bleached", "Dead", "CCA", "Echinometra", "Macroalgae",
                   "Other", "Other Invertebrate", 
                   "Turf", "Urchin", "Zoanthid")))


#calculate proportion of total coral
coral_prop <- coral %>%
  mutate(
    total = rowSums(across(Acanthastrea:Stylophora))  # Add a total column
  ) %>%
  mutate(
    across(Acanthastrea:Stylophora, ~ (.x / total), .names = NULL)  # Overwrite columns with proportions
  )




# Replace NaN with 0 in the dataset
coral_prop <- coral_prop %>% replace(is.na(.), 0)



# Prepare the data
prop_data <- coral_prop %>%
  pivot_longer(
    cols = Acanthastrea:Stylophora,  # Adjust to include all coral genus columns
    names_to = "Genus",
    values_to = "Proportion"
  )

# get mean of each site (remove site)
prop_data <- prop_data %>% group_by(Region, Site, time_period, Genus) %>%
  summarise(mean_prop = mean(Proportion), se_prop = se(Proportion))

#get difference in pre and post 
epsilon <- 1e-6 
differences <- prop_data %>%
  filter(time_period %in% c("Pre bleaching", "Post bleaching")) %>%  # Filter relevant time periods
  group_by(Region, Site, Genus, time_period) %>%
  summarise(mean_prop = mean(mean_prop, na.rm = TRUE), .groups = "drop") %>%  # Calculate means if not already grouped
  pivot_wider(
    names_from = time_period,
    values_from = mean_prop
  ) %>%
  mutate(
    difference = (`Post bleaching` - `Pre bleaching`)  # Calculate the difference
  ) %>%
  select(Region, Site, Genus, `Pre bleaching`, `Post bleaching`, difference)  # Select desired columns



#east coast
differences_east <- filter(differences, Region == "East coast")
differences_east$difference_log_signed <- sign(differences_east$difference) * log1p(abs(differences_east$difference))
differences_east$difference_sqrt <- sign(differences_east$difference) * sqrt(abs(differences_east$difference))


priors <- c(
  prior(student_t(3, 0, 2.5), class = Intercept),  # Centered at 0, allowing some flexibility
  prior(normal(0, 1), class = b),  # Shrinks coefficients but allows variation
  prior(student_t(3, 0, 2.5), class = sigma)  # Ensures positive values for residual standard deviation
)


m1_east_p <- brm(
  formula = bf(difference ~ Site + Genus),  # Fixed effects model
  data = differences_east,  # Dataset
  family = gaussian(),
  control = list(adapt_delta = 0.95),  # For convergence
  iter = 5000,  # Increase iterations for better estimation
  warmup = 2000,  # Warmup period
  chains = 4,  # Parallel chains
  cores = 4,  # Parallel computation
  seed = 10,
  prior = priors
)
summary(m1_east_p)
#tab_model(m1_east_p)
plot(m1_east_p)
bayes_R2(m1_east_p)
pp_check(m1_east_p)
pp_check(m1_east_p, type = "boxplot")
pp_check(m1_east_p, type = "stat_2d")
neff_ratio(m1_east_p)
coef <- summary(m1_east_p)$fixed
east_cond <- conditional_effects(m1_east_p)$Genus




# Generate posterior predictions for ploting predictive checks
differences_east$Site <- droplevels(differences_east$Site)
posterior_data <- add_predicted_draws(m1_east_p, newdata = differences_east, ndraws= 1000, seed = 10)

# Extract posterior predictive intervals per genus
posterior_intervals <- posterior_data %>%
  group_by(Genus) %>%  # Retain Genus before summarizing
  median_qi(.prediction, .width = c(0.9, 0.5)) 

# Extract observed data intervals per genus
observed_intervals <- differences_east %>%
  group_by(Genus) %>%
  median_qi(difference, .width = c(0.9, 0.5)) 

# Plot
pp_east <- ggplot() +
  # Posterior predictive intervals (blue)
  geom_pointinterval(data = posterior_intervals, 
                     aes(x = .prediction, xmin = .lower, xmax = .upper, y = Genus), 
                     color = "blue", alpha = 0.7) +
  
  # Observed data intervals (black, shifted downward)
  geom_pointinterval(data = observed_intervals, 
                     aes(x = difference, xmin = .lower, xmax = .upper, y = Genus), 
                     color = "black", alpha = 0.7, position = position_nudge(y = -0.5)) +
  
  facet_wrap(~Genus, scales = "free_y") +
  labs(title = "Posterior Predictive Checks by Genus", 
       x = "Difference (%)") +
  
  theme_minimal() + 
  theme(legend.position = "none", 
        axis.text.y = element_blank(),  
        axis.title.y = element_blank())

jpeg("Output/Figures/east_propcoral_pp.jpg", width = 1200, height = 800, res = 150)
# Generate the trace plots
pp_east
# Close the graphics device
dev.off()

jpeg("Output/Figures/east_propcoral_traceplots_%d.jpg", width = 1200, height = 800, res = 150)
# Generate the trace plots
plot(m1_east_p)
# Close the graphics device
dev.off()


#sbn
differences_sbn <- filter(differences, Region == "SBN")

m1_sbn_p <- brm(
  formula = bf(difference ~ Site + Genus ),  # Fixed effects model
  data = differences_sbn,  # Dataset
  family = gaussian(),  # 
  control = list(adapt_delta = 0.95),  # For convergence
  iter = 5000,  # Increase iterations for better estimation
  warmup = 2000,  # Warmup period
  chains = 4,  # Parallel chains
  cores = 4, # Parallel computation
  seed = 10,
  prior = priors
)
summary(m1_sbn_p)
plot(m1_sbn_p)
sbn_cond <- conditional_effects(m1_sbn_p, "Genus")$Genus
bayes_R2(m1_sbn_p)
coef_sbn <- summary(m1_sbn_p)$fixed

# Generate posterior predictions for plotting predictive checks
differences_sbn$Site <- droplevels(differences_sbn$Site)
posterior_data <- add_predicted_draws(m1_sbn_p, newdata = differences_sbn, ndraws = 1000, seed = 10)

# Extract posterior predictive intervals per genus
posterior_intervals <- posterior_data %>%
  group_by(Genus) %>%  # Retain Genus before summarizing
  median_qi(.prediction, .width = c(0.9, 0.5)) 

# Extract observed data intervals per genus
observed_intervals <- differences_sbn %>%
  group_by(Genus) %>%
  median_qi(difference, .width = c(0.9, 0.5)) 

# Plot
pp_sbn <- ggplot() +
  # Posterior predictive intervals (blue)
  geom_pointinterval(data = posterior_intervals, 
                     aes(x = .prediction, xmin = .lower, xmax = .upper, y = Genus), 
                     color = "blue", alpha = 0.7) +
  
  # Observed data intervals (black, shifted downward)
  geom_pointinterval(data = observed_intervals, 
                     aes(x = difference, xmin = .lower, xmax = .upper, y = Genus), 
                     color = "black", alpha = 0.7, position = position_nudge(y = -0.5)) +
  
  facet_wrap(~Genus, scales = "free_y") +
  labs(title = "Posterior Predictive Checks by Genus", 
       x = "Difference (%)") +
  
  theme_minimal() + 
  theme(legend.position = "none", 
        axis.text.y = element_blank(),  
        axis.title.y = element_blank())



jpeg("Output/Figures/sbn_propcoral_pp.jpg", width = 1200, height = 800, res = 150)
# Generate the trace plots
pp_sbn 
# Close the graphics device
dev.off()

jpeg("Output/Figures/sbn_propcoral_traceplots_%d.jpg", width = 1200, height = 800, res = 150)
# Generate the trace plots
plot(m1_sbn_p)
# Close the graphics device
dev.off()



#conditional plots
east_cond$color <- ifelse(
  east_cond$lower__ > 0, "positive",                # Positive effect (does not overlap 0)
  ifelse(east_cond$upper__ < 0, "negative", "zero") # Negative effect or overlaps 0
)

east_cond_plot <- ggplot(east_cond, aes(x = factor(Genus, levels = rev(unique(Genus))), y = estimate__)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +  # Add dashed line at 0
  geom_point(aes(color = color),size = 5) +
  scale_color_manual(
    values = c("positive" = "darkgreen", "negative" = "darkred", "zero" = "grey"),
  ) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__), width = 0.2) +
  theme_minimal() +
  labs(
    x = "",
    y = "Conditional effect"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")+
  theme+
  coord_flip()

sbn_cond$color <- ifelse(
  sbn_cond$lower__ > 0, "positive",                # Positive effect (does not overlap 0)
  ifelse(sbn_cond$upper__ < 0, "negative", "zero") # Negative effect or overlaps 0
)

sbn_cond_plot <- ggplot(sbn_cond, aes(x = factor(Genus, levels = rev(unique(Genus))), y = estimate__)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +  # Add dashed line at 0
  geom_point(aes(color = color),size = 5) +
  scale_color_manual(
    values = c("positive" = "darkgreen", "negative" = "darkred", "zero" = "grey"),
  ) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__), width = 0.2) +
  theme_minimal() +
  labs(
    x = "",
    y = "Conditional effect"
  ) +
  ylim(-2,2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_blank(), legend.position = "none")+
  theme+
  coord_flip()

#save figures

pdf(file = "Output/Figures/difference_genus_cond_plot.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 10) # The height of the plot in inches
plot_grid(
  plot_grid(east_cond_plot, sbn_cond_plot,  ncol = 2, labels = c("A", "B"),  rel_widths = c(1.2,1), label_size = 24)
)
# Step 3: Run dev.off() to create the file!
dev.off()



