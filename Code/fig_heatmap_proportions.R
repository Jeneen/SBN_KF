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

# Helper function for standard error
se <- function(x) sd(x) / sqrt(length(x))

# Custom function to format median with two decimal places
format_median <- function(x) {
  sprintf("%.2f", median(x, na.rm = TRUE))
}

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

# Load data and clean
coral <- read.csv("Data/SBN_KF_surveys_averages.csv") %>% na.omit()
coralGroup <- read.csv("Data/SBN_KF_surveys_MajorGroups.csv") %>% na.omit()

# Define regions and bleaching periods
East <- c("MartiniBay", "MartiniWall", "SharkIsland")
SBN <- c("SBN-N", "SBN-S", "SBN-W")

pre2021 <- c("11/05/2021", "18/05/2021")
post2021 <- c("04/09/2021", "19/09/2021")
recovery2022 <- c("04/05/2022", "09/06/2022")

# Add region and bleaching stage to data
coralGroup <- coralGroup %>%
  mutate(
    Region = case_when(Site %in% East ~ "East coast", Site %in% SBN ~ "SBN"),
    Bleaching = case_when(
      Date %in% pre2021 ~ "Pre bleaching",
      Date %in% post2021 ~ "During bleaching",
      Date %in% recovery2022 ~ "Post bleaching"
    )
  )

coralGroup$Bleaching <- factor(coralGroup$Bleaching, levels = c("Pre bleaching", "During bleaching", "Post bleaching"))

# Set site factor levels
coralGroup$Site <- factor(
  coralGroup$Site,
  levels = c("MartiniBay", "MartiniWall", "SharkIsland", "SBN-N", "SBN-W", "SBN-S"))




# Heatmap
#just use genus level
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

# get mean of each site 
prop_data_mean <- prop_data %>% group_by(Region, Site, time_period, Genus) %>%
  summarise(mean_prop = mean(Proportion), se_prop = se(Proportion))



# set levels
prop_data_mean$Site <-factor(prop_data_mean$Site, levels = c("MartiniBay", "MartiniWall",
                                                        "SharkIsland", "SBN-N",
                                                        "SBN-W", "SBN-S"))
prop_data_mean$time_period <- factor(prop_data_mean$time_period, levels = c("Pre bleaching",
                                                                  "During bleaching",
                                                                  "Post bleaching"))

#only pre bleaching
prop_data_mean <- filter(prop_data_mean, time_period == "Pre bleaching")


#plot heatmap

prop_heatmap <- ggplot(prop_data_mean, aes(x = Site, y = factor(Genus, levels = rev(unique(Genus))), fill = mean_prop)) +
  geom_tile(color = "white") +  # Add a white border to tiles
  scale_fill_gradient(low = "white", high = "darkblue", name = "Proportion",
                      limits = c(0, 1)) +  # Gradient from white to blue
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "right"  # Place legend on the right
  )



#save figures

pdf(file = "Output/Figures/heatmap_proportions.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 10) # The height of the plot in inches
prop_heatmap 
dev.off()


