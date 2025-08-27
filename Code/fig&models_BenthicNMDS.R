library(ggplot2)
library(vegan)
library(tidyr)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
library(tibble)
library(broom)
library(tidyselect)
library(dplyr)

cols_region <- c("darkgrey", "lightgrey")
cols_bleaching <- brewer.pal(3, "Dark2")

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




#coral nmds
coral <- read.csv("Data/SBN_KF_surveys_averages.csv")
coralGroup <- read.csv("Data/SBN_KF_surveys_MajorGroups.csv")
coralGroup <- na.omit(coralGroup)

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

coralGroup <- coralGroup %>% mutate(time_period = case_when(
  Date %in% pre2021 ~ "Pre bleaching",
  Date %in% post2021 ~ "During bleaching",
  Date %in% recovery2022 ~ "Post bleaching"))

# Filter data for the "SBN" region
sbn_coral <- coral[coral$Region == "SBN", ]

# Run nMDS ordination with Bray Curtis metric for SBN region
sbnMDS <- metaMDS(sbn_coral[-c(1:5, 43:44)], metric = "bray", trymax = 1000)

# MDS summary for SBN region
sbnMDS # Stress: 0.15

#save nmds
distance <- sbnMDS$points
sbn_coral$NMDS1 <- distance[,1]
sbn_coral$NMDS2 <- distance[,2]
write.csv(sbn_coral, "Output/Tables/SBN_NMDS.csv")

# Check stress for SBN region
stressplot(sbnMDS)

# Run PERMANOVA to test for the effect of time period for SBN region
#perm_sbn <- adonis2(sbn_coral[-c(1:5, 43:44)] ~ time_period, sbn_coral, distance = "Bray", permutations = 999)
perm_sbn <- adonis2(sbn_coral[-c(1:5, 43:46)] ~ Site + time_period, sbn_coral, distance = "Bray", permutations = 999)
perm_sbn

# Test dispersion using PERMDISP
sbn_coral_disp <- betadisper(vegdist(sbn_coral[-c(1:5, 43:46)], method = "bray"), sbn_coral$time_period)
sbn_coral_disp

# Test for overall differences
anova(sbn_coral_disp) #need to transform, homogeneity not met


## Permutation test for pairwise comparisons
permutest(sbn_coral_disp, pairwise = TRUE)
# (note: permutation test, so slightly different p-values each time)

# Extract mean dispersion values (distance to centroid), along with se's and 95% CI
mod_sbn_coral.mean <- tapply(sbn_coral_disp$distances, sbn_coral$time_period, mean)
mod_sbn_coral.sd <- tapply(sbn_coral_disp$distances, sbn_coral$time_period, sd)
mod_sbn_coral.length <- tapply(sbn_coral_disp$distances, sbn_coral$time_period, length)
mod_sbn_coral.se <- mod_sbn_coral.sd / sqrt(mod_sbn_coral.length)
mod_sbn_coral.ci_low <- mod_sbn_coral.mean - (1.96 * mod_sbn_coral.se)
mod_sbn_coral.ci_high <- mod_sbn_coral.mean + (1.96 * mod_sbn_coral.se)

# Combine into a dataframe
mod_sbn_coral.out <- data.frame(mod_sbn_coral.mean, mod_sbn_coral.se, mod_sbn_coral.ci_low, mod_sbn_coral.ci_high)
mod_sbn_coral.out <- cbind(Treat_Year = rownames(mod_sbn_coral.out), mod_sbn_coral.out)
mod_sbn_coral.out


# Store MDS values for SBN region
site_scores_sbn <- as_tibble(scores(sbnMDS)$sites) %>%
  mutate(Region = sbn_coral$Region,
         Site = sbn_coral$Site,
         Transect = sbn_coral$Transect,
         Bleaching = sbn_coral$time_period)

# Create convex hulls function
convex_hull <- function(site_scores) site_scores[chull(site_scores$NMDS1, site_scores$NMDS2), ]

# Run convex hulls on time_period for SBN region
bleaching_hulls_sbn <- site_scores_sbn %>%
  group_by(Bleaching) %>%
  do(data.frame(convex_hull(.))) %>%
  ungroup()

# Now do SIMPER analysis to identify species with the highest contribution for SBN region
comp_sim_sbn <- as.list(simper(sbn_coral[-c(1:5, 43:44)], sbn_coral$time_period)) # removing metadata
simsum_sbn <- summary(comp_sim_sbn)
simsum_sbn 

# extract species with >0.025 contribution to cumsum
simspecies_sbn <- as.data.frame(rbind(simsum_sbn$`Post bleaching_Pre bleaching`, 
                                      simsum_sbn$`Post bleaching_During bleaching`,
                                      simsum_sbn$`Pre bleaching_During bleaching`)) %>%
  mutate(genus_species = rownames(.)) %>%
  filter(average >0.025)
simspecies_sbn
#write.csv(simspecies_sbn, "Output/Tables/simper_sbn.csv")

#remove number from genus_sp label for merge
simspecies_sbn$genus_species <- gsub("[[:digit:]]", "", simspecies_sbn$genus_species) 



# store point estimates of species and combine with simper results
mds.points_sbn <- as.data.frame(scores(sbnMDS, "species")) %>%
  add_column(genus_species = rownames(.)) %>%
  right_join(simspecies_sbn) %>%
  dplyr::select(NMDS1, NMDS2, genus_species) %>%
  unique(.)


# Update levels for Site and Bleaching factors for consistent plotting
site_scores_sbn$Site <- factor(site_scores_sbn$Site, levels = c("MartiniBay", "MartiniWall", "SharkIsland", "SBN-N","SBN-W", "SBN-S"))
bleaching_hulls_sbn$Bleaching <- factor(bleaching_hulls_sbn$Bleaching, levels = c("Pre bleaching", "During bleaching", "Post bleaching"))


# Create plot for SBN region
sbn_benthic_fig <- ggplot(site_scores_sbn, aes(x = NMDS1, y = NMDS2)) +
  geom_polygon(data = bleaching_hulls_sbn, aes(x = NMDS1, y = NMDS2, group = Bleaching, color = Bleaching, fill = Bleaching), color = "black", alpha = 0.3, lwd = 0.1) +
  geom_point(data = site_scores_sbn, aes(x = NMDS1, y = NMDS2, shape = Site)) +
  geom_segment(data = mds.points_sbn, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = mds.points_sbn, aes(x = NMDS1, NMDS2, label = genus_species), size = 5, inherit.aes = FALSE) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title.align = 0,
    legend.box = "horizontal"
  ) +
  scale_color_manual(values = cols_bleaching) +
  scale_fill_manual(values = cols_bleaching) +
  scale_shape_manual(
    name = "Site",
    values = c(
      "MartiniBay" = 15,
      "MartiniWall" = 16,
      "SharkIsland" = 17,
      "SBN-N" = 0,
      "SBN-W" = 1,
      "SBN-S" = 2
    )
  ) +
  guides(
    fill = guide_legend(nrow = 3, order =1),
    shape = guide_legend(ncol = 1)
  ) +
  scale_x_continuous(limits = c(-0.8, 0.6), breaks = seq(-0.6, 0.6, 0.4))+
  scale_y_continuous(breaks = seq(-0.5, 0.25, 0.25))+
  xlab("NMDS1") +
  ylab("NMDS2") +ggtitle("SBN")+
  theme +
  theme(legend.position = "none")


# Plot for SBN region
sbn_benthic_fig






# Filter data for the "East coast" region
east_coast_coral <- coral[coral$Region == "East coast", ]

# Run nMDS ordination with Bray Curtis metric for East coast region
east_coast_MDS <- metaMDS(east_coast_coral[-c(1:5, 43:44)], metric = "bray", trymax = 1000)

# MDS summary for East coast region
east_coast_MDS # Stress: 0.2

#save nmds of each site
distance <- east_coast_MDS$points
east_coast_coral$NMDS1 <- distance[,1]
east_coast_coral$NMDS2 <- distance[,2]
write.csv(east_coast_coral, "Output/Tables/EastCoast_NMDS.csv")

# Check stress for East coast region
stressplot(east_coast_MDS)

# Run PERMANOVA to test for the effect of time period  for East coast region

#perm_east_coast <- adonis2(east_coast_coral[-c(1:5, 43:44)] ~ time_period, east_coast_coral, distance = "Bray", permutations = 999)
perm_east_coast <- adonis2(east_coast_coral[-c(1:5, 43:46)] ~ Site + time_period, east_coast_coral, distance = "Bray", permutations = 999)
perm_east_coast


# Test dispersion using PERMDISP
east_coast_coral_disp <- betadisper(vegdist(east_coast_coral[-c(1:5, 43:46)], method = "bray"), 
                                    east_coast_coral$time_period)
east_coast_coral_disp 


# Test for overall differences
anova(east_coast_coral_disp) #assumption of homogeneity is not met --> transform


## Permutation test for pairwise comparisons - post hoc comparison
permutest(east_coast_coral_disp_trans, pairwise = TRUE)
# (note: permutation test, so slightly different p-values each time)



# Extract mean dispersion values (distance to centroid), along with se's and 95% CI
east_coast_coral_disp <- betadisper(vegdist(east_coast_coral[-c(1:5, 43:46)], method = "bray"), 
                                    east_coast_coral$time_period)
mod_east_coast_coral.mean <- tapply(east_coast_coral_disp$distances, east_coast_coral$time_period, mean)
mod_east_coast_coral.sd <- tapply(east_coast_coral_disp$distances, east_coast_coral$time_period, sd)
mod_east_coast_coral.length <- tapply(east_coast_coral_disp$distances, east_coast_coral$time_period, length)
mod_east_coast_coral.se <- mod_east_coast_coral.sd / sqrt(mod_east_coast_coral.length)
mod_east_coast_coral.ci_low <- mod_east_coast_coral.mean - (1.96 * mod_east_coast_coral.se)
mod_east_coast_coral.ci_high <- mod_east_coast_coral.mean + (1.96 * mod_east_coast_coral.se)

# Combine into a dataframe
mod_east_coast_coral.out <- data.frame(mod_east_coast_coral.mean, mod_east_coast_coral.se, mod_east_coast_coral.ci_low, mod_east_coast_coral.ci_high)
mod_east_coast_coral.out <- cbind(Treat_Year = rownames(mod_east_coast_coral.out), mod_east_coast_coral.out)
mod_east_coast_coral.out

# Store MDS values for East coast region
site_scores_east_coast <- as_tibble(scores(east_coast_MDS)$sites) %>%
  mutate(Region = east_coast_coral$Region,
         Site = east_coast_coral$Site,
         Transect = east_coast_coral$Transect,
         Bleaching = east_coast_coral$time_period)

# Run convex hulls on time_period for East coast region
bleaching_hulls_east_coast <- site_scores_east_coast %>%
  group_by(Bleaching) %>%
  do(data.frame(convex_hull(.))) %>%
  ungroup()

# Now do SIMPER analysis to identify species with the highest contribution for East coast region
comp_sim_east_coast <- as.list(simper(east_coast_coral[-c(1:5, 43:46)], east_coast_coral$time_period))
simsum_east_coast <- summary(comp_sim_east_coast)
#write.csv(simspecies_east_coast, "Output/Tables/simper_east.csv")


# extract species with >0.025 contribution to cumsum
simspecies_east_coast <- as.data.frame(rbind(simsum_east_coast$`Post bleaching_During bleaching`, 
                                             simsum_east_coast$`Post bleaching_Pre bleaching`,
                                             simsum_east_coast$`During bleaching_Pre bleaching`)) %>%
  mutate(genus_species = rownames(.)) %>%
  filter(average >0.025)
simspecies_east_coast

#remove number from genus_sp label for merge
simspecies_east_coast$genus_species <- gsub("[[:digit:]]", "", simspecies_east_coast$genus_species) 

# store point estimates of species and combine with simper results
mds.points_east_coast <- as.data.frame(scores(east_coast_MDS, "species")) %>%
  add_column(genus_species = rownames(.)) %>%
  right_join(simspecies_east_coast) %>%
  dplyr::select(NMDS1, NMDS2, genus_species) %>%
  unique(.)


# Update levels for Site and Bleaching factors for consistent plotting
site_scores_east_coast$Site <- factor(site_scores_east_coast$Site, levels = c("MartiniBay", "MartiniWall", "SharkIsland", "SBN-N","SBN-W", "SBN-S"))
bleaching_hulls_east_coast$Bleaching <- factor(bleaching_hulls_east_coast$Bleaching, levels = c("Pre bleaching", "During bleaching", "Post bleaching"))

# Create plot for East coast region
east_coast_benthic_fig <- ggplot(site_scores_east_coast, aes(x = NMDS1, y = NMDS2)) +
  geom_polygon(data = bleaching_hulls_east_coast, aes(x = NMDS1, y = NMDS2, group = Bleaching, color = Bleaching, fill = Bleaching), color = "black", alpha = 0.3, lwd = 0.1) +
  geom_point(data = site_scores_east_coast, aes(x = NMDS1, y = NMDS2, shape = Site)) +
  geom_segment(data = mds.points_east_coast, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = mds.points_east_coast, aes(x = NMDS1, y = NMDS2, label = genus_species), size = 5, inherit.aes = FALSE) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title.align = 0,
    legend.box = "horizontal"
  ) +
  scale_color_manual(values = cols_bleaching) +
  scale_fill_manual(values = cols_bleaching) +
  scale_shape_manual(
    name = "Site",
    values = c(
      "MartiniBay" = 15,
      "MartiniWall" = 16,
      "SharkIsland" = 17,
      "SBN-N" = 0,
      "SBN-W" = 1,
      "SBN-S" = 2
    )
  ) +
  guides(
    fill = guide_legend(nrow = 3, order =1),
    shape = guide_legend(ncol = 1)
  ) +
  scale_x_continuous(limits = c(-0.9, 1), breaks = seq(-1.0, 1, 0.5))+
  scale_y_continuous(breaks = seq(-1.0, 1, 0.5))+
  xlab("NMDS1") +
  ylab("NMDS2") + ggtitle("East coast")+
  theme +
  theme(legend.position = "none")

# Plot for East coast region
east_coast_benthic_fig


#save
pdf(file = "Output/Figures/benthic_nmds_east_sbn_comb_fig.pdf",   # The directory you want to save the file in
        width = 15, # The width of the plot in inches
        height = 10) # The height of the plot in inches

plot_grid(plot_grid(east_coast_benthic_fig, sbn_benthic_fig, ncol = 2, labels = c("A", "B")),nrow = 1)

dev.off()