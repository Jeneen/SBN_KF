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

# Pre-compute the medians grouped by Bleaching and Region
medians <- coralGroup %>%
  group_by(Bleaching, Region) %>%
  summarize(median_value = round(median(Live.coral, na.rm = TRUE), 2), .groups = "drop")



# Boxplot all coral -------------------------------------------------------


#east alone
east_alone <- filter(coralGroup, Region == "East coast")
median_east <- filter(medians, Region == "East coast")
bleaching_east_fig <- ggplot(east_alone, aes(x = Bleaching, 
                                             y = Live.coral, 
                                             fill = Bleaching)) + 
  geom_boxplot(aes(alpha = Region)) +
  stat_summary(
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.75),
    shape = 4,
    fill = "black",
    size = 3
  ) +
  geom_text(
    data = median_east,
    aes(x = Bleaching, y = median_value, label = sprintf("%.2f", median_value), alpha = Region),
    position = position_dodge(width = 0.75),
    vjust = -20,
    color = "black",
    size = 5
  ) +
  geom_jitter(aes(alpha = Region, fill = Bleaching, shape = Site), 
              position = position_jitterdodge(), size = 3) +
  scale_fill_manual(values = c(cols_bleaching[1], cols_bleaching[2], cols_bleaching[3])) +
  scale_alpha_manual(values = c(1,  1, 1)) +
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
  theme_minimal() +
  guides(fill = guide_legend(nrow = 3),
         shape = guide_legend(ncol = 2)) +
  ylab("Live coral (%)") +
  xlab("") +
  ylim(0,100)+
  theme+
  theme(legend.position = "none")

bleaching_east_fig



#sbn alone
sbn_alone <- filter(coralGroup, Region == "SBN")
median_sbn <- filter(medians, Region == "SBN")
bleaching_sbn_fig <- ggplot(sbn_alone, aes(x = Bleaching, 
                                           y = Live.coral, 
                                           fill = Bleaching)) + 
  geom_boxplot(aes(alpha = Region)) +
  geom_jitter(aes(alpha = Region, fill = Bleaching, shape = Site), 
              position = position_jitterdodge(), size = 3) +
  stat_summary(
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.75),
    shape = 4,
    fill = "black",
    size = 3
  ) +
  geom_text(
    data = median_sbn,
    aes(x = Bleaching, y = median_value, label = sprintf("%.2f", median_value), alpha = Region),
    position = position_dodge(width = 0.75),
    vjust = -20,
    color = "black",
    size = 5
  ) +
  scale_fill_manual(values = c(cols_bleaching[1], cols_bleaching[2], cols_bleaching[3])) +
  scale_alpha_manual(values = c(0.5, 1, 0.5)) +
  scale_shape_manual(
    name = "Site",
    values = c(
      "SBN-N" = 0,
      "SBN-W" = 1,
      "SBN-S" = 2
    )
  ) +
  theme_minimal() +
  guides(fill = guide_legend(nrow = 3),
         shape = guide_legend(ncol = 2)) +
  ylab("Live coral (%)") +
  xlab("") +
  ylim(0,100)+
  theme+
  theme(legend.position = "none",   axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

bleaching_sbn_fig


coralGroup$Site <- factor(coralGroup$Site, c("MartiniBay", "MartiniWall", "SharkIsland",
                                         "SBN-N", "SBN-W", "SBN-S"))
legend_fig <- ggplot(coralGroup, aes(x = Bleaching, y = Live.coral)) + 
  geom_jitter(aes(shape = Site), position = position_jitterdodge(),
              size =3) +
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
  guides(fill = guide_legend(nrow = 3),
         shape = guide_legend(ncol = 2)) +
  ylab("Percentage live coral") +
  xlab("") +
  theme


legend <- get_legend(legend_fig)
legend <- as_ggplot(legend)
legend 



#Save PDF 
pdf(file = "Output/Figures/live_coral_boxplot.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 10) # The height of the plot in inches
plot_grid(
  plot_grid(bleaching_east_fig, bleaching_sbn_fig,  ncol = 2, labels = c("A", "B"), label_size = 24)
)
dev.off()



#save summary table
summary_table <- coralGroup %>% select(Site:Live.coral)
summary_table <- summary_table %>% select(c(Site,Transect,Date,Live.coral))
write.csv(summary_table, "Output/Tables/summary_live_coral.csv")
