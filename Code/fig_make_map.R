library(rstudioapi)
library(tidyverse)
library(ggplot2)
library(ggmap)
library(lubridate)
library(viridis)
library(patchwork)
library(cowplot)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)



#register_google(key = #WRITE YOUR CONFIDENTIAL GOOGLE KEY HERE#, write = TRUE)

theme <-  theme(panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.text = element_text(size =25),
                axis.title = element_text(size =25),
                axis.line = element_line(size = 0.5, linetype = "solid",
                                         colour = "black"),
                legend.title=element_text(size=25), 
                legend.text=element_text(size=25),
                plot.title = element_text(size=30))

# Coordinates of the points
points <- data.frame(
  Site = c("MartiniWall", "MartiniBay", "SharkIsland", "SBN-S", "SBN-N", "SBN-W"),
  Latitude = c(25.33904, 25.33492, 25.35151, 25.21504, 25.25874, 25.24414),
  Longitude = c(56.37819, 56.37917, 56.37575, 54.22265, 54.20809, 54.19723)
)

scale_shape_manual(name = "Site", values = c("MartiniBay"=15,
                                             "MartiniWall" = 16,
                                             "SharkIsland" = 17,
                                             "SBN-N" = 0,
                                             "SBN-W" = 1,
                                             "SBN-S" = 2))

base_map <- get_googlemap(center = c(lon = 55, lat = 25),
                          zoom = 7,
                          maptype = "satellite")


#as google map
overview_map <-
  ggmap(get_googlemap(center = c(lon = 55, lat = 25),
                                    zoom = 7,
                                    maptype = "satellite")) +
  xlab("Longitude") + ylab("Latitude") +
  geom_point(data = points[c(1,6),],
             aes(x = Longitude, y = Latitude),
             size = 8, shape = 0, color = "red", stroke = 2) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "tl", width_hint = 0.2) +
  coord_sf(crs = st_crs(4326))

overview_map


sbn_map <- ggmap(get_googlemap(center = c(long = 54.22, lat = 25.23), zoom = 13, maptype = "satellite",
                               color = "bw")) +
  xlab("Longitude") + ylab("Latitude") +
  geom_point(data = points, aes(x = Longitude, y = Latitude, shape = Site), size = 8,stroke = 2 , color = "white") +
  geom_text(data = points, aes(x = Longitude, y = Latitude + 0.006, label = Site), size = 8, color = "white") +
  scale_shape_manual(name = "Site", values = c(
    "SBN-N" = 0,
    "SBN-W" = 1,
    "SBN-S" = 2
  )) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(54.18, 54.27, length.out = 4), labels = sprintf("%.2f", seq(54.18, 54.27, length.out = 4))) +
  scale_y_continuous(breaks = seq(25.19, 25.275, length.out = 4), labels = sprintf("%.2f", seq(25.19, 25.275, length.out = 4))) +
  coord_fixed(1 / cos(25.23 * pi / 180), xlim = c(54.175, 54.27), ylim = c(25.19, 25.275))+
  xlab("")+ylab("")+
  theme

sbn_map



east_map <-ggmap(get_googlemap(center = c(long = 56.39, lat = 25.345), zoom = 14, maptype = "satellite", 
                               color = "bw")) +
  xlab("Longitude") + ylab("Latitude") +
  geom_point(data = points, aes(x = Longitude, y = Latitude, shape = Site), size = 8, color = "white") +
  geom_text(data = points, aes(x = Longitude, y = Latitude - 0.002, label = Site), size = 8, color = "white") +
  scale_shape_manual(name = "Site", values = c(
    "MartiniBay" = 15,
    "MartiniWall" = 16,
    "SharkIsland" = 17
  )) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(56.365, 56.39, by = 0.01), labels = sprintf("%.2f", seq(56.365, 56.39, by = 0.01))) +
  scale_y_continuous(breaks = seq(25.325, 25.365, by = 0.01), labels = sprintf("%.2f", seq(25.325, 25.365, by = 0.01))) +
  coord_fixed(1 / cos(25.345 * pi / 180), xlim = c(56.365, 56.39), ylim = c(25.325, 25.365))+
  xlab("")+ylab("")+
  theme

east_map










############ SAVE ########

pdf(file = "Output/Figures/overview_map_fig.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 10) # The height of the plot in inches
overview_map
# Step 3: Run dev.off() to create the file!
dev.off()


pdf(file = "Output/Figures/east_coast_map_fig.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 10) # The height of the plot in inches
east_map
# Step 3: Run dev.off() to create the file!
dev.off()

pdf(file = "Output/Figures/sbn_map_fig.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 10) # The height of the plot in inches
sbn_map
# Step 3: Run dev.off() to create the file!
dev.off()


