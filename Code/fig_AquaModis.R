library(raster)
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)

# Define the latitude and longitude points
latitude <- c(25.33904, 25.33492, 25.35151, 25.21504, 25.25874, 25.24414)
longitude <- c(56.37819, 56.37917, 56.37575, 54.22265, 54.20809, 54.19723)

points <- data.frame(
  Site = c("MartiniWall", "MartiniBay", "SharkIsland", "SBN-S", "SBN-N", "SBN-W"),
  Region = c("EastCoast", "EastCoast", "EastCoast", "SBN", "SBN", "SBN"),
  latitude = c(25.33904, 25.33492, 25.35151, 25.21504, 25.25874, 25.24414),
  longitude = c(56.37819, 56.37917, 56.37575, 54.22265, 54.20809, 54.19723)
)

# Get the list of .nc files in the folder
files <- list.files(path = "Data/AquaModis", pattern = "sst.nc", full.names = TRUE)

# Create an empty dataframe to store the results
results <- data.frame()

# Loop through each file
for (file in files) {
  # Load the .nc file as a RasterStack
  stack <- stack(file, varname = "sst")
  
  # Extract the data for the specified latitude and longitude points
  sst_values <- extract(stack, data.frame(lon = longitude, lat = latitude))
  
  # Get the date from the filename (format "YYYYMMDD")
  date <- as.Date(str_extract(basename(file), "(?<=\\.)\\d{8}"), format = "%Y%m%d")
  
  # For each missing point, find the nearest non-missing value
  for (i in which(is.na(sst_values))) {
    # Calculate the distance to all other points
    distances <- pointDistance(matrix(c(longitude[i], latitude[i]), nrow = 1),
                               cbind(longitude, latitude), lonlat = TRUE)
    
    # Identify the nearest point with a non-missing value
    ordered <- order(distances)
    for (j in ordered) {
      if (!is.na(sst_values[j])) {
        sst_values[i] <- sst_values[j]
        break
      }
    }
  }
  
  # Create a temporary dataframe with the extracted data
  temp_df <- data.frame(date = rep(date, length(latitude)), 
                        latitude = latitude, 
                        longitude = longitude, 
                        sst = sst_values)
  
  # Append the temporary dataframe to the result dataframe
  results <- rbind(results, temp_df)
}

# Return the result dataframe
results

#merge with site
all_sst <- inner_join(results, points)


#number of days above bleaching thresholds
thresh <- filter(all_sst, date <="2021-12-01")
sbn_thresh <- filter(thresh, Region == "SBN")
sbn_thresh <- sbn_thresh %>% group_by(Region, date) %>% summarise(mean = mean(Sea.Surface.Temperature))
sbn_thresh <- filter(sbn_thresh, mean >= 34)
sbn_thresh2 <- filter(sbn_thresh, mean >= 35)
east_thresh <- filter(thresh, Region == "EastCoast")
east_thresh <- east_thresh %>% group_by(Region, date) %>% summarise(mean = mean(Sea.Surface.Temperature))
east_thresh <- filter(east_thresh, mean >= 32)





#plot
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
cols_bleaching <- brewer.pal(3, "Dark2")
cols_bleaching <- c(cols_bleaching[1], cols_bleaching[1], cols_bleaching[2], cols_bleaching[2],
                    cols_bleaching[3], cols_bleaching[3])
#average across region
all_sst <- all_sst %>% group_by(Region, date) %>% summarise(sst = mean(Sea.Surface.Temperature))

# Define the dates for the vertical lines
vlines <- as.Date(c("2021-05-11", "2021-05-18", "2021-09-04", "2021-09-19", "2022-05-04", "2022-06-09"))

sst_no_na <- filter(all_sst, !is.na(sst))
# Plot the data
all_sst_fig <- ggplot(sst_no_na, aes(x=date, y=sst, color=Region)) +
  geom_line() +
  scale_color_manual(values = c("EastCoast" = "black", "SBN" = "grey")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  geom_vline(all_sst, aes(x=time, y=sst, color=Region), xintercept = vlines, linetype="solid", 
             color = cols_bleaching, size =1)+
  geom_hline(yintercept = c(32,34), linetype = 2)+
  labs(x="",
       y="SST") +
  theme+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
all_sst_fig 

#just sbn
sst_sbn <- filter(sst_no_na, Region == "SBN")
# Define the dates for the vertical lines
vlines_sbn <- as.Date(c("2021-05-18",  "2021-09-19",  "2022-06-09"))
sbn_sst_fig <- ggplot(sst_sbn, aes(x=date, y=sst, color=Region)) +
  geom_line() +
  scale_color_manual(values = c("SBN" = "black")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  geom_vline(sst_sbn, 
             xintercept = vlines_sbn, linetype="solid", 
             color = cols_bleaching[c(1,3,6)], size =1)+
  geom_hline(yintercept = c(34), linetype = 2)+
  labs(x="",
       y="SST") +
  theme+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))
sbn_sst_fig



#just kf
sst_east <- filter(sst_no_na, Region == "EastCoast")
# Define the dates for the vertical lines
vlines_east <- as.Date(c("2021-05-11", "2021-09-04", "2022-05-04"))
east_sst_fig <- ggplot(sst_east, aes(x=date, y=sst, color=Region)) +
  geom_line() +
  scale_color_manual(values = c("EastCoast" = "black")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  geom_vline(sst_sbn, 
             xintercept = vlines_east, linetype="solid", 
             color = cols_bleaching[c(1,3,6)], size =1)+
  geom_hline(yintercept = c(34), linetype = 2)+
  labs(x="",
       y="SST") +
  theme+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))
east_sst_fig


pdf(file = "Output/Figures/aquamodis_fig.pdf",   
    width = 15, # The width of the plot in inches
    height = 10) # The height of the plot in inches
all_sst_fig 
dev.off()

pdf(file = "Output/Figures/aquamodis_sbn_fig.pdf",   
    width = 15, 
    height = 10) 
sbn_sst_fig
dev.off()

pdf(file = "Output/Figures/aquamodis_east_fig.pdf",  
    width = 15, 
    height = 10) 
east_sst_fig
dev.off()



