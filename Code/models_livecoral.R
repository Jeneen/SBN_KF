library(car)
library(brms)
library(nlme)
library(sjPlot)
library(easystats)
library(dplyr)
library(gtools)
options(mc.cores = parallel::detectCores())

#load data
dat_all_bl <- readRDS("Data/Processed/df_for_model.rds")

#clean data
dat_all_bl$Live.coral <- dat_all_bl$Live.coral/100
colnames(dat_all_bl) <- gsub("\\[,1\\]", "", colnames(dat_all_bl))
dat_all_bl$Bleaching <- factor(dat_all_bl$Bleaching, levels = c("Pre bleaching", "During bleaching",
                                                                "Post bleaching"))
dat_all_bl$logitLiveCoral <- logit(dat_all_bl$Live.coral)


#drop extra
dat_all_bl <- na.omit(dat_all_bl)

#summary
summary <- dat_all_bl %>% filter(Bleaching == "Pre bleaching") %>% group_by(Site) %>% 
  summarise(mean = mean(live_acrop))


# split into separate models for each region
east_data_all_bl <- dat_all_bl %>% filter(Region == "East coast")
sbn_data_all_bl <- dat_all_bl %>% filter(Region == "SBN")


#East 
# Define the priors
priors <- c(
  prior(student_t(3, 0, 2.5), class = Intercept),  
  prior(normal(0, 10), class = b),      
  prior(gamma(2, 0.1), class = phi) 
)

east_m_livecoral_p <- brm(Live.coral~ Site + Bleaching, data = east_data_all_bl,
                        control = list(adapt_delta = 0.9),  family="beta", prior = priors,
                        iter=5000,  warmup=2000, chains=4, cores = 4, seed = 10)
summary(east_m_livecoral_p)
plot(east_m_livecoral_p)
bayes_R2(east_m_livecoral_p)


# Open a PDF device to save the plot
pdf("Output/Figures/east_livecoral_traceplots.pdf", width = 8, height = 6)
# Plot the trace plots
plot(east_m_livecoral_p)
# Close the PDF device
dev.off()

jpeg("Output/Figures/east_livecoral_traceplots_%d.jpg", width = 1200, height = 800, res = 150)
# Generate the trace plots
plot(east_m_livecoral_p)
# Close the graphics device
dev.off()

jpeg("Output/Figures/east_livecoral_pp.jpg", width = 1200, height = 800, res = 150)
# Generate the trace plots
pp_check(east_m_livecoral_p)
# Close the graphics device
dev.off()


#SBN
# Define the priors
priors <- c(
  prior(student_t(3, 0, 2.5), class = Intercept),  
  prior(normal(0, 10), class = b),      
  prior(gamma(2, 0.1), class = phi) 
)

sbn_m_livecoral_p <- brm(Live.coral ~ Site + Bleaching, data = sbn_data_all_bl,
                       control = list(adapt_delta = 0.9),  family="beta", prior = priors,
                       iter=5000,  warmup=2000, chains=4, cores = 4, seed = 10)
summary(sbn_m_livecoral_p)
plot(sbn_m_livecoral_p)
bayes_R2(sbn_m_livecoral_p)
pp_check(sbn_m_livecoral_p)

# Open a PDF device to save the plot
pdf("Output/Figures/sbn_livecoral_traceplots.pdf", width = 8, height = 6)
# Plot the trace plots
plot(sbn_m_livecoral_p)
# Close the PDF device
dev.off()

jpeg("Output/Figures/sbn_livecoral_traceplots_%d.jpg", width = 1200, height = 800, res = 150)
# Generate the trace plots
plot(sbn_m_livecoral_p)
# Close the graphics device
dev.off()

jpeg("Output/Figures/sbn_livecoral_pp.jpg", width = 1200, height = 800, res = 150)
# Generate the trace plots
pp_check(sbn_m_livecoral_p)
# Close the graphics device
dev.off()
