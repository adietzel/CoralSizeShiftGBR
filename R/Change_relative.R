### Long-term shifts in the mean, standard deviation, and 10th and 90th percentile of colony size
### Author: Andreas Dietzel
### Last modified: June 16, 2020

# This script
# - fits multi-level models to line-intercept transect data of communities and taxa
# - calculates changes in the mean, standard deviation, 10th and 90th percentile of colony size

rm(list = ls())

# load packages and data
library(tidyverse)
library(ggthemes)
library(ggridges)
library(brms)
library(tidybayes)
library(modelr)

LIT <- read.csv("data/Intercept_data.csv", strip.white = T) %>%
  mutate(SecHab = paste0(Sector, Habitat),
         logIntercept = log(Intercept),
         Time = ifelse(Year %in% c(1995, 1996), "Historic", "Recent"),
         Sector = as.factor(Sector))

# custom plot theme
theme_simple <- theme_classic() +
  theme(text       = element_text(color = "black"),
        strip.text = element_text(color = "black"),
        axis.text  = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        line       = element_line(color = "black"),
        plot.background   = element_rect(fill = "white", color = "transparent"),
        panel.background  = element_rect(fill = "white", color = "black"),
        strip.background  = element_rect(fill = "white", color = "transparent"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.background = element_rect(fill = "white", color = "transparent"),
        legend.key        = element_rect(fill = "white", color = "transparent"))





### Changes in the size structure of communities in different sectors and habitats ###

# Fit model

Model_Sector <- brm(bf(logIntercept ~ Time * Sector * Habitat,
               sigma ~ Time *  Sector * Habitat), data = LIT,
            iter = 2000, warmup = 200, thin = 5, chains = 3)

save(Model_Sector, file = "RData/Model_Sector.RData")
load("RData/Model_Sector.RData")

marginal_effects(Model_Sector)

# Model diagnostics
summary(Model_Sector)
plot(Model_Sector)

# Model fit
y_pred <- posterior_predict(Model_Sector)

newdata <- LIT %>%
  dplyr::select(Time, Habitat, Sector, logIntercept) %>%
  cbind(t(y_pred)) %>%
  gather(key = "Rep",
         value = "logIntercept", -Time,-Habitat,-Sector,-logIntercept)

newdata[newdata$Rep %in% seq(1, 360, by = 20), ] %>%
  ggplot() +
  geom_density_ridges(aes(x = exp(logIntercept), y = Time, colour = factor(Rep)),
                      fill = NA, scale = .5, size = .2) +
  geom_density_ridges(data = LIT, aes(x = exp(logIntercept), y = Time),
                      colour = "red", fill = NA, scale = .5) +
  facet_grid(Sector ~ Habitat) +
  scale_x_continuous(trans = "log10", breaks = c(1,10,100,1000)) +
  scale_colour_manual(breaks = unique(newdata2$Rep),
                      values = rep("#0000001A", length(unique(newdata2$Rep)))) +
  labs(y = "", x = "log intercept length (cm)") +
  theme_classic() +
  theme(legend.position = "")

# Residual plot
resid = resid(Model_Sector)[,'Estimate']
fit = fitted(Model_Sector)[,'Estimate']
ggplot() + geom_point(data = NULL, aes(y = resid, x = fit))

# Extract posterior draws
grid <- LIT %>%
  data_grid(Time, Sector, Habitat) %>%
  add_fitted_draws(Model_Sector, dpar = TRUE) %>%
  mutate(p90 = qnorm(p = .9, mean = mu, sd = sigma),
         p10 = qnorm(p = .1, mean = mu, sd = sigma),
         CV = sigma/mu * 100)

# Calculate percent change in mu, sig, 10th and 90th percentile
mu.long <- grid %>%
  ungroup() %>%
  mutate(Spread = c(1:(0.5*nrow(grid)), 1:(0.5*nrow(grid)))) %>%
  dplyr::select(Time, Habitat, Sector, mu, Spread, .draw) %>%
  spread(key = Time, value = mu) %>%
  mutate(PC.change.mu = 100 * (exp(Recent)-exp(Historic))/exp(Historic)) %>%
  dplyr::select(-(Spread:Recent))

sig.long <- grid %>%
  ungroup() %>%
  mutate(Spread = c(1:(0.5*nrow(grid)), 1:(0.5*nrow(grid)))) %>%
  dplyr::select(Time, Habitat, Sector, sigma, Spread, .draw) %>%
  spread(key = Time, value = sigma) %>%
  mutate(PC.change.sig = 100 * ((Recent)-(Historic))/(Historic)) %>%
  dplyr::select(-(Spread:Recent))

CV.long <- grid %>%
  ungroup() %>%
  mutate(Spread = c(1:(0.5*nrow(grid)), 1:(0.5*nrow(grid)))) %>%
  dplyr::select(Time, Habitat, Sector, CV, Spread, .draw) %>%
  spread(key = Time, value = CV) %>%
  mutate(PC.change.CV = 100 * ((Recent)-(Historic))/(Historic)) %>%
  dplyr::select(-(Spread:Recent))

p90.long <- grid %>%
  ungroup() %>%
  mutate(Spread = c(1:(0.5*nrow(grid)), 1:(0.5*nrow(grid)))) %>%
  dplyr::select(Time, Habitat, Sector, p90, Spread, .draw) %>%
  spread(key = Time, value = p90) %>%
  mutate(PC.change.p90 = 100 * (exp(Recent)-exp(Historic))/exp(Historic)) %>%
  dplyr::select(-(Spread:Recent))

p10.long <- grid %>%
  ungroup() %>%
  mutate(Spread = c(1:(0.5*nrow(grid)), 1:(0.5*nrow(grid)))) %>%
  dplyr::select(Time, Habitat, Sector, p10, Spread, .draw) %>%
  spread(key = Time, value = p10) %>%
  mutate(PC.change.p10 = 100 * (exp(Recent)-exp(Historic))/exp(Historic)) %>%
  dplyr::select(-(Spread:Recent))

# Combine all metrics
PCchange.all <- cbind(mu.long, sig.long[,3], CV.long[,3],
                      p10.long[,3], p90.long[,3]) %>%
  gather(key, value, 3:7) %>%
  mutate(key = substr(key, 11, 15))

# Change order and names of factor levels
PCchange.all$Sector <- factor(PCchange.all$Sector,
                              levels = c(5,4,3,2,1))

PCchange.all$key <- factor(PCchange.all$key,
                              levels = c("mu","sig","CV","p10","p90"))

levels(PCchange.all$key) <- c("mean","sigma","CV",
                              "10th percentile","90th percentile")

# Summarise and save as table

Summary.HDI.SecHab <- PCchange.all %>%
  group_by(key, Habitat, Sector) %>%
  median_hdi(value, .width = .95) %>%
  arrange(key, Habitat, Sector) %>%
  write.csv("figures/Summary.HDI.SecHab.csv", row.names = F)





### Changes in size structure of major taxa in different habitats ###

# Fit model

Model_Taxa <- brm(bf(logIntercept ~ Time * Habitat * Taxa,
                     sigma ~ Time * Habitat * Taxa),
                  data = LIT, iter = 2000, warmup = 200,
                  thin = 5, chains = 3)

save(Model_Taxa, file = "RData/Model_Taxa.RData")

marginal_effects(Model_Taxa)

# Model diagnostics
summary(Model_Taxa)
plot(Model_Taxa)

# Examine model fit
y_pred <- posterior_predict(Model_Taxa)

newdata <- LIT %>%
  dplyr::select(Time, Habitat, Taxa, logIntercept) %>%
  cbind(t(y_pred)) %>%
  gather(key = "Rep",
         value = "logIntercept", -Time,-Habitat,-Taxa,-logIntercept)

newdata2 <- newdata[newdata$Rep %in% seq(1, 360, by = 20), ] %>%
  mutate(TaxaID = ifelse(Taxa %in% c("Isopora","Montipora","Poritidae","Faviidae",
                                     "Mussidae","Seriatopora"),1,2))
LIT2 <- LIT %>%
  mutate(TaxaID = ifelse(Taxa %in% c("Isopora","Montipora","Poritidae","Faviidae",
                                     "Mussidae","Seriatopora"),1,2))

# Model fit
ggplot(newdata2[newdata2$TaxaID == 1,]) +
  geom_density_ridges(aes(x = exp(logIntercept), y = Time, colour = factor(Rep)),
                      fill = NA, scale = .5, size = .2) +
  geom_density_ridges(data = LIT2[LIT2$TaxaID==1,],
                      aes(x = exp(logIntercept), y = Time),
                      colour = "red", fill = NA, scale = .5) +
  facet_grid(Taxa ~ Habitat) +
  scale_x_continuous(trans = "log10", breaks = c(1,10,100,1000)) +
  scale_colour_manual(breaks = unique(newdata2$Rep[newdata2$TaxaID==1]),
                      values = rep("#0000001A",
                                   length(unique(newdata2$Rep[newdata2$TaxaID==1])))) +
  labs(y = "", x = "log intercept length (cm)") +
  theme_classic() +
  theme(legend.position = "")

ggplot(newdata2[newdata2$TaxaID == 2,]) +
  geom_density_ridges(aes(x = exp(logIntercept), y = Time, colour = factor(Rep)),
                      fill = NA, scale = .5, size = .2) +
  geom_density_ridges(data = LIT2[LIT2$TaxaID==2,],
                      aes(x = exp(logIntercept), y = Time),
                      colour = "red", fill = NA, scale = .5) +
  facet_grid(Taxa ~ Habitat) +
  scale_x_continuous(trans = "log10", breaks = c(1,10,100,1000)) +
  scale_colour_manual(breaks = unique(newdata2$Rep[newdata2$TaxaID==2]),
                      values = rep("#0000001A",
                                   length(unique(newdata2$Rep[newdata2$TaxaID==2])))) +
  labs(y = "", x = "log intercept length (cm)") +
  theme_classic() +
  theme(legend.position = "")

ggsave("figures/ModelFit.Ridges.Taxa.pdf", width = 4, height = 7)

# Residual plot
resid = resid(Model_Taxa)[,'Estimate']
fit = fitted(Model_Taxa)[,'Estimate']
ggplot() + geom_point(data = NULL, aes(y = resid, x = fit))

ggsave("figures/ResidualPlot_brm.Taxa.pdf", width = 4, height = 4)

# Extract posterior draws
grid.Taxa <- LIT %>%
  data_grid(Time, Habitat, Taxa) %>%
  add_fitted_draws(Model_Taxa, dpar = TRUE) %>%
  mutate(p90 = qnorm(p = .9, mean = mu, sd = sigma),
         p10 = qnorm(p = .1, mean = mu, sd = sigma),
         CV = sigma/mu * 100)

grid.Taxa.summ <- grid.Taxa %>%
  group_by(Time, Habitat, Taxa) %>%
  summarise(mean = median(mu), sigma = median(sigma)) %>%
  write.csv("figures/Summary.TaxSec.MuSig.csv", row.names = F)

# Calculate percent change in mu, sig, 10th and 90th percentiles

mu.long <- grid.Taxa %>%
  ungroup() %>%
  mutate(Spread = c(1:(0.5*nrow(grid.Taxa)), 1:(0.5*nrow(grid.Taxa)))) %>%
  dplyr::select(Time, Habitat, Taxa, mu, Spread, .draw) %>%
  spread(key = Time, value = mu) %>%
  mutate(PC.change.mu = 100 * (exp(Recent)-exp(Historic))/exp(Historic)) %>%
  dplyr::select(-(Spread:Recent))

sig.long <- grid.Taxa %>%
  ungroup() %>%
  mutate(Spread = c(1:(0.5*nrow(grid.Taxa)), 1:(0.5*nrow(grid.Taxa)))) %>%
  dplyr::select(Time, Habitat, Taxa, sigma, Spread, .draw) %>%
  spread(key = Time, value = sigma) %>%
  mutate(PC.change.sig = 100 * ((Recent)-(Historic))/(Historic)) %>%
  dplyr::select(-(Spread:Recent))

p90.long <- grid.Taxa %>%
  ungroup() %>%
  mutate(Spread = c(1:(0.5*nrow(grid.Taxa)), 1:(0.5*nrow(grid.Taxa)))) %>%
  dplyr::select(Time, Habitat, Taxa, p90, Spread, .draw) %>%
  spread(key = Time, value = p90) %>%
  mutate(PC.change.p90 = 100 * (exp(Recent)-exp(Historic))/exp(Historic)) %>%
  dplyr::select(-(Spread:Recent))

p10.long <- grid.Taxa %>%
  ungroup() %>%
  mutate(Spread = c(1:(0.5*nrow(grid.Taxa)), 1:(0.5*nrow(grid.Taxa)))) %>%
  dplyr::select(Time, Habitat, Taxa, p10, Spread, .draw) %>%
  spread(key = Time, value = p10) %>%
  mutate(PC.change.p10 = 100 * (exp(Recent)-exp(Historic))/exp(Historic)) %>%
  dplyr::select(-(Spread:Recent))

CV.long <- grid.Taxa %>%
  ungroup() %>%
  mutate(Spread = c(1:(0.5*nrow(grid.Taxa)), 1:(0.5*nrow(grid.Taxa)))) %>%
  dplyr::select(Time, Habitat, Taxa, CV, Spread, .draw) %>%
  spread(key = Time, value = CV) %>%
  mutate(PC.change.CV = 100 * ((Recent)-(Historic))/(Historic)) %>%
  dplyr::select(-(Spread:Recent))

# Combine all metrics
PCchange.all.Tax <- cbind(mu.long, sig.long[,3],CV.long[,3],
                          p10.long[,3],p90.long[,3]) %>%
  gather(key, value, 3:7) %>%
  mutate(key = substr(key, 11, 15))

# Change order and names of factor levels
PCchange.all.Tax$key <- factor(PCchange.all.Tax$key,
                           levels = c("mu","sig","CV","p10","p90"))

levels(PCchange.all.Tax$key) <- c("mean","sigma","CV","10th percentile","90th percentile")

# Summarise and save as table
Summary.PCchange.all.Tax <- PCchange.all.Tax %>%
  group_by(key, Taxa, Habitat) %>%
  median_hdi(value, .width = c(.66,.95)) %>%
  arrange(key, Taxa, Habitat) %>%
  write.csv("figures/Summary.HDI.TaxHab.csv", row.names = F)



#############################################################
### COMBINE COMMUNITY AND TAXA AND PRODUCE SUMMARY FIGURE ###

ChangeTaxSecHab <- PCchange.all %>%
  mutate(Sector = paste("Sector", Sector)) %>%
  rename(Taxa = Sector) %>%
  bind_rows(., PCchange.all.Tax) %>%
  mutate(Taxa = as.factor(gsub("_", " ", Taxa)))

ChangeTaxSecHab$Taxa <- factor(ChangeTaxSecHab$Taxa,
                               levels(ChangeTaxSecHab$Taxa)[c(17:15, 9:1, 14:10)])

Summ <- ChangeTaxSecHab %>%
  group_by(Habitat, Taxa, key) %>%
  median_hdi(value, .width = c(.66, .95))

Summ1 <- Summ %>% filter(key %in% c("mean","sigma"))

MeanSigma <- Summ1 %>%
  ggplot() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_hline(yintercept = -50, linetype = 2, colour = "lightgrey") +
  geom_hline(yintercept = 50, linetype = 2, colour = "lightgrey") +
  geom_point(aes(y = value, x = Taxa), size = 1.5) +
  geom_segment(data = Summ1[Summ1$.width==.95,],
               aes(x = Taxa, xend = Taxa, y = .lower, yend = .upper), size = .5) +
  geom_segment(data = Summ1[Summ1$.width==.66,],
               aes(x = Taxa, xend = Taxa, y = .lower, yend = .upper), size = 1) +
  facet_grid(key ~ Habitat) +
  theme_simple +
  scale_y_continuous(breaks = c(-100,-50,0,50,100)) +
  labs(y = "change (%)", x = "") +
  coord_flip() + theme(legend.position = "")

Summ2 <- Summ %>% filter(key %in% c("10th percentile","90th percentile"))

Percentiles <- Summ2 %>%
  ggplot() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_hline(yintercept = -50, linetype = 2, colour = "lightgrey") +
  geom_hline(yintercept = 50, linetype = 2, colour = "lightgrey") +
  geom_point(aes(y = value, x = Taxa), size = 1.5) +
  geom_segment(data = Summ2[Summ2$.width==.95,],
               aes(x = Taxa, xend = Taxa, y = .lower, yend = .upper), size = .5) +
  geom_segment(data = Summ2[Summ2$.width==.66,],
               aes(x = Taxa, xend = Taxa, y = .lower, yend = .upper), size = 1) +
  facet_grid(key ~ Habitat) +
  theme_simple +
  scale_y_continuous(breaks = c(-100,-50,0,50,100)) +
  labs(y = "change (%)", x = "") +
  coord_flip() + theme(legend.position = "")

cowplot::plot_grid(MeanSigma, Percentiles, ncol = 2)

ggsave("figures/Figure4_RelativeAbundances.pdf", width = 7, height = 5)

