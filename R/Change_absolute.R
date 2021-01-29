### Long-term changes in size class abundances on the Great Barrier Reef
### Long-term shifts in size-class abundances of communities and taxa
### Author: Andreas Dietzel
### Last modified: June 16, 2020

# This script calculates:
# - Changes in the size frequency distributions of coral communities and taxa
# - Changes in size-class abundances

rm(list = ls())

# load packages and data
library(tidyverse)
library(moments)
library(tidybayes)

LIT <- read.csv("data/Intercept_data.csv", strip.white = T) %>%
  mutate(SecHab = paste0(Sector, Habitat),
         logIntercept = log(Intercept),
         Time = ifelse(Year %in% c(1995, 1996), "Historic", "Recent"),
         Sector = as.factor(Sector))

NrTransects <- read.csv("data/Number_of_transects.csv", strip.white = T) %>%
  mutate(Time = ifelse(Year %in% c(1995, 1996), "Historic", "Recent"),
         SecHab = paste0(Sector, Habitat))

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



###############################################
### CHANGES IN SIZE-FREQUENCY DISTRIBUTIONS ###

# Calculate and plot the size-frequency distributions of communities in different sectors

plotList <- list()

# For loop to produce individual plots, necessary to account for different number of transects
for (i in 1:length(unique(LIT$SecHab))) {

  print(i)

  SecHab <- unique(LIT$SecHab)[i]
  print(SecHab)
  data <- LIT[LIT$SecHab == SecHab, ]
  NrT.Hist <- NrTransects$NrTrans[NrTransects$SecHab == SecHab &
                                    NrTransects$Time == "Historic"]
  print(NrT.Hist)
  NrT.Rec <- NrTransects$NrTrans[NrTransects$SecHab == SecHab &
                                 NrTransects$Time == "Recent"]
  print(NrT.Rec)

  plotList[[i]] <- ggplot() +
    geom_density(data = data[data$Time == "Historic",], fill = "#F8766D",
                 aes(x = Intercept, y = (..count..)/NrT.Hist),
                 alpha = .5, colour = NA) +
    geom_density(data = data[data$Time == "Recent",], fill = "#00BFC4",
                 aes(x = Intercept, y = (..count..)/NrT.Rec),
                 alpha = .5, colour = NA) +
    geom_line(data = data[data$Time == "Historic",],
                 aes(x = Intercept, y = (..count..)/NrT.Hist),
                 colour = "#F8766D", stat = "density") +
    geom_line(data = data[data$Time == "Recent",],
                 aes(x = Intercept, y = (..count..)/NrT.Rec),
                 colour = "#00BFC4", stat = "density") +
    scale_x_continuous(trans = "log10", limits = c(0.9,300)) +
    coord_cartesian(ylim = c(0,40)) +
    theme_minimal() +
    theme(text = element_text(colour = "black"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(),
          axis.ticks = element_line()) +
    labs(x = "", y = "")
}

names(plotList) <- unique(data$SecHab)

# combine individual plots
cowplot::plot_grid(plotList[[1]], plotList[[6]] + facet_grid(Sector ~.),
                   plotList[[2]], plotList[[7]] + facet_grid(Sector ~.),
                   plotList[[3]] + labs(y = "Intercepts per transect"),
                   plotList[[8]] + facet_grid(Sector ~.),
                   plotList[[4]], plotList[[9]] + facet_grid(Sector ~.),
                   plotList[[5]], plotList[[10]] + facet_grid(Sector ~.),
                   ncol = 2, rel_widths = c(.9,1))

# save plot
ggsave("figures/CountDensity_Sectors.pdf", width = 5, height = 7)

# Change in size-frequency distributions pooled across all sectors

data <- LIT[LIT$Habitat == "Crest",] %>% mutate(Sector = "GBR wide")

# Calculate number of transects on crest for historic and recent surveys
NrT.Hist <- sum(NrTransects$NrTrans[NrTransects$Habitat == "Crest" &
                                  NrTransects$Time == "Historic"])
NrT.Rec <- sum(NrTransects$NrTrans[NrTransects$Habitat == "Crest" &
                               NrTransects$Time == "Recent"])

GBR_crest <- ggplot() +
  geom_density(data = data[data$Time == "Historic",], fill = "#F8766D",
               aes(x = Intercept, y = (..count..)/NrT.Hist),
               alpha = .5, colour = NA) +
  geom_density(data = data[data$Time == "Recent",], fill = "#00BFC4",
               aes(x = Intercept, y = (..count..)/NrT.Rec),
               alpha = .5, colour = NA) +
  geom_line(data = data[data$Time == "Historic",],
            aes(x = Intercept, y = (..count..)/NrT.Hist),
            colour = "#F8766D", stat = "density") +
  geom_line(data = data[data$Time == "Recent",],
            aes(x = Intercept, y = (..count..)/NrT.Rec),
            colour = "#00BFC4", stat = "density") +
  scale_x_continuous(trans = "log10", limits = c(0.9,300)) +
  coord_cartesian(ylim = c(0,40)) +
  facet_grid(. ~ Habitat) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  theme_minimal() +
  theme(text = element_text(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line()) +
  labs(x = "log colony size", y = "")

# Same for slope pooled across all sectors

data <- LIT[LIT$Habitat == "Slope",] %>% mutate(Sector = "GBR wide")

NrT.Hist <- sum(NrTransects$NrTrans[NrTransects$Habitat == "Slope" &
                                  NrTransects$Time == "Historic"])
NrT.Rec <- sum(NrTransects$NrTrans[NrTransects$Habitat == "Slope" &
                               NrTransects$Time == "Recent"])

GBR_slope <- ggplot() +
  geom_density(data = data[data$Time == "Historic",], fill = "#F8766D",
               aes(x = Intercept, y = (..count..)/NrT.Hist),
               alpha = .5, colour = NA) +
  geom_density(data = data[data$Time == "Recent",], fill = "#00BFC4",
               aes(x = Intercept, y = (..count..)/NrT.Rec),
               alpha = .5, colour = NA) +
  geom_line(data = data[data$Time == "Historic",],
            aes(x = Intercept, y = (..count..)/NrT.Hist),
            colour = "#F8766D", stat = "density") +
  geom_line(data = data[data$Time == "Recent",],
            aes(x = Intercept, y = (..count..)/NrT.Rec),
            colour = "#00BFC4", stat = "density") +
  scale_x_continuous(trans = "log10", limits = c(0.9,300)) +
  coord_cartesian(ylim = c(0,40)) +
  facet_grid(Sector ~ Habitat) +
  theme_minimal() +
  theme(text = element_text(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line()) +
  labs(x = "log colony size (cm)", y = "")

# combine individual plots
cowplot::plot_grid(GBR_crest + labs(y = "Intercepts per transect"),
                   GBR_slope, nrow = 1, rel_widths = c(.9,1))

# save plot
ggsave("figures/CountDensity_GBR.pdf", width = 5, height = 1.7)



# Calculate and plot the size-frequency distributions of major taxa

data <- LIT[LIT$Habitat == "Crest", ] %>%
  mutate(Taxa2 = gsub("Acropora", "Acr", Taxa),
         Taxa2 = gsub("_", " ", Taxa2),
         Taxa2 = gsub("Pocillopora", "Poc", Taxa2),
         Taxa2 = gsub("scleractinians", "sclera", Taxa2))

NrT.Hist <- sum(NrTransects$NrTrans[NrTransects$Habitat == "Crest" &
                                      NrTransects$Time == "Historic"])
NrT.Rec <- sum(NrTransects$NrTrans[NrTransects$Habitat == "Crest" &
                                     NrTransects$Time == "Recent"])

(Taxa_crest <- ggplot() +
  geom_density(data = data[data$Time == "Historic",], fill = "#F8766D",
               aes(x = Intercept, y = (..count..)/NrT.Hist),
               alpha = .5, colour = NA) +
  geom_density(data = data[data$Time == "Recent",], fill = "#00BFC4",
               aes(x = Intercept, y = (..count..)/NrT.Rec),
               alpha = .5, colour = NA) +
  geom_line(data = data[data$Time == "Historic",],
            aes(x = Intercept, y = (..count..)/NrT.Hist),
            colour = "#F8766D", stat = "density") +
  geom_line(data = data[data$Time == "Recent",],
            aes(x = Intercept, y = (..count..)/NrT.Rec),
            colour = "#00BFC4", stat = "density") +
  scale_x_continuous(trans = "log10", limits = c(0.9,300)) +
  facet_grid(Taxa2 ~ Habitat, scales = "free_y") +
  theme_minimal() +
  theme(text = element_text(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line()) +
  labs(x = "log colony size (cm)", y = "intercepts per transect"))

# Same for slope

data <- LIT[LIT$Habitat == "Slope", ] %>%
  mutate(Taxa2 = gsub("Acropora", "Acr", Taxa),
         Taxa2 = gsub("_", " ", Taxa2),
         Taxa2 = gsub("Pocillopora", "Poc", Taxa2),
         Taxa2 = gsub("scleractinians", "sclera", Taxa2))

NrT.Hist <- sum(NrTransects$NrTrans[NrTransects$Habitat == "Slope" &
                                      NrTransects$Time == "Historic"])
NrT.Rec <- sum(NrTransects$NrTrans[NrTransects$Habitat == "Slope" &
                                     NrTransects$Time == "Recent"])

(Taxa_slope <- ggplot() +
  geom_density(data = data[data$Time == "Historic",], fill = "#F8766D",
               aes(x = Intercept, y = (..count..)/NrT.Hist),
               alpha = .5, colour = NA) +
  geom_density(data = data[data$Time == "Recent",], fill = "#00BFC4",
               aes(x = Intercept, y = (..count..)/NrT.Rec),
               alpha = .5, colour = NA) +
  geom_line(data = data[data$Time == "Historic",],
            aes(x = Intercept, y = (..count..)/NrT.Hist),
            colour = "#F8766D", stat = "density") +
  geom_line(data = data[data$Time == "Recent",],
            aes(x = Intercept, y = (..count..)/NrT.Rec),
            colour = "#00BFC4", stat = "density") +
  scale_x_continuous(trans = "log10", limits = c(0.9,300)) +
  facet_grid(Taxa2 ~ Habitat, scales = "free_y") +
  theme_minimal() +
  theme(text = element_text(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line()) +
  labs(x = "log colony size (cm)", y = ""))

cowplot::plot_grid(Taxa_crest, Taxa_slope, nrow = 1)

ggsave("figures/CountDensity_Taxa.pdf", width = 6, height = 11)



##################################################
### PERCENTAGE CHANGE IN SIZE-CLASS ABUNDANCES ###

# Percentage change in size-class abundances by taxa
PCchangeAbsolute.TaxHab <- data.frame()

GBRhist.temp <- LIT %>%
  group_by(Taxa, Habitat) %>%
  mutate(Bin = cut_number(Intercept, n = 5,
                   labels = c("small","medium","medium","medium","large"))) %>%
  ungroup()

# Bootstrap uncertainties
for (i in 1:1000) {
  a <- GBRhist.temp %>%
    sample_n(size = nrow(GBRhist.temp), replace = T) %>%
    group_by(Taxa, Habitat, Time, Bin) %>%
    summarise(Freq = n()) %>%
    ungroup() %>%
    complete(Taxa, Habitat, Time, Bin, fill = list(Freq = 0)) %>%
    left_join(NrTransects[NrTransects$Sector==0,], by =
                c("Time"="Time","Habitat"="Habitat")) %>%
    mutate(Freq = Freq/NrTrans) %>%
    dplyr::select(Taxa:Freq) %>%
    spread(key = Time, value = Freq) %>%
    mutate(PCchange = 100*(Recent - Historic)/Historic,
           draw = i)

  PCchangeAbsolute.TaxHab <- bind_rows(PCchangeAbsolute.TaxHab, a)
}

# Calculate highest density intervals
Summ.Tax <- PCchangeAbsolute.TaxHab %>%
  group_by(Taxa, Habitat, Bin) %>%
  median_hdi(PCchange, .width = c(.66,.95)) %>%
  arrange(Taxa, Habitat, Bin) %>%
  ungroup() %>%
  mutate(Taxa = gsub("_"," ", Taxa))

# Create plot
Summ.Tax %>%
  ggplot() +
  facet_grid(Bin ~ Habitat) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_hline(yintercept = -50, linetype = 2, colour = "grey") +
  geom_hline(yintercept = 50, linetype = 2, colour = "grey") +
  geom_point(aes(y = PCchange, x = Taxa), size = 1.5) +
  geom_segment(data = Summ.Tax[Summ.Tax$.width==.95,],
               aes(x = Taxa, xend = Taxa, y = .lower, yend = .upper), size = .5) +
  geom_segment(data = Summ.Tax[Summ.Tax$.width==.66,],
               aes(x = Taxa, xend = Taxa, y = .lower, yend = .upper), size = 1) +
  coord_flip() +
  theme_simple +
  labs(y = "change (%)", x = "")

# Save plot
ggsave("figures/PCchangeAbsolute.TaxHab.pdf", height = 6, width = 5)




# Percentage change in size-class abundances of communities in different sectors

PCchangeAbsolute.SecHab <- data.frame()

GBRhist.temp <- LIT %>%
  group_by(Sector, Habitat) %>%
  mutate(Bin = cut_number(Intercept, n = 5,
                          labels = c("small","medium","medium","medium","large"))) %>%
  ungroup()

Bins <- GBRhist.temp %>%
  group_by(Sector, Habitat, Bin) %>%
  summarise(n = n())

# Bootstrap uncertainties
for (i in 1:1000) {
  a <- GBRhist.temp %>%
    sample_n(size = nrow(GBRhist.temp), replace = T) %>%
    group_by(Sector, Habitat, Time, Bin) %>%
    summarise(Freq = n()) %>%
    ungroup() %>%
    complete(Sector, Habitat, Time, Bin, fill = list(Freq = 0)) %>%
    left_join(NrTransects, by = c("Time"="Time","Habitat"="Habitat",
                                "Sector"="Sector")) %>%
    mutate(Freq = Freq/NrTrans) %>%
    dplyr::select(Sector:Freq) %>%
    spread(key = Time, value = Freq) %>%
    mutate(PCchange = 100*(Recent - Historic)/Historic,
           draw = i)

  PCchangeAbsolute.SecHab <- bind_rows(PCchangeAbsolute.SecHab, a)
}

# Calculate highest density intervals
Summ <- PCchangeAbsolute.SecHab %>%
  group_by(Sector, Habitat, Bin) %>%
  median_hdi(PCchange, .width = c(.66,.95)) %>%
  arrange(Habitat, Sector, Bin) %>%
  ungroup() %>%
  mutate(Sector = ifelse(Sector == 0, "GBR wide",
                         paste("Sector",Sector)))

# Create plot
Summ %>%
  ggplot() +
  facet_grid(Bin ~ Habitat) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_hline(yintercept = -50, linetype = 2, colour = "grey") +
  geom_hline(yintercept = 50, linetype = 2, colour = "grey") +
  geom_point(aes(y = PCchange, x = reorder(Sector, desc(Sector))), size = 1.5) +
  geom_segment(data = Summ[Summ$.width==.95,],
               aes(x = Sector, xend = Sector, y = .lower, yend = .upper), size = .5) +
  geom_segment(data = Summ[Summ$.width==.66,],
               aes(x = Sector, xend = Sector, y = .lower, yend = .upper), size = 1) +
  coord_flip() +
  theme_classic() +
  theme_simple +
  labs(y= "change (%)", x = "")

# Save plot
ggsave("figures/PCchangeAbsolute.SecHab.pdf", height = 3.5, width = 3.5)
