# Script to plot sample map

library(rgdal)
library(ggplot2)
library(rgeos)
library(gpclib)
library(maptools)
library(dplyr)
library(stringr)
options(scipen=999)
library(ggmap)
library(ggthemr)
source("scripts/theme_emily.R")
library(ggtext)

# Read in map data 

gpclibPermit()
world.map <- readOGR(dsn="data/map/", layer="TM_WORLD_BORDERS_SIMPL-0.3")
world.ggmap <- fortify(world.map, region = "NAME")

n <- length(unique(world.ggmap$id))

sample_sites <- read.csv("data/map/map_info.csv") %>%
  mutate(Species = gsub(". ", "_", Species)) %>%
  filter(Species != "R_bonasus" & Species != "A_narinari")

colnames(sample_sites) <- c("sample", "species", "broad_location",
                            "country", "site", "lat", "long")

# Summarise number of samples per site
sample_sites <- sample_sites %>%
  dplyr::group_by(species, country, lat, long) %>%
  dplyr::summarise(n = n()) %>%
  mutate(lon2 = long,
         lat2 = lat)


# Remove Antarctica
world.ggmap <- world.ggmap %>%
  filter(id != "Antarctica")


# Combine map and sample data
df <- full_join(world.ggmap, sample_sites, by = c("long", "lat")) %>%
  mutate(long = ifelse(!is.na(lon2), NA, long),
         lat = ifelse(!is.na(lat2), NA, lat),
         species = as.factor(species))


levels(df$species)

# Specify names and colours

sp <- c("M_alfredi", "M_birostris", "M_eregoodootenkee","M_hypostoma",
        "M_japanica", "M_kuhlii", "M_mobular", "M_munkiana", "M_rochebrunei",
        "M_tarapacana", "M_thurstoni")

pal <- c("#e41a1c", # alfredi red
         "#0c2c84", # birostris blue
         "#ae017e", # eregoo pink
         "goldenrod",   # hypostoma yellow
         "#984ea3", # japanica purple
         "#1d91c0", # kuhlii light blue
         "#35978f", # mobular green
         "#a65628", # munkiana brown
         "black", # rochebrunei black
         "lightseagreen", # tarapacana turquoise
         "#1b7837") # thurstoni dark green

# Plot map

m <- ggplot(df, aes(map_id = id)) +
  geom_map(colour = "grey99", fill = "grey80", size = 0.1, map = df) +
  expand_limits(x = df$long, y = df$lat) +
  #geom_point(aes(size = n, color = species), 
  #           x = df$lon2, y = df$lat2, alpha = 0.6) +
  geom_point(aes(x = lon2, y = lat2, color = species, size = n), 
             data = df, position = position_jitter(w = 3, h = 2, seed = 100), alpha = 0.7) + # 3, 2
  #scale_size_continuous(range = c(1, 10)) +
  scale_size_area(max_size = 8,
                  breaks=c(1,5,10,15,20),
                  name = "Number of samples") +
  scale_color_manual(breaks = sp, 
                     values = pal,
                     name = "Species",
                     labels = c("*M. alfredi* (*n* = 18)",
                                "*M. birostris* (*n* = 23)",
                                "*M. eregoodootenkee* (*n* = 5)",
                                "*M. hypostoma* (*n* = 14)",
                                "*M. japanica* (*n* = 17)",
                                "*M. kuhlii* (*n* = 8)",
                                "*M. mobular* (*n* = 5)",
                                "*M. munkiana* (*n* = 12)",
                                "*M. rochebrunei* (*n* = 1)",
                                "*M. tarapacana* (*n* = 4)",
                                "*M. thurstoni* (*n* = 9)")) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.line = element_line(colour = "white")) +
  theme_emily() +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.text = element_markdown(size = 11)) +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.8)))


m

ggsave("figs/Figure_1.png", m, width = 11, height = 5)
ggsave("figs/Figure_1.tiff", m, width = 11, height = 5)

# Sample numbers

sample_sites %>%
  group_by(species) %>%
  summarise(n = sum(n))


