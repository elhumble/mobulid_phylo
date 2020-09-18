# Code to plot individual based phylogeny and model results

library(phangorn)
library(ape)
#source("https://bioconductor.org/biocLite.R")
#biocLite("ggtree")
#library(ggtree)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggtext)
# bioconductor ggtree stopped working
library(devtools)
#install_github("GuangchuangYu/ggtree")
library(ggtree)
library(treeio)
library(forcats)

#~~~~~~~~~~~~~~~~~~~~~~~~~#
#        7K SNPs          #
#~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Read in tree

tree <- read.newick("data/p10_individual_long.tree")
tree$tip.label

#~~ Read in BFD model matrix

model <- read.csv("data/BFD_raw_wide.csv", header = T, colClasses = "character")
model$individual

rownames(model) <- model$individual
model <- model %>%
  select(-individual)


#~~ Plot basic tree

p <- ggtree(tree) + geom_text(aes(label=node), size = 3, hjust = 1, vjust = 1)
p

#~~ Colour code dummy individuals

node_df <- data.frame(1:nrow(tree$edge)) %>%
  mutate(colour = case_when(X1.nrow.tree.edge. == 106 ~ "white",
                            X1.nrow.tree.edge. == 107 ~ "white",
                            X1.nrow.tree.edge. == 108 ~ "white",
                            X1.nrow.tree.edge. == 109 ~ "white",
                            X1.nrow.tree.edge. == 110 ~ "white",
                            X1.nrow.tree.edge. == 111 ~ "white",
                            X1.nrow.tree.edge. == 82:87 ~ "white",
                            X1.nrow.tree.edge. == 31:36 ~ "white",
                            TRUE ~ as.character("gray55")))

#colour <- c(rep("red", 1), rep("gray55", nrow(node_df)-1))
#node_df$colour <- colour
colnames(node_df) <- c("node", "colour")


#~~ Read in location metadata

meta <- read.csv("data/metadata_original.csv", header = T)

sp <- levels(meta$species)

levels(meta$region)
meta$region <- factor(meta$region, levels(meta$region)[c(1,2,5,3,4)])
loc <- levels(meta$region)

#~~ Specific colour palette

pal <- c("#e41a1c", # alfredi
         "#046C9A", # birostris
         "#ae017e", # eregoo pink
         "goldenrod",   # hypostoma yellow
         "#984ea3", # japanica purple
         "#1d91c0", # kuhlii light blue
         "#35978f", # mobular dark turquoise
         "#a65628", # munkiana brown
         "lightseagreen", # tarapacana turquoise
         "#1b7837", # thurstoni dark green
         "#ff7f00", # new species orange
         "#999999") # white

#~~ Plot tree with colours and regional info

p <- ggtree(tree, color = "black", size = 0.4) %<+% meta + 
  #geom_tiplab(size = 2) +
  geom_tippoint(aes(x= x + 0.00001, fill = I(species), shape = region),  # colour species
                alpha = 0.6, size = 2.5, colour = "black") +
  #scale_colour_manual(breaks = sp, values = pal) +
  scale_fill_manual(breaks = sp, 
                    values = pal,
                    name = "phylo", guide = "none") +
  scale_shape_manual(values = c(21, 23, 22, 24, 25),
                     breaks = loc,
                     labels = c("Atlantic", "Indian", "Pacific", 
                                "Eastern Indian Ocean", "Western Indian Ocean")) +
  theme(legend.position="left") +
  guides(fill = guide_legend(override.aes=list(shape=21)))

#~~ Flip groupings

p <- flip(p, 150, 146) # 141 and 142
p <- flip(p, 151, 152)
p <- flip(p, 155, 156) # 141 and 142

#p <- flip(p, 152, 153) # 1

#~~ Identify nodes with >95% support

d <- p$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 95,]

#~~ Bootstrap support

#p <- p + geom_text(data = d, aes(x = branch, label=label), size = 3, vjust=-.5, color='grey10')
p <- p + geom_nodepoint(data = d, size = 3, shape = 21, fill = "white")

p
#p <- p + geom_tiplab(size=2)

#~~ Add BFD models

# heat map palette:
heat <- c("#B40F20", # alfredi *
          "#9ecae1", # light blue birostris
          "goldenrod", # gold hypostoma 1 *
          "#79402E", # brown 6*
          "#fb9a99", # peach alfredi 12
          "#046C9A",  # birostris blue, 11,*
          "#DD8D29", # orange new species 5*
          "#984ea3", # light blue 2*
          "#35978f", # mobular light green 8*
          "#1b7837", # green thurstoni 10*
          "#ae017e", # eregoo pink 9*
          "#1d91c0", # kuhlii blue 7 *
          "#B40F20", # red*
          "#046C9A",  # birostris blue,*
          "#ae017e", # eregoo pink*
          "goldenrod", # gold*
          "#984ea3", # purple*
          "#1d91c0", # kuhlii light blue*
          "#35978f", # mobular light green *
          "#79402E", # brown*
          "lightseagreen", # peach*
          "#1b7837", # green*
          "#DD8D29", # orange*
          "#999999") # grey *

colnames(model) <- c("1", "2", "3", "4", "5")

#~~ Column 1
q <- gheatmap(p, model[1], width = 0.05, 
              colnames_position = "top", colnames_offset_y = 1, 
              font.size = 3) +
  scale_fill_manual(values=heat, name = "model") +
  theme(legend.position = "none")

#~~ Column 2
r <- gheatmap(q, model[2], width = 0.05, offset = 0.00008, 
              colnames_position = "top", colnames_offset_y = 1, 
              font.size = 3) +
  scale_fill_manual(values=heat)
#theme(legend.position = "none")

#~~ Column 3
s <- gheatmap(r, model[3], width = 0.05, offset = 0.00016, 
              colnames_position = "top", colnames_offset_y = 1, 
              font.size = 3) +
  scale_fill_manual(values=heat)
#theme(legend.position = "none")

#~~ Column 4
t <- gheatmap(s, model[4], width = 0.05, offset = 0.00024, 
              colnames_position = "top", colnames_offset_y = 1, 
              font.size = 3) +
  scale_fill_manual(values=heat)
#theme(legend.position = "none")

#~~ Column 5
u <- gheatmap(t, model[5], width = 0.05, offset = 0.00032, 
              colnames_position = "top", colnames_offset_y = 1, 
              font.size = 3) +
  scale_fill_manual(values=heat) +
  theme(legend.position = "none")

u


ggsave("figs/Figure_2.tiff", height=9, width=9)


#~~~~~~~~~~~~~~~~~~~~~~~~~#
#        1K SNPs          #
#~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Read in tree
tree <- read.newick("data/p90_individual_long.tree")
tree$tip.label

#~~ Read in BFD model matrix

model <- read.csv("data/BFD_raw_wide.csv", header = T, colClasses = "character")
model$individual
rownames(model) <- model$individual
model <- model %>%
  select(-individual)

#~~ Plot basic tree
p <- ggtree(tree) + geom_text(aes(label=node), size = 3, hjust = 1, vjust = 1)
p

#~~ Colour code dummy individuals
node_df <- data.frame(1:nrow(tree$edge)) %>%
  mutate(colour = case_when(X1.nrow.tree.edge. == 106 ~ "white",
                            X1.nrow.tree.edge. == 107 ~ "white",
                            X1.nrow.tree.edge. == 108 ~ "white",
                            X1.nrow.tree.edge. == 109 ~ "white",
                            X1.nrow.tree.edge. == 110 ~ "white",
                            X1.nrow.tree.edge. == 111 ~ "white",
                            X1.nrow.tree.edge. == 82:87 ~ "white",
                            X1.nrow.tree.edge. == 31:36 ~ "white",
                            TRUE ~ as.character("gray55")))

#colour <- c(rep("red", 1), rep("gray55", nrow(node_df)-1))
#node_df$colour <- colour
colnames(node_df) <- c("node", "colour")

#~~ Read in location metadata

meta <- read.csv("data/metadata_original.csv", header = T)
levels(meta$species)
meta$species <- fct_relevel(meta$species, "Manta_birostris", "Manta_alfredi",
                            "Putative_manta", "Mobula_japanica",
                            "Mobula_mobular", "Mobula_kuhlii",
                            "Mobula_eregoodootenkee", "Mobula_thurstoni",
                            "Mobula_hypostoma", "Mobula_munkiana",
                            "Mobula_tarapacana", "Rhinopterabonasus")
levels(meta$species)
sp <- levels(meta$species)

levels(sp)

levels(meta$region)
meta$region <- factor(meta$region, levels(meta$region)[c(1,2,5,3,4)])
loc <- levels(meta$region)

#~~ Specify colour palette

pal <- c("#046C9A", # birostris
         "#e41a1c", # alfredi
         "#ff7f00", # new species orange
         "#984ea3", # japanica purple
         "#35978f", # mobular dark turquoise
         "#1d91c0", # kuhlii light blue
         "#ae017e", # eregoo pink
         "#1b7837", # thurstoni dark green
         "goldenrod",   # hypostoma yellow
         "#a65628", # munkiana brown
         "lightseagreen", # tarapacana turquoise
         "#999999") # grey

ggtree(tree) + geom_tiplab(size = 3)

#~~ Plot tree with colours and regional info

p <- ggtree(tree, color = "black", size = 0.4) %<+% meta + 
  #geom_tiplab(size = 2) +
  geom_tippoint(aes(x= x + 0.00001, fill = I(species), shape = region),  # colour species
                alpha = 0.6, size = 2.5, colour = "black") +
  #scale_colour_manual(breaks = sp, values = pal) +
  scale_fill_manual(breaks = sp, 
                    values = pal,
                    labels = c("*M. birostris*",
                               "*M. alfredi*",
                               "Putative new species",
                               "*M. mobular cf. japanica*",
                               "*M. mobular*",
                               "*M. kuhlii*",
                               "*M. eregoodoo*",
                               "*M. thurstoni*",
                               "*M. hypostoma*",
                               "*M. munkiana*",
                               "*M. tarapacana*",
                               "*R. bonasus*"),
                    name = "Species") +
  scale_shape_manual(values = c(21, 23, 22, 24, 25),
                     breaks = loc,
                     labels = c("Atlantic", "Indian", "Pacific", 
                                "Eastern Indian Ocean", "Western Indian Ocean"),
                     name = "Location") +
  theme(legend.position="left",
        legend.text = element_markdown(size=13),
        legend.title = element_text(size=13)) +
  scale_y_discrete(expand = expand_scale(add = c(2, 2))) +
  guides(fill = guide_legend(override.aes=list(shape=21, size = 3, alpha = 0.8),
                             keyheight = 0.2,
                             default.unit = "inch"),
         shape = guide_legend(override.aes = list(size = 3),
                              keyheight = 0.2,
                              default.unit = "inch"))

p

#~~ Identify nodes with >95% support

d <- p$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 95,]

#~~ Bootstrap support
#p <- p + geom_text(data = d, aes(x = branch, label=label), size = 3, vjust=-.5, color='grey10')
p <- p + geom_nodepoint(data = d, size = 3, shape = 21, fill = "white")

p


ggsave("figs/Figure_S1.png", height = 7, width = 9)
ggsave("figs/Figure_S1.tiff", height = 7, width = 9)

