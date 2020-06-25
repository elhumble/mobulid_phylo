# COI tree

library(phangorn)
library(ape)
library(ggtree)
library(ggplot2)
library(tidyr)
library(dplyr)
library(devtools)
#install_github("GuangchuangYu/ggtree")
library(ggtree)
library(treeio)
library(ggtext)

#~~~~~~~~~~~~~~~~~~~~~~~~~#
#        1K SNPs          #
#~~~~~~~~~~~~~~~~~~~~~~~~~#

tree <- read.newick("data/COI.tree")
tree$tip.label

p <- ggtree(tree) + geom_text(aes(label=node), size = 3, hjust = 1, vjust = 1)
p

meta <- read.csv("data/metadata_original.csv", header = T)

rownames(meta) <- meta$taxa

levels(meta$species)

sp <- levels(meta$species)

levels(meta$region)
meta$region <- factor(meta$region, levels(meta$region)[c(1,2,5,3,4)])
loc <- levels(meta$region)

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


p <- ggtree(tree, color = "black", size = 0.4) %<+% meta + 
  #geom_tiplab(size = 2) +
  geom_tippoint(aes(x = x + 0.0001, fill = I(species), shape = region),  # colour species
                alpha = 0.6, size = 2.5, colour = "black") +
  #geom_tippoint(aes(x= x + 0.00001, fill = I(species), shape = region),  # colour species
  #              alpha = 0.6, size = 2.5, colour = "black") +
  #scale_colour_manual(breaks = sp, values = pal) +
  scale_fill_manual(breaks = sp, 
                    values = pal,
                    labels = c("*M. alfredi*", 
                               "*M. birostris*", 
                               "*M. eregoodoo*",
                               "*M. hypostoma*", 
                               "*M. mobular cf. japanica*",
                               "*M. kuhlii*",
                               "*M. mobular*", 
                               "*M. munkiana*", 
                               "*M. tarapacana*", 
                               "*M. thurstoni*", 
                               "*Undescribed species*", 
                               "*R. bonasus*"),
                    name = "Species", guide = "none") +
  scale_shape_manual(values = c(21, 23, 22, 24, 25),
                     breaks = loc,
                     labels = c("Atlantic", "Indian", "Eastern Indian Ocean", 
                                "Western Indian Ocean","Pacific"),
                     name = "Location") +
  theme(legend.position="right",
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

d <- p$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 85,]

p <- p + geom_text(data = d, aes(x = branch, label=label), 
                   size = 3, vjust=-.5, color='grey10')
#p <- p + geom_nodepoint(data = d, size = 3, shape = 21, fill = "white")

p

ggsave("figs/Figure_S3.png", p, height = 6, width = 9)



