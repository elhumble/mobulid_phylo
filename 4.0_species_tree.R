# Species trees for Supp Figure 6

library(ggtree)
library(ggplot2)
library(ape)
library(treeio)
library(ggtext)

#~~~~~~~~~~~~~~~~~~~~~~~~~#
#        7K SNPs          #
#~~~~~~~~~~~~~~~~~~~~~~~~~#

tr <- read.newick("data/p10_species.tree")
ggtree(tr) + geom_text(aes(label=node), size = 3, hjust = 1, vjust = 1)
# flip 7, 6
tr$tip.label

#~~ Change tip labels

genus <- c("M.", "M.", "M.",
           "M.", "M.", "M.",
           "Undescribed", "M.", "M.",
           "M.", "R.")

species <- c("tarapacana", "thurstoni", "kuhlii",
             "kuhlii", "alfredi", "birostris",
             "species", "mobular", "hypostoma",
             "munkiana", "bonasus")

cf <- c("", "", "",
        "cf.", "", "",
        "", "", "",
        "", "")

synonym <- c("", "", "",
             "eregoodootenkee", "", "",
             "", "", "",
             "", "")

d <- data.frame(label = tr$tip.label, genus = genus,
                species = species)
d


p <- ggtree(tr, color = "black", size = 0.4) %<+% d + xlim(NA, 0.002) +
  geom_tiplab(aes(label=paste0('italic(', genus, ')~italic(', species, ')~italic(', cf, ')~italic(', synonym,')')), parse=T, color = "black")
p

p <- flip(p, 7, 6)


d <- p$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 0,]

p <- p + geom_text(data = d, aes(x = branch, label=label), size = 3, vjust=-.5, color='black')
#p <- p + geom_nodepoint(data = d, size = 3, shape = 21, fill = "white")
p
              
ggsave("figs/Figure_S6A.tiff", p, height = 7, width = 9)


#~~~~~~~~~~~~~~~~~~~~~~~~~#
#        1K SNPs          #
#~~~~~~~~~~~~~~~~~~~~~~~~~#

tr <- read.newick("data/p90_species.tree")
ggtree(tr) + geom_text(aes(label=node), size = 4, hjust = 1, vjust = 1)
# flip 7, 6
tr$tip.label

#~~ Change tip labels

genus <- c("M.", "Putative.", "M.",
           "M.", "M.", "M.",
           "M.", "M.", "M.",
           "M.", "R.")

species <- c("birostris", "species", "alfredi",
             "mobular", "eregoo", "kuhlii",
             "thurstoni", "hypostoma", "munkiana",
             "tarapacana", "bonasus")


d <- data.frame(label = tr$tip.label, genus = genus,
                species = species)
d


p <- ggtree(tr, color = "black", size = 0.4) %<+% d + xlim(NA, 0.002) +
  geom_tiplab(aes(label=paste0('italic(', genus, ')~italic(', species, ')')), parse=T, color = "black")
p


d <- p$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 0,]

p <- p + geom_text(data = d, aes(x = branch, label=label), size = 3, vjust=-.5, color='black')
#p <- p + geom_nodepoint(data = d, size = 3, shape = 21, fill = "white")
p

ggsave("figs/Figure_S6B.tiff", p, height = 6, width = 10)
