#~~ Script to plot PCAs mobulid clades for Dataset B

library(ggplot2)
library(adegenet)
library(dplyr)
library(patchwork)
library(ggtext)

# Function to read plink data and get sp names and geographic info

prep_genlight <- function(plink_raw, meta) {
  
  # read plink raw file
  x <- read.PLINK(plink_raw, n.cores = 1)
  
  # remove unsampled loci
  toRemove <- is.na(glMean(x, alleleAsUnit = F)) # TRUE where NA
  data <- x[, !toRemove]
  
  # assign species info
  
  ids <- read.table(meta)
  sp <- as.factor(ids$V1)
  
  data@pop <- sp
  
  return(data)
  
}

# Run function on each clade

mantas_data <- prep_genlight(plink_raw = "data/p90_mantas.raw",
                             meta = "data/IDs_p10_mantas.txt")

mobjap_data <- prep_genlight(plink_raw = "data/p90_mobjap.raw",
                             meta = "data/IDs_p10_mobjap.txt")

thurerekuh_data <- prep_genlight(plink_raw = "data/p90_thurerekuh.raw",
                                 meta = "data/IDs_p10_thurerekuh.txt")

hypmunk_data <- prep_genlight(plink_raw = "data/p90_hypmunk.raw",
                              meta = "data/IDs_p10_hypmunk.txt")

# Get geographic data

ocean_mantas <- read.table("data/IDs_p10_mantas.txt")
ocean_mantas <- ocean_mantas$V2

ocean_mobjap <- read.table("data/IDs_p10_mobjap_geog.txt")
ocean_mobjap <- as.factor(ocean_mobjap$V1)

ocean_thurerekuh <- read.table("data/IDs_p10_thurerekuh_geog.txt")
ocean_thurerekuh <- as.factor(ocean_thurerekuh$V1)

ocean_hypmunk <- read.table("data/IDs_p10_hypmunk_geog.txt")
ocean_hypmunk <- as.factor(ocean_hypmunk$V1)


#~~~~~~~~~~~~~~~~~~~~#
#      Run PCAs      #
#~~~~~~~~~~~~~~~~~~~~#

#~~ Mantas
pca_mantas <- glPca(mantas_data)
3

#~~ mobular/japanica
pca_mobjap <- glPca(mobjap_data)
3

#~~ thurstoni/kuhlii/eregoodootenkee
pca_thurerekuh <- glPca(thurerekuh_data)
3

#~~ hypostoma/munkiana
pca_hypmunk <- glPca(hypmunk_data)
3

#~~~~~~~~~~~~~~~~~~~~~#
#      Plot PCAs      #
#~~~~~~~~~~~~~~~~~~~~~#

#~~ Mantas

#~~ Assign PCA axes

mantas_pc1 <- pca_mantas$scores[,1]
mantas_pc2 <- pca_mantas$scores[,2]
mantas_pc3 <- pca_mantas$scores[,3]
ind_names_mantas <- mantas_data@pop

ggplot_mantas_1_2 <- as.data.frame(cbind(x = mantas_pc1, y = mantas_pc2)) %>%
  mutate(ind_names = ind_names_mantas,
         ocean = ocean_mantas,
         plot = "a",
         pc = 12)

ggplot_mantas_1_3 <- as.data.frame(cbind(x = mantas_pc1, y = mantas_pc3)) %>%
  mutate(ind_names = ind_names_mantas,
         ocean = ocean_mantas,
         plot = "b",
         pc = 13)

ggplot_mantas_all <- as.data.frame(cbind(PC1 = mantas_pc1, 
                                         PC2 = mantas_pc2,
                                         PC3 = mantas_pc3)) %>%
  mutate(ind_names = ind_names_mantas) %>%
  mutate(ocean = ocean_mantas)

# eig

mantas_eig <- data.frame(pca_mantas$eig)
mantas_eig$percentage = (mantas_eig[, 1]/sum(mantas_eig$pca_mantas.eig))*100
sum(mantas_eig$percentage)
sum(mantas_eig$percentage[1:2])

mantas_eig$percentage <- round(mantas_eig$percentage, digits = 1)
mantas_eig$percentage[1]
mantas_eig$percentage[2]
mantas_eig$percentage[3]

#~~ mobular/japanica

## Assign PCA axes
mobjap_pc1 <- pca_mobjap$scores[,1]
mobjap_pc2 <- pca_mobjap$scores[,2]
mobjap_pc3 <- pca_mobjap$scores[,3]
ind_names_mobjap <- mobjap_data@pop

ggplot_mobjap_1_2 <- as.data.frame(cbind(x = mobjap_pc1, y = mobjap_pc2)) %>%
  mutate(ind_names = ind_names_mobjap,
         ocean = ocean_mobjap,
         plot = "c",
         pc = 12)

ggplot_mobjap_1_3 <- as.data.frame(cbind(x = mobjap_pc1, y = mobjap_pc3)) %>%
  mutate(ind_names = ind_names_mobjap,
         ocean = ocean_mobjap,
         plot = "d",
         pc = 13)

ggplot_mobjap_all <- as.data.frame(cbind(PC1 = mobjap_pc1, 
                                         PC2 = mobjap_pc2,
                                         PC3 = mobjap_pc3)) %>%
  mutate(ind_names = ind_names_mobjap) %>%
  mutate(ocean = ocean_mobjap) 

# eig

mobjap_eig <- data.frame(pca_mobjap$eig)
mobjap_eig$percentage = (mobjap_eig[, 1]/sum(mobjap_eig$pca_mobjap.eig))*100
sum(mobjap_eig$percentage)
sum(mobjap_eig$percentage[1:2])

mobjap_eig$percentage <- round(mobjap_eig$percentage, digits = 1)
mobjap_eig$percentage[1]
mobjap_eig$percentage[2]
mobjap_eig$percentage[3]

## thurstoni/kuhlii/eregoodootenkee

## Assign PCA axes
thurerekuh_pc1 <- pca_thurerekuh$scores[,1]
thurerekuh_pc2 <- pca_thurerekuh$scores[,2]
thurerekuh_pc3 <- pca_thurerekuh$scores[,3]
ind_names_thurerekuh <- thurerekuh_data@pop


ggplot_thurerekuh_1_2 <- as.data.frame(cbind(x = thurerekuh_pc1, y = thurerekuh_pc2)) %>%
  mutate(ind_names = ind_names_thurerekuh,
         ocean = ocean_thurerekuh,
         plot = "e",
         pc = 12)

ggplot_thurerekuh_1_3 <- as.data.frame(cbind(x = thurerekuh_pc1, y = thurerekuh_pc3)) %>%
  mutate(ind_names = ind_names_thurerekuh,
         ocean = ocean_thurerekuh,
         plot = "f",
         pc = 13)

ggplot_thurerekuh_all <- as.data.frame(cbind(PC1 = thurerekuh_pc1, 
                                             PC2 = thurerekuh_pc2,
                                             PC3 = thurerekuh_pc3)) %>%
  mutate(ind_names = ind_names_thurerekuh) %>%
  mutate(ocean = ocean_thurerekuh) 

# eig

thurerekuh_eig <- data.frame(pca_thurerekuh$eig)
thurerekuh_eig$percentage = (thurerekuh_eig[, 1]/sum(thurerekuh_eig$pca_thurerekuh.eig))*100
sum(thurerekuh_eig$percentage)
sum(thurerekuh_eig$percentage[1:2])

thurerekuh_eig$percentage <- round(thurerekuh_eig$percentage, digits = 1)
thurerekuh_eig$percentage[1]
thurerekuh_eig$percentage[2]
thurerekuh_eig$percentage[3]


## hypostoma/munkiana

## Assign PCA axes
hypmunk_pc1 <- pca_hypmunk$scores[,1]
hypmunk_pc2 <- pca_hypmunk$scores[,2]
hypmunk_pc3 <- pca_hypmunk$scores[,3]
ind_names_hypmunk <- hypmunk_data@pop


ggplot_hypmunk_1_2 <- as.data.frame(cbind(x = hypmunk_pc1, y = hypmunk_pc2)) %>%
  mutate(ind_names = ind_names_hypmunk,
         ocean = ocean_hypmunk,
         plot = "g",
         pc = 12)


ggplot_hypmunk_1_3 <- as.data.frame(cbind(x = hypmunk_pc1, y = hypmunk_pc3)) %>%
  mutate(ind_names = ind_names_hypmunk,
         ocean = ocean_hypmunk,
         plot = "h",
         pc = 13)

ggplot_hypmunk_all <- as.data.frame(cbind(PC1 = hypmunk_pc1, 
                                          PC2 = hypmunk_pc2,
                                          PC3 = hypmunk_pc3)) %>%
  mutate(ind_names = ind_names_hypmunk) %>%
  mutate(ocean = ocean_hypmunk) 


# eig

hypmunk_eig <- data.frame(pca_hypmunk$eig)
hypmunk_eig$percentage = (hypmunk_eig[, 1]/sum(hypmunk_eig$pca_hypmunk.eig))*100
sum(hypmunk_eig$percentage)
sum(hypmunk_eig$percentage[1:2])

hypmunk_eig$percentage <- round(hypmunk_eig$percentage, digits = 1)
hypmunk_eig$percentage[1]
hypmunk_eig$percentage[2]
hypmunk_eig$percentage[3]


#~~ Generate global dataframe for plotting

ggplot_mob <- rbind(ggplot_mantas_1_2, ggplot_mantas_1_3,
                    ggplot_mobjap_1_2, ggplot_mobjap_1_3,
                    ggplot_thurerekuh_1_2, ggplot_thurerekuh_1_3,
                    ggplot_hypmunk_1_2, ggplot_hypmunk_1_3) %>%
  mutate(plot = toupper(plot))


ggplot_mob_full <- rbind(ggplot_mantas_all,
                         ggplot_mobjap_all,
                         ggplot_thurerekuh_all,
                         ggplot_hypmunk_all)

#~~ Plot

source("scripts/theme_emily.R")

pal <- c("#0c2c84", # birostris
         "#e41a1c", # alfredi
         "#ff7f00", # new species
         "#984ea3", # japanica
         "#35978f", # mobular
         "#ae017e", # eregoo
         "#1d91c0", # kuhlii
         "#1b7837", # thurstoni 
         "gold2",   # hypostoma
         "#a65628") # munkiana


ggplot_mob$ind_names <- factor(ggplot_mob$ind_names,
                               levels = c("M.birostris", "M.alfredi",
                                          "New.manta", "M.japanica",
                                          "M.mobular", "M.eregoodootenkee",
                                          "M.kuhlii", "M.thurstoni",
                                          "M.hypostoma","M.munkiana"))

ggplot_mob_full$ind_names <- factor(ggplot_mob_full$ind_names,
                                    levels = c("M.birostris", "M.alfredi",
                                               "New.manta", "M.japanica",
                                               "M.mobular", "M.eregoodootenkee",
                                               "M.kuhlii", "M.thurstoni",
                                               "M.hypostoma","M.munkiana"))

sp <- levels(ggplot_mob$ind_names)
loc <- levels(ggplot_mob$ocean)

ggplot_mob <- ggplot_mob %>%
  mutate(plot = as.factor(plot))

#~~ Mantas

a <- ggplot(filter(ggplot_mob_full, ind_names == "M.birostris" |
                     ind_names == "M.alfredi" | ind_names == "New.manta"), 
            aes(PC1, PC2, colour = factor(ind_names), shape = factor(ocean), fill = ind_names)) +
  scale_fill_manual(values = pal,
                    breaks = sp,
                    labels = c("*M. birostris*", 
                               "*M. alfredi*",  
                               "Putative new species", 
                               "*M. mobular cf. japanica*",
                               "*M. mobular*", 
                               "*M. eregoodoo*",
                               "*M. kuhlii*",
                               "*M. thurstoni*", 
                               "*M. hypostoma*", 
                               "*M. munkiana*"),
                    name = "Species") +
  scale_shape_manual(values = c(21, 23, 22, 24, 25),
                     breaks = loc,
                     labels = c("Atlantic", "Indian", "Pacific", "Eastern Indian Ocean", "Western Indian Ocean"),
                     guide = "none") +
  geom_point(size = 3.5, alpha = 0.7, colour = "black") +
  ylab("PC2") +
  xlab("PC1") +
  theme_emily() +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.text = element_markdown(size=12),
        legend.box.margin=margin(-10,-10,-15,-10),
        legend.position = "top", 
        plot.margin = unit(c(0,0,0,2),"cm"),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        strip.text.x = element_text(face = "plain"),
        text = element_text(color='black')) +
  labs(color="Species") +
  labs(shape="Ocean") + 
  guides(fill = guide_legend(override.aes=list(shape=21)))

b <- ggplot(filter(ggplot_mob_full, ind_names == "M.birostris" |
                     ind_names == "M.alfredi" | ind_names == "New.manta"), 
            aes(PC1, PC3, colour = factor(ind_names), shape = factor(ocean), fill = ind_names)) +
  scale_fill_manual(values = pal,
                    breaks = sp,
                    labels = c("*M. birostris*", 
                               "*M. alfredi*",  
                               "Putative new species", 
                               "*M. mobular cf. japanica*",
                               "*M. mobular*", 
                               "*M. eregoodoo*",
                               "*M. kuhlii*",
                               "*M. thurstoni*", 
                               "*M. hypostoma*", 
                               "*M. munkiana*"),
                    name = "Species") +
  scale_shape_manual(values = c(21, 23, 22, 24, 25),
                     breaks = loc,
                     labels = c("Atlantic", "Indian", "Pacific", "Eastern Indian Ocean", "Western Indian Ocean"),
                     guide = "none") +
  geom_point(size = 3.5, alpha = 0.7, colour = "black") +
  ylab("PC3") +
  xlab("PC1") +
  theme_emily() +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "none", 
        text = element_text(color='black')) +
  labs(color="Species") +
  labs(shape="Ocean") + 
  guides(fill = guide_legend(override.aes=list(shape=21)))


a + b 


#~~ Mob / japanica

c <- ggplot(filter(ggplot_mob_full, ind_names == "M.japanica" |
                     ind_names == "M.mobular"), 
            aes(PC1, PC2, colour = factor(ind_names), shape = factor(ocean), fill = ind_names)) +
  scale_fill_manual(values = pal,
                    breaks = sp,
                    labels = c("*M. birostris*", 
                               "*M. alfredi*",  
                               "Putative new species", 
                               "*M. mobular cf. japanica*",
                               "*M. mobular*", 
                               "*M. eregoodoo*",
                               "*M. kuhlii*",
                               "*M. thurstoni*", 
                               "*M. hypostoma*", 
                               "*M. munkiana*"),
                    name = "Species") +
  scale_shape_manual(values = c(21, 23, 22, 24, 25),
                     breaks = loc,
                     labels = c("Atlantic", "Indian", "Pacific", "Eastern Indian Ocean", "Western Indian Ocean"),
                     guide = "none") +
  geom_point(size = 3.5, alpha = 0.7, colour = "black") +
  ylab("PC2") +
  xlab("PC1") +
  theme_emily() +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.placement =  "inside",
        legend.text = element_markdown(size=12),
        legend.box.margin=margin(-10,-10,-15,-10),
        legend.position = "top", 
        plot.margin = unit(c(0,0,0,2),"cm"),
        legend.title = element_blank(),
        strip.text.x = element_text(face = "plain"),
        text = element_text(color='black')) +
  labs(color="Species") +
  labs(shape="Ocean") + 
  guides(fill = guide_legend(override.aes=list(shape=21)))

d <- ggplot(filter(ggplot_mob_full, ind_names == "M.japanica" |
                     ind_names == "M.mobular"), 
            aes(PC1, PC3, colour = factor(ind_names), shape = factor(ocean), fill = ind_names)) +
  scale_fill_manual(values = pal,
                    breaks = sp,
                    labels = c("*M. birostris*", 
                               "*M. alfredi*",  
                               "Putative new species", 
                               "*M. mobular cf. japanica*",
                               "*M. mobular*", 
                               "*M. eregoodoo*",
                               "*M. kuhlii*",
                               "*M. thurstoni*", 
                               "*M. hypostoma*", 
                               "*M. munkiana*"),
                    name = "Species") +
  scale_shape_manual(values = c(21, 23, 22, 24, 25),
                     breaks = loc,
                     labels = c("Atlantic", "Indian", "Pacific", "Eastern Indian Ocean", "Western Indian Ocean"),
                     guide = "none") +
  geom_point(size = 3.5, alpha = 0.7, colour = "black") +
  ylab("PC3") +
  xlab("PC1") +
  theme_emily() +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "none", 
        text = element_text(color='black')) +
  labs(color="Species") +
  labs(shape="Ocean") + 
  guides(fill = guide_legend(override.aes=list(shape=21)))

c + d

# eregoodoo, kuhlii, thurstoni

e <- ggplot(filter(ggplot_mob_full, ind_names == "M.eregoodootenkee" |
                     ind_names == "M.kuhlii" | ind_names == "M.thurstoni"), 
            aes(PC1, PC2, colour = factor(ind_names), shape = factor(ocean), fill = ind_names)) +
  scale_fill_manual(values = pal,
                    breaks = sp,
                    labels = c("*M. birostris*", 
                               "*M. alfredi*",  
                               "Putative new species", 
                               "*M. mobular cf. japanica*",
                               "*M. mobular*", 
                               "*M. eregoodoo*",
                               "*M. kuhlii*",
                               "*M. thurstoni*", 
                               "*M. hypostoma*", 
                               "*M. munkiana*"),
                    name = "Species") +
  scale_shape_manual(values = c(21, 23, 22, 24, 25),
                     breaks = loc,
                     labels = c("Atlantic", "Indian", "Pacific", "Eastern Indian Ocean", "Western Indian Ocean"),
                     guide = "none") +
  geom_point(size = 3.5, alpha = 0.7, colour = "black") +
  ylab("PC2") +
  xlab("PC1") +
  theme_emily() +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.placement =  "inside",
        legend.text = element_markdown(size=12),
        legend.box.margin=margin(-10,-10,-15,-10),
        legend.position = "top", 
        plot.margin = unit(c(0,0,0,2),"cm"),
        legend.title = element_blank(),
        strip.text.x = element_text(face = "plain"),
        text = element_text(color='black')) +
  labs(color="Species") +
  labs(shape="Ocean") + 
  guides(fill = guide_legend(override.aes=list(shape=21)))

f <- ggplot(filter(ggplot_mob_full, ind_names == "M.eregoodootenkee" |
                     ind_names == "M.kuhlii" | ind_names == "M.thurstoni"), 
            aes(PC1, PC3, colour = factor(ind_names), shape = factor(ocean), fill = ind_names)) +
  scale_fill_manual(values = pal,
                    breaks = sp,
                    labels = c("*M. birostris*", 
                               "*M. alfredi*",  
                               "Putative new species", 
                               "*M. mobular cf. japanica*",
                               "*M. mobular*", 
                               "*M. eregoodoo*",
                               "*M. kuhlii*",
                               "*M. thurstoni*", 
                               "*M. hypostoma*", 
                               "*M. munkiana*"),
                    name = "Species") +
  scale_shape_manual(values = c(21, 23, 22, 24, 25),
                     breaks = loc,
                     labels = c("Atlantic", "Indian", "Pacific", "Eastern Indian Ocean", "Western Indian Ocean"),
                     guide = "none") +
  geom_point(size = 3.5, alpha = 0.7, colour = "black") +
  ylab("PC3") +
  xlab("PC1") +
  theme_emily() +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "none", 
        text = element_text(color='black')) +
  labs(color="Species") +
  labs(shape="Ocean") + 
  guides(fill = guide_legend(override.aes=list(shape=21)))


e + f

# hypostoma + munkiana

g <- ggplot(filter(ggplot_mob_full, ind_names == "M.hypostoma" |
                     ind_names == "M.munkiana"), 
            aes(PC1, PC2, colour = factor(ind_names), shape = factor(ocean), fill = ind_names)) +
  scale_fill_manual(values = pal,
                    breaks = sp,
                    labels = c("*M. birostris*", 
                               "*M. alfredi*",  
                               "Putative new species", 
                               "*M. mobular cf. japanica*",
                               "*M. mobular*", 
                               "*M. eregoodoo*",
                               "*M. kuhlii*",
                               "*M. thurstoni*", 
                               "*M. hypostoma*", 
                               "*M. munkiana*"),
                    name = "Species") +
  scale_shape_manual(values = c(21, 23, 22, 24, 25),
                     breaks = loc,
                     labels = c("Atlantic", "Indian", "Pacific", "Eastern Indian Ocean", "Western Indian Ocean"),
                     guide = "none") +
  geom_point(size = 3.5, alpha = 0.7, colour = "black") +
  ylab("PC2") +
  xlab("PC1") +
  theme_emily() +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.text = element_markdown(size=12),
        legend.box.margin=margin(-10,-10,-15,-10),
        legend.position = "top", 
        plot.margin = unit(c(0,1,0,2),"cm"),
        legend.title = element_blank(),
        text = element_text(color='black')) +
  labs(color="Species") +
  labs(shape="Ocean") + 
  guides(fill = guide_legend(override.aes=list(shape=21)))


h <- ggplot(filter(ggplot_mob_full, ind_names == "M.hypostoma" |
                     ind_names == "M.munkiana"), 
            aes(PC1, PC3, colour = factor(ind_names), shape = factor(ocean), fill = ind_names)) +
  scale_fill_manual(values = pal,
                    breaks = sp,
                    labels = c("*M. birostris*", 
                               "*M. alfredi*",  
                               "Putative new species", 
                               "*M. mobular cf. japanica*",
                               "*M. mobular*", 
                               "*M. eregoodoo*",
                               "*M. kuhlii*",
                               "*M. thurstoni*", 
                               "*M. hypostoma*", 
                               "*M. munkiana*"),
                    name = "Species") +
  scale_shape_manual(values = c(21, 23, 22, 24, 25),
                     breaks = loc,
                     labels = c("Atlantic", "Indian", "Pacific", "Eastern Indian Ocean", "Western Indian Ocean"),
                     guide = "none") +
  geom_point(size = 3.5, alpha = 0.7, colour = "black") +
  ylab("PC3") +
  xlab("PC1") +
  theme_emily() +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "none", 
        text = element_text(color='black')) +
  labs(color="Species") +
  labs(shape="Ocean") + 
  guides(fill = guide_legend(override.aes=list(shape=21)))

g + h



png(file="figs/Figure_S2.png", units = "in", res = 300, height=11, width=7)
a + b + c + d + e + f + g + h + plot_layout(ncol = 2)
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      All in one Plot       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

pal <- c("#0c2c84", # birostris
         "#e41a1c", # alfredi
         "#ff7f00", # new species
         "#984ea3", # japanica
         "#35978f", # mobular
         "#ae017e", # eregoo
         "#1d91c0", # kuhlii
         "#1b7837", # thurstoni 
         "gold2",   # hypostoma
         "#a65628") # munkiana

x <- ggplot(filter(ggplot_mob, pc == 12), 
            aes(x, y, shape = factor(ocean), fill = ind_names)) +
  scale_fill_manual(values = pal,
                    breaks = sp) +
  scale_shape_manual(values = c(21, 23, 22, 24, 25),
                     breaks = loc,
                     labels = c("Atlantic", "Indian", "Pacific", 
                                "Eastern Indian Ocean", "Western Indian Ocean")) +
  geom_point(size = 3.5, alpha = 0.7, colour = "black") +
  ylab("PC2") +
  xlab("PC1") +
  theme_emily() +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        legend.position = "none",
        strip.text = element_text(size=13, hjust = 0),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(face = "plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  labs(color="Species", shape = "Ocean", tag = "A")


y <- ggplot(filter(ggplot_mob, pc == 13), 
            aes(x, y, colour = factor(ind_names), shape = factor(ocean), fill = ind_names)) +
  scale_fill_manual(values = pal,
                    breaks = sp,
                    labels = c("*M. birostris*", 
                               "*M. alfredi*",  
                               "Putative new species", 
                               "*M. mobular cf. japanica*",
                               "*M. mobular*", 
                               "*M. kuhlii cf. eregoodootenkee*",
                               "*M. kuhlii*",
                               "*M. thurstoni*", 
                               "*M. hypostoma*", 
                               "*M. munkiana*"),
                    name = "Species") +
  scale_shape_manual(values = c(21, 23, 22, 24, 25),
                     breaks = loc,
                     labels = c("Atlantic", "Indian", "Pacific", "Eastern Indian Ocean", "Western Indian Ocean")) +
  geom_point(size = 3.5, alpha = 0.7, colour = "black") +
  ylab("PC3") +
  xlab("PC1") +
  theme_emily() +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.text = element_text(size=13, hjust = 0),
        legend.text = element_markdown(size=13),
        legend.title = element_text(size=13),
        strip.text.x = element_text(face = "plain"),
        text = element_text(color='black')) +
  labs(color="Species", shape = "Ocean", tag = "B") +
  guides(fill = guide_legend(override.aes=list(shape=21)))

x + y


png(file="figs/Full_PCA_Dataset_B.png", units = "in", res = 300, height=5, width=13)
x + y
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Supplementary Plot       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

mantas_eig <- pca_mantas$eig
mobjap_eig <- pca_mobjap$eig
thurerekuh_eig <- pca_thurerekuh$eig
hypmunk_eig <- pca_hypmunk$eig

ggplot_mantas_eig <- as.data.frame(mantas_eig) %>%
  mutate(eig = mantas_eig, 
         plot = "a")

ggplot_mobjap_eig <- as.data.frame(mobjap_eig) %>%
  mutate(eig = mobjap_eig, plot = "b")

ggplot_thurerekuh_eig <- as.data.frame(thurerekuh_eig) %>%
  mutate(eig = thurerekuh_eig, plot = "c")

ggplot_hypmunk_eig <- as.data.frame(hypmunk_eig) %>%
  mutate(eig = hypmunk_eig, plot = "d")

all_eig <- read.csv("data/all_eig.csv", header=TRUE) %>%
  mutate(plot = toupper(plot)) %>%
  mutate(plot = case_when(plot == "B" ~ "A",
                          plot == "C" ~ "B",
                          plot == "D" ~ "B",
                          plot == "E" ~ "C",
                          plot == "F" ~ "C",
                          plot == "G" ~ "D",
                          plot == "H" ~ "D",
                          TRUE ~ plot))


#~~ Plot

eig_plot_a <- ggplot(filter(all_eig, pc == 12), aes(x=X, y=eig, fill = factor(color))) + 
  scale_fill_manual(values = c("gray70", "black", "gray40")) +
  geom_col() +
  theme_emily() +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.text = element_text(size=13, hjust = 0),
        legend.text = element_markdown(size=13),
        legend.title = element_text(size=13),
        strip.text.x = element_text(face = "plain"),
        legend.position = "none",
        text = element_text(color='black')) +
  xlab("PC") + 
  ylab("Eigenvector") +
  facet_wrap(~ plot, scales = "free", ncol = 1)


eig_plot_b <- ggplot(filter(all_eig, pc == 13), 
                     aes(x=X, y=eig, fill = factor(color))) + 
  scale_fill_manual(values = c("gray70", "black", "gray40")) +
  geom_col() +
  theme_emily() +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.text = element_blank(),
        strip.text.x = element_text(face = "plain"),
        legend.position = "none",
        text = element_text(color='black')) +
  xlab("PC") +
  ylab("Eigenvector") +
  facet_wrap(~ plot, scales = "free", ncol = 1)

eig_plot_a + eig_plot_b

png(file="figs/Figure_S4.png", units = "in", res = 300, height=10, width=9)
eig_plot_a + eig_plot_b
dev.off()

