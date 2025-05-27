library(ape) # Tree manipulation
library(adephylo) # distRoot function
library(tidyverse) # Pipe and dataframe operations
library(boot) # Bootstrap
library(ggplot2) # Plots
library(ggridges) # Distribution plots
library(dunn.test) # Dunn's Test
library(patchwork) # Plot grids

setwd("ERC")

#### Load phylogenetic trees ####

mt.tree <- read.tree(file = "best_tree_mt.nwk") # Mitochondrial genes
nmt.tree <- read.tree(file = "best_tree_N-mt.nwk") # N-mt genes
gly.tree <- read.tree(file= "best_tree_glyco.nwk") # Glycolysis genes
rand.tree <- read.tree(file = "best_tree_rand.nwk") # Random orthologs
no2.tree <- read.tree(file = "best_tree_N-mt_noCII.nwk") # N-mt genes, excluding CII
con.tree <- read.tree(file = "best_tree_contact.nwk") # N-mt genes that contact mt
noncon.tree <- read.tree(file = "best_tree_noncontact.nwk") # N-mt genes that do not contact mt
uce.tree <- read.tree(file = "best_tree_uce_astral_rescaled.nwk")

#### Prune the outgroup (mellifera) from trees ####

mt.tree <- drop.tip(mt.tree, "mellifera")
nmt.tree <- drop.tip(nmt.tree, "mellifera")
gly.tree <- drop.tip(gly.tree, "mellifera")
rand.tree <- drop.tip(rand.tree, "mellifera")
no2.tree <- drop.tip(no2.tree, "mellifera")
con.tree <- drop.tip(con.tree, "mellifera")
noncon.tree <- drop.tip(noncon.tree, "mellifera")
uce.tree <- drop.tip(uce.tree, "mellifera")

#### Extract tip labels ####

mt.spp <- mt.tree$tip.label
nmt.spp <- nmt.tree$tip.label
gly.spp <- gly.tree$tip.label
rand.spp <- rand.tree$tip.label
no2.spp <- no2.tree$tip.label
con.spp <- con.tree$tip.label
noncon.spp <- noncon.tree$tip.label
uce.spp <- uce.tree$tip.label

#### Read species table ####

spp.table <- read.csv("species.csv")

#### Extract branch lengths ####

mt.dist <- distRoot(mt.tree)
nmt.dist <- distRoot(nmt.tree)
gly.dist <- distRoot(gly.tree)
rand.dist <- distRoot(rand.tree)
no2.dist <- distRoot(no2.tree)
con.dist <- distRoot(con.tree)
noncon.dist <- distRoot(noncon.tree)

#### Convert branch lengths to a data frame ####

mt.dist <- data.frame(as.list(mt.dist))
nmt.dist <- data.frame(as.list(nmt.dist))
gly.dist <- data.frame(as.list(gly.dist))
rand.dist <- data.frame(as.list(rand.dist))
no2.dist <- data.frame(as.list(no2.dist))
con.dist <- data.frame(as.list(con.dist))
noncon.dist <- data.frame(as.list(noncon.dist))

mt.dist %>% 
  mutate(species = "mt.dist") %>% 
  gather("species", "mt.dist") -> mt.dist
nmt.dist %>% 
  mutate(species = "nmt.dist") %>% 
  gather("species", "nmt.dist") -> nmt.dist
gly.dist %>% 
  mutate(species = "gly.dist") %>% 
  gather("species", "gly.dist") -> gly.dist
rand.dist %>% 
  mutate(species = "rand.dist") %>% 
  gather("species", "rand.dist") -> rand.dist
no2.dist %>% 
  mutate(species = "no2.dist") %>% 
  gather("species", "no2.dist") -> no2.dist
con.dist %>% 
  mutate(species = "con.dist") %>% 
  gather("species", "con.dist") -> con.dist
noncon.dist %>% 
  mutate(species = "noncon.dist") %>% 
  gather("species", "noncon.dist") -> noncon.dist

spp.table %>% 
  full_join(mt.dist)->dist.table
dist.table %>% 
  full_join(nmt.dist) ->dist.table
dist.table %>% 
  full_join(gly.dist)->dist.table
dist.table %>% 
  full_join(rand.dist)->dist.table
dist.table %>% 
  full_join(no2.dist)->dist.table
dist.table %>%
  full_join(con.dist)->dist.table
dist.table %>%
  full_join(noncon.dist)->dist.table

## Add a column with the mt / nmt rates
## And also normalize rates against the random dataset
dist.table %>%
  mutate(mt.norm = (mt.dist/rand.dist))->dist.table
dist.table %>%
  mutate(nmt.norm = (nmt.dist/rand.dist))->dist.table
dist.table %>%
  mutate(gly.norm = (gly.dist/rand.dist))->dist.table
dist.table %>%
  mutate(no2.norm = (no2.dist/rand.dist))->dist.table
dist.table %>%
  mutate(con.norm = (con.dist/rand.dist))->dist.table
dist.table %>%
  mutate(noncon.norm = (noncon.dist/rand.dist))->dist.table
dist.table %>%
  mutate(mt.nmt.rate = (mt.dist/nmt.dist))->dist.table
dist.table %>%
  mutate(mt.rand.rate = (mt.dist/rand.dist))->dist.table
dist.table %>%
  mutate(nmt.rand.rate = (nmt.dist/rand.dist))->dist.table

#### Summarize branch length results ####

dist.table %>% 
  summarise(mt.avg = mean(mt.dist),
            mt.min = min(mt.dist),
            mt.max = max(mt.dist),
            nmt.avg = mean(nmt.dist),
            nmt.min = min(nmt.dist),
            nmt.max = max(nmt.dist),
            low.mt.nmt.rate = min((mt.dist/nmt.dist)),
            high.mt.nmt.rate = max((mt.dist/nmt.dist)),
            low.mt.rand.rate = min((mt.dist/rand.dist)),
            high.mt.rand.rate = max((mt.dist/rand.dist)),
            fastest.slowest = high.mt.nmt.rate/low.mt.nmt.rate,
            )

#### Convert into phylogenetic independent contrasts ####

## Select the tree which topology will be used

pic.tree <- uce.tree

## Calculate phylogenetic independent contrasts (PICs)

pic <- cbind(mt.pic = pic(dist.table[, "mt.norm"], pic.tree),
             nmt.pic = pic(dist.table[, "nmt.norm"], pic.tree),
             gly.pic = pic(dist.table[, "gly.norm"], pic.tree),
             rand.pic = pic(dist.table[, "rand.dist"], pic.tree),
             no2.pic = pic(dist.table[, "no2.norm"], pic.tree),
             con.pic = pic(dist.table[, "con.norm"], pic.tree),
             noncon.pic = pic(dist.table[, "noncon.norm"], pic.tree)
             )

#### Correlation tests ####

## First some basic correlation tests

## mt x N-mt

cor.test(~ mt.pic + nmt.pic, method = "pearson", data = pic)

## mt x glycolysis

cor.test(~ mt.pic + gly.pic, method = "pearson", data = pic)

## mt x random

cor.test(~ mt.pic + rand.pic, method = "pearson", data = pic)

## mt x contact

cor.test(~ mt.pic + con.pic, method = "pearson", data = pic)

## mt x non-contact

cor.test(~ mt.pic + noncon.pic, method = "pearson", data = pic)

## Now bootstrapping the correlations

## mt x N-mt

set.seed(1)

# Define the correlation function
corr.fun <- function(data, i) {
  # Subset the data using the indices
  subset.data <- data[i, ]
    # Calculate the correlation coefficient
  cor.value <- cor.test(~ mt.pic + nmt.pic, method = "pearson", data = subset.data)$estimate
  return(cor.value)
}

mt.nmt.boot <- boot(data = pic, 
                    statistic = corr.fun,
                    R = 10000
)

mt.nmt.ci <- boot.ci(mt.nmt.boot, type = "bca")

sum.mt.nmt <- do.call(rbind, Map(data.frame,
                                 rho = mt.nmt.ci$t0,
                                 lower = mt.nmt.ci$bca[4],
                                 upper = mt.nmt.ci$bca[5],
                                 comparison ="mt-Nmt"))

## mt x glycolysis

set.seed(1)

# Define the correlation function
corr.fun <- function(data, i) {
  # Subset the data using the indices
  subset.data <- data[i, ]
  # Calculate the correlation coefficient
  cor.value <- cor.test(~ mt.pic + gly.pic, method = "pearson", data = subset.data)$estimate
  return(cor.value)
}

mt.gly.boot <- boot(data = pic, 
                    statistic = corr.fun,
                    R = 10000
)

mt.gly.ci <- boot.ci(mt.gly.boot, type = "bca")

sum.mt.gly <- do.call(rbind, Map(data.frame,
                                 rho = mt.gly.ci$t0,
                                 lower = mt.gly.ci$bca[4],
                                 upper = mt.gly.ci$bca[5],
                                 comparison ="mt-gly"))

## mt x random

set.seed(1)

# Define the correlation function
corr.fun <- function(data, i) {
  # Subset the data using the indices
  subset.data <- data[i, ]
  # Calculate the correlation coefficient
  cor.value <- cor.test(~ mt.pic + rand.pic, method = "pearson", data = subset.data)$estimate
  return(cor.value)
}

mt.rand.boot <- boot(data = pic, 
                    statistic = corr.fun,
                    R = 10000
)

mt.rand.ci <- boot.ci(mt.rand.boot, type = "bca")

sum.mt.rand <- do.call(rbind, Map(data.frame,
                                 rho = mt.rand.ci$t0,
                                 lower = mt.rand.ci$bca[4],
                                 upper = mt.rand.ci$bca[5],
                                 comparison ="mt-random"))

## mt x contact

set.seed(1)

# Define the correlation function
corr.fun <- function(data, i) {
  # Subset the data using the indices
  subset.data <- data[i, ]
  # Calculate the correlation coefficient
  cor.value <- cor.test(~ mt.pic + con.pic, method = "pearson", data = subset.data)$estimate
  return(cor.value)
}

mt.con.boot <- boot(data = pic, 
                     statistic = corr.fun,
                     R = 10000
)

mt.con.ci <- boot.ci(mt.con.boot, type = "bca")

sum.mt.con <- do.call(rbind, Map(data.frame,
                                  rho = mt.con.ci$t0,
                                  lower = mt.con.ci$bca[4],
                                  upper = mt.con.ci$bca[5],
                                  comparison ="mt-contact"))

## mt x non-contact

set.seed(1)

# Define the correlation function
corr.fun <- function(data, i) {
  # Subset the data using the indices
  subset.data <- data[i, ]
  # Calculate the correlation coefficient
  cor.value <- cor.test(~ mt.pic + noncon.pic, method = "pearson", data = subset.data)$estimate
  return(cor.value)
}

mt.noncon.boot <- boot(data = pic, 
                    statistic = corr.fun,
                    R = 10000
)

mt.noncon.ci <- boot.ci(mt.noncon.boot, type = "bca")

sum.mt.noncon <- do.call(rbind, Map(data.frame,
                                   rho = mt.noncon.ci$t0,
                                   lower = mt.noncon.ci$bca[4],
                                   upper = mt.noncon.ci$bca[5],
                                   comparison ="mt-noncontact"))

## Long-faced

lf.table <- dist.table %>% filter(subgenus=="Thoracobombus" | 
                                      subgenus=="Psithyrus" | 
                                      subgenus=="Subterraneobombus" | 
                                      subgenus=="Orientalibombus" |
                                      subgenus=="Megabombus")

lf.tree <- keep.tip(pic.tree, lf.table[,"species"])

lf.pic <- cbind(mt.pic = pic(lf.table[, "mt.dist"], lf.tree),
                  nmt.pic = pic(lf.table[, "nmt.dist"], lf.tree))

set.seed(1)

# Define the correlation function
corr.fun <- function(data, i) {
  # Subset the data using the indices
  subset.data <- data[i, ]
  # Calculate the correlation coefficient
  cor.value <- cor.test(~ mt.pic + nmt.pic, method = "pearson", data = subset.data)$estimate
  return(cor.value)
}

mt.nmt.lf.boot <- boot(data = lf.pic, 
                       statistic = corr.fun,
                       R = 10000
)

mt.nmt.lf.ci <- boot.ci(mt.nmt.lf.boot, type = "bca")

sum.mt.nmt.lf <- do.call(rbind, Map(data.frame,
                                    rho = mt.nmt.lf.ci$t0,
                                    lower = mt.nmt.lf.ci$bca[4],
                                    upper = mt.nmt.lf.ci$bca[5],
                                    comparison ="long-faced"))

## Short-faced

sf.table <- dist.table %>% filter(subgenus=="Bombus" | 
                                    subgenus=="Alpinobombus" | 
                                    subgenus=="Pyrobombus" | 
                                    subgenus=="Subterraneobombus" |
                                    subgenus=="Cullumanobombus"|
                                    subgenus=="Alpigenobombus"|
                                    subgenus=="Melanobombus")

sf.tree <- keep.tip(pic.tree, sf.table[,"species"])

sf.pic <- cbind(mt.pic = pic(sf.table[, "mt.dist"], sf.tree),
                nmt.pic = pic(sf.table[, "nmt.dist"], sf.tree))

set.seed(1)

# Define the correlation function
corr.fun <- function(data, i) {
  # Subset the data using the indices
  subset.data <- data[i, ]
  # Calculate the correlation coefficient
  cor.value <- cor.test(~ mt.pic + nmt.pic, method = "pearson", data = subset.data)$estimate
  return(cor.value)
}

mt.nmt.sf.boot <- boot(data = sf.pic, 
                       statistic = corr.fun,
                       R = 10000
)

mt.nmt.sf.ci <- boot.ci(mt.nmt.sf.boot, type = "bca")

sum.mt.nmt.sf <- do.call(rbind, Map(data.frame,
                                    rho = mt.nmt.sf.ci$t0,
                                    lower = mt.nmt.sf.ci$bca[4],
                                    upper = mt.nmt.sf.ci$bca[5],
                                    comparison ="short-faced"))


## Compile results in a data frame

sum.mt.nmt %>% 
  full_join(sum.mt.gly) -> summary.stats
summary.stats %>% 
  full_join(sum.mt.rand) -> summary.stats
summary.stats %>% 
  full_join(sum.mt.con) -> summary.stats
summary.stats %>% 
  full_join(sum.mt.noncon) -> summary.stats

#### Extract data from lists and make dataframe

boot.dist.mt.nmt <- do.call(rbind, Map(data.frame, val = mt.nmt.boot$t, comp = "N-mt OXPHOS"))
boot.dist.mt.gly <- do.call(rbind, Map(data.frame, val = mt.gly.boot$t, comp="Glycolysis"))
boot.dist.mt.rand <- do.call(rbind, Map(data.frame, val = mt.rand.boot$t, comp="Random"))
boot.dist.mt.con <- do.call(rbind, Map(data.frame, val = mt.con.boot$t, comp = "Contact"))
boot.dist.mt.noncon <- do.call(rbind, Map(data.frame, val = mt.noncon.boot$t, comp = "Non-contact"))

boot.dist.mt.nmt %>% 
  full_join(boot.dist.mt.gly) -> boot.table
boot.table %>%
  full_join(boot.dist.mt.rand) -> boot.table

boot.dist.mt.con %>% 
  full_join(boot.dist.mt.noncon) -> contact.boot.table

#### Correlation tests - clades/subsets ####

# Subset tree to clade Short-Faced

sf.table <- dist.table %>% filter(subgenus %in% c("Melanobombus","Alpigenobombus",
                                                  "Cullumanobombus", "Sibiricobombus",
                                                  "Alpinobombus", "Bombus", "Pyrobombus"))

sf.tree <- keep.tip(pic.tree, sf.table[,"species"])

# Subset tree to clade Long-Faced

lf.table <- dist.table %>% filter(subgenus %in% c("Subterraneobombus","Megabombus",
                                                  "Orientalibombus", "Psithyrus",
                                                  "Thoracobombus"))

lf.tree <- keep.tip(pic.tree, lf.table[,"species"])

# Subset tree to Psithyrus

ps.table <- dist.table %>% filter(subgenus=="Psithyrus")

ps.tree <- keep.tip(pic.tree, ps.table[,"species"])

# Subset tree to Megabombus

mg.table <- dist.table %>% filter(subgenus=="Megabombus")

mg.tree <- keep.tip(pic.tree, mg.table[,"species"])

# Subset tree to Thoracobombus

th.table <- dist.table %>% filter(subgenus=="Thoracobombus")

th.tree <- keep.tip(pic.tree, th.table[,"species"])

# Subset tree to Pyrobombus

pr.table <- dist.table %>% filter(subgenus=="Pyrobombus")

pr.tree <- keep.tip(pic.tree, pr.table[,"species"])

# Subset tree to Alpinobombus

al.table <- dist.table %>% filter(subgenus=="Alpinobombus")

al.tree <- keep.tip(pic.tree, al.table[,"species"])

# Subset tree to Bombus

bo.table <- dist.table %>% filter(subgenus=="Bombus")

bo.tree <- keep.tip(pic.tree, bo.table[,"species"])

# Calculate phylogenetic independent contrasts (PICs)

male.pic <- cbind(mt.pic = pic(male.table[, "mt.norm"], male.tree),
                  nmt.pic = pic(male.table[, "nmt.norm"], male.tree))

female.pic <- cbind(mt.pic = pic(female.table[, "mt.norm"], female.tree),
                    nmt.pic = pic(female.table[, "nmt.norm"], female.tree))

parasite.pic <- cbind(mt.pic = pic(parasite.table[, "mt.norm"], parasite.tree),
                      nmt.pic = pic(parasite.table[, "nmt.norm"], parasite.tree))

lf.pic <- cbind(mt.pic = pic(lf.table[, "mt.norm"], lf.tree),
                nmt.pic = pic(lf.table[, "nmt.norm"], lf.tree))

sf.pic <- cbind(mt.pic = pic(sf.table[, "mt.norm"], sf.tree),
                nmt.pic = pic(sf.table[, "nmt.norm"], sf.tree))

ps.pic <- cbind(mt.pic = pic(ps.table[, "mt.norm"], ps.tree),
                nmt.pic = pic(ps.table[, "nmt.norm"], ps.tree))

mg.pic <- cbind(mt.pic = pic(mg.table[, "mt.norm"], mg.tree),
                nmt.pic = pic(mg.table[, "nmt.norm"], mg.tree))

th.pic <- cbind(mt.pic = pic(th.table[, "mt.norm"], th.tree),
                nmt.pic = pic(th.table[, "nmt.norm"], th.tree))

pr.pic  <- cbind(mt.pic = pic(pr.table[, "mt.norm"], pr.tree),
                 nmt.pic = pic(pr.table[, "nmt.norm"], pr.tree))

al.pic  <- cbind(mt.pic = pic(al.table[, "mt.norm"], al.tree),
                 nmt.pic = pic(al.table[, "nmt.norm"], al.tree))

bo.pic  <- cbind(mt.pic = pic(bo.table[, "mt.norm"], bo.tree),
                 nmt.pic = pic(bo.table[, "nmt.norm"], bo.tree))

# Conduct some simple correlation tests

## Test for Long-Faced

cor.test(~ mt.pic + nmt.pic, method = "pearson", data = lf.pic)

## Test for Short-Faced

cor.test(~ mt.pic + nmt.pic, method = "pearson", data = sf.pic)

## Test for Psithyrus

cor.test(~ mt.pic + nmt.pic, method = "pearson", data = ps.pic)

## Test for Megabombus

cor.test(~ mt.pic + nmt.pic, method = "pearson", data = mg.pic)

## Test for Thoracobombus

cor.test(~ mt.pic + nmt.pic, method = "pearson", data = th.pic)

## Test for Pyrobombus

cor.test(~ mt.pic + nmt.pic, method = "pearson", data = pr.pic)

## Test for Alpinobombus

cor.test(~ mt.pic + nmt.pic, method = "pearson", data = al.pic)

## Test for Bombus

cor.test(~ mt.pic + nmt.pic, method = "pearson", data = bo.pic)

# Bootstrap correlations

# Define the correlation function
corr.fun <- function(data, i) {
  # Subset the data using the indices
  subset.data <- data[i, ]
  # Calculate the correlation coefficient
  cor.value <- cor.test(~ mt.pic + nmt.pic, method = "pearson", data = subset.data)$estimate
  return(cor.value)
}

## Long-Faced

set.seed(1)

lf.boot <- boot(data = lf.pic, 
                      statistic = corr.fun,
                      R = 10000
)

lf.ci <- boot.ci(lf.boot, type = "bca")

sum.lf <- do.call(rbind, Map(data.frame,
                                   rho = lf.ci$t0,
                                   lower = lf.ci$bca[4],
                                   upper = lf.ci$bca[5],
                                   comparison ="Long-Faced"))

## Short-Faced

set.seed(1)

sf.boot <- boot(data = sf.pic, 
                statistic = corr.fun,
                R = 10000
)

sf.ci <- boot.ci(sf.boot, type = "bca")

sum.sf <- do.call(rbind, Map(data.frame,
                             rho = sf.ci$t0,
                             lower = sf.ci$bca[4],
                             upper = sf.ci$bca[5],
                             comparison ="Short-Faced"))

#### Plot the data ####

# Reorder labels
boot.table$comp <- factor(boot.table$comp, levels = c("N-mt OXPHOS", "Glycolysis", "Random"))

fig4a <- ggplot(boot.table, aes(y = val, x = comp, fill = comp)) +
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) +
  stat_summary(fun=mean, geom="point", shape=19, size=2) +
  ylab(expression(paste("ERC strength  ",(italic(ρ))))) +
  xlab(("ERC between mt and")) +
  scale_x_discrete(labels = c("N-mt", "Glycolysis", "Random")) +
  scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits=c(-0.5, 1.0)) +
  geom_hline(yintercept= 0, linetype="dashed") +
  scale_fill_manual(values=c("#C5E5FB", "#ECECEC", "#ECECEC")) +
  theme_bw() + theme(legend.position = "none")

fig4a

# Reorder labels
contact.boot.table$comp <- factor(contact.boot.table$comp, levels = c("Contact", "Non-contact"))

fig4c <- ggplot(contact.boot.table, aes(y = val, x = comp, fill = comp)) +
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) +
  stat_summary(fun=mean, geom="point", shape=19, size=2) +
  ylab(element_blank()) +
  xlab(("mt N-mt ERC of subunits in")) +
  scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits=c(-0.5, 1.0)) +
  geom_hline(yintercept= 0, linetype="dashed")+
  scale_fill_manual(values=c("#C5E5FB", "#C5E5FB"))+
  theme_bw() + theme(legend.position = "none")

fig4c

#### Complex data #####
#### Load phylogenetic trees ####

nmt.ci.tree <- read.tree(file = "complex/best_tree_N-mt_CI.nwk") # N-mt, CI
nmt.cii.tree <- read.tree(file = "complex/best_tree_N-mt_CII_var2.nwk") # N-mt, CII
nmt.ciii.tree <- read.tree(file = "complex/best_tree_N-mt_CIII.nwk") # N-mt, CIII
nmt.civ.tree <- read.tree(file = "complex/best_tree_N-mt_CIV.nwk") # N-mt, CIV
nmt.cv.tree <- read.tree(file = "complex/best_tree_N-mt_CV.nwk") # N-mt, CV
mt.ci.tree <- read.tree(file = "complex/best_tree_mt_CI.nwk") # mt, CI
mt.ciii.tree <- read.tree(file = "complex/best_tree_mt_CIII.nwk") # mt, CIII
mt.civ.tree <- read.tree(file = "complex/best_tree_mt_CIV.nwk") # mt, CIV
mt.cv.tree <- read.tree(file = "complex/best_tree_mt_CV.nwk") # mt, CV

#### Prune the outgroup (mellifera) from trees ####

nmt.ci.tree <- drop.tip(nmt.ci.tree, "mellifera")
nmt.cii.tree <- drop.tip(nmt.cii.tree, "mellifera")
nmt.ciii.tree <- drop.tip(nmt.ciii.tree, "mellifera")
nmt.civ.tree <- drop.tip(nmt.civ.tree, "mellifera")
nmt.cv.tree <- drop.tip(nmt.cv.tree, "mellifera")
mt.ci.tree <- drop.tip(mt.ci.tree, "mellifera")
mt.ciii.tree <- drop.tip(mt.ciii.tree, "mellifera")
mt.civ.tree <- drop.tip(mt.civ.tree, "mellifera")
mt.cv.tree <- drop.tip(mt.cv.tree, "mellifera")

#### Extract tip labels ####

nmt.ci.spp <- nmt.ci.tree$tip.label
nmt.cii.spp <- nmt.cii.tree$tip.label
nmt.ciii.spp <- nmt.ciii.tree$tip.label
nmt.civ.spp <- nmt.civ.tree$tip.label
nmt.cv.spp <- nmt.cv.tree$tip.label
mt.ci.spp <- mt.ci.tree$tip.label
mt.ciii.spp <- mt.ciii.tree$tip.label
mt.civ.spp <- mt.civ.tree$tip.label
mt.cv.spp <- mt.cv.tree$tip.label

#### Extract branch lengths ####

nmt.ci.dist <- distRoot(nmt.ci.tree)
nmt.cii.dist <- distRoot(nmt.cii.tree)
nmt.ciii.dist <- distRoot(nmt.ciii.tree)
nmt.civ.dist <- distRoot(nmt.civ.tree)
nmt.cv.dist <- distRoot(nmt.cv.tree)
mt.ci.dist <- distRoot(mt.ci.tree)
mt.ciii.dist <- distRoot(mt.ciii.tree)
mt.civ.dist <- distRoot(mt.civ.tree)
mt.cv.dist <- distRoot(mt.cv.tree)

#### Convert branch lengths to a data frame ####

nmt.ci.dist <- data.frame(as.list(nmt.ci.dist))
nmt.cii.dist <- data.frame(as.list(nmt.cii.dist))
nmt.ciii.dist <- data.frame(as.list(nmt.ciii.dist))
nmt.civ.dist <- data.frame(as.list(nmt.civ.dist))
nmt.cv.dist <- data.frame(as.list(nmt.cv.dist))
mt.ci.dist <- data.frame(as.list(mt.ci.dist))
mt.ciii.dist <- data.frame(as.list(mt.ciii.dist))
mt.civ.dist <- data.frame(as.list(mt.civ.dist))
mt.cv.dist <- data.frame(as.list(mt.cv.dist))

nmt.ci.dist %>% 
  mutate(species = "nmt.ci.dist") %>% 
  gather("species", "nmt.ci.dist") -> nmt.ci.dist
nmt.cii.dist %>% 
  mutate(species = "nmt.cii.dist") %>% 
  gather("species", "nmt.cii.dist") -> nmt.cii.dist
nmt.ciii.dist %>% 
  mutate(species = "nmt.ciii.dist") %>% 
  gather("species", "nmt.ciii.dist") -> nmt.ciii.dist
nmt.civ.dist %>% 
  mutate(species = "nmt.civ.dist") %>% 
  gather("species", "nmt.civ.dist") -> nmt.civ.dist
nmt.cv.dist %>% 
  mutate(species = "nmt.cv.dist") %>% 
  gather("species", "nmt.cv.dist") -> nmt.cv.dist
mt.ci.dist %>% 
  mutate(species = "mt.ci.dist") %>% 
  gather("species", "mt.ci.dist") -> mt.ci.dist
mt.ciii.dist %>% 
  mutate(species = "mt.ciii.dist") %>% 
  gather("species", "mt.ciii.dist") -> mt.ciii.dist
mt.civ.dist %>% 
  mutate(species = "mt.civ.dist") %>% 
  gather("species", "mt.civ.dist") -> mt.civ.dist
mt.cv.dist %>% 
  mutate(species = "mt.cv.dist") %>% 
  gather("species", "mt.cv.dist") -> mt.cv.dist

spp.table %>% 
  full_join(nmt.ci.dist) -> complex.table
complex.table %>% 
  full_join(nmt.cii.dist) -> complex.table
complex.table %>% 
  full_join(nmt.ciii.dist) -> complex.table
complex.table %>% 
  full_join(nmt.civ.dist) -> complex.table
complex.table %>% 
  full_join(nmt.cv.dist) -> complex.table
complex.table %>%
  full_join(mt.ci.dist) -> complex.table
complex.table %>%
  full_join(mt.ciii.dist) -> complex.table
complex.table %>% 
  full_join(mt.civ.dist) -> complex.table
complex.table %>% 
  full_join(mt.cv.dist) -> complex.table
complex.table %>% 
  full_join(rand.dist) -> complex.table

# Normalize rates against the random genes
complex.table %>%
  mutate(nmt.ci.dist = (nmt.ci.dist/rand.dist))->complex.table
complex.table %>%
  mutate(nmt.cii.dist = (nmt.cii.dist/rand.dist))->complex.table
complex.table %>%
  mutate(nmt.ciii.dist = (nmt.ciii.dist/rand.dist))->complex.table
complex.table %>%
  mutate(nmt.civ.dist = (nmt.civ.dist/rand.dist))->complex.table
complex.table %>%
  mutate(nmt.cv.dist = (nmt.cv.dist/rand.dist))->complex.table
complex.table %>%
  mutate(mt.ci.dist = (mt.ci.dist/rand.dist))->complex.table
complex.table %>%
  mutate(mt.ciii.dist = (mt.ciii.dist/rand.dist))->complex.table
complex.table %>%
  mutate(mt.civ.dist = (mt.civ.dist/rand.dist))->complex.table
complex.table %>%
  mutate(mt.cv.dist = (mt.cv.dist/rand.dist))->complex.table



#### Convert into phylogenetic independent contrasts ####

pic.tree <- uce.tree

pic.complex <- cbind(nmt.ci.pic = pic(complex.table[, "nmt.ci.dist"], pic.tree),
             nmt.cii.pic = pic(complex.table[, "nmt.cii.dist"], pic.tree),
             nmt.ciii.pic = pic(complex.table[, "nmt.ciii.dist"], pic.tree),
             nmt.civ.pic = pic(complex.table[, "nmt.civ.dist"], pic.tree),
             nmt.cv.pic = pic(complex.table[, "nmt.cv.dist"], pic.tree),
             mt.ci.pic = pic(complex.table[, "mt.ci.dist"], pic.tree),
             mt.ciii.pic = pic(complex.table[, "mt.ciii.dist"], pic.tree),
             mt.civ.pic = pic(complex.table[, "mt.civ.dist"], pic.tree),
             mt.cv.pic = pic(complex.table[, "mt.cv.dist"], pic.tree),
             gly.pic = pic(dist.table[, "gly.norm"], pic.tree),
             rand.pic = pic(dist.table[, "rand.dist"], pic.tree),
             mt.pic = pic(dist.table[, "mt.norm"], pic.tree),
             nmt.pic = pic(dist.table[, "nmt.norm"], pic.tree))

#### Correlation tests ####

## First some basic correlation tests

## mt CI x N-mt CI

cor.test(~ nmt.ci.pic + mt.ci.pic, method = "pearson", data = pic.complex)

## all mt x N-mt CII

cor.test(~ mt.pic + nmt.cii.pic, method = "pearson", data = pic.complex)

## mt CIII x N-mt CIII

cor.test(~ nmt.ciii.pic + mt.ciii.pic, method = "pearson", data = pic.complex)

## mt x contact

cor.test(~ nmt.civ.pic + mt.civ.pic, method = "pearson", data = pic.complex)

## mt x non-contact

cor.test(~ nmt.cv.pic + mt.cv.pic, method = "pearson", data = pic.complex)

## all N-mt x N-mt CII
## CII could act as a control dataset

cor.test(~ mt.pic + nmt.cii.pic, method = "pearson", data = pic.complex)

## Now bootstrapping the correlations

## mt CI x N-mt CI

set.seed(1)

# Define the correlation function
corr.fun <- function(data, i) {
  # Subset the data using the indices
  subset.data <- data[i, ]
  # Calculate the correlation coefficient
  cor.value <- cor.test(~ nmt.ci.pic + mt.ci.pic, method = "pearson", data = subset.data)$estimate
  return(cor.value)
}

ci.boot <- boot(data = pic.complex, 
                    statistic = corr.fun,
                    R = 10000
)

ci.ci <- boot.ci(ci.boot, type = "bca")

sum.ci <- do.call(rbind, Map(data.frame,
                                 rho = ci.ci$t0,
                                 lower = ci.ci$bca[4],
                                 upper = ci.ci$bca[5],
                                 comparison ="CI"))

## all mt x N-mt CII

set.seed(1)

# Define the correlation function
corr.fun <- function(data, i) {
  # Subset the data using the indices
  subset.data <- data[i, ]
  # Calculate the correlation coefficient
  cor.value <- cor.test(~ nmt.cii.pic + mt.pic, method = "pearson", data = subset.data)$estimate
  return(cor.value)
}

cii.boot <- boot(data = pic.complex, 
                statistic = corr.fun,
                R = 10000
)

cii.ci <- boot.ci(cii.boot, type = "bca")

sum.cii <- do.call(rbind, Map(data.frame,
                             rho = cii.ci$t0,
                             lower = cii.ci$bca[4],
                             upper = cii.ci$bca[5],
                             comparison ="CII"))

## mt CIII x N-mt CIII

set.seed(1)

# Define the correlation function
corr.fun <- function(data, i) {
  # Subset the data using the indices
  subset.data <- data[i, ]
  # Calculate the correlation coefficient
  cor.value <- cor.test(~ nmt.ciii.pic + mt.ciii.pic, method = "pearson", data = subset.data)$estimate
  return(cor.value)
}

ciii.boot <- boot(data = pic.complex, 
                statistic = corr.fun,
                R = 10000
)

ciii.ci <- boot.ci(ciii.boot, type = "bca")

sum.ciii <- do.call(rbind, Map(data.frame,
                             rho = ciii.ci$t0,
                             lower = ciii.ci$bca[4],
                             upper = ciii.ci$bca[5],
                             comparison ="CIII"))

## all mt x N-mt CIII

set.seed(1)

# Define the correlation function
corr.fun <- function(data, i) {
  # Subset the data using the indices
  subset.data <- data[i, ]
  # Calculate the correlation coefficient
  cor.value <- cor.test(~ nmt.ciii.pic + mt.pic, method = "pearson", data = subset.data)$estimate
  return(cor.value)
}

ciii.all.boot <- boot(data = pic.complex, 
                  statistic = corr.fun,
                  R = 10000
)

ciii.all.ci <- boot.ci(ciii.all.boot, type = "bca")

sum.ciii.all <- do.call(rbind, Map(data.frame,
                               rho = ciii.all.ci$t0,
                               lower = ciii.all.ci$bca[4],
                               upper = ciii.all.ci$bca[5],
                               comparison ="CIII (vs all mt)"))

## mt CIV x N-mt CIV

set.seed(1)

# Define the correlation function
corr.fun <- function(data, i) {
  # Subset the data using the indices
  subset.data <- data[i, ]
  # Calculate the correlation coefficient
  cor.value <- cor.test(~ nmt.civ.pic + mt.civ.pic, method = "pearson", data = subset.data)$estimate
  return(cor.value)
}

civ.boot <- boot(data = pic.complex, 
                  statistic = corr.fun,
                  R = 10000
)

civ.ci <- boot.ci(civ.boot, type = "bca")

sum.civ <- do.call(rbind, Map(data.frame,
                               rho = civ.ci$t0,
                               lower = civ.ci$bca[4],
                               upper = civ.ci$bca[5],
                               comparison ="CIV"))

## mt CV x N-mt CV

set.seed(1)

# Define the correlation function
corr.fun <- function(data, i) {
  # Subset the data using the indices
  subset.data <- data[i, ]
  # Calculate the correlation coefficient
  cor.value <- cor.test(~ nmt.cv.pic + mt.cv.pic, method = "pearson", data = subset.data)$estimate
  return(cor.value)
}

cv.boot <- boot(data = pic.complex, 
                 statistic = corr.fun,
                 R = 10000
)

cv.ci <- boot.ci(cv.boot, type = "bca")

sum.cv <- do.call(rbind, Map(data.frame,
                              rho = cv.ci$t0,
                              lower = cv.ci$bca[4],
                              upper = cv.ci$bca[5],
                              comparison ="CV"))

## Compile results in a data frame

sum.ci %>% 
  full_join(sum.cii) -> c.summary.stats
c.summary.stats %>% 
  full_join(sum.ciii) -> c.summary.stats
c.summary.stats %>%
  full_join(sum.ciii.all) -> c.summary.stats
c.summary.stats %>% 
  full_join(sum.civ) -> c.summary.stats
c.summary.stats %>% 
  full_join(sum.cv) -> c.summary.stats

#### Extract data from lists and make dataframe

boot.dist.ci <- do.call(rbind, Map(data.frame, val = ci.boot$t, comp = "Complex I"))
boot.dist.cii <- do.call(rbind, Map(data.frame, val = cii.boot$t, comp = "Complex II"))
boot.dist.ciii <- do.call(rbind, Map(data.frame, val = ciii.boot$t, comp = "Complex III"))
boot.dist.civ <- do.call(rbind, Map(data.frame, val = civ.boot$t, comp = "Complex IV"))
boot.dist.cv <- do.call(rbind, Map(data.frame, val = cv.boot$t, comp = "Complex V"))

boot.dist.ci %>% 
  full_join(boot.dist.cii) -> complex.boot.table
complex.boot.table %>%
  full_join(boot.dist.ciii) -> complex.boot.table
complex.boot.table %>%
  full_join(boot.dist.civ) -> complex.boot.table
complex.boot.table %>%
  full_join(boot.dist.cv) -> complex.boot.table

#### Plot the data ####

# Reorder labels
fig4b <- ggplot(complex.boot.table, aes(y = val, x = comp, fill = comp)) +
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) +
  stat_summary(fun=mean, geom="point", shape=19, size=2) +
  ylab(element_blank()) +
  xlab(("mt N-mt ERC within")) +
  scale_y_continuous(breaks = seq(-0.5, 1.0, by = 0.5), limits=c(-0.5, 1.0)) +
  geom_hline(yintercept= 0, linetype="dashed")+
  scale_fill_manual(values=c("#C5E5FB","#ECECEC","#C5E5FB","#C5E5FB","#C5E5FB"))+
  theme_bw() + theme(legend.position = "none")

#### dN/dS plots ####

dnds.dat <- read.csv(file = "dNdS_results.csv")

dnds.dat$group <- factor(dnds.dat$group, c("OXPHOS","Nmt", "Glycolysis", "Random"))
dnds.dat$contact <- factor(dnds.dat$contact, c("Contact", "Non-contact"))

dnds.dat %>% 
  mutate(mt_interact = case_when(group == "Nmt" ~ "yes",
                                 TRUE ~ ("no"))) -> dnds.dat
dnds.dat %>% 
  group_by(group) %>% 
  summarise(avg = round(mean(omega),digits = 3),
            se = round(sd(omega)/sqrt(length(omega)),digits = 3)) -> dnds.plot.dat

# Plot violin plots for omega values
omega.plot <- ggplot(dnds.dat, aes(x=group, y=omega)) + 
  geom_violin(draw_quantiles = c(0.5), aes(fill = mt_interact)) + 
  geom_jitter(shape=21, position=position_jitter(0.2), alpha = 0.6, fill="#084594", size=2) +
  labs(x=NULL, y = expression(paste(italic("d"["N"]), " / ", italic("d"["S"])))) +
  scale_x_discrete(labels = c("mt", "N-mt", "Glycolysis", "Random")) +
  scale_y_continuous(breaks = seq(0, 1.0, by = 0.5), limits=c(0.0, 1.05)) +
  theme_bw() + theme(panel.background = element_blank(), legend.position="none") +
  scale_fill_manual(values = c("#ECECEC", "#C5E5FB")) +
  stat_summary(geom = "text", label = c("a", "b", "ac", "c"), fun.y = max, vjust = -1)

omega.plot

# Run Kruskal-Wallis and Dunn's Test
omega.test <- kruskal.test(omega~group, data=dnds.dat)
omega.test <- dunn.test(dnds.dat$omega, dnds.dat$group, method = "bonferroni", list = T)

# Create a dataset for the complexes
dnds.dat %>% 
  filter(group == "Nmt") %>%
  mutate(complex = case_when(str_detect(gene, "CI_") ~ "CI",
                             str_detect(gene, "CII_") ~ "CII",
                             str_detect(gene, "CIII_") ~ "CIII",
                             str_detect(gene, "CIV_") ~ "CIV",
                             str_detect(gene, "CV_") ~ "CV"),
         mt_interact = case_when(complex == "CII"~ "no",
                                 TRUE ~ ("yes"))) -> dnds.complex.dat

omega.plot.complexes <- ggplot(dnds.complex.dat, aes(x=complex, y=omega)) + 
  geom_violin(draw_quantiles = c(0.5), aes(fill = mt_interact)) + 
  geom_jitter(shape=21, position=position_jitter(0.2), alpha = 0.6, fill="#084594", size=2) +
  ylab(element_blank()) + xlab("OXPHOS complex") +
  scale_y_continuous(breaks = seq(0, 1.0, by = 0.5), limits=c(0.0, 1.05)) +
  theme_bw() + theme(panel.background = element_blank(), legend.position="none") +
  scale_fill_manual(values = c("#ECECEC", "#C5E5FB")) +
  stat_summary(geom = "text", label = "a", fun.y = max, vjust = -1)

omega.plot.complexes

# Run Kruskal-Wallis and Dunn's Test
omega.complex.test <- kruskal.test(omega~complex, data=dnds.complex.dat)
omega.complex.test <- dunn.test(dnds.complex.dat$omega, dnds.complex.dat$complex, method = "bonferroni", list = T)

# Contact and non-contact

omega.plot.contact <- ggplot(na.omit(dnds.dat), aes(x=contact, y=omega)) + 
  geom_violin(draw_quantiles = c(0.5), aes(fill = mt_interact), width = 0.5) + 
  geom_jitter(shape=21, position=position_jitter(0.2), alpha = 0.6, fill="#084594", size=2) +
  ylab(element_blank()) + xlab(element_blank()) +
  scale_y_continuous(breaks = seq(0, 1.0, by = 0.5), limits=c(0.0, 1.05)) +
  theme_bw() + theme(panel.background = element_blank(), legend.position="none") +
  scale_fill_manual(values = c("#C5E5FB")) +
  stat_summary(geom = "text", label = c("a", "b"), fun.y = max, vjust = -1)

omega.plot.contact 

# Run Kruskal-Wallis and Dunn's Test
omega.contact.test <- kruskal.test(omega~contact, data=dnds.dat)
omega.contact.test <- dunn.test(dnds.dat$omega, dnds.dat$contact, method = "bonferroni", list = T)

#### Plot the dn/ds data ####

omega.plot.complexes <- omega.plot.complexes + guides(fill = "none")
omega.plot.contact <- omega.plot.contact + guides(fill = "none")
fig4a <- fig4a + guides(fill = "none")
fig4b <- fig4b + guides(fill = "none")
fig4c <- fig4c + guides(fill = "none")

panel <- (omega.plot + omega.plot.complexes + omega.plot.contact + plot_layout(widths = c(4, 5, 2))) / 
          (fig4a + fig4b + fig4c + plot_layout(widths = c(4, 5, 2)))


svglite::svglite(filename = "panel.svg", width = 10.3, height = 6.4)

panel + plot_layout(guides = "collect") & theme(legend.position = "bottom") &
  labs(fill = "Mitonuclear interaction") &
  plot_annotation(tag_levels = "a")

dev.off()

#### Multivariate linear model (MLM) for environmental data

library(nlme)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ape)

clim.var <- read.csv("variables.csv")

clim.table <- dist.table %>% filter(species %in% clim.var$species)

mlm.data <- merge(clim.table,clim.var,by="species")
mlm.data <- na.omit(mlm.data)

# PCA to reduce multicollinearity
env.vars <- mlm.data[, c("elevation", "long", "lat", "bio1", "bio4", "bio5", "bio6", "bio7")]

pca.env <- prcomp(env.vars, scale. = TRUE)

# View PCA loadings
# Helps to interpret how much each variable contributed to PCs
pca.env$rotation

# View importance of components
# Important to check cumulative Proportion
summary(pca.env)

# Add first principal components as predictors
mlm.data$PC1 <- pca.env$x[,1]
mlm.data$PC2 <- pca.env$x[,2] # Multiplied by -1 to improve visualization

# Ensure species names match between the tree and the dataset
mlm.tree <- keep.tip(pic.tree, mlm.data[,"species"])

# Create PICs
mlm.pic <- cbind(mt.pic = pic(mlm.data[, "mt.dist"], mlm.tree),
                 nmt.pic = pic(mlm.data[, "nmt.dist"], mlm.tree),
                 mt.nmt.rate.pic = pic(mlm.data[, "mt.nmt.rate"], mlm.tree),
                 mt.rand.rate.pic = pic(mlm.data[, "mt.rand.rate"], mlm.tree),
                 nmt.rand.rate.pic = pic(mlm.data[, "nmt.rand.rate"], mlm.tree),
                 PC1.pic = pic(mlm.data[, "PC1"], mlm.tree),
                 PC2.pic = pic(mlm.data[, "PC2"], mlm.tree))

# Make a multivariate linear model
# Dependent variables: mt, nmt, mt/nmt, mt/rand, nmt/rand
# Independent variables: PC1 and PC2 (environmental variables)
# Since PCs are already centered on zero, we use -1 to ommit the intercept
model <- lm(cbind(mt.pic, nmt.pic, mt.nmt.rate.pic, mt.rand.rate.pic, nmt.rand.rate.pic) ~ PC1.pic + PC2.pic -1, data = as.data.frame(mlm.pic))

# MANOVA test with Pillai's Trace
manova <- manova(model)
summary(manova, test = "Pillai")

# Analyzing univariate results
summary(model)

# Comparing species at the extremes
# Sort dataframe by PC1
test_by_group <- as.data.frame(mlm.pic)
test_by_group <- test_by_group[order(test_by_group$PC1.pic), ]

# Calculate 33% of samples
n_33 <- floor(42 * 0.33)

# Assign group labels (coldest 33% = "cold", hottest 33% = "hot")
test_by_group$group <- NA
test_by_group$group[1:n_33] <- "cold"
test_by_group$group[(42 - n_33 + 1):42] <- "hot"

# Subset to only include cold and hot species
df_extremes <- subset(test_by_group, group %in% c("cold", "hot"))

# Compare evolutionary rates
wilcox.test(mt.pic ~ group, data = df_extremes)
wilcox.test(nmt.pic ~ group, data = df_extremes)
wilcox.test(mt.nmt.rate.pic ~ group, data = df_extremes)
wilcox.test(mt.rand.rate.pic ~ group, data = df_extremes)
wilcox.test(nmt.rand.rate.pic ~ group, data = df_extremes)

# Convert data into long format for ggplot
data_long <- as.data.frame(mlm.pic) %>%
  select(PC1.pic, PC2.pic, mt.pic, nmt.pic, mt.nmt.rate.pic, mt.rand.rate.pic, nmt.rand.rate.pic) %>%
  pivot_longer(cols = c(mt.pic, nmt.pic, mt.nmt.rate.pic, mt.rand.rate.pic, nmt.rand.rate.pic),
               names_to = "Response_Variable", values_to = "Value")

# Plot PC1 vs response variables
pc1_plot<- ggplot(transform(data_long, Response_Variable=factor(Response_Variable,levels=c("mt.pic","nmt.pic","mt.nmt.rate.pic","mt.rand.rate.pic", "nmt.rand.rate.pic"))),
       aes(x = PC1.pic, y = Value)) +
  geom_point(size = 2, color = "firebrick4", alpha = 0.6) +  # Add points
  geom_smooth(method = "lm", se = TRUE, color = "firebrick4", fill = "firebrick4", alpha = 1/5) +  # Regression line
  facet_wrap(~Response_Variable, scales = "free_y", ncol=5) +  # Separate plots for each response
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(x = "PC1 (Environmental Gradient)", y = "Evolutionary Rate")

pc2_plot<- ggplot(transform(data_long, Response_Variable=factor(Response_Variable,levels=c("mt.pic","nmt.pic","mt.nmt.rate.pic","mt.rand.rate.pic", "nmt.rand.rate.pic"))),
                  aes(x = PC2.pic*(-1), y = Value)) +
  geom_point(size = 2, color = "firebrick4", alpha = 0.6) +  # Add points
  geom_smooth(method = "lm", se = TRUE, color = "firebrick4", fill = "firebrick4", alpha = 1/5) +  # Regression line
  facet_wrap(~Response_Variable, scales = "free_y", ncol = 5) +  # Separate plots for each response
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(x = "PC2 (Environmental Gradient)", y = "Evolutionary Rate")

svglite::svglite(filename = "environment_panel.svg", width = 10.5, height = 4.5)

pc1_plot/pc2_plot

dev.off()
