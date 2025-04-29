## 1. Obtaining a UCE species tree
### 1.1. UCE extraction, processing, and alignment
We have followed the exact steps of the tutorials available at the [Phyluce documentation page](https://phyluce.readthedocs.io/en/latest/tutorials/index.html).
### 1.2. Inferring UCE gene trees with RAxML
```bash
for i in uce-75p/*.fasta; do

file=$(basename "$i")

raxmlHPC-PTHREADS-SSE3 -T 32 -${i} -n ${file%.fasta} -m GTRGAMMAI -# 1000 -p ${RANDOM} -f a -x ${RANDOM} 

done
```
### 1.3. Inferring a coalescent-based species tree with ASTRAL
```bash
java -jar astral.5.7.8.jar -i all_uce_trees.nwk -o uce_astral.nwk
```
### 1.4. Rescaling the ASTRAL tree to proportional branch lengths
Obtains a tree from the concatenated UCE alignment, setting the ASTRAL tree as a topological constraint.
```bash
raxmlHPC-PTHREADS-SSE3 -T 8 -s mafft-nexus-edge-trimmed-gblocks-clean-75p-raxml.phylip \
  -n uce_rescaled_astral -m GTRGAMMAI -o mellifera \
  -# 1000 -p ${RANDOM} -f a -x ${RANDOM} -g uce_astral.nwk 
```
# 2. Unconstrained phylogenies
Obtained for the comparison of phylogenetic signal among datasets.
```bash
# Using the random nuclear gene dataset as an example

raxmlHPC-PTHREADS-SSE3 -T 8 -s random_AA.fasta -n random_unconstrained  \
  -m PROTGAMMAAUTO -o mellifera -# 1000 -p ${RANDOM} -f a -x ${RANDOM}
```
# 3. Tree topology comparison
```R
setwd("Mitonuclear/discordance_test")

library(phangorn)
library(TreeDist)

Nmt <- read.tree("N-mt_unconstrained.nwk")
mt <- read.tree("mt_unconstrained.nwk")
rand <- read.tree("random_unconstrained.nwk")
gly <- read.tree("glyco_unconstrained.nwk")

# Comparing tree topologies using
# Information-based generalized Robinson-Foulds distances
# as proposed in Smith (2020)
# (package TreeDist)
# Full concordance at value of 0, full discordance at a value of 1

Nmt_mt <- TreeDistance(Nmt, mt)
Nmt_rand <- TreeDistance(Nmt, rand)
mt_rand <- TreeDistance(mt, rand)
Nmt_gly <- TreeDistance(Nmt, gly)
mt_gly <- TreeDistance(mt, gly)
gly_rand <- TreeDistance(gly, rand)

Nmt_mt
Nmt_rand
mt_rand
Nmt_gly
mt_gly
gly_rand

# Plot this into a heatmap
library(ggplot2)

df <- read.delim("tree_distance_heatmap.txt")
df$TreeDistance <- round(df$TreeDistance, digits = 2)

svglite::svglite(filename = "distance_heatmap.svg", width = 6, height = 4)

ggplot(df, aes(x = dataset1, y = dataset2, fill = TreeDistance)) +
  geom_tile(color = "white", lwd = 1, linetype = 1) + geom_text(aes(label = TreeDistance)) +
  theme_minimal() + 
  theme(panel.grid = element_blank()) + #theme(legend.position="none") +
  scale_fill_distiller(direction = +1) +
  labs(x = NULL, y = NULL)

dev.off()
```
