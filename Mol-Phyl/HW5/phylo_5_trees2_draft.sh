iqtree2 -s SUP35_aln_prank.trim.fas -m TIM3+F+G4 -pre SUP35_TIM3

time iqtree2 -s SUP35_aln_prank.trim.fas -m TIM3+F+G4 -redo -pre SUP35_TIM3_b -b 10

time iqtree2 -s SUP35_aln_prank.trim.fas -m TIM3+F+G4 -redo -pre SUP35_TIM3_ufb -bb 1000

iqtree2 -s SUP35_aln_prank.trim.fas -m TIM3+F+G4 -pre SUP35_TIM3_B_alrt_abayes -bb 1000 -alrt 1000 -abayes

#https://itol.embl.de/ 

```{r}
library(ggtree)
tree_alrt_abayes_ufb <- read.tree("SUP35_TIM3_B_alrt_abayes.treefile")
ggtree(tree_alrt_abayes_ufb) + 
  geom_tiplab() + geom_nodelab() +
  geom_treescale() + xlim(0, 0.7)
# funny labels
label <- tree_alrt_abayes_ufb$node.label
alrt <- as.numeric(sapply(strsplit(label, "/"), FUN = "[", 1)) #sun
abayes <- as.numeric(sapply(strsplit(label, "/"), FUN = "[", 2)) #yin yang
ufb <- as.numeric(sapply(strsplit(label, "/"), FUN = "[", 3)) #star
large_alrt <- ifelse(alrt > 70, intToUtf8(9728), "")
large_abayes <- ifelse(abayes > 0.7, intToUtf8(9775), "")
large_ufb <- ifelse(ufb > 95, intToUtf8(9733), "")
newlabel <- paste0(large_alrt, large_abayes, large_ufb)
tree_alrt_abayes_ufb$node.label <- newlabel
ggtree(tree_alrt_abayes_ufb) + 
  geom_tiplab() + geom_nodelab(nudge_x = -.01, nudge_y = .1) +
  geom_treescale() + xlim(0, 0.7)
``` 

iqtree2 -s SUP35_aln_prank.trim.fas -m TIM3+F+G4 -pre SUP35_TIM3_root_outgroup -bb 1000 -alrt 1000 -abayes  -o SUP35_Kla_AB039749,SUP35_Agos_ATCC_10895_NM_211584


#https://github.com/mooreryan/midpoint-root
```{r}
#install.packages("phytools")
library(phytools)
midpoint.root(tree_alrt_abayes_ufb)
```

iqtree2 -s SUP35_aln_prank.trim.fas -m TIM3+F+G4 -pre SUP35_TIM3_root_auto --model-joint 12.12 -B 1000

iqtree2 -s SUP35_aln_prank.trim.fas -m JC -pre SUP35_JC -bb 1000 -alrt 1000 -abayes -o SUP35_Kla_AB039749,SUP35_Agos_ATCC_10895_NM_211584

#http://phylo.io/ 
#https://beta.phylo.io/viewer/#

```{r}
library(ggtree)
treeTIM3 <- read.tree("SUP35_TIM3_root_outgroup.treefile")
treeJC <- read.tree("SUP35_JC.treefile")
library(ggplot2)
tim3tree <- 
  ggtree(treeTIM3) + geom_tiplab() +
  geom_nodelab(color = "blue", nudge_x = -.05) + 
  theme_tree2() + 
  xlim(0,1) + ggtitle("TIM3")
jctree <- 
  ggtree(treeJC) + geom_tiplab() +
  geom_nodelab(color = "red", nudge_x = -.05) + 
  theme_tree2() + 
  xlim(0,1) + ggtitle("JC")
library(ggpubr)
ggarrange(tim3tree, jctree)
association <- cbind(treeTIM3$tip.label, treeJC$tip.label)
cophyloplot(treeTIM3, treeJC, assoc=association, length.line=4, space=28, gap=3)
library(phytools)
trees.cophylo<-cophylo(treeTIM3, treeJC, rotate = TRUE)
png("cophylo.png", width = 1200, height = 800)
plot(trees.cophylo, link.type="curved",link.lwd=4,
     link.lty="solid",link.col=make.transparent("red", 0.25), size = 1)
dev.off()

