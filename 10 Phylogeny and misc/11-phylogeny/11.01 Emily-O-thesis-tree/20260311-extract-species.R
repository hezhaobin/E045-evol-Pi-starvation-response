# title: extract species for plotting a species tree
# author: Bin He
# date: 2026-03-11

# Load libraries
require(tidyverse)
require(ggtree)
require(treeio)

# Use TreeHouse
## running the app
shiny::runGitHub("treehouse", "JLSteenwyk")

# Read tree
spsInfo <- read_tsv("species-info.tsv", col_types = cols())
sps.tree <- read.tree("20260311-extracted-species-tree.nwk") %>% 
  as_tibble() %>% 
  #mutate(label = gsub("_", " ", label)) %>% 
  left_join(spsInfo, by = c("label" = "treeName")) %>% 
  as.treedata()

# Plot the tree
p.tree <- ggtree(sps.tree, ladderize = TRUE, branch.length = "branch.length") +
  theme_tree2() +
  scale_y_discrete(limits = rev) +
  geom_tiplab(aes(label = species), size = 10, face = 3, as_ylab = TRUE) +
  #geom_tiplab(size = 3.2, fontface = "italic", align = TRUE, linesize = 0.1, offset = 0.05) +
  #geom_treescale() +
  #geom_tippoint(aes(color = pathogen), size = 2) +
  theme(legend.position = c(0.227, 0.80))
p.tree
ggsave(p.tree, file = "20260311-Emily-thesis-species-tree.png", width = 4, height = 4)
