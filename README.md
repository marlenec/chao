# chao
####################################################################################
# FUNCTION to compute Chao's alpha and beta diversity, expressed as Hill Numbers,
# for a set of communities
#
# This function requires the installation of packages entropart (REF PACKAGE  + papier MARCON)
#
# Arguments :
# - matrix : a matrix with raw or relative abundances of species (columns) in
#            communities (rows)
# -tree_phylog (optional) : a rooted ultrametric phylogenetic tree in phylog format
#                           see ?newick2phylog and nota bene below for explanation on phylog format
#                           The tip labels and the column names of the matrix must be identical ;
#                           but not necessarely in the same order.
#
# - q : value(s) for the q parameter, among 0,1,2. Default is q=c(0,1,2)
#       See Chao et al. 2014 for explanations about the q parameter.
# - run_example : if TRUE, a simple example provided at the end of this script is runned.
#
#
# NB : The conversion of a tree from the phylo format into the phylog format is easy :
# write.tree(tree, file="tree.txt") # 1°/ you need to make a .txt file in a newick format
# tree<-paste(readLines('tree.txt')) # 2°/ you need to read it as a character string
# tree_phylog<-newick2phylog(tree) # 3°/ newick2phylog{ade4} converts the character string into a phylog object
#
#
# Details :
# - The taxonomic and phylogenetic diversities are calculated following Chao et al (2014). Please see the
#   entropart documentation for more details.
# - If tree_phylog is provided, the phylogenetic diversity is calculated.
# - Chao's phylogenetic diversity can be calculated on non-ultrametric trees, but the meaning is not clear.
#   If a non-ultrametric is provided, phylogenetic diversity will be calculated an returned with a warning.
# - If raw abundances are provided, unequal number of indivdi is allowed. In this case, the total abundance
#   in each community will be taken into account when calculating pairwise beta diversity.
#
#
# OUTPUTS : a list with for each biodiversity facet (i.e. taxo and phylo):
# - $alpha_taxo: matrix of taxonomic alpha diversity values for each value of q (column) for all communities
#   (rows);
# - $alpha_phylo: matrix of phylogenetic alpha diversity values for each value of q (column) for all communities
#   (rows);
# - $beta_taxo: up to three matrices containing pairwise taxonomic beta diversity for each value of q;
# - $beta_phylo: up to three matrices containing pairwise phylogenetic beta diversity for each value of q;
#
#
# NB : The conversion of a tree from the phylo format into the phylog format is easy :
# write.tree(tree, file="tree.txt") # 1°/ you need to make a .txt file in a newick format
# tree<-paste(readLines('tree.txt')) # 2°/ you need to read it as a character string
# tree_phylog<-newick2phylog(tree) # 3°/ newick2phylog{ade4} converts the character string into a phylog object
#
####################################################################################
####################################################################################
