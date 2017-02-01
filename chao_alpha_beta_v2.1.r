####################################################################################
# FUNCTION to compute Chao's alpha and beta diversity, expressed as Hill Numbers,
# for a set of communities
#
# This function requires the installation of packages entropart (Marcon & Hérault, 2015 ; Marcon et al.,
# PloS One, 2014)
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
#
# - run_example : if TRUE, a simple example provided at the end of this script is runned.
#
#
# Details :
# - The taxonomic and phylogenetic diversities are calculated following Chao et al (2014). Please see the
#   entropart documentation for more details.
# - If tree_phylog is provided, the phylogenetic diversity is calculated.
# - Chao's phylogenetic diversity can be calculated on non-ultrametric trees, but the meaning is not clear.
#   If a non-ultrametric is provided, phylogenetic diversity will be calculated an returned with a warning.
# - If raw abundances are provided, unequal numbers of individuals between communities are allowed. In this case,
#   the total abundance in each community will be taken into account when calculating pairwise beta diversity.
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
rm(list=ls())

chao_alpha_beta <- function(matrix, tree_phylog=NULL, q=c(0,1,2), beta=TRUE, run_example=F)
{

# check function arguments
if(missing(matrix)) stop("A matrix with relative abundance data must be provided")
    if(is.null(tree_phylog)==F & class(tree_phylog) != "phylog") stop("Please convert the tree into a 'phylog' format")
    if(any(is.na(matrix))) stop("The abundances matrix contains NA(s). Please check")
  
# load libraries
require(ape)
require(ade4)
require(entropart)
  
  ########################## Alpha and beta phylo
  if(is.null(tree_phylog)==F){
      
      if(is.ultrametric(tree_phylog==F)) warning("Your tree is not ultrametric")

alpha_matrix<-matrix(NA,nrow=length(rownames(matrix)), ncol=length(q), dimnames=list(rownames(matrix),paste0("q",q)))
    
    for(i in 1:length(rownames(matrix)))
    {

      rel_abundances_i<-matrix[i,]/sum(matrix[i,])
      
      for(j in 1:length(q))
      {
        alpha_matrix[i,j]<-ChaoPD(rel_abundances_i, q=q[j],tree_phylog,  CheckArguments=T)
      } # end of j
    } # end of i
    
    if(beta==T){ # if.beta
        
    combin<-combn(rownames(matrix),2)
    
    res_beta<-lapply(1:length(q), function(x) matrix(0,nrow=length(rownames(matrix)), ncol=length(rownames(matrix)), dimnames=list(rownames(matrix),rownames(matrix))))
    names(res_beta)<-paste0("q",q)
    
    for(k in 1:(length(combin)/2))  
    {
      com_1<-combin[1,k]
      com_2<-combin[2,k]
      
      for(l in 1:length(q))
      {
      alpha_com_1<-alpha_matrix[com_1,paste0("q",q[l])]
      alpha_com_2<-alpha_matrix[com_2,paste0("q",q[l])]
      
      weight_com_1<-sum(matrix[com_1,])
      weight_com_2<-sum(matrix[com_2,])
      
      weighted_abundances_com_1<-(matrix[com_1,]*weight_com_1)/(weight_com_1+weight_com_2)
      weighted_abundances_com_2<-(matrix[com_2,]*weight_com_2)/(weight_com_1+weight_com_2)

      rel_abundances_2_com<-rbind(weighted_abundances_com_1, weighted_abundances_com_2)
      rel_abundances_2_com<-apply(rel_abundances_2_com, 2, sum)
      rel_abundances_2_com<-rel_abundances_2_com/sum(rel_abundances_2_com)
      
      gamma<-ChaoPD(rel_abundances_2_com, q=q[l], tree_phylog, CheckArguments=F)
      
      if(q[l]==0) #2
      {
        mean_alpha<-alpha_com_1*(weight_com_1/(weight_com_1+weight_com_2))+alpha_com_2*(weight_com_2/(weight_com_1+weight_com_2))
        
        res_beta$q0[com_1, com_2]<-(gamma/mean_alpha)-1
        res_beta$q0[com_2, com_1]<-res_beta$q0[com_1, com_2]
      } # end of if() #2
      
      if(q[l]==1) #3
      {
        alpha_com_1<-log(alpha_com_1)
        alpha_com_2<-log(alpha_com_2)
        
        mean_alpha<-alpha_com_1*(weight_com_1/(weight_com_1+weight_com_2))+alpha_com_2*(weight_com_2/(weight_com_1+weight_com_2))
        
        res_beta$q1[com_1, com_2]<-(gamma/exp(mean_alpha))-1
        res_beta$q1[com_2, com_1]<-res_beta$q1[com_1, com_2]
      } # end of if() #3
      
      if(q[l]==2) #4
      {
        alpha_com_1<-1-1/alpha_com_1
        alpha_com_2<-1-1/alpha_com_2
        
        mean_alpha<-alpha_com_1*(weight_com_1/(weight_com_1+weight_com_2))+alpha_com_2*(weight_com_2/(weight_com_1+weight_com_2))

        res_beta$q2[com_1, com_2]<-(gamma/(1/(1-mean_alpha)))-1
        res_beta$q2[com_2, com_1]<-res_beta$q2[com_1, com_2]
      } # end of if() #4      
      
      } # end of l
    } # end of k
    } # end of if.beta
  } # end of if.null(tree_phylog)
  
  ################################# Alpha and beta taxo
  alpha_matrix_taxo<-matrix(NA,nrow=length(rownames(matrix)), ncol=length(q), dimnames=list(rownames(matrix),paste0("q",q)))
  
  for(i in 1:length(rownames(matrix)))
  {
    
    rel_abundances_i<-matrix[i,]/sum(matrix[i,])
    
    for(j in 1:length(q))
    {
      if(q[j]==0) #if1
      {
        alpha_matrix_taxo[i,j]<-Tsallis(rel_abundances_i, q = 0, CheckArguments=T)
      } # end of if1
      
      if(q[j]==1) # if2
      {
        alpha_matrix_taxo[i,j]<-exp(Tsallis(rel_abundances_i, q = 1, CheckArguments=T))
      } # end of if2
      
      if(q[j]==2) # if3
      {
        alpha_matrix_taxo[i,j]<-1/(1-Tsallis(rel_abundances_i, q = 2, CheckArguments=F))
      } # if3
    } # end of j
  } # end of i
  
  combin<-combn(rownames(matrix),2)
  
  if(beta==T){ # if.beta
  res_beta_taxo<-lapply(1:length(q), function(x) matrix(0,nrow=length(rownames(matrix)), ncol=length(rownames(matrix)), dimnames=list(rownames(matrix),rownames(matrix))))
  names(res_beta_taxo)<-paste0("q",q)
  
  for(k in 1:(length(combin)/2))  
  {
    com_1<-combin[1,k]
    com_2<-combin[2,k]
    
    for(l in 1:length(q))
    {
      alpha_com_1<-alpha_matrix_taxo[com_1,paste0("q",q[l])]
      alpha_com_2<-alpha_matrix_taxo[com_2,paste0("q",q[l])]
      
      weight_com_1<-sum(matrix[com_1,])
      weight_com_2<-sum(matrix[com_2,])
      
      weighted_abundances_com_1<-(matrix[com_1,]*weight_com_1)/(weight_com_1+weight_com_2)
      weighted_abundances_com_2<-(matrix[com_2,]*weight_com_2)/(weight_com_1+weight_com_2)
      
      rel_abundances_2_com<-rbind(weighted_abundances_com_1, weighted_abundances_com_2)
      rel_abundances_2_com<-apply(rel_abundances_2_com, 2, sum)
      rel_abundances_2_com<-rel_abundances_2_com/sum(rel_abundances_2_com)
      
      if(q[l]==0) #2
      {
        mean_alpha<-alpha_com_1*(weight_com_1/(weight_com_1+weight_com_2))+alpha_com_2*(weight_com_2/(weight_com_1+weight_com_2))
          
        gamma<-Tsallis(rel_abundances_2_com, q = 0, CheckArguments=T)
        res_beta_taxo$q0[com_1, com_2]<-(gamma/mean_alpha)-1
        res_beta_taxo$q0[com_2, com_1]<-res_beta_taxo$q0[com_1, com_2]
      } # end of if() #2
      
      if(q[l]==1) #3
      {
        alpha_com_1<-log(alpha_com_1)
        alpha_com_2<-log(alpha_com_2)
          
          mean_alpha<-alpha_com_1*(weight_com_1/(weight_com_1+weight_com_2))+alpha_com_2*(weight_com_2/(weight_com_1+weight_com_2))
          
        gamma<-exp(Tsallis(rel_abundances_2_com, q = 1, CheckArguments=F))
        res_beta_taxo$q1[com_1, com_2]<-(gamma/exp(mean_alpha))-1
        res_beta_taxo$q1[com_2, com_1]<-res_beta_taxo$q1[com_1, com_2]
      } # end of if() #3
      
      if(q[l]==2) #4
      {
        alpha_com_1<-1-1/alpha_com_1
        alpha_com_2<-1-1/alpha_com_2
          
        mean_alpha<-alpha_com_1*(weight_com_1/(weight_com_1+weight_com_2))+alpha_com_2*(weight_com_2/(weight_com_1+weight_com_2))
          
        gamma<-1/(1-Tsallis(rel_abundances_2_com, q = 2, CheckArguments=F))
        res_beta_taxo$q2[com_1, com_2]<-(gamma/(1/(1-mean_alpha)))-1
        res_beta_taxo$q2[com_2, com_1]<-res_beta_taxo$q2[com_1, com_2]
      } # end of if() #4      
      
    } # end of l
  } # end of k
  } # end of if.beta
  
  if(is.null(tree_phylog)==F) {
	if(beta==T){
  res_chao_alpha_beta<-list(alpha_matrix_taxo,res_beta_taxo, alpha_matrix, res_beta)
  names(res_chao_alpha_beta)<-c("alpha_taxo","beta_taxo", "alpha_phylo", "beta_phylo")
  return(res_chao_alpha_beta)
	} # end of if.beta
	else {
  res_chao_alpha_beta<-list(alpha_matrix_taxo, alpha_matrix)
  names(res_chao_alpha_beta)<-c("alpha_taxo", "alpha_phylo")
  return(res_chao_alpha_beta)
	} # end of else
  } # end of is.null(tree_phylog)==F
  
  if(is.null(tree_phylog)) {
	if(beta==T)
	{
  res_chao_alpha_beta<-list(alpha_matrix_taxo,res_beta_taxo)
  names(res_chao_alpha_beta)<-c("alpha_taxo","beta_taxo")
  return(res_chao_alpha_beta)
	} # end of if.beta
	else {
  res_chao_alpha_beta<-list(alpha_matrix_taxo)
  names(res_chao_alpha_beta)<-c("alpha_taxo")
  return(res_chao_alpha_beta)
	} # end of else
  } # end of if.taxo
  #####################################

if( run_example==TRUE) {
  

# INPUTS
matrix<-matrix(NA, nrow=7, ncol=5)
matrix[1,]<-c(1,2,3,4,5)
matrix[2,]<-c(5,4,3,2,0)
matrix[3,]<-c(3,2,1,4,5)
matrix[4,]<-c(2,0,5,4,1)
matrix[5,]<-c(4,3,2,0,5)
matrix[6,]<-c(8,3,2,0,5)
matrix[7,]<-c(4,3,2,0,10)


rownames(matrix)<-c("Community_A", "Community_B", "Community_C", "Community_D", "Community_E", "Community_F", "Community_G")
colnames(matrix)<-paste0("OTU_", 1:5)

 library(ape)
 library(ade4)
tree<-rcoal(n=5, rooted=T, tip.label=colnames(matrix))
write.tree(tree, file="tree.txt")
tree<-paste(readLines('tree.txt'))
tree_phylog<-newick2phylog(tree)

} #end of run_example
  
} # end of chao_alpha_beta

####################################################################################
# Example

# EXECUTE THE FUNCTION
#test<-chao_alpha_beta(matrix, tree_phylog)
#test$alpha_phylo # alpha diversit values in a matrix format
#test$alpha_taxo
#test$beta_phylo$q0 # beta values for q=0
#test2<-chao_alpha_beta(matrix)
#test2

####################################################################################