library(ggplot2)
library(tibble)
library(vegan)
library(dplyr)
library(rio)
library(zip)
library(car)
library(MASS)
library(ggeffects)
library(viridis)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(swamp)
library(dendsort)
library(cluster)
library(ClusterR)
library(factoextra)
library(ggsignif)
library(rstatix)
library(wordspace)
library(ggforce)
library(ggalt)
library(chisq.posthoc.test)
library(usedist)
library(dendextend)
library(flextable)
library(ggpubr)
library(XNomial)
library(emstreeR)
library(fpc)
library(tidyverse)
library(plyr)
library(pscl)
library(boot)
library(sjPlot)

############################# cluster compounds from all data - more data = better clustering.####################################
# Identify optimal cluster number 
fviz_nbclust(fing_diss_matx,cluster::pam,method="silhouette",k.max=20) # suggests 2
fviz_nbclust(fing_diss_matx,cluster::pam,method="wss",k.max=20) # suggests around 9
fviz_nbclust(fing_diss_matx,cluster::pam,method="gap_stat") # suggests 10
# k = 9

chem_clust = pam(fing_diss_matx,k=9) # cluster compounds according to k-medoids clustering w. optim. cluster no.  
# plot clustering 
clust_plot = fviz_cluster(chem_clust,geom="point",ggtheme = theme_bw(),main="Cluster Plot of Chemical Distances of Compounds") # looks good
chemical_coordinates = data.frame(clust_plot$data[c("name","x","y","cluster")],row.names=NULL) #extract coordinates generated from the distance matrix for future plotting 
# save file
saveRDS(chemical_coordinates,"chemical_coordinates_gita.rds")
chemical_coordinates = readRDS("chemical_coordinates_gita.rds")


#####################################################################################################################################

################################################## Creation of matrices to be analysed ########################################################
# matrix of non-native species
nnat_spec = compiled_pa_matrix[which(compiled_pa_matrix$Species %in% c("Lantana camara","Melaleuca quinquenervia","Psidium cattleianum",
                                                                       "Ageratum conyzoides","Achillea millefolium","Rosmarinus officinalis",
                                                                       "Tagetes erecta")),]
nnat_spec = nnat_spec[,c(1,2,c((ncol(nnat_spec)-3):ncol(nnat_spec)),
                         (which(colSums(nnat_spec[3:(ncol(nnat_spec)-4)]) > 0)+2))] # remove compounds absent from any profile

# matrix of non-native families
nnat_fam = compiled_pa_matrix[which(compiled_pa_matrix$Family %in% c("Verbenaceae","Myrtaceae","Asteraceae","Compositae","Lamiaceae")),]
nnat_fam = nnat_fam[,c(1,2,c((ncol(nnat_fam)-3):ncol(nnat_fam)),
                        (which(colSums(nnat_fam[3:(ncol(nnat_fam)-4)]) > 0)+2))] # remove compounds absent from any profile


# Do each step of the analysis 1st by species, 2nd by family
# Do invasive vs non-invasive for L. camara, M. quinquenervia and P. cattleainum (+ fam?), non-native vs native for other species
# compare results between
unique(nnat_fam$Species[which(nnat_fam$Family %in% c("Verbenaceae","Myrtaceae"))]) # this is fine, other family members are not in the non-native species set

# rename compositae to asteraceae 
nnat_fam$Family[which(nnat_fam$Family %in% "Compositae")] = "Asteraceae"

# get important vectors for loops and functions
species = unique(nnat_spec$Species)
family = unique(nnat_fam$Family)
compound_s = colnames(nnat_spec[7:ncol(nnat_spec)])
compound_col_s = 7:ncol(nnat_spec)
compound_f = colnames(nnat_fam[7:ncol(nnat_fam)])
compound_col_f = 7:ncol(nnat_fam)
profile_s = nnat_spec$Profile
profile_f = nnat_fam$Profile


#############################################################################################################################################

####################### Profile plots (presence and absence)  ############################
# plots are generated as heatmaps to look like a BLAST plot
# colouring indicates presence of a compounds
# profiles are coloured according to whether they are invasive (purple) /non-invasive (green) for all species but A. conyzoides, which are coloured by native and non-native status


colbreak_proplot = c(0,0.99,1.99,2) # split colouring so that values 0 = black, 1 = green and 2 = purple
cols = c("black","#66CC00","#CC66CC")

## for species
for(spec in species) {
  
  if(spec == species[1] | spec == species[2] | spec == species[3]) {
    var_newcolor = "invasive_status" # colour according to invasivity
  } else {
    var_newcolor = "native_status" # colour according to native status
  }
  submat = nnat_spec[which(nnat_spec$Species %in% spec),compound_col_s] # gets compound data for each species
  rownames(submat) = nnat_spec$Profile[which(nnat_spec$Species %in% spec)] # rownames = profile no. for the df
  submat = submat[,which(colSums(submat) != 0)] # removes any absent compounds from the df
  
  # current presence/absence data is binary
  # non-invasive/native no. stays 1, presence of compounds in invasive/non-native profiles is converted to 2 for plotting
  recolour = nnat_spec$Profile[which(nnat_spec[[var_newcolor]] %in% 1 & nnat_spec$Species %in% spec)] # identify all invasive/non-native profiles
  for(pro in recolour) {
    row = which(rownames(submat) == pro)
    submat[row,which(submat[row,] %in% 1)] = 2 # set all present compounds to equal 2 
  }
  
  # reorder data so that all invasive/non-native profiles are grouped together and all native/non-invasive profiles are together
  rowinds = which(rownames(submat) %in% recolour) # identifies invasive/non-native profile indices 
  submat = as.data.frame(rbind(submat[rowinds,],submat[-c(rowinds),])) # invasive goes first, non-invasive after
  # plot heatmap
  pheatmap(t(submat),cluster_rows=FALSE,cluster_cols=FALSE,color=cols,breaks=colbreak_proplot,border_color = "black",show_rownames=FALSE,show_colnames=FALSE)
}
# from inspection, it's not clear if the new invasive species can be included - balance or native/non-native is poor

## for family
for(spec in family) {

  if(spec == family[1] | spec == family[2]) {
    var_newcolor = "invasive_status" # colour according to invasivity
  } else {
    var_newcolor = "native_status" # colour according to native status
  }
  submat = nnat_fam[which(nnat_fam$Family %in% spec),compound_col_f] # gets compound data for each species
  rownames(submat) = nnat_fam$Profile[which(nnat_fam$Family %in% spec)] # rownames = profile no. for the df
  submat = submat[,which(colSums(submat) != 0)] # removes any absent compounds from the df
  
  # current presence/absence data is binary
  # non-invasive/native no. stays 1, presence of compounds in invasive/non-native profiles is converted to 2 for plotting
  recolour = nnat_fam$Profile[which(nnat_fam[[var_newcolor]] %in% 1 & nnat_fam$Family %in% spec)] # identify all invasive/non-native profiles
  for(pro in recolour) {
    row = which(rownames(submat) == pro)
    submat[row,which(submat[row,] %in% 1)] = 2 # set all present compounds to equal 2 
  }
  
  # reorder data so that all invasive/non-native profiles are grouped together and all native/non-invasive profiles are together
  rowinds = which(rownames(submat) %in% recolour) # identifies invasive/non-native profile indices 
  submat = as.data.frame(rbind(submat[rowinds,],submat[-c(rowinds),])) # invasive goes first, non-invasive after
  # plot heatmap
  pheatmap(t(submat),cluster_rows=FALSE,cluster_cols=FALSE,color=cols,breaks=colbreak_proplot,border_color = "black",show_rownames=FALSE,show_colnames=FALSE,
           main=bquote(italic(.(spec))))
}
# family data looks a lot better - can definitely be included except for the last 

########################################### Updating Matrices #####################################################
species = species[c(1,2,3,4,6)]
nnat_spec = nnat_spec[which(nnat_spec$Species %in% species),]
nnat_spec = nnat_spec[,c(1:6,which(colSums(nnat_spec[7:ncol(nnat_spec)]) > 0)+6)]
nnat_spec$invasive_status = as.character(nnat_spec$invasive_status)
nnat_spec$native_status = as.character(nnat_spec$native_status)

family = family[c(1,2,3,4)]
nnat_fam = nnat_fam[which(nnat_fam$Family %in% family),]
nnat_fam = nnat_fam[,c(1:6,which(colSums(nnat_fam[7:ncol(nnat_fam)]) > 0)+6)]
nnat_fam$invasive_status = as.character(nnat_fam$invasive_status)
nnat_fam$native_status = as.character(nnat_fam$native_status)

compound_s = colnames(nnat_spec[7:ncol(nnat_spec)])
compound_col_s = 7:ncol(nnat_spec)
compound_f = colnames(nnat_fam[7:ncol(nnat_fam)])
compound_col_f = 7:ncol(nnat_fam)
profile_s = nnat_spec$Profile
profile_f = nnat_fam$Profile

saveRDS(nnat_spec,"nnat_spec.rds")
saveRDS(nnat_fam,"nnat_fam.rds")

# I think from now on focus on family matrix
#######################################################################################################################

############################# cluster compounds from family matrix ####################################
# Identify optimal cluster number 
fing_diss_matx_fam = fing_diss_matx[which(rownames(fing_diss_matx) %in% compound_f),
                                    which(colnames(fing_diss_matx) %in% compound_f)]

fviz_nbclust(fing_diss_matx_fam,cluster::pam,method="silhouette") # 2
fviz_nbclust(fing_diss_matx_fam,cluster::pam,method="wss") # 6
fviz_nbclust(fing_diss_matx_fam,cluster::pam,method="gap_stat") # 10
# k = 6? middle ground

chem_clust_fam = pam(fing_diss_matx_fam,k=6) # cluster compounds according to k-medoids clustering w. optim. cluster no.  
# plot clustering 
clust_plot_fam = fviz_cluster(chem_clust_fam,geom="point",ggtheme = theme_bw(),main="Clustering of Compounds Present in Profiles by Family") # looks good
chemical_coordinates_fam = data.frame(clust_plot_fam$data[c("name","x","y","cluster")],row.names=NULL) #extract coordinates generated from the distance matrix for future plotting 
# save file
saveRDS(chemical_coordinates_fam,"chemical_coordinates_fam.rds")
chemical_coordinates_fam = readRDS("chemical_coordinates_fam.rds")

# looks quite good overall, perhaps worth using
#####################################################################################################################################

################################ First Data ######################################################################

# Sample size per family

for(fam in family) {
  
  if(fam == family[1] | fam == family[2]) {
    
    sam_1 = nrow(nnat_fam[which(nnat_fam$Family %in% fam & nnat_fam$invasive_status %in% "1"),]) # get sample size
    sam_0 = nrow(nnat_fam[which(nnat_fam$Family %in% fam & nnat_fam$invasive_status %in% "0"),])
    
  } else {
    
    sam_1 = nrow(nnat_fam[which(nnat_fam$Family %in% fam & nnat_fam$native_status %in% "1"),])
    sam_0 = nrow(nnat_fam[which(nnat_fam$Family %in% fam & nnat_fam$native_status %in% "0"),])
    
  }
  
  print(c(sam_1,sam_0))

}

# Total chemodiversity per family

for(fam in family) {
  
  if(fam == family[1] | fam == family[2]) {
    
    sub_data = nnat_fam[which(nnat_fam$Family %in% fam & nnat_fam$invasive_status %in% "1"),compound_col_f]
    sub_data2 = nnat_fam[which(nnat_fam$Family %in% fam & nnat_fam$invasive_status %in% "0"),compound_col_f]
    
  } else {
    
    sub_data = nnat_fam[which(nnat_fam$Family %in% fam & nnat_fam$native_status %in% "1"),compound_col_f]
    sub_data2 = nnat_fam[which(nnat_fam$Family %in% fam & nnat_fam$native_status %in% "0"),compound_col_f]
      
  }
  
  chem_1 = length(which(colSums(sub_data) > 0)) # get total no. unique compounds per fam per population
  chem_0 = length(which(colSums(sub_data2) > 0))
  
  print(c(chem_1,chem_0))
  
}


##########################################################################################################

########################################### Profile Size #############################################

profile_sizes = NULL

for(row in 1:nrow(nnat_fam)) {
  
  pro_size = sum(nnat_fam[row,compound_col_f])
  
  data_add = c(nnat_fam$Profile[row],nnat_fam$Family[row],nnat_fam$invasive_status[row],
               nnat_fam$native_status[row],pro_size)
  
  profile_sizes = rbind(profile_sizes,data_add)
  
}

profile_sizes = as.data.frame(profile_sizes,row.names=FALSE)
colnames(profile_sizes) = c("Profile","Family","invasive_status","native_status","Profile_Size")
profile_sizes$Profile_Size = as.numeric(profile_sizes$Profile_Size)

for(fam in family) {
  
  
  if(fam == family[1] | fam == family[2]) {
    
    print(tapply(profile_sizes$Profile_Size[which(profile_sizes$Family %in% fam)],
           profile_sizes$invasive_status[which(profile_sizes$Family %in% fam)],
           summary)) 
    
    print(kruskal.test(Profile_Size ~ invasive_status,profile_sizes[which(profile_sizes$Family %in% fam),]))
    
  } else {
    
    print(tapply(profile_sizes$Profile_Size[which(profile_sizes$Family %in% fam)],
           profile_sizes$native_status[which(profile_sizes$Family %in% fam)],
           summary)) 
    
    print(kruskal.test(Profile_Size ~ native_status,profile_sizes[which(profile_sizes$Family %in% fam),]))
    
  }
}


# K-W tests for comparing emission profile size between and within groups
kruskal.test(Profile_Size~invasive_status,profile_sizes) # not signif
kruskal.test(Profile_Size~native_status,profile_sizes) # not signif
kruskal.test(Profile_Size~Family,profile_sizes) # almost signif (0.06) - could be difference between families

inv_data_kw = c(nnat_fam$invasive_status[which(nnat_fam$Family %in% family[1] | nnat_fam$Family %in% family[2])],
                nnat_fam$native_status[which(nnat_fam$Family %in% family[3])])
profile_sizes = cbind(profile_sizes,inv_data_kw)
kruskal.test(Profile_Size~inv_data_kw,data=profile_sizes)

##################################################### Cluster Profiles by P/A ############################################

pa_clusters_f = list()

for(fam in family) {
  
  # create dendograms
  
  submat = nnat_fam[which(nnat_fam$Family %in% fam),compound_col_f] # matrix of all profiles of the specified species, data = compounds only
  rownames(submat) = nnat_fam$Profile[which(nnat_fam$Family %in% fam)] # set rownames as profile number
  submat = submat[,which(colSums(submat) != 0)] # remove any compounds that are not present in any chemical profile of this subset
  
  hc_cut = hcut(submat,hc_func = "hclust",k=2,hc_metric = "euclidean", hc_method = "ward.D2") # cluster matrix by hierachical clustering using Ward Clustering and Euclidean distances, cut into 2 groups  
  
  if(fam == family[1] | fam == family[2]) {
    
    annotate = nnat_fam$invasive_status[which(nnat_fam$Family %in% fam)] 
    
  } else {
    
    annotate = nnat_fam$native_status[which(nnat_fam$Family %in% fam)] 
    
  }
  # convert annotated groups, which can take binary values 0 and 1, to 1 and 2 so that it can be directly compared with cluster number
  #   note this means if the first cluster is set as 1, but mainly contains invasive/non-native profiles, this will result in -ve clustering
  #   this is ok as we only have 2 cluster groups, but note the modulus of the clustering should be taken into account, rather than whether it is +ve/-ve
  ind1 = which(annotate %in% 1) 
  ind0 = which(annotate %in% 0)
  annotate[ind1] = 1
  annotate[ind0] = 2
  
  # create vector of annotation colours 
  anno_cols = vector(length=length(annotate))
  anno_cols[ind1] = "#CC66CC"  # invasive/non-native = purple
  anno_cols[ind0] = "#66CC00"  # non-invasive/native = green
  
  hc_plot = as.dendrogram(hc_cut) # produce dendrogram object from hierachal clustering
  hc_plot = colour_branches(hc_plot,k=2) # colour branches according to true clustering into 2 groups
  labels_colors(hc_plot) = "white" # removes profile no. labels
  plot(hc_plot, main = bquote(italic(.(fam))), cex.main = 3) # plots the dendrogram 
  colored_bars(anno_cols,hc_plot,rowLabels = "") # adds annotation bar to the plot 
  
  # cluster validation
  
  expect_cluster = as.numeric(annotate)
  
  clust_dif = cluster.stats(dist(submat),expect_cluster,hc_cut$cluster) 
  clust_sim = clust_dif$corrected.rand 
  print(clust_sim)
  
  pa_clusters_f[[length(pa_clusters_f)+1]] = list(hc_cut$cluster,clust_sim)
 
}

names(pa_clusters_f) = family

########################################################################################################

############################################## Cluster by chemical similarity ################################################

###################### functions ##########################################

# calculate centroid distance from tanimoto's distance matrix of compounds, where each profile = group of compounds belonging to it

dist_profiles_chemsim = function(df,distmat,profilesub) {
  # df = dataframe of each profile with p/a data
  # distmat = distance matrix containing items that will be in the profile
  # profilesub = list of profile numbers for which the average distance is going to be calculated
  
  out_distmat = as.data.frame(matrix(NA,nrow=length(profilesub),ncol=length(profilesub))) # initialise the outputted distance matrix
  colnames(out_distmat) = profilesub
  rownames(out_distmat) = profilesub
  
  for(pro1 in profilesub) {
    
    for(pro2 in profilesub) {
      
      cell_val = out_distmat[which(rownames(out_distmat) == pro2),which(colnames(out_distmat) == pro1)]
      
      if(is.na(cell_val) == TRUE) {
        right_row1 = which(rownames(df) == pro1)
        comps_1 = colnames(df)[which(df[right_row1,] %in% 1)]
        right_row2 = which(rownames(df) == pro2)
        comps_2 = colnames(df)[which(df[right_row2,] %in% 1)]
        
        dist_bet_group = dist_between_centroids(distmat,comps_1,comps_2)
        
        out_distmat[which(rownames(out_distmat) == pro2),which(colnames(out_distmat) == pro1)] = dist_bet_group
        out_distmat[which(rownames(out_distmat) == pro1),which(colnames(out_distmat) == pro2)] = dist_bet_group
        
      }
    }
  }
  return(out_distmat)
}

#################################################################

chemsim_clusters_f = list()
list_ind_c_f = length(chemsim_clusters_f)+1

for(fam in family) {
  
  fam_distcomp = dist_profiles_chemsim(as.data.frame(nnat_fam[,compound_col_f],row.names=nnat_fam$Profile),
                                        fing_diss_matx_fam,nnat_fam$Profile[which(nnat_fam$Family %in% fam)])
  
  for(i in 1:nrow(fam_distcomp)) {
    
    
    fam_distcomp[i,i] = 0
  }
  
  fam_distcomp = as.dist(fam_distcomp)
  hc_cut = hcut(fam_distcomp,k=2,hc_func = "hclust",hc_method = "ward.D2",hc_metric = "euclidean")
  
  if(fam == family[1] | fam == family[2]) {
    
    annotate = nnat_fam$invasive_status[which(nnat_fam$Family %in% fam)]
    
  } else {
    annotate = nnat_fam$native_status[which(nnat_fam$Family %in% fam)]
    
      
  }
  ind1 = which(annotate %in% 1)
  ind0 = which(annotate %in% 0)
  annotate[ind1] = 1
  annotate[ind0] = 2
  
  anno_cols = vector(length=length(annotate))
  anno_cols[ind1] = "#CC66CC" 
  anno_cols[ind0] = "#66CC00"
  
  hc_plot = as.dendrogram(hc_cut)
  hc_plot = colour_branches(hc_plot,k=2)
  labels_colors(hc_plot) = "white"
  plot(hc_plot, main = bquote(italic(.(fam))), cex.main = 3)
  colored_bars(anno_cols,hc_plot,rowLabels = "")
  
  clusters = hc_cut$cluster
  clu_stat = cluster.stats(fam_distcomp,clustering = clusters, alt.clustering = annotate)
  print(clu_stat$corrected.rand)
  
  chemsim_clusters_f[[list_ind_c_f]] = list(clusters,clu_stat$corrected.rand)
  list_ind_c_f = length(chemsim_clusters_f) + 1
  
}

names(chemsim_clusters_f) = family

######################################################################################################################

# compare clusters produced (cluster no. = 2) to see if clustering agrees with each other, if not invasive/non-invasive clustering
# corrected.rand just compares cluster but cluster.stats requires a distance matrix - a dummy matrix can be supplies and will not affect results
#  if "compareonly" parameter is set to TRUE, meaning the distance matrix is not used.

clusters_compared = matrix(NA,ncol=3,nrow=0)

for(fam_i in 1:length(family)) {

  list_compare = list(pa_clusters_f[[fam_i]][1],chemsim_clusters_f[[fam_i]][1])
  names(list_compare) = c("p/a","chem_sim")
  
  c_stat = cluster.stats(dist(bin_cluster),clustering = unlist(list_compare[[1]]), 
                         alt.clustering = unlist(list_compare[[2]]), compareonly = TRUE)
  
  data_add = c(family[fam_i],paste(names(list_compare)[1],names(list_compare)[2]),c_stat$corrected.rand)
  clusters_compared = rbind(clusters_compared,data_add)
  
}

# total agreement w. middle 2 families, very weak agreement between the 1st and last (not sure this counts for much, could be biased by collection)
c_stat_n = cluster.stats(dist(bin_cluster),clustering = unlist(pa_clusters_f[[3]][[1]]), 
                       alt.clustering = unlist(chemsim_clusters_f[[7]][[1]]), compareonly = TRUE)

c_stat_n$corrected.rand
##############################################################################################################

################################################# Identifying Unique Compounds #############################################

# Getting the unique-to-population compounds
uniq_comp_f = NULL

for(fam in family) {
  
  if(fam == family[1] | fam == family[2]) {
    
    sub_data = nnat_fam[which(nnat_fam$Family %in% fam & nnat_fam$invasive_status %in% "1"),compound_col_f] # get data for invasive population  + family
    sub_data2 = nnat_fam[which(nnat_fam$Family %in% fam & nnat_fam$invasive_status %in% "0"),compound_col_f] # non-invasive
    
  } else {
    
    sub_data = nnat_fam[which(nnat_fam$Family %in% fam & nnat_fam$native_status %in% "1"),compound_col_f] # non-native
    sub_data2 = nnat_fam[which(nnat_fam$Family %in% fam & nnat_fam$native_status %in% "0"),compound_col_f] # native
    
  }
  
  chem_1 = colnames(sub_data)[which(colSums(sub_data) > 0)] # get compounds found in the invasive/non-native population
  chem_0 = colnames(sub_data2)[which(colSums(sub_data2) > 0)] # get compounds found in the non-invasive/native population
  
  uniq_1 = chem_1[-c(which(chem_1 %in% chem_0))] # get unique-to-invasive/non-native by removing compounds that are also found in the other population
  uniq_0 = chem_0[-c(which(chem_0 %in% chem_1))] # same as above but the other way round
  
  clu_comp1 = chemical_coordinates_fam$cluster[which(chemical_coordinates_fam$name %in% uniq_1)] # get the cluster data for the unique compounds
  clu_count1 = table(clu_comp1) # get the no. of counts/compounds found per cluster (given as named vector)
  clu_comp0 = chemical_coordinates_fam$cluster[which(chemical_coordinates_fam$name %in% uniq_0)] 
  clu_count0 = table(clu_comp0)
    
  data_add1 = c(fam,"1",unname(clu_count1),sum(unname(clu_count1))) # prepare data to add to the df
  data_add0 = c(fam,"0",unname(clu_count0),sum(unname(clu_count0)))
  
  uniq_comp_f = rbind(uniq_comp_f,data_add1,data_add0) # add to the df
  
}

uniq_comp_f = as.data.frame(uniq_comp_f,row.names = FALSE)
colnames(uniq_comp_f) = c("Family","Population",unique(chemical_coordinates_fam$cluster),"Total")
uniq_comp_f[,3:9] = lapply(uniq_comp_f[,3:9],as.numeric)

# statistical tests

## binomial tests total
tot_uniq = matrix(c(260,84,57,46,13,46,14,295,179,48,27,17,
                    17,7),ncol=7,nrow=2,byrow=TRUE)

chisq.test(tot_uniq[,2:ncol(tot_uniq)])
binom.test(tot_uniq[1,1],sum(tot_uniq[1,1],tot_uniq[2,1]),
           0.5,alternative="two.sided")

for(c in 2:ncol(tot_uniq)) {
 
   bin = binom.test(tot_uniq[1,c],sum(tot_uniq[1,c],tot_uniq[2,c]),
                   0.5,alternative="two.sided")
  print(bin$p.value)
  
}



## binomial tests between clusters and totals per family
uniq_test_res = NULL

for(i in seq(1,5,2)) {
  
  test_res = NULL
  
  chi_test = chisq.test(matrix(rbind(as.numeric(uniq_comp_f[i,3:8]),as.numeric(uniq_comp_f[i+1,3:8]))))
  bin_test_tot = binom.test(uniq_comp_f[i,9],sum(uniq_comp_f[i,9],uniq_comp_f[i+1,9]),0.5,alternative="two.sided")
  
  test_res = append(test_res,c(chi_test$p.value,bin_test_tot$p.value))
  
  for(j in 3:8) {
    
    bin_test = binom.test(uniq_comp_f[i,j],sum(uniq_comp_f[i,j],uniq_comp_f[i+1,j]),0.5,alternative="two.sided")
    
    test_res = append(test_res,bin_test$p.value)
    
  }
  
  uniq_test_res = rbind(uniq_test_res,test_res)
  
} 

colnames(uniq_test_res) = c("chi","binom",1:6)
rownames(uniq_test_res) = family

# absurdly low p-values - is this correct?

# solution - iterate over random samples of equal sample size where new_A = B, get average or range of average and p-vals

uniq_bootstrapped = function(fam) {
  
  if(fam == family[1] | fam == family[2]) {
    
    sub_data = nnat_fam[which(nnat_fam$Family %in% fam & nnat_fam$invasive_status %in% "1"),c(4,compound_col_f)] # get data for invasive population  + family
    sub_data2 = nnat_fam[which(nnat_fam$Family %in% fam & nnat_fam$invasive_status %in% "0"),c(4,compound_col_f)] # non-invasive
    
  } else {
    
    sub_data = nnat_fam[which(nnat_fam$Family %in% fam & nnat_fam$native_status %in% "1"),c(5,compound_col_f)] # non-native
    sub_data2 = nnat_fam[which(nnat_fam$Family %in% fam & nnat_fam$native_status %in% "0"),c(5,compound_col_f)] # native
    
  }
  
  data_list = list(sub_data,sub_data2)
  samp_size = c(nrow(sub_data),nrow(sub_data2))
  small_size = which.min(samp_size)
  data_used = data_list[[-small_size]]
  data_based = data_list[[small_size]]
  
  if(fam == family[1] | fam == family[2]) {
    
    invstat = c(unique(data_used$invasive_status),unique(data_based$invasive_status))
    
  } else {
    
    invstat = c(unique(data_used$native_status),unique(data_based$native_status))
    
  }
  
  data_used = data_used[,-1]
  data_based = data_based[,-1]
  
  multi_comp = NULL
  chi_bin_p = NULL
  
  for(itr in 1:10) {
    
    rand_samp = sample(1:samp_size[-small_size],samp_size[small_size])
    data_usedf = data_used[rand_samp,]
    
   if(invstat[1] == "1") {
     
     chem_1 = colnames(data_usedf)[which(colSums(data_usedf) > 0)] 
     chem_0 = colnames(data_based)[which(colSums(data_based) > 0)] 
     
   
    } else if(invstat[1] == "0") {
     
       chem_1 = colnames(data_based)[which(colSums(data_based) > 0)] 
       chem_0 = colnames(data_usedf)[which(colSums(data_usedf) > 0)] 
      
   }
    
    uniq_1 = chem_1[-c(which(chem_1 %in% chem_0))] 
    uniq_0 = chem_0[-c(which(chem_0 %in% chem_1))] 
    
    clu_comp1 = chemical_coordinates_fam$cluster[which(chemical_coordinates_fam$name %in% uniq_1)] 
    clu_count1 = table(clu_comp1) 
    clu_comp0 = chemical_coordinates_fam$cluster[which(chemical_coordinates_fam$name %in% uniq_0)] 
    clu_count0 = table(clu_comp0)
    
    data_add1 = c(sum(unname(clu_count1)),unname(clu_count1)) 
    data_add0 = c(sum(unname(clu_count0)),unname(clu_count0))
    multi_comp = rbind(multi_comp,data_add1,data_add0)

    chi_test = chisq.test(clu_count1,clu_count0)
    bin_tot = binom.test(sum(unname(clu_count1)),sum(sum(unname(clu_count1)),sum(unname(clu_count0))),
                         0.5,alternative="two.sided")
      
    test_res = NULL
    
    for(i in 1:6) {
        
        bin_test = binom.test(unname(clu_count1[i]),sum(unname(clu_count1[i]),unname(clu_count0[i])),
                              0.5,alternative="two.sided")
        test_res = append(test_res,bin_test$p.value)
        
      }
      
    data_add_p = c(chi_test$p.value,bin_tot$p.value,test_res)
    chi_bin_p = rbind(chi_bin_p,data_add_p)  
      
  } 
    
  analysis_rep = list()
   
  for(j in 1:7) {
    
    av_size1 = median(multi_comp[seq(1,19,2),j],na.rm=TRUE)
    av_size0 = median(multi_comp[seq(2,20,2),j],na.rm=TRUE)
    range_1 = c(min(multi_comp[seq(1,19,2),j],na.rm=TRUE),max(multi_comp[seq(1,9,2),j],na.rm=TRUE))
    range_0 = c(min(multi_comp[seq(2,20,2),j],na.rm=TRUE),max(multi_comp[seq(2,10,2),j],na.rm=TRUE))
    
    mean_ps = median(chi_bin_p[,j+1])
    minmax_ps = c(min(chi_bin_p[,j+1]),max(chi_bin_p[,j+1]))
    
    
    
    if(j == 1) {
      
      mean_chi = median(chi_bin_p[,j])
      minmax_chi = c(min(chi_bin_p[,j]),max(chi_bin_p[,j]))
      
      list_p = list(mean_chi,minmax_chi,mean_ps,minmax_ps)
      
    } else {
      
      list_p = list(mean_ps,minmax_ps)
    
    }
    
    analysis_rep[[length(analysis_rep)+1]] = c(list(av_size1,range_1),list(av_size0,range_0),list_p)
  
  } 

  return(analysis_rep)
  
}
  
v_ana = uniq_bootstrapped(family[1])
m_ana = uniq_bootstrapped(family[2])
a_ana = uniq_bootstrapped(family[3]) 
# = uniq_bootstrapped(family[4])

#multi_test_uniq = list(v_ana,m_ana,c_ana)
#names(multi_test_uniq) = family[-3]
multi_test_uniq = list(v_ana,m_ana,a_ana)
names(multi_test_uniq) = family
  
################### see where clusters align with our original clustering ###########################

chem_coords_og = readRDS("chemical_coordinates.rds")

smi_clu_new = chemical_coordinates_fam$cluster[which(chemical_coordinates_fam$name %in% chem_coords_og$name)]
comp_clust_compare = as.data.frame(cbind(chem_coords_og$cluster,smi_clu_new))
plot(comp_clust_compare[,1],comp_clust_compare[,2])

clu_counting = NULL

for(i in 1:3) {
  
  clu_count = NULL
  
  for(j in 1:6) {
    
    count_clu = length(which(comp_clust_compare[,1] %in% i & comp_clust_compare[,2] %in% j))
    clu_count = append(clu_count,count_clu)
    
  }
  
  clu_counting = rbind(clu_counting,clu_count)
  
}

clu_counting = as.data.frame(clu_counting)
rownames(clu_counting) = c(1:3)
colnames(clu_counting) = c(1:6)
clu_counting[,1:ncol(clu_counting)] = lapply(clu_counting[,1:ncol(clu_counting)],as.numeric)
clu_counting = as.matrix(clu_counting)
clu_counting[which(clu_counting %in% 0)] = NA
heatmap(clu_counting,Rowv=NA,Colv=NA,col=colorRampPalette(brewer.pal(7, "YlOrRd"))(25),cexCol=2,cexRow=2)

################################################ Analysing cluster representation compared to abiotic and biotic factors ##################

# invasives only 
# clean data first - ideally ones with complete datasets 
# get clu count per profile 

factors_cplt = readRDS("factors_cplt.rds")

inv_data = which(nnat_fam$invasive_status[1:51] %in% 1)
inv_data = c(inv_data,(which(nnat_fam$native_status[52:nrow(nnat_fam)] %in% 1)+51))
fact_inv = factors_cplt[inv_data,] # only invasive/non-native profiles
fact_inv = fact_inv[which(is.na(fact_inv$Tree_Loss_.)==FALSE),] # clear NAs, data with NAs are not filled across any factor columns
# I was wrong 
rownames(fact_inv) = 1:nrow(fact_inv)
fact_inv = fact_inv[-c(23,35),] # remove data with NAs
rownames(fact_inv) = 1:nrow(fact_inv)
fact_inv[24,8] = 1 # correct false NA data

# get count data 
count_dat = NULL

for(pro in fact_inv$Profile) {
  
  sub_data = nnat_fam[which(nnat_fam$Profile %in% pro),compound_col_f]
  comps = colnames(sub_data)[which(colSums(sub_data) > 0)]
  clu_comp = chemical_coordinates_fam$cluster[which(chemical_coordinates_fam$name %in% comps)] 
  clu_count = table(clu_comp) 
  data_add = c(pro,unname(clu_count)) # prepare data to add to the df
  count_dat = rbind(count_dat,data_add)
  
}

count_dat = as.data.frame(count_dat)
rownames(count_dat) = 1:nrow(count_dat)
colnames(count_dat) = c("Profile",1:6)

model_data = as.data.frame(cbind(count_dat,fact_inv[,-1]))
rownames(model_data) = 1:nrow(model_data)
model_data[19,18] = 1355
model_data[20,18] = 5
model_data[,2:7] = lapply(model_data[,2:7] ,as.numeric)
model_data[,2:7] = lapply(model_data[,2:7] ,as.integer)
model_data[,c(11,12,15:ncol(model_data))] = lapply(model_data[,c(11,12,15:ncol(model_data))] ,as.numeric)
model_data[,c(13,14)] = lapply(model_data[,c(13,14)], as.character)
model_data[,8] = as.character(model_data[,8])
colnames(model_data)[c(11,14,15,16)] =c("Tree_Loss","Anthropological_Activity",
                                        "Temperature_High","Temperature_Low")
colnames(model_data)[2:7] = c("Cluster_1","Cluster_2","Cluster_3","Cluster_4","Cluster_5","Cluster_6")

glm_generator = function(data,indep_var,dep_var) {
  

  glm_res = list()
  
  indep_var_string = indep_var[1]
  
  for(var in indep_var[2:length(indep_var)]) {
    
    indep_var_string = paste0(indep_var_string,paste0("+",var))
  
    }
  
  for(vard in dep_var) {
    
    formulap = as.formula(paste0(vard,"~",indep_var_string))
    p_glm = glm(formulap,family=poisson(),data=data)
    glm_res[[length(glm_res)+1]] = p_glm
    
  }
  
  return(glm_res)
  
}

glm_tree = glm_generator(model_data,colnames(model_data)[c(8,11,13:18)],colnames(model_data)[2:7])
glm_fire = glm_generator(model_data,colnames(model_data)[c(8,12:18)],colnames(model_data)[2:7])
glm_sr = glm_generator(model_data,colnames(model_data)[c(8,13:19)],colnames(model_data)[2:7])

glm1 = plot_models(glm_tree[[1]],glm_fire[[1]],glm_sr[[1]],transform=NULL,
            axis.title=c("Coefficient","Variable"),p.shape=TRUE,title="Cluster 1") + 
  theme_sjplot2(base_size = 20, base_family = "") + 
  scale_color_discrete(name = "Model",labels = c("With Plant SR","With Fire Risk","With Tree Loss"))
glm_legend = get_legend(glm1)

glm1 = plot_models(glm_tree[[1]],glm_fire[[1]],glm_sr[[1]],transform=NULL,
                   axis.title=c("Coefficient","Variable"),p.shape=TRUE,title="Cluster 1",
                   show.legend=FALSE) + 
  theme_sjplot2(base_size = 17, base_family = "") + 
  scale_color_discrete(name = "Model",labels = c("With Plant SR","With Fire Risk","With Tree Loss"))

glm2 = plot_models(glm_tree[[2]],glm_fire[[2]],glm_sr[[2]],transform=NULL,
            axis.title=c("Coefficient","Variable"),p.shape=TRUE,title="Cluster 2",
            show.legend=FALSE) + 
  theme_sjplot2(base_size = 17, base_family = "") + 
  scale_color_discrete(name = "Model",labels = c("With Plant SR","With Fire Risk","With Tree Loss"))

glm3 = plot_models(glm_tree[[3]],glm_fire[[3]],glm_sr[[3]],transform=NULL,
            axis.title=c("Coefficient","Variable"),p.shape=TRUE,title="Cluster 3",
            show.legend=FALSE) + 
  theme_sjplot2(base_size = 17, base_family = "") + 
  scale_color_discrete(name = "Model",labels = c("With Plant SR","With Fire Risk","With Tree Loss"))

glm4 = plot_models(glm_tree[[4]],glm_fire[[4]],glm_sr[[4]],transform=NULL,
            axis.title=c("Coefficient","Variable"),p.shape=TRUE,title="Cluster 4",
            show.legend=FALSE) + 
  theme_sjplot2(base_size = 17, base_family = "") + 
  scale_color_discrete(name = "Model",labels = c("With Plant SR","With Fire Risk","With Tree Loss"))

glm5 = plot_models(glm_tree[[5]],glm_fire[[5]],glm_sr[[5]],transform=NULL,
            axis.title=c("Coefficient","Variable"),p.shape=TRUE,title="Cluster 5",
            show.legend=FALSE) + 
  theme_sjplot2(base_size = 17, base_family = "") + 
  scale_color_discrete(name = "Model",labels = c("With Plant SR","With Fire Risk","With Tree Loss"))

grid.arrange(arrangeGrob(glm1,glm2,glm3,glm4,glm5,glm_legend,ncol=3,nrow=2))


