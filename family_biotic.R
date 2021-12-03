library(dplyr)
library(rgbif)
library(dplyr)
library(countrycode)
library(geosphere)

#########################################################
plantae_key = name_backbone("Plantae")
plantae_key = plantae_key %>%
  dplyr::select(kingdomKey)
plantae_key = plantae_key[[1]]
taxonp = plantae_key

animalia_key = name_backbone("Animalia")
animalia_key = animalia_key %>%
  dplyr::select(kingdomKey)
animalia_key = animalia_key[[1]]
taxona = animalia_key

overall_sample = c(1970:2020)

iso3cs = unique(nnat_fam$iso3c[68:nrow(nnat_fam)])
iso3cs = iso3cs[!is.na(iso3cs)]

# convert iso3c to iso2c
iso2cs = NULL

for(iso in iso3cs) {
  
  iso2c = countrycode(iso,origin="iso3c",destination="iso2c")
  iso2cs = append(iso2cs,iso2c)
  
}

gbif_data = list()

for(iso in iso2cs) {
  
  plant = occ_data(country=iso,year=overall_sample,hasCoordinate = TRUE,kingdomKey=taxonp)
  animal = occ_data(country=iso,year=overall_sample,hasCoordinate = TRUE,kingdomKey=taxona)
  
  format_country_p = NULL
  
  for(i in 1:length(plant)) {
    
    country_data2 = plant[i]
    year = names(plant)[[i]]
    country_data2 = country_data2[[year]]$data
    
    if(is.null(country_data2) == FALSE & is.null(country_data2$species) == FALSE) {
      
      country_data2 = as.data.frame(country_data2) # creates a dataframe of downloaded data
      country_subset = subset(country_data2, select=c(species,year))
      format_country_p = rbind(format_country_p,country_subset)
      
    }
    
  }
  
  format_country_a = NULL
  
  for(i in 1:length(animal)) {
    
    country_data2 = animal[i]
    year = names(animal)[[i]]
    country_data2 = country_data2[[year]]$data
    
    if(is.null(country_data2) == FALSE & is.null(country_data2$species) == FALSE) {
      
      country_data2 = as.data.frame(country_data2) # creates a dataframe of downloaded data
      country_subset = subset(country_data2, select=c(species,year))
      format_country_p = rbind(format_country_a,country_subset)
      
    }
    
  }
  
  list_add = list(format_country_p,format_country_a)
  names(list_add) = c("plant","animal")
  gbif_data[[length(gbif_data)+1]] = list_add
  
}

names(gbif_data) = iso2cs # worked for plants not animals, do plants stuff

#richness_table = as.data.frame(matrix(NA,ncol=4,nrow=length(68:nrow(nnat_fam))))
richness_table = as.data.frame(matrix(NA,ncol=3,nrow=length(68:nrow(nnat_fam))))

collection_date = species_data[which(species_data$family %in% family[-3]),c(1,5)] # get species and sample year
collection_date = collection_date[-c(which(collection_date$Plant.Name %in% "Lantana camara  ")),2] # get sample years 
# tidy data, not totally numeric
unique(collection_date)
collection_date[which(collection_date %in% "Jan-dec 2014")] = 2014
collection_date = as.numeric(collection_date)
# check data matches
length(68:nrow(nnat_fam)) == length(collection_date) # does match, good to go
                               
for(row in 68:nrow(nnat_fam)) {
  
  sample_y = collection_date[row]
  
  if(is.na(sample_y) == TRUE) {
    
    sampling = c(1980:2020)
    
  } else {
    
    sampling_low = sample_y-10
    
    if((sample_y+10) > 2020) {
      
      sampling_high = 2020
      
    } else {
      
      sampling_high = sample_y+10
      
    }
    
    sampling = c(sampling_low:sampling_high)
    
  }
  
  country = nnat_fam$iso3c[row]
  country = iso2cs[which(iso3cs %in% country)]
  
  if(is.na(country) == FALSE) {
    
    data_p = gbif_data[[which(names(gbif_data) %in% country)]][[1]]
    #  data_a = gbif_data[[which(names(gbif_data) %in% country)]][[2]]
    
    species_p = unique(data_p$species[which(data_p$year %in% sampling)])
    sr_p = length(species_p)
    #  species_a = unique(data_a$species[which(data_p$year %in% sampling)])
    #  sr_a = length(species_a)
    
    #  richness_table[(row-67),] = c(nnat_fam$Profile[row],sample_y,sr_p,sr_a)
    richness_table[(row-67),] = c(nnat_fam$Profile[row],sample_y,sr_p)
    
  }
}
                               
names(richness_table) = c("Profile","Collection_Date","Plantae_SR")

############################################ Abiotic data ###############################################

locations = species_data[which(species_data$family %in% family[-3]),c(1,3)] # get species and location
locations = locations[-c(which(locations$Plant.Name %in% "Lantana camara  ")),2] # get locations
unique(locations)
unique(matrix_data_feat$Location)
factors_to_collect = as.data.frame(unique(locations))
write.csv(factors_to_collect,"~/factors_gita.csv")
factors_to_collect = read.csv("~/factors_gita.csv")

# compiling factors 
factors_cplt = as.data.frame(matrix(NA,nrow=length(locations),ncol=length(factors_to_collect))) # expand factors_to_collect

for(row in 1:nrow(factors_to_collect)) {
  
  loc = factors_to_collect$Location[row]
  loc_ind = which(locations %in% loc)
  factors_cplt[loc_ind,] = factors_to_collect[row,]
  
}

factors_cplt = as.data.frame(cbind(collection_date,factors_cplt,richness_table$Plantae_SR))
factors_cplt = factors_cplt[,-7]
colnames(factors_cplt) = c("Collection_Date",colnames(factors_to_collect)[-6],"Plant_SR")

mat_add = subset(matrix_data_feat,select=c(Sample_year,Location,Deforestation..total.tree.cover.loss.over.18yrs.,
                                                      Fire.Risk..VIIRS.alerts.2020.,Flood_risk..1.yes..0.no.,Agronomic_act_close_prox..1.yes.0.no.,
                                                      Average.high.Temperature..annual.or.month...celcius.,Average.low.Temperature..celcius.,
                                                      Preciptation..mm...annual,Elevation,plantae_20ySR_country))
colnames(mat_add) = colnames(factors_cplt)
factors_cplt = as.data.frame(rbind(mat_add,factors_cplt))
factors_cplt = as.data.frame(cbind(nnat_fam$Profile,nnat_fam$Family,factors_cplt))
colnames(factors_cplt)[1:2] = c("Profile","Family")

saveRDS(factors_cplt,"factors_cplt.rds")


