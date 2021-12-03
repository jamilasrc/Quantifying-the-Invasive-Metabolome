library(countrycode)
library(stringr)
library(webchem)
library(rJava)
library(rcdk)


######## Combining data matrices ############################
# Aims: move gitas data into the p/a matrix
#  From gita's matrix extract country and sampling date as well
#  From GISD for the new species establish if a species is invasive and invasivity by distribution

# GISD Data
taxon_gisd = read.delim("~/dwca-gisd-v1.4/taxon.txt",header=TRUE) # taxon data
occurance_gisd = read.delim("~/dwca-gisd-v1.4/distribution.txt",header=TRUE) # distribution - country, native/alien, invasive/non-invasive/unspecified

# Gita's Data
species_data = read.delim("~/gita_databyspecies.tsv",header=TRUE) 
family_data = read.delim("~/gita_databyfamilies.tsv",header=TRUE) # this is just a list of compounds with no other info. Not sure how helpful this is

# Issue 1 - GISD invasivity status given as a hyperlink - convert using str_remove
# Issue 2 - Gita's location data is not given according to country, so country will have to be extracted. 

# Issue 2 
## get list of countries of the world
## go through dataframe location column, see if one of the words matches, replact with matched column

gitaspec_newctry = as.data.frame(matrix(NA,ncol=1,nrow=nrow(species_data))) # initiate new dataframe

all_ctry = read.csv("~/list-countries-world.csv") # downloaded list of countries of the world
all_ctry = c(all_ctry[,1]) # remove excess columns

old_loc = unique(species_data$Location)

for(ctry in all_ctry) {
  
  ctry_extr = grep(ctry,old_loc,ignore.case=TRUE,value=TRUE)
  loc_ind = which(species_data$Location %in% ctry_extr)
  gitaspec_newctry[loc_ind,1] = ctry
  
}

# some location data did not include countries - add this manually
manual_loc = species_data$Location[which(is.na(gitaspec_newctry[,1]))]
unique(manual_loc)
# Table of values
## "Dakal-e-Derazno, Golestan" = "Iran"
## "Karaj, Tehran" = "Iran"
## "Orlando, Ethipio" = unknown (remove)
## "Jalisco" = "Mexico"
## "Adriatic coast" = "Italy"
## "Jadran coast" = "Croatia"
## "Mugla, Marmaris" = "Turkey"
## "Reunion Island" = technically France, but is a pacific island. Keep (might have a country code)
## NA = REMOVE DATA

gitaspec_newctry[which(species_data$Location %in% "Dakal-e-Derazno, Golestan"),1] = "Iran"
gitaspec_newctry[which(species_data$Location %in% "Karaj, Tehran"),1] = "Iran"
gitaspec_newctry[which(species_data$Location %in% "Orlando, Ethipio"),1] = NA
gitaspec_newctry[which(species_data$Location %in% "Jalisco"),1] = "Mexico"
gitaspec_newctry[which(species_data$Location %in% "Adriatic coast"),1] = "Italy"
gitaspec_newctry[which(species_data$Location %in% "Jadran coast"),1] = "Croatia"
gitaspec_newctry[which(species_data$Location %in% "Mugla, Marmaris"),1] = "Turkey"
gitaspec_newctry[which(species_data$Location %in% "Reunion Island"),1] = "Reunion Island"

gitaspec_newctry = as.data.frame(cbind(gitaspec_newctry,rep(NA,length=nrow(gitaspec_newctry))))
colnames(gitaspec_newctry) = c("Country_Name","Country_Code")

# next step is to convert countries to country code
# call library(countrycode) at the top
# note - gisd is 3c vs my original 2c matrices 
uniq_cntry = unique(gitaspec_newctry[,1])

for(cntry in uniq_cntry) {
  
  code3c = countrycode(cntry,origin="country.name",destination="iso3c") # gisd country code is iso3c
  cntry_ind = which(gitaspec_newctry$Country_Name %in% cntry)
  gitaspec_newctry$Country_Code[cntry_ind] = code3c
  
}

# create matrices 
# 1 matrix for all data, the other for invasives only.

# identify which species from the matrix are invasive
# firstly because species names are not given we are going to have to go by genus and compare manually
gita_genus = NULL

for(name in species_data$Plant.Name) {
  
  genus = strsplit(name," ")[[1]][1]
  gita_genus = append(gita_genus,genus)
  
}

invasives = unique(gita_genus[which(gita_genus %in% taxon_gisd$genus)]) # suggests 3 are invasive
# inspect species
unique(species_data$Plant.Name[which(gita_genus %in% invasives)]) # "Lantana camara " "Ocimum basilicum  " "Piper nigrum  " 
unique(taxon_gisd$scientificName[which(taxon_gisd$genus %in% invasives)]) # "Lantana camara L.1753" "Piper aduncum L." "Ocimum gratissimum  L."
# only L. camara 

# will need to check L. camara sources before adding to the invasive matrix
# not in the dataframe - ask gita for sources 

which(species_data$Plant.Name %in% "Melaleuca quinquenervia")
which(species_data$Plant.Name %in% "Psidium cattleainum")
which(species_data$Plant.Name %in% "Ageratum conyzoides")
# no data on any other my other invasive species

# Start creating matrices without L. camara for all data 

matrix_data_feat = readRDS("matrix_format_features.rds") # reload original p/a matrix by me 

species_add1st = species_data[-c(which(species_data$Plant.Name %in% "Lantana camara  ")),]

# inspect if plant names need normalising 
unique(species_add1st$Plant.Name) # all good 

# BEFORE CONTINUING, CHECK FAMILIES OF MY DATA + GITA'S DATA TO SEE IF HER DATA CAN BE USED WITH MINE
unique(species_data$family)
# L. camara = Verbenaceae
# M. quinquenervia = Myrtaceae
# P. cattleainum = Myrtaceae
# A. conyzoides = Asteraceae
## all have matches except A. conyzoides
unique(species_data$Plant.Name[which(species_data$family %in% "Verbenaceae")])
## only L. camara from this family, this dataset cannot be bulked out 
# CONCLUSION: A COMPARISON OF THE MYRTACEAE FAMILY CAN BE DONE

# inspect compounds
## load all smiles 
SMILES = colnames(matrix_data_feat[3:381])
## get common names from other databases
comp_feat_csv_cleaned = read.csv("~/University Work/3rd Year/Research Project Paper Publication/data_tableform.csv")

smiles_comp = as.data.frame(cbind(family_data,rep(NA,length=nrow(family_data))))
colnames(smiles_comp) = c("common_name","SMILES")

for(comp in smiles_comp$common_name) {
  
  comp_match = which(tolower(comp_feat_csv_cleaned$Compound) %in% tolower(comp))[1]
  SMILE = comp_feat_csv_cleaned$SMILES[comp_match]
  smiles_comp$SMILES[which(smiles_comp$common_name %in% comp)] = SMILE
  
}

length(which(!is.na(smiles_comp$SMILES))) # 345 out of 379, seems like any shared compounds haven't been missed out
# put the rest through the functions 
# create a cleaner for loop now the internet is better
# select 1st SMILE, this was generally the unique 1 

for(comp_i in which(is.na(smiles_comp$SMILES))) { # index for all NA values for SMILES
  
   # remove when completed, just done to continue loop from where it stopped
    
    comp = smiles_comp$common_name[comp_i] # get corresponding common name
    smi_opt = cir_query(identifier = comp, representation = "smiles") # convert to smiles
    smi_opt = smi_opt[[1]][1] # get 1st SMILE (tends to be the unique 1)
    
    # Can't remember if all SMILES are given as canonical or isomeric
    # Add rcdk functions to return canonical smiles
    if(is.na(smi_opt) == FALSE) {
      smi_final = parse.smiles(smi_opt)[[1]] # convert any smiles to a new data type
      smi_final = get.smiles(smi_final,smiles.flavors("Canonical")) # convert data type to canonical smiles
      smiles_comp$SMILES[comp_i] = smi_final # add to df
      
    }
}

# 1/3 of the data had no smiles found - lots of cleaning to do

# option 1 - get CID from pubchem then convert to canonical smiles - seems more forgiving of spelling mistakes
for(comp_i in which(is.na(smiles_comp$SMILES))) { # index for all NA values for SMILES
  
  comp = smiles_comp$common_name[comp_i] # get corresponding common name
  
  cid_pc = get_cid(comp) # get pubchem cid tibble
  cid_pc = cid_pc[[2]] # get cid string
  
  if(is.na(cid_pc) == FALSE) {
    
    smiles_pc = pc_prop(cid_pc,properties="CanonicalSMILES") # get canonical SMILES vector from the pubchem database for the specific cid
    smiles_pc = smiles_pc[[2]] # get SMILES string
    
    smiles_comp$SMILES[which(smiles_comp$common_name %in% comp)] = smiles_pc # add to dataframe
    
  }
}

# successful but the other errors are non-trival aka specific to how a compound was entered to pubchem database
# process this manually
na_smiles = smiles_comp[which(is.na(smiles_comp$SMILES)),]

write.csv(na_smiles[1:150,],"~/comp_for_jb.csv",row.names=FALSE) # compounds for jonathan to go through
write.csv(na_smiles[151:nrow(na_smiles),],"~/no_smiles.csv",row.names=FALSE) # compounds for me to edit

# Note - for some of the data I am not sure if they are canonical or isomerical SMILES - run through parse.smiles + get.smiles

# Got the manually collected smiles, clean them
upd_smiles1 = read.csv("~/no_smiles.csv")

for(comp in upd_smiles1$SMILES) {
  
  if(is.na(comp) == FALSE) {
    
    smi_in = parse.smiles(upd_smiles1$SMILES[comp_i])[[1]]
    smi_out = get.smiles(smi_in,smiles.flavors("Canonical"))
    
    upd_smiles1$SMILES[comp_i] = smi_out
    
  }
  
}

# SMILES from Jonathan
 upd_smiles2 = read.csv("~/comp_for_jb.csv")

for(comp in upd_smiles2$SMILES) {
  
  if(is.na(comp) == FALSE) {
    
    smi_in = parse.smiles(upd_smiles2$SMILES[comp_i])[[1]]
    smi_out = get.smiles(smi_in,smiles.flavors("Canonical"))
    
    upd_smiles2$SMILES[comp_i] = smi_out
    
  }
  
}

# Compiling Dataframes

 ## first get all SMILES
SMILES_complete = unique(c(SMILES,upd_smiles1$SMILE,upd_smiles2$SMILE,smiles_comp$SMILES)) # total 2866 cleaned compounds
SMILES_complete = SMILES_complete[which(is.na(SMILES_complete) == FALSE)] # remove NA value

## Start the dataframes
compiled_pa_matrix = as.data.frame(matrix(NA,nrow=sum(nrow(matrix_data_feat),nrow(species_add1st)),
                                          ncol=(length(SMILES_complete)+2))) # initialise completed data frame
colnames(compiled_pa_matrix) = c("Profile","Species",SMILES_complete)

## add my data to the matrix
sub_matrix_data_feat = matrix_data_feat[,3:381]

for(row in 1:nrow(matrix_data_feat)) {
  
  present = colnames(sub_matrix_data_feat[which(sub_matrix_data_feat[row,] %in% 1)]) # get SMILES present in the profile
  
  compiled_pa_matrix[row,which(colnames(compiled_pa_matrix) %in% present)] = 1 # add in present compounds
  compiled_pa_matrix[row,1:2] = matrix_data_feat[row,1:2] # add in profile no. and species
  compiled_pa_matrix[row,which(is.na(compiled_pa_matrix[row,]))] = 0 # all NA-filled columns will now be absent compounds, fill in with 0s
  
}

## add gita's data
compiled_gita_smiles = as.data.frame(rbind(smiles_comp,upd_smiles1,upd_smiles2))
colnames(compiled_gita_smiles) = c("Name","SMILE")

colnames(species_add1st)[8:ncol(species_add1st)] = family_data[,1] # compound columns should now exactly match the names of the compounds
sub_specdat = species_add1st[,8:ncol(species_add1st)]

for(row in 1:nrow(sub_specdat)) {
  
  pres_comp = colnames(sub_specdat[which(sub_specdat[row,] %in% 1)]) # get common names of present compounds
  
  smi_pres = NULL # initialise vector of present SMILES
  
  for(p in pres_comp) {
    
    smi_opti = compiled_gita_smiles$SMILE[which(compiled_gita_smiles$Name %in% p)] # get SMILES that match with the name (should be 1 SMILE + NA in some circumstances)
    smi_opti = smi_opti[which(is.na(smi_opti) == FALSE)] # make sure the SMILE and not the NA value is obtained
    
    if((length(smi_opti) > 0) == TRUE) { # if there is a SMILE left in the vector (so removes any compound where the name could not be converted)
      
      smi_pres = append(smi_pres,smi_opti) # add to list of present SMILES
      
    }
  }
  
  new_row = nrow(matrix_data_feat) + row
  compiled_pa_matrix[new_row,which(colnames(compiled_pa_matrix) %in% smi_pres)] = 1 # add in present compounds
  pro_no = matrix_data_feat[nrow(matrix_data_feat),1] + row # get new profile no. 
  compiled_pa_matrix[new_row,1:2] = c(pro_no,species_add1st[row,1]) # add in profile no. and species
  compiled_pa_matrix[new_row,which(is.na(compiled_pa_matrix[new_row,]))] = 0 # all NA-filled columns will now be absent compounds, fill in with 0s
  
}

empty_cols = NULL # since colSums check isn't working

for(col in 3:ncol(compiled_pa_matrix)) {
  
  summed = sum(compiled_pa_matrix[,col])
  
  if(summed == 0) {
    
    empty_cols = append(empty_cols,col)
    
  }
}
# turns out colsums was working and the data was just not cleaned fully to remove absent compounds 

compiled_pa_matrix = compiled_pa_matrix[,-empty_cols] # remove any columns containing 0s
# loads were missing 
saveRDS(compiled_pa_matrix,"compiled_pa_matrix.rds")
compiled_pa_matrix = readRDS("compiled_pa_matrix.rds")
# next, add needed columns - family data, invasive status and native status

compiled_pa_matrix = cbind(compiled_pa_matrix,rep(NA,length=nrow(compiled_pa_matrix)),rep(NA,length=nrow(compiled_pa_matrix)),
                            rep(NA,length=nrow(compiled_pa_matrix)))

colnames(compiled_pa_matrix)[(ncol(compiled_pa_matrix)-2):ncol(compiled_pa_matrix)] = c("Family","invasive_status","native_status")
# note invasive_status: 1 = invasive, 0 = non_invasive; native_status: 1 = non-native, 0 = native


# first country data will need to be added
## correct my matrix's country data
for(cntry in unique(matrix_data_feat$iso2c)) {
  
  code3c = countrycode(cntry,origin="iso2c",destination="iso3c") # gisd country code is iso3c
  cntry_ind = which(matrix_data_feat$iso2c %in% cntry)
  matrix_data_feat$iso2c[cntry_ind] = code3c
  
}

compiled_pa_matrix = cbind(compiled_pa_matrix,rep(NA,length=nrow(compiled_pa_matrix)))
colnames(compiled_pa_matrix)[ncol(compiled_pa_matrix)] = c("iso3c")
compiled_pa_matrix$iso3c[1:nrow(matrix_data_feat)] = matrix_data_feat$iso2c # add my data
gitaspec_newctry_add = gitaspec_newctry[-which(species_data$Plant.Name %in% "Lantana camara  "),] # get country codes for Gita's data
compiled_pa_matrix$iso3c[(nrow(matrix_data_feat)+1):nrow(compiled_pa_matrix)] = gitaspec_newctry_add[,2] # add gita's data
compiled_pa_matrix$iso3c = as.character(unlist(compiled_pa_matrix$iso3c)) # convert countries from list to characters

## add my data for invasive and native status
compiled_pa_matrix$native_status[1:nrow(matrix_data_feat)] = matrix_data_feat$native_status
compiled_pa_matrix$invasive_status[1:nrow(matrix_data_feat)] = matrix_data_feat$invasive_status

# note: some of the new species are non-native but naturalised - this should maybe be looked at. 
#  I am not entirely sure how to get this info - might have to do it manually 
# also some are reported as invasive according to CABI
# From quick web search, where invasiveness isn't reported, set as 0

unique(compiled_pa_matrix$Species) # go through manually

# Invasive species: Rosmarinus officinalis, Achillea millefolium, Tagetes erecta
unique(species_data$family[which(species_data$Plant.Name %in% c("Achillea millefolium  ","Rosmarinus officinalis  ","Tagetes erecta  "))])
# Families with invasive species: Verbenaceae, Myrtaceae, Asteraceae, Compositae, Lamiaceae
unique(compiled_pa_matrix$Species[which(compiled_pa_matrix$Family %in% c("Verbenaceae","Myrtaceae","Asteraceae","Compositae","Lamiaceae"))])
# lots of species can be included in the family analysis

# go through country data to assign invasive status
unique(compiled_pa_matrix$iso3c[which(compiled_pa_matrix$Species %in% "Rosmarinus officinalis  ")])   # R. officianlis countries
# invasive in: 
# non-native in: BRA, CHN
unique(compiled_pa_matrix$iso3c[which(compiled_pa_matrix$Species %in% "Achillea millefolium  ")])
# invasive in:
# non-native in: CUB
unique(compiled_pa_matrix$iso3c[which(compiled_pa_matrix$Species %in% "Tagetes erecta  ")])
# invasive in:
# non-native in: IND, ITA, BGR

# add data
compiled_pa_matrix$invasive_status[which(is.na(compiled_pa_matrix$invasive_status))] = 0 # none are reported as invasive in our countries
# add non-native data
compiled_pa_matrix$native_status[which(compiled_pa_matrix$Species %in% "Rosmarinus officinalis  " &
                                         compiled_pa_matrix$iso3c %in% c("BRA","CHN"))] = 1 
compiled_pa_matrix$native_status[which(compiled_pa_matrix$Species %in% "Achillea millefolium  " &
                                         compiled_pa_matrix$iso3c %in% "CUB")] = 1 
compiled_pa_matrix$native_status[which(compiled_pa_matrix$Species %in% "Tagetes erecta  " &
                                         compiled_pa_matrix$iso3c %in% c("IND","ITA","BGR"))] = 1 
compiled_pa_matrix$native_status[which(is.na(compiled_pa_matrix$native_status))] = 0 # add in native sites


# add family data
compiled_pa_matrix$Family[which(compiled_pa_matrix$Species %in% "Lantana camara")] = "Verbenaceae"
compiled_pa_matrix$Family[which(compiled_pa_matrix$Species %in% c("Melaleuca quinquenervia","Psidium cattleianum"))] = "Myrtaceae"
compiled_pa_matrix$Family[which(compiled_pa_matrix$Species %in% "Ageratum conyzoides")] = "Asteraceae"
compiled_pa_matrix$Family[(nrow(matrix_data_feat)+1):nrow(compiled_pa_matrix)] = species_add1st$family

# tidy species names to remove spaces from the end of names
compiled_pa_matrix$Species[(nrow(matrix_data_feat)+1):nrow(compiled_pa_matrix)] = substr(compiled_pa_matrix$Species[(nrow(matrix_data_feat)+1):nrow(compiled_pa_matrix)],
                                                                                        1,nchar(compiled_pa_matrix$Species[(nrow(matrix_data_feat)+1):nrow(compiled_pa_matrix)])-2)

# final save
saveRDS(compiled_pa_matrix,"compiled_pa_matrix.rds")

# Make similarity and dissimilarity matrix for compound morgan circular fingerprints

compounds = colnames(compiled_pa_matrix)[3:(ncol(compiled_pa_matrix)-4)]
mols = parse.smiles(compounds) # parse smiles to get a molecule object
circ_fing = lapply(mols, get.fingerprint, type="circular") # fingerpints

# similarity matrix w. tanimoto distance metric 
fing_sim_matx = fingerprint::fp.sim.matrix(circ_fing,method="tanimoto")
fing_sim_matx = as.data.frame(fing_sim_matx)
colnames(fing_sim_matx) = compounds
rownames(fing_sim_matx) = compounds

# distance matrix
dis_sim_converter = function(measure) {
  new_measure = 1 - measure
  return(new_measure)
}

fing_diss_matx = as.data.frame(dis_sim_converter(fing_sim_matx))

# save files 
saveRDS(fing_sim_matx,"fing_sim_matrx_gita.rds")
saveRDS(fing_diss_matx,"fing_diss_matrx_gita.rds")

fing_diss_matx = readRDS("fing_diss_matrx_gita.rds")






