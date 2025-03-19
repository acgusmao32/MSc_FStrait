#Data Cleaning - 18S - EUKARYOTES
#MSc FStrait - AGusmao

#Modified DataLoad - MWietz

library(gtools)
library(dplyr)
library(tibble)
library(tidyr)

#------------------- LOADING FILES -----
#DADA2 files
ASV.euk <- read.csv("Input/All_EUK_asv.txt",h = T, sep = "\t",check.names=F, row.names = 1)
TAX.euk <- read.table("Input/All_EUK_taxa.txt",h = T, sep = "\t", row.names =1)

#For Alpha
##OBS: The ones that includes singletons
ASV.euk.alph <- read.csv("Input/ASV_EUK_FALPHA.txt",h = T, sep = "\t",check.names=F, row.names = 1)
TAX.euk.alph <- read.table("Input/TAX_EUK_FALPHA.txt",h = T, sep = "\t", row.names =1)

#Metadata
MET.euk <- read.csv("Input/sequencing_metadata.txt", h=T, sep="\t", stringsAsFactors=F, skipNul=T) %>%
  filter(!locus_tag %in% c("16S")) %>%
  dplyr::select(c(
    "locus_tag","Mooring","Type","date1","date2","cycles","aftertrim_id"))

#------------------- EDIT TAX TABLES ----

# Rename PR2 taxranks: "Division" to "Phylum" 
# OK since taxnames are consistent with Silva
# enables cross-compatibility with 16S patterns

colnames(TAX.euk)<- c("Domain","Kingdom","Supergroup","Phylum","Class",
                      "Order","Family","Genus","Species")
TAX.euk$Supergroup <- NULL #remove supergroup

#How much ASVs at Kingdom level are NA
kingdom_na_count <- sum(is.na(TAX.euk$Kingdom)) 
kingdom_na_count #1588

TAX.euk.good <- TAX.euk[!is.na(TAX.euk[, 2]), ] #Remove ASVs with NA at Kingdom level
ASV.euk.good <- ASV.euk[rownames(TAX.euk.good),] #Match ASV and TAX tables

#For Alpha
colnames(TAX.euk.alph)<- c("Domain","Kingdom","Supergroup","Phylum","Class",
                           "Order","Family","Genus","Species")
TAX.euk.alph$Supergroup <- NULL #remove Supergroup

#------------------- REMOVE NEGATIVE CONTROLS ----

#There is no negative controls for the trap samples
#So it will only be removed from the RAS samples

# Define samples PCR-amplified with 25/30/35 cycles
NK25 <- MET.euk %>% filter(cycles=="25" & locus_tag=="18S"& Type =="RAS") %>% pull(aftertrim_id)
NK30 <- MET.euk %>% filter(cycles=="30" & locus_tag=="18S"& Type =="RAS") %>% pull(aftertrim_id)
NK35 <- MET.euk %>% filter(cycles=="35" & locus_tag=="18S"& Type =="RAS") %>% pull(aftertrim_id)

#Duplicate the table (transpose if necessary)
#ASVs needs to be rows
ASV.euk.good2 <- ASV.euk.good

# Define negative controls
ASV.euk.good2$NK25 <- rowMeans(ASV.euk.good2[,c(
  "NK_2880_S11","NK_82126_S35",
  "NK_PS126_25c","NK_negcontrol_25c")])
ASV.euk.good2$NK30 <- rowMeans(ASV.euk.good2[,c(
  "NK_PS126_30c",
  "NK_negcontrol_30c")])
ASV.euk.good2$NK35 <- rowMeans(ASV.euk.good2[,c(
  "NK_PS126_35c",
  "NK_negcontrol_35c")])

# Round the mean negative control values to the nearest whole number.
ASV.euk.good2$NK25 <- round(ASV.euk.good2$NK25, 0)
ASV.euk.good2$NK30 <- round(ASV.euk.good2$NK30, 0)
ASV.euk.good2$NK35 <- round(ASV.euk.good2$NK35, 0)

# Subtract negative counts 
asv1 = ASV.euk.good2[NK25] - ASV.euk.good2$NK25
asv2 = ASV.euk.good2[NK30] - ASV.euk.good2$NK30
asv3 = ASV.euk.good2[NK35] - ASV.euk.good2$NK35

#Subset the ASV with only the trap samples
trap_samples <- MET.euk %>% filter(Type == "Trap") %>% pull(aftertrim_id)
asv4 = ASV.euk.good2[, colnames(ASV.euk.good2) %in% trap_samples]

# Rejoin everything/ new ASV table
ASV.euk.clean <- cbind(asv1, asv2, asv3, asv4)

# Set negative values to zero
ASV.euk.clean[ASV.euk.clean < 0] <- 0

# Remove NegCtr columns 
ASV.euk.clean <- ASV.euk.clean[, !grepl('NK', names(ASV.euk.clean))]

#remove ASVs with zero abundance in all samples, caused by the negative controls
zero.ab.asvs <- rownames(ASV.euk.clean)[apply(ASV.euk.clean == 0, 1, all)]
cat("Number of ASVs with zero abundance across all samples:", length(zero.ab.asvs), "\n")

ASV.euk.clean <- ASV.euk.clean[!rownames(ASV.euk.clean) %in% zero.ab.asvs, ]
zero.ab.asvs2 <- rownames(ASV.euk.clean)[apply(ASV.euk.clean == 0, 1, all)]
cat("Number of ASVs with zero abundance across all samples:", length(zero.ab.asvs2), "\n")

# Match TAX and ASV table
TAX.euk.good <- TAX.euk.good[row.names(ASV.euk.clean),]

#------------------- REMOVE NEGATIVE CONTROLS - ALPHA ----

#Duplicate the table (transpose if necessary)
#ASVs needs to be rows

ASV.euk.alph.good <- ASV.euk.alph

# Define negative controls
ASV.euk.alph.good$NK25 <- rowMeans(ASV.euk.alph.good[,c(
  "NK_2880_S11","NK_82126_S35",
  "NK_PS126_25c","NK_negcontrol_25c")])
ASV.euk.alph.good$NK30 <- rowMeans(ASV.euk.alph.good[,c(
  "NK_PS126_30c",
  "NK_negcontrol_30c")])
ASV.euk.alph.good$NK35 <- rowMeans(ASV.euk.alph.good[,c(
  "NK_PS126_35c",
  "NK_negcontrol_35c")])

# Round the mean negative control values to the nearest whole number.
ASV.euk.alph.good$NK25 <- round(ASV.euk.alph.good$NK25, 0)
ASV.euk.alph.good$NK30 <- round(ASV.euk.alph.good$NK30, 0)
ASV.euk.alph.good$NK35 <- round(ASV.euk.alph.good$NK35, 0)

# Subtract negative counts 
asv5 = ASV.euk.alph.good[NK25] - ASV.euk.alph.good$NK25
asv6 = ASV.euk.alph.good[NK30] - ASV.euk.alph.good$NK30
asv7 = ASV.euk.alph.good[NK35] - ASV.euk.alph.good$NK35

#Subset the ASV with only the trap samples
trap_samples <- MET.euk %>% filter(Type == "Trap") %>% pull(aftertrim_id)
asv8 = ASV.euk.alph.good[, colnames(ASV.euk.alph.good) %in% trap_samples]

# Rejoin everything/ new ASV table
ASV.euk.alph.clean <- cbind(asv5, asv6, asv7, asv8)

# Set negative values to zero
ASV.euk.alph.clean[ASV.euk.alph.clean < 0] <- 0

# Remove NegCtr columns 
ASV.euk.alph.clean <- ASV.euk.alph.clean[, !grepl('NK', names(ASV.euk.alph.clean))]

#remove ASVs with zero abundance in all samples, caused by the negative controls
zero.ab.asvs <- rownames(ASV.euk.alph.clean)[apply(ASV.euk.alph.clean == 0, 1, all)]
cat("Number of ASVs with zero abundance across all samples:", length(zero.ab.asvs), "\n")

ASV.euk.alph.clean <- ASV.euk.alph.clean[!rownames(ASV.euk.alph.clean) %in% zero.ab.asvs, ]
zero.ab.asvs2 <- rownames(ASV.euk.alph.clean)[apply(ASV.euk.alph.clean == 0, 1, all)]
cat("Number of ASVs with zero abundance across all samples:", length(zero.ab.asvs2), "\n")

# Match TAX and ASV table
TAX.euk.alph <- TAX.euk.alph[row.names(ASV.euk.alph.clean),]

#------------------- REMOVE METAZOA ----

# Export Animalia/Metazoa to new DF
TAX.meta <- TAX.euk.good[grep('Metazoa', TAX.euk.good$Phylum),]
ASV.meta <- ASV.euk.clean[row.names(TAX.meta),]

# Remove Animalia/Metazoa from org table
TAX.euk.good <- TAX.euk.good[-grep('Metazoa', TAX.euk.good$Phylum),]

# Match TAX
ASV.euk.clean <- ASV.euk.clean[row.names(TAX.euk.good),]

#For Alpha
# Export Animalia/Metazoa to new DF
TAX.meta.alph <- TAX.euk.alph[grep('Metazoa', TAX.euk.alph$Phylum),]
ASV.meta.alph <- ASV.euk.alph.clean[row.names(TAX.meta.alph),]

# Remove Animalia/Metazoa from org table
TAX.euk.alph.good <- TAX.euk.alph[-grep('Metazoa', TAX.euk.alph$Phylum),]

# Match TAX after contaminant removal
ASV.euk.alph.clean <- ASV.euk.alph.clean[row.names(TAX.euk.alph.good),]

#---------------------- FORMAT TAXONOMY ----

# Rename NAs with last known taxrank + "UC"
k <- ncol(TAX.euk.good)-1
for (i in 2:k) {
  if (sum(is.na(TAX.euk.good[, i])) >1) {
    temp <- TAX.euk.good[is.na(TAX.euk.good[, i]), ]
    for (j in 1:nrow(temp)) {
      if (sum(is.na(
        temp[j, i:(k+1)])) == length(temp[j, i:(k+1)])) {
        temp[j, i] <- paste(temp[j, (i-1)], " UC", sep = "")
        temp[j, (i+1):(k+1)] <- temp[j, i]
      }
    }
    TAX.euk.good[is.na(TAX.euk.good[, i]), ] <- temp}
  if (sum(is.na(TAX.euk.good[, i]))==1) {
    temp <- TAX.euk.good[is.na(TAX.euk.good[, i]), ]
    if (sum(is.na(temp[i:(k+1)])) == length(temp[i:(k+1)])) {
      temp[i] <- paste(temp[(i-1)], " UC", sep="")
      temp[(i+1):(k+1)] <- temp[i]
    }
    TAX.euk.good[is.na(TAX.euk.good[, i]),] <- temp
  }
}
TAX.euk.good[is.na(TAX.euk.good[, (k+1)]), (k+1)] <- paste(
  TAX.euk.good[is.na(TAX.euk.good[, (k+1)]), k], " UC", sep="")

## Shorten/modify names 
TAX.euk.good <- TAX.euk.good %>%
  mutate(across(everything(),~gsub("Dino-Group-I-Clade-","Dino-I-", .))) %>%
  mutate(across(everything(),~gsub("Dino-Group-II-Clade-","Dino-II-", .))) %>%
  mutate(across(everything(),~gsub("Polar-centric-","", .))) %>%
  mutate(across(everything(),~gsub("Radial-centric-basal-","", .))) %>%
  mutate(across(everything(),~gsub("Chrysophyceae_Clade-","Chrysophyceae ", .))) %>%
  mutate(across(everything(),~gsub("Stephanoecidae_Group_","Stephanoecidae ", .))) %>%
  mutate(across(everything(),~gsub("Pirsonia_Clade_","Pirsonia", .))) %>%
  mutate(across(everything(),~gsub("_X|_XX|_XXX|_XXXX"," uc", .))) 

#---------------------- FORMAT TAXONOMY - ALPHA ----
# Rename NAs with last known taxrank + "UC"
k <- ncol(TAX.euk.alph.good)-1
for (i in 2:k) {
  if (sum(is.na(TAX.euk.alph.good[, i])) >1) {
    temp <- TAX.euk.alph.good[is.na(TAX.euk.alph.good[, i]), ]
    for (j in 1:nrow(temp)) {
      if (sum(is.na(
        temp[j, i:(k+1)])) == length(temp[j, i:(k+1)])) {
        temp[j, i] <- paste(temp[j, (i-1)], " UC", sep = "")
        temp[j, (i+1):(k+1)] <- temp[j, i]
      }
    }
    TAX.euk.alph.good[is.na(TAX.euk.alph.good[, i]), ] <- temp}
  if (sum(is.na(TAX.euk.alph.good[, i]))==1) {
    temp <- TAX.euk.alph.good[is.na(TAX.euk.alph.good[, i]), ]
    if (sum(is.na(temp[i:(k+1)])) == length(temp[i:(k+1)])) {
      temp[i] <- paste(temp[(i-1)], " UC", sep="")
      temp[(i+1):(k+1)] <- temp[i]
    }
    TAX.euk.alph.good[is.na(TAX.euk.alph.good[, i]),] <- temp
  }
}
TAX.euk.alph.good[is.na(TAX.euk.alph.good[, (k+1)]), (k+1)] <- paste(
  TAX.euk.alph.good[is.na(TAX.euk.alph.good[, (k+1)]), k], " UC", sep="")

## Shorten/modify names 
TAX.euk.alph.good <- TAX.euk.alph.good %>%
  mutate(across(everything(),~gsub("Dino-Group-I-Clade-","Dino-I-", .))) %>%
  mutate(across(everything(),~gsub("Dino-Group-II-Clade-","Dino-II-", .))) %>%
  mutate(across(everything(),~gsub("Polar-centric-","", .))) %>%
  mutate(across(everything(),~gsub("Radial-centric-basal-","", .))) %>%
  mutate(across(everything(),~gsub("Chrysophyceae_Clade-","Chrysophyceae ", .))) %>%
  mutate(across(everything(),~gsub("Stephanoecidae_Group_","Stephanoecidae ", .))) %>%
  mutate(across(everything(),~gsub("Pirsonia_Clade_","Pirsonia", .))) %>%
  mutate(across(everything(),~gsub("_X|_XX|_XXX|_XXXX"," uc", .))) 

#-------------------- EXPORT ----

write.table(ASV.euk.clean, file="ASV_EUK_AfC.txt", sep="\t", row.names=T, quote=F)
write.table(TAX.euk.good, file="TAX_EUK_AfC.txt", sep="\t", row.names=T, quote=F)

write.table(ASV.meta, file="ASV_MetaEUK_AfC.txt", sep="\t", row.names=T, quote=F)
write.table(TAX.meta, file="TAX_MetaEUK_AfC.txt", sep="\t", row.names=T, quote=F)

save.image("RDATA_cleaning_EUK_All.Rdata")

#[2]
write.table(ASV.euk.alph.clean, file="ASV_EUK_FALPHA_AfC.txt", sep="\t", row.names=T, quote=F)
write.table(TAX.euk.alph.good, file="TAX_EUK_FALPHA_AfC.txt", sep="\t", row.names=T, quote=F)

save.image("RDATA_Cleaning_EUK_FORALPH.Rdata")
