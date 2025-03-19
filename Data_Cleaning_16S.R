#DATA CLEANING - 16S/PROKARYOTES
#MSc FStrait - AnaGusmao

#Modified from DataLoad - M.Wietz

library(dplyr)
library(gtools)
library(tibble)
library(tidyr)

#------------------- LOADING FILES -----
#DADA2 Files
ASV.bac <- read.table("Input/seqtab_PROK_All.txt",h = T, sep = "\t", check.names=F, row.names = 1)
TAX.bac <- read.table("Input/tax_PROK_All.csv", h = T, sep = "\t", row.names = 1)

#For Alpha
#OBS: The ones that includes singletons
ASV.bac.alph <- read.table("Input/seqtab_PROK_FORALPHA.txt",h = T, sep = "\t", check.names=F, row.names = 1)
TAX.bac.alph <- read.table("Input/tax_PROK_FORALPHA.txt", h = T, sep = "\t", row.names = 1)

#Metadata
MET.bac <- read.csv("Input/sequencing_metadata.txt", h=T, sep="\t", stringsAsFactors=F, skipNul=T) %>%
  filter(!locus_tag %in% c("18S")) %>%
  dplyr::select(c(
    "locus_tag","Mooring","Type","date1","date2","cycles","aftertrim_id"))

#------------------- SOME CORRECTIONS ----

colnames(TAX.bac)<- c("Domain","Phylum","Class","Order","Family","Genus","Species")
ASV.bac.transp <- as.data.frame(t(ASV.bac))
TAX.match <- TAX.bac[!is.na(TAX.bac[, 2]), ] #Filter ASVs with unclassified phylum
TAX.match <- TAX.bac[row.names(ASV.bac.transp),] #Match TAX and ASV table

#For Alpha
colnames(TAX.bac.alph)<- c("Domain","Phylum","Class","Order","Family","Genus","Species")
ASV.alph.transp <- as.data.frame(t(ASV.bac.alph))
TAX.match2 <- TAX.bac.alph[!is.na(TAX.bac.alph[, 1]), ] #Filter ASVs with unclassified Domain
ASV.alph.transp <- ASV.alph.transp[rownames(TAX.match2), ] #Match ASV with TAX table

#------------------- REMOVE NEGATIVE CONTROLS ----

#There is no negative controls for the trap samples
#So it will only be removed from the RAS samples

# Define samples PCR-amplified with 25/30/35 cycles
NK25 <- MET.bac %>% filter(cycles == "25" & locus_tag == "16S" & Type == "RAS") %>% pull(aftertrim_id)
NK30 <- MET.bac %>% filter(cycles == "30" & locus_tag == "16S" & Type == "RAS") %>% pull(aftertrim_id)
NK35 <- MET.bac %>% filter(cycles == "35" & locus_tag == "16S" & Type == "RAS") %>% pull(aftertrim_id)

# Define corresponding negative controls
ASV.bac.transp$NK25 <- rowMeans(ASV.bac.transp[, c(
  "NK_RAS_16S_1", "NK_RAS_16S_2",
  "NK_RAS_bac_S39", "NK_RAS_bac_S49",
  "NK_RAS_negativ_controlbac",
  "NK_RAS_PS126_bac25c",
  "NK_RAS_negcontrol_bac25c")], na.rm = TRUE)
ASV.bac.transp$NK30 <- rowMeans(ASV.bac.transp[, c(
  "NK_RAS_negcontrol_bac30c",
  "NK_RAS_PS126_bac30c")], na.rm = TRUE)
ASV.bac.transp$NK35 <- rowMeans(ASV.bac.transp[, c(
  "NK_RAS_negcontrol_bac35c",
  "NK_RAS_PS126_bac35c")], na.rm = TRUE)

# Round
ASV.bac.transp$NK25 <- round(ASV.bac.transp$NK25, 0)
ASV.bac.transp$NK30 <- round(ASV.bac.transp$NK30, 0)
ASV.bac.transp$NK35 <- round(ASV.bac.transp$NK35, 0)

# Subtract negative controls only for RAS samples
asv1 = ASV.bac.transp[NK25] - ASV.bac.transp$NK25
asv2 = ASV.bac.transp[NK30] - ASV.bac.transp$NK30
asv3 = ASV.bac.transp[NK35] - ASV.bac.transp$NK35

#Subset the ASV with only the trap samples
trap_samples <- MET.bac %>% filter(Type == "Trap") %>% pull(aftertrim_id)
asv4 = ASV.bac.transp[, colnames(ASV.bac.transp) %in% trap_samples]

ASV.bac.new <- cbind(asv1, asv2, asv3, asv4) #Combine clean data

ASV.bac.new[ASV.bac.new < 0] <- 0 #set negative values to zero

# Remove NegCtr columns 
ASV.bac.new <- ASV.bac.new[, !grepl(
  'NK', names(ASV.bac.new))]

#remove ASVs with zero abundance in all samples, caused by the negative controls
zero.ab.asvs <- rownames(ASV.bac.new)[apply(ASV.bac.new == 0, 1, all)]
cat("Number of ASVs with zero abundance across all samples:", length(zero.ab.asvs), "\n")

ASV.bac.good <- ASV.bac.new[!rownames(ASV.bac.new) %in% zero.ab.asvs, ]
zero.ab.asvs2 <- rownames(ASV.bac.good)[apply(ASV.bac.good == 0, 1, all)]
#to confirm it was removed
cat("Number of ASVs with zero abundance across all samples:", length(zero.ab.asvs2), "\n")

# Match TAX and ASV table
TAX.match <- TAX.match[row.names(ASV.bac.good),]

#------------------- REMOVE NEGATIVE CONTROLS - ALPHA ----
#For Alpha 

# Define samples PCR-amplified with 25/30/35 cycles
NK25 <- MET.bac %>% filter(cycles == "25" & locus_tag == "16S" & Type == "RAS") %>% pull(aftertrim_id)
NK30 <- MET.bac %>% filter(cycles == "30" & locus_tag == "16S" & Type == "RAS") %>% pull(aftertrim_id)
NK35 <- MET.bac %>% filter(cycles == "35" & locus_tag == "16S" & Type == "RAS") %>% pull(aftertrim_id)

# Define corresponding negative controls
ASV.alph.transp$NK25 <- rowMeans(ASV.alph.transp[, c(
  "NK_RAS_16S_1", "NK_RAS_16S_2",
  "NK_RAS_bac_S39", "NK_RAS_bac_S49",
  "NK_RAS_negativ_controlbac",
  "NK_RAS_PS126_bac25c",
  "NK_RAS_negcontrol_bac25c")], na.rm = TRUE)
ASV.alph.transp$NK30 <- rowMeans(ASV.alph.transp[, c(
  "NK_RAS_negcontrol_bac30c",
  "NK_RAS_PS126_bac30c")], na.rm = TRUE)
ASV.alph.transp$NK35 <- rowMeans(ASV.alph.transp[, c(
  "NK_RAS_negcontrol_bac35c",
  "NK_RAS_PS126_bac35c")], na.rm = TRUE)

# Round
ASV.alph.transp$NK25 <- round(ASV.alph.transp$NK25, 0)
ASV.alph.transp$NK30 <- round(ASV.alph.transp$NK30, 0)
ASV.alph.transp$NK35 <- round(ASV.alph.transp$NK35, 0)

# Subtract negative controls only for RAS samples
asv5 = ASV.alph.transp[NK25] - ASV.alph.transp$NK25
asv6 = ASV.alph.transp[NK30] - ASV.alph.transp$NK30
asv7 = ASV.alph.transp[NK35] - ASV.alph.transp$NK35

#Subset the ASV with only the trap samples
trap_samples <- MET.bac %>% filter(Type == "Trap") %>% pull(aftertrim_id)
asv8 = ASV.alph.transp[, colnames(ASV.alph.transp) %in% trap_samples]

# Combine Cleaned Data
ASV.alph.new <- cbind(asv5, asv6, asv7, asv8)

# Set negative values to zero
ASV.alph.new[ASV.alph.new < 0] <- 0

# Remove NegCtr columns 
ASV.alph.new <- ASV.alph.new[, !grepl(
  'NK', names(ASV.alph.new))]

#remove ASVs with zero abundance in all samples, caused by the negative controls
zero.ab.asvs.alph <- rownames(ASV.alph.new)[apply(ASV.alph.new == 0, 1, all)]
cat("Number of ASVs with zero abundance across all samples:", length(zero.ab.asvs.alph), "\n")

ASV.alph.good <- ASV.alph.new[!rownames(ASV.alph.new) %in% zero.ab.asvs.alph, ]
zero.ab.asvs.alph2 <- rownames(ASV.alph.good)[apply(ASV.alph.good == 0, 1, all)]
cat("Number of ASVs with zero abundance across all samples:", length(zero.ab.asvs.alph2), "\n")

# Match TAX and ASV table
TAX.match2 <- TAX.match2[row.names(ASV.alph.good),]


#------------------- REMOVE MITOCHONDRIA/CHLOROPLAST ----

# Remove EUK
sum(grepl('Eukaryota', TAX.match$Domain, ignore.case = TRUE)) #Check how many Eukaryotas has
#TAX.match <- TAX.match[-grep('Eukaryota', TAX.match$Domain),] #Remove, if necessary

#Remove Mitochondria and chloroplast
TAX.good <- TAX.match[-grep('Mitochondria', TAX.match$Family),]

#Remove Chloroplast
TAX.good <- TAX.good[-grep('Chloroplast', TAX.good$Order),]

#Remove more potential contaminants
TAX.good <- TAX.good %>% filter(!Family %in% c(
    "Corynebacteriaceae","Bacillaceae",
    "Weeksellaceae","Enterococcaceae",
    "Streptococcaceae","Propionibacteriaceae",
    "Xanthobacteraceae","Staphylococcaceae", "Burkholderiaceae"))

# Match TAX after contaminant removal
ASV.bac.good <- ASV.bac.good[row.names(TAX.good),]

#------------------- REMOVE MITOCHONDRIA/CHLOROPLAST - ALPHA ----
# Remove EUK
sum(grepl('Eukaryota', TAX.match2$Domain, ignore.case = TRUE)) #Check how many Eukaryotas has
TAX.match2 <- TAX.match2[-grep('Eukaryota', TAX.match2$Domain),] #Remove, if necessary

#Remove Mitochondria and chloroplast
TAX.good2 <- TAX.match2[-grep('Mitochondria', TAX.match2$Family),]

#Remove Chloroplast
TAX.good2 <- TAX.good2[-grep('Chloroplast', TAX.good2$Order),]

#Remove more potential contaminants
TAX.good2 <- TAX.good2 %>% filter(!Family %in% c(
  "Corynebacteriaceae","Bacillaceae",
  "Weeksellaceae","Enterococcaceae",
  "Streptococcaceae","Propionibacteriaceae",
  "Xanthobacteraceae","Staphylococcaceae", "Burkholderiaceae"))

# Match TAX after contaminant removal
ASV.alph.good <- ASV.alph.good[row.names(TAX.good2),]

#------------------- FORMAT TAXONOMY TABLE ----

#Rename NAs with last known tax rank + "UC"
k <- ncol(TAX.good)-1
for (i in 2:k) {
  if (sum(is.na(TAX.good[, i])) >1) {
    temp <- TAX.good[is.na(TAX.good[, i]), ]
    for (j in 1:nrow(temp)) {
      if (sum(is.na(
        temp[j, i:(k+1)])) == length(temp[j, i:(k+1)])) {
        temp[j, i] <- paste(temp[j, (i-1)], " UC", sep = "")
        temp[j, (i+1):(k+1)] <- temp[j, i]
      }
    }
    TAX.good[is.na(TAX.good[, i]), ] <- temp}
  if (sum(is.na(TAX.good[, i]))==1) {
    temp <- TAX.good[is.na(TAX.good[, i]), ]
    if (sum(is.na(temp[i:(k+1)])) == length(temp[i:(k+1)])) {
      temp[i] <- paste(temp[(i-1)], " UC", sep="")
      temp[(i+1):(k+1)] <- temp[i]
    }
    TAX.good[is.na(TAX.good[, i]),] <- temp
  }
}
  TAX.good[is.na(TAX.good[, (k+1)]), (k+1)] <- paste(
  TAX.good[is.na(TAX.good[, (k+1)]), k], " UC", sep="")

# Shorten/modify names
TAX.good <- TAX.good %>%
  mutate(across(everything(),~gsub("Clade ","SAR11 Clade ", .))) %>%
  mutate(across(everything(),~gsub("_clade","", .))) %>%
  mutate(across(everything(),~gsub("Candidatus","Cand", .))) %>%
  mutate(across(everything(),~gsub("Roseobacter NAC11-7 lineage","NAC11-7", .))) %>%
  mutate(across(everything(),~gsub("_marine_group","", .))) %>%
  mutate(across(everything(),~gsub("_terrestrial_group","", .))) %>%
  mutate(across(everything(),~gsub("_CC9902","", .))) %>%
  mutate(across(everything(),~gsub("(Marine_group_B)","", ., fixed=T))) %>%
  mutate(across(everything(),~gsub("(SAR406)","SAR406", ., fixed=T)))


#------------------- FORMAT TAXONOMY TABLE - ALPHA ----
#Rename NAs with last known tax rank + "UC"
k <- ncol(TAX.good2)-1
for (i in 2:k) {
  if (sum(is.na(TAX.good2[, i])) >1) {
    temp <- TAX.good2[is.na(TAX.good2[, i]), ]
    for (j in 1:nrow(temp)) {
      if (sum(is.na(
        temp[j, i:(k+1)])) == length(temp[j, i:(k+1)])) {
        temp[j, i] <- paste(temp[j, (i-1)], " UC", sep = "")
        temp[j, (i+1):(k+1)] <- temp[j, i]
      }
    }
    TAX.good2[is.na(TAX.good2[, i]), ] <- temp}
  if (sum(is.na(TAX.good2[, i]))==1) {
    temp <- TAX.good2[is.na(TAX.good2[, i]), ]
    if (sum(is.na(temp[i:(k+1)])) == length(temp[i:(k+1)])) {
      temp[i] <- paste(temp[(i-1)], " UC", sep="")
      temp[(i+1):(k+1)] <- temp[i]
    }
    TAX.good2[is.na(TAX.good2[, i]),] <- temp
  }
}
TAX.good2[is.na(TAX.good2[, (k+1)]), (k+1)] <- paste(
  TAX.good2[is.na(TAX.good2[, (k+1)]), k], " UC", sep="")

# Shorten/modify names
TAX.good2 <- TAX.good2 %>%
  mutate(across(everything(),~gsub("Clade ","SAR11 Clade ", .))) %>%
  mutate(across(everything(),~gsub("_clade","", .))) %>%
  mutate(across(everything(),~gsub("Candidatus","Cand", .))) %>%
  mutate(across(everything(),~gsub("Roseobacter NAC11-7 lineage","NAC11-7", .))) %>%
  mutate(across(everything(),~gsub("_marine_group","", .))) %>%
  mutate(across(everything(),~gsub("_terrestrial_group","", .))) %>%
  mutate(across(everything(),~gsub("_CC9902","", .))) %>%
  mutate(across(everything(),~gsub("(Marine_group_B)","", ., fixed=T))) %>%
  mutate(across(everything(),~gsub("(SAR406)","SAR406", ., fixed=T)))



#------------------- EXPORT FILES ----

write.table(ASV.bac.good, "ASV_PROK_AfC.txt", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
write.table(TAX.good, "TAX_PROK_AfC.txt", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
save.image("RDATA_cleaning_PROK_All.Rdata")

#For Alpha
write.table(ASV.alph.good, file = "ASV_PROK_FALPHA_AfC.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(TAX.good2,file = "TAX_PROK_FALPHA_AfC.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
save.image("RDATA_cleaning_PROK_FORALPH.Rdata")


