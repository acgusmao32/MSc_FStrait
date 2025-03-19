#Rarefaction
#Rarefaction once by sequencing depth

# -------- Loading R DATA ----
load("Input/phy_prok.RData")
load("Input/phy_euk.RData")

#Including singletons and unclassified taxa at phylum or kingdom level
load("Input/phy_prok_alph.RData")
load("Input/phy_euk_alph.RData")

# -------- Subsets ----
phy_prok_trap <- subset_samples(phy_prok_fl, Habitat == "Particle-Associated")
phy_prok_ras <- subset_samples(phy_prok_fl, Habitat == "Free-living Surface Water")

phy_euk_trap <- subset_samples(phy_euk, Habitat == "Particle-Associated")
phy_euk_ras <- subset_samples(phy_euk, Habitat == "Free-living Surface Water")

#----------------------- NORMALIZATION - RAREFACTION ----
#Lower value > 10.000 reads

#Prok All
head(sort(sample_sums(phy_prok_fl)), 5) #13614
size_prok <- head(sort(sample_sums(phy_prok_fl)), 5)[2] #RAS_HS1_P3365

p.rar <- rarefy_even_depth(phy_prok_fl, trimOTUs = T, sample.size = size_prok,
                           rngseed = 42)

#Euk All
head(sort(sample_sums(phy_euk)), 5) #11066
size_euk <- head(sort(sample_sums(phy_euk)), 5)[4] 
#Removed: Fevi36 Up15, Fevi 36 Up09, HS3_PS121_E35
e.rar <- rarefy_even_depth(phy_euk, trimOTUs = T, sample.size = size_euk,
                           rngseed = 42)

#ProkTrap
head(sort(sample_sums(phy_prok_trap)), 5) #29577
size_prok1 <- head(sort(sample_sums(phy_prok_trap)), 5)[1] #no sample removed
pt.rar <- rarefy_even_depth(phy_prok_trap, trimOTUs = T, sample.size = size_prok1,
                            rngseed = 42)

#ProkRAS
head(sort(sample_sums(phy_prok_ras)), 5) #13614
size_prok2 <- head(sort(sample_sums(phy_prok_ras)), 5)[2] #RAS_HS1_P3365
pr.rar <- rarefy_even_depth(phy_prok_ras, trimOTUs = T, sample.size = size_prok2,
                            rngseed = 42)

#EukTrap
head(sort(sample_sums(phy_euk_trap)), 5) #21058
size_euk1 <- head(sort(sample_sums(phy_euk_trap)), 5)[3] #Fevi36 Up15 and Fevi36 Up09
et.rar <- rarefy_even_depth(phy_euk_trap, trimOTUs = T, sample.size = size_euk1, 
                            rngseed = 42)

#EukRAS
head(sort(sample_sums(phy_euk_ras)), 5) #11066
size_euk2 <- head(sort(sample_sums(phy_euk_ras)), 5)[2] #HS3 Ps121 E35
er.rar <- rarefy_even_depth(phy_euk_ras, trimOTUs = T, sample.size = size_euk2, 
                            rngseed = 42)

#For Alpha Div

#Prokaryotes
head(sort(sample_sums(phy_prok_fl_alph)), 5) #First 5 lower numbers
size_prokal <- head(sort(sample_sums(phy_prok_fl_alph)), 5)[2] #14175
#Sample removed: RAS_HS1_P3365 (7566)

# Rarefy
p.rar.alph <- rarefy_even_depth(phy_prok_fl_alph, trimOTUs = T, 
                           sample.size = size_prokal, rngseed = 42)

#Eukaryotes 

head(sort(sample_sums(phy_euk_alph)), 5) #First 5 lower number
size_eukal <- head(sort(sample_sums(phy_euk_alph)), 5)[4] #11074

#Samples removed: Fevi36_Up15_euk (1041), HS3_PS121_E35 (1992), Fevi36_Up09_euk (4674)

#Rarefy
e.rar.alph <- rarefy_even_depth(phy_euk_alph, trimOTUs = T, 
                                sample.size = size_eukal, rngseed = 42)

# -------- Export ----

#Alpha Div
save(p.rar.alph, file = "p_rar_alph.RData") #[1]
save(e.rar.alph, file = "e_rar_alph.RData") #[2]

# -------- Stats ----

#Number of ASVs
#Obs after cleaning, but without removing singletons
#Also it was not "filtered", so it was not removed ASVs with NA in the kingdom level

phy_prok_alph #15084 ASVs in 148 samples
phy_euk_alph #12348 ASVs in 139 samples

#Sequence Counts per sample (RAW ASV TABLE)

# Prok
mean(sample_sums(phy_prok_alph)) #Average: 127599
sd(sample_sums(phy_prok_alph)) # SD: 94707.79
sum(sample_sums(phy_prok_alph)) #Total number of sequences: 18884724

sample_counts <- sample_sums(phy_prok_alph) #Sequences per sample
sample_summary_prok <- data.frame(SampleID = names(sample_counts), Total_Sequences = sample_counts) #dataframe
rm(sample_counts)
write.csv(sample_summary_prok, "Sequence_Summary_PROK.csv", row.names = FALSE) #Export Table

# Euk
mean(sample_sums(phy_euk_alph)) #Average: 113741.7
sd(sample_sums(phy_euk_alph)) # SD: 69154.13
sum(sample_sums(phy_euk_alph)) #Total number of sequences: 15810102

sample_counts <- sample_sums(phy_euk_alph) #Sequences per sample
sample_summary_euk <- data.frame(SampleID = names(sample_counts), Total_Sequences = sample_counts) #dataframe
rm(sample_counts)
write.csv(sample_summary_euk, "Sequence_Summary_EUK.csv", row.names = FALSE) #Export Table



#Hellinger transformation and Ordination ----

#Hellinger
prok.hel <- decostand(as(otu_table(p.rar), "matrix"), method = "hellinger")
proktrap.hel <- decostand(as(otu_table(pt.rar), "matrix"), method = "hellinger")
prokras.hel <- decostand(as(otu_table(pr.rar), "matrix"), method = "hellinger")

euk.hel <- decostand(as(otu_table(e.rar), "matrix"), method = "hellinger")
euktrap.hel <- decostand(as(otu_table(et.rar), "matrix"), method = "hellinger")
eukras.hel <- decostand(as(otu_table(er.rar), "matrix"), method = "hellinger")

#recreating phy object
nphy_pall <- phyloseq(otu_table(prok.hel, taxa_are_rows = TRUE), tax_table(phy_prok_fl), sample_data(phy_prok_fl))
nphy_pktrap <- phyloseq(otu_table(proktrap.hel, taxa_are_rows = TRUE), tax_table(phy_prok_trap), sample_data(phy_prok_trap))
nphy_pkras <- phyloseq(otu_table(prokras.hel, taxa_are_rows = TRUE), tax_table(phy_prok_ras), sample_data(phy_prok_ras))

nphy_eall <- phyloseq(otu_table(euk.hel, taxa_are_rows = TRUE), tax_table(phy_euk), sample_data(phy_euk))
nphy_ektrap <- phyloseq(otu_table(euktrap.hel, taxa_are_rows = TRUE), tax_table(phy_euk_trap), sample_data(phy_euk_trap))
nphy_ekras <- phyloseq(otu_table(eukras.hel, taxa_are_rows = TRUE), tax_table(phy_euk_ras), sample_data(phy_euk_ras))

#Export phy objects
save(nphy_pall, file = "nphy_pall.RData") #[1]
save(nphy_eall, file = "nphy_eall.RData") #[2]

save(nphy_pktrap, file = "nphy_pktrap.RData") #[1]
save(nphy_pkras, file = "nphy_pkras.RData") #[2]
save(nphy_ektrap, file = "nphy_ektrap.RData") #[3]
save(nphy_ekras, file = "nphy_ekras.RData") #[3]