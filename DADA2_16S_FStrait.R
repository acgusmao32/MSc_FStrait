#DADA2 Workflow (16S)
#MSc FStrait - Ana Gusmao
#Modified from Wietz.M

# [1] PROK Sed Trap
# [2] PROK Ras

install.packages('gridExtra')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dada2")
BiocManager::install("ShortReads")

library(dada2)
library(ShortRead)
library(ggplot2)
library(gridExtra)
library(Rcpp)

setwd("/isibhv/projects/FRAMdata/MolObs/SedTrap_AtlantECO/Traps_Amplicon/DADA2/PROK_Trap") 

#------------------ List files

#SED TRAP [1]
path1 <- "/isibhv/projects/FRAMdata/MolObs/SedTrap_AtlantECO/Clipped/Clipped_16S_Trap"
fns1 <- list.files(path)
fns1


#RAS [2]
path2 <- "/isibhv/projects/FRAMdata/MolObs/SedTrap_AtlantECO/Clipped/Clipped_16S_RAS"
fns2 <- list.files(path)
fns2

#------------------ Separating Forward and Reverse of 16S

#[1]
Fs_16S_t <- sort(list.files(path1, pattern="*.-16S-.*_R1_clean_clipped.fastq.gz"))  # Forward reads for 16S
Rv_16S_t <- sort(list.files(path1, pattern="*.-16S-.*_R2_clean_clipped.fastq.gz"))  # Reverse reads for 16S

# Define sample names
sampleNames1 <- sort(read.table("dada2_names_PROK_TRAP_fstrait.txt", h=F, stringsAsFactors=F)$V1)

# Specify the full path to Fwd & Rvs
Fs_16S_t <- file.path(pat1h, Fs_16S_t)
Rv_16S_t <- file.path(path1, Rv_16S_t)

#[2]
Fs_16S_r <- sort(list.files(path2, pattern="*._R1_.*.fastq.gz"))  # Forward reads for 16S
Rv_16S_r <- sort(list.files(path2, pattern="*._R2_.*.fastq.gz"))  # Reverse reads for 16S

# Define sample names
sampleNames2 <- sort(read.table("d2_names_PROK_RAS_fstrait.txt", h=F, stringsAsFactors=F)$V1)

# Specify the full path to Fwd & Rvs
Fs_16S_r <- file.path(path2, Fs_16S_r)
Rv_16S_r <- file.path(path2, Rv_16S_r)


#####----------------------- Quality Plots

#Visualize the Quality Profile to set the "truncLen" in the next step (often low-Q rev-reads)
#Since it takes a lot of time to generate them, it is better to just save everything in .pdf file

#[1]
#Forward
QualityProfileFs1 <- list()

for(i in 1:length(Fs_16S_t)) { #for each file on Fs_16S
  QualityProfileFs1[[i]] <- list() #create a new empty list
  QualityProfileFs1[[i]][[1]] <- plotQualityProfile(Fs_16S_t[i])} #and plot a quality profile
pdf("QualityProfileForward_TRAP.pdf") #saving into .pdf

for(i in 1:length(Fs_16S_t)) #organize the plots
{do.call("grid.arrange", QualityProfileFs1[[i]])}
dev.off()
rm(QualityProfileFs1)

#Reverse
QualityProfileRv1 <- list()

for(i in 1:length(Rv_16S_t)) {
  QualityProfileRv1[[i]] <- list()
  QualityProfileRv1[[i]][[1]] <- plotQualityProfile(Rv_16S_t[i])}
pdf("QualityProfileReverse_TRAP.pdf")

for(i in 1:length(Rv_16S_t)) 
{do.call("grid.arrange", QualityProfileRv1[[i]])}
dev.off()
rm(QualityProfileRv1)

#[2]
#Forward
QualityProfileFs2 <- list()

for(i in 1:length(Fs_16S_r)) { #for each file on Fs_16S
  QualityProfileFs2[[i]] <- list() #create a new empty list
  QualityProfileFs2[[i]][[1]] <- plotQualityProfile(Fs_16S_r[i])} #and plot a quality profile
pdf("QualityProfileForward_RAS.pdf") #saving into .pdf

for(i in 1:length(Fs_16S_r)) #organize the plots
{do.call("grid.arrange", QualityProfileFs2[[i]])}
dev.off()
rm(QualityProfileFs2)

#Reverse
QualityProfileRv2 <- list()

for(i in 1:length(Rv_16S_r)) {
  QualityProfileRv2[[i]] <- list()
  QualityProfileRv2[[i]][[1]] <- plotQualityProfile(Rv_16S_r[i])}
pdf("QualityProfileReverse_RAS.pdf")

for(i in 1:length(Rv_16S_r)) 
{do.call("grid.arrange", QualityProfileRv2[[i]])}
dev.off()
rm(QualityProfileRv2)

#####----------------------- Preparation for filtering

#define the directory path for the filtered/trimmed files
filt_path1 <- file.path("/isibhv/projects/FRAMdata/MolObs/SedTrap_AtlantECO/DADA2/PROK/Filtered_T")
filt_path2 <- file.path("/isibhv/projects/FRAMdata/MolObs/SedTrap_AtlantECO/DADA2/PROK/Filtered_R") 

#check if the directory exist, if not it create one
if(!file_test("-d", filt_path1)) dir.create(filt_path1)
if(!file_test("-d", filt_path2)) dir.create(filt_path2)

#Assign the file names and the directory path
filtFs1 <- file.path1(filt_path, paste0(sampleNames1, "_F_filt.fastq"))
filtRv1 <- file.path1(filt_path, paste0(sampleNames1, "_R_filt.fastq"))

filtFs2 <- file.path2(filt_path, paste0(sampleNames2, "_F_filt.fastq"))
filtRv2 <- file.path2(filt_path, paste0(sampleNames2, "_R_filt.fastq"))

#####-----------------------Calculate the overlap

#The minimum necessary is 20 bp of overlap + biological variation.
#Info of the samples: Primer (515F - 926R), Sequencing 2x300
#Amplicon ~400bp

#[1]
#truncLen = c(250, 190): 515F+250bp: until RNA position 765. 926R-190bp: until RNA position 736
#765 - 736 = min overlap up to 39 in merging

#[2]
#truncLen values = F250 + R190 = 440bp
#Overlap = 440bp - 400 = 40bp

#####----------------------- Filter

# We played around with minQ, maxEE and truncQ, without big differences (might differ between datasets though)

# Standard values: 
# Maximum Number of ambiguous nucleotides (maxN) = 0
# Minimum Quality Score for any base (minQ)
# Maximum number of expected errors per read (maxEE) = 2
# Truncation Quality Score (trunQ) = 2

#Error message? often something about "cores": try multithread=F

#[1]
out1 <- filterAndTrim(
  Fs_16S_t, filtFs1, 
  Rv_16S_t, filtRv1,
  truncLen = c(250, 190),
  maxN = 0, 
  minQ = 2,
  maxEE = c(3, 3), 
  truncQ = 0, 
  rm.phix = T, #Remove reads that match the phix (bacteriophage) genome
  compress = F,
  multithread = F)  

#Check the results:
#Should retain >70% -- OK!

head(out1)
summary(out1[, 2]/out1[, 1])

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9157  0.9385  0.9534  0.9483  0.9598  0.9655 
 
#Median: Half of the samples retained more than 95%


#[2]
out2 <- filterAndTrim(
  Fs_16S_r, filtFs2, 
  Rv_16S_r, filtRv2,
  truncLen = c(250, 190),
  maxN = 0, 
  minQ = 2,
  maxEE = c(3, 3), 
  truncQ = 0, 
  rm.phix = T, #Remove reads that match the phix (bacteriophage) genome
  compress = F,
  multithread = F)  

print(filtFs2)
#Check the results:
#Should retain >70% -- OK!

head(out2)
summary(out2[, 2]/out2[, 1])

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.4888  0.8262  0.8828  0.8565  0.9011  0.9585

#Median: Half of the samples retained more than 88%

#####----------------------- Quality Check

#Visualize the Quality Profiles of the filtered/trimmed files

#[1]
#Forward
QualityProfileFs.filt1 <- list()

for(i in 1:length(filtFs1)) {
  QualityProfileFs.filt1[[i]] <- list()
  QualityProfileFs.filt1[[i]][[1]] <- plotQualityProfile(filtFs1[i])}
pdf("QualityProfileForwardFiltered_TRAP.pdf")

for(i in 1:length(filtFs1)) {do.call(
  "grid.arrange", QualityProfileFs.filt1[[i]])}
dev.off()
rm(QualityProfileFs.filt1)

#Reverse
QualityProfileRv.filt1 <- list()

for(i in 1:length(filtRv1)) {
  QualityProfileRv.filt1[[i]] <- list()
  QualityProfileRv.filt1[[i]][[1]] <- plotQualityProfile(filtRv1[i])}
pdf("QualityProfileReverseFiltered_TRAP.pdf")

for(i in 1:length(filtRv1)) {do.call(
  "grid.arrange", QualityProfileRv.filt1[[i]])}
dev.off()
rm(QualityProfileRv.filt1)

#[2]
#Forward
QualityProfileFs.filt2 <- list()

for(i in 1:length(filtFs2)) {
  QualityProfileFs.filt2[[i]] <- list()
  QualityProfileFs.filt2[[i]][[1]] <- plotQualityProfile(filtFs2[i])}
pdf("QualityProfileForwardFiltered_RAS.pdf")

for(i in 1:length(filtFs2)) {do.call(
  "grid.arrange", QualityProfileFs.filt2[[i]])}
dev.off()
rm(QualityProfileFs.filt2)

#Reverse
QualityProfileRv.filt2 <- list()

for(i in 1:length(filtRv2)) {
  QualityProfileRv.filt2[[i]] <- list()
  QualityProfileRv.filt2[[i]][[1]] <- plotQualityProfile(filtRv2[i])}
pdf("QualityProfileReverseFiltered_RAS.pdf")

for(i in 1:length(filtRv2)) {do.call(
  "grid.arrange", QualityProfileRv.filt2[[i]])}
dev.off()
rm(QualityProfileRv.filt2)

#####----------------------- Errors Rates

#Learn errors 
#err: parametric error model
#verbose: how much details will show during processing
#max_consist: number of iteraction to refine the error

#[1]
errF1 <- learnErrors(filtFs1, multithread=T,randomize=T, verbose=1, MAX_CONSIST=20)
errR1 <- learnErrors(filtRv1, multithread=T, randomize=T, verbose=1, MAX_CONSIST=20)

#[2]
errF2 <- learnErrors(filtFs2, multithread=T,randomize=T, verbose=1, MAX_CONSIST=20)
errR2 <- learnErrors(filtRv2, multithread=T, randomize=T, verbose=1, MAX_CONSIST=20)

# Plot error profiles
#nominalQ: the expected error rates based on the raw quality scores assigned by the sequencing machine.

#[1]
pdf("ErrorProfiles_TRAP.pdf")
plotErrors(errF1, nominalQ=T)
plotErrors(errR1, nominalQ=T)
dev.off()

#[2]
pdf("ErrorProfiles_RAS.pdf")
plotErrors(errF2, nominalQ=T)
plotErrors(errR2, nominalQ=T)
dev.off()

# Check the results:
# convergence after 6/6 rounds - ok!
# few outliers outside black line - ok!

#####----------------------- Dereplication
#Step removed in newer versions, but can be useful
#See https://github.com/benjjneb/dada2/issues/1095#issuecomment-671333174

#[1]
derepFs1 <- derepFastq(filtFs1, verbose=T)
derepRv1 <- derepFastq(filtRv1, verbose=T)

# Rename by clip-filenames
names(derepFs1) <- sampleNames1
names(derepRv1) <- sampleNames1

#[2]
derepFs2 <- derepFastq(filtFs2, verbose=T)
derepRv2 <- derepFastq(filtRv2, verbose=T)

# Rename by clip-filenames
names(derepFs2) <- sampleNames2
names(derepRv2) <- sampleNames2

#####----------------------- Denoising
#Remove sequence errors
#many samples / memory error: consider pool="pseudo"

dadaFs1 <- dada(derepFs1, err=errF, multithread=T, pool=T)     
dadaRv1 <- dada(derepRv1, err=errR, multithread=T, pool=T)

dadaFs2 <- dada(derepFs2, err=errF, multithread=T, pool=T)     
dadaRv2 <- dada(derepRv2, err=errR, multithread=T, pool=T)

#####----------------------- Read Merging

#[1]
mergers1 <- mergePairs(
  dadaFs1, 
  derepFs1, 
  dadaRv1,
  derepRv1,
  minOverlap=20, #default now 12
  verbose=T)

# Create sequence table
seqtab1 <- makeSequenceTable(mergers1)

# 50 samples -- 29911 sequences
dim(seqtab1) 

# Save
saveRDS(seqtab1,"seqtab_PROK_Trap.rds")

#[2]
mergers2 <- mergePairs(
  dadaFs2, 
  derepFs2, 
  dadaRv2,
  derepRv2,
  minOverlap=20, #default now 12
  verbose=T)

# Create sequence table
seqtab2 <- makeSequenceTable(mergers2)

#109 samples were pooled: 12798 reads
dim(seqtab2) 

# Save
saveRDS(seqtab2,"seqtab_RAS_PROK.rds")


####----------------------- Summary stats

#[1]
getN1 <- function(x) sum(getUniques(x))
track1 <- cbind(out, sapply(
  dadaFs1, getN1), sapply(mergers, getN1), 
  rowSums(seqtab1))
colnames(track1) <- c(
  "input","filtered","denoised",
  "merged","tabled")
rownames(track1) <- sampleNames1
track1 <- data.frame(track1)
head(track1)

write.table(track1, "dadastats_PROK_Trap.txt", quote=F, sep="\t")

save.image("PROK_Trap.Rdata")

#[2]
getN2 <- function(x) sum(getUniques(x))
track2 <- cbind(out, sapply(
  dadaFs2, getN2), sapply(mergers, getN2), 
  rowSums(seqtab2))
colnames(track2) <- c(
  "input","filtered","denoised",
  "merged","tabled")
rownames(track2) <- sampleNames2
track <- data.frame(track2)
head(track2)

write.table(track2, "dadastats_PROK_RAS_Amplicon.txt", quote=F, sep="\t")

save.image("RAS_PROK.Rdata")

#######
#---------------MERGE THE DATASETS
# To create one ASV table from multiple MiSeq runs

#We now merge the runs to create a joint ASV table, remove chimeras, and assign taxonomy
#If only 1 run was processed: start directly with removeBimeraDenovo

# Load seqtab from each MiSeq run
sq1 <- readRDS("seqtab_PROK_Trap.rds")
sq2 <- readRDS("seqtab_RAS_PROK.rds")

## Merge
seqtab <- mergeSequenceTables(
  sq1,sq2, repeats="error")

#------------------Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread=T, verbose=T)

#Merged dataset -> 23500 bimeras out of 40122 sequences

dim(seqtab.nochim)  # 159 samples -- 16622 seqs

summary(rowSums(seqtab.nochim)/rowSums(seqtab))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.8910  0.9801  0.9850  0.9791  0.9896  1.0000


#------------------ Remove singletons and "junk" sequences

#Determine amplicon length/size range
table(rep(nchar(
  colnames(seqtab.nochim)), 
  colSums(seqtab.nochim)))

# Remove singletons and "junk" sequences
# "c" adjusted to size range of amplicons, have a look in the output of the last command
# For alpha-div metrics, do NOT remove singletons

seqtab.nochim.all <- seqtab.nochim[, nchar( #removing singletons
  colnames(seqtab.nochim)) %in% c(354:412) & 
    colSums(seqtab.nochim) > 1]

seqtab.nochim.all.foralph <- seqtab.nochim[, nchar( #not removing singletons, only junk sequences. For alpha-div.
  colnames(seqtab.nochim)) %in% c(354:412)]

dim(seqtab.nochim.all) #Result: 159 samples - 15024 seqs
dim(seqtab.nochim.all.foralph) #159 samples - 16345 seqs


summary(rowSums(seqtab.nochim.all))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#6   64364  100930  124615  161834  645142

summary(rowSums(seqtab.nochim.all)/rowSums(seqtab.nochim))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9784  0.9995  0.9999  0.9993  1.0000  1.0000 

summary(rowSums(seqtab.nochim.all.foralph)/rowSums(seqtab.nochim))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9784  0.9996  1.0000  0.9994  1.0000  1.0000  

#------------------ TAXONOMY -- Silva v138.2 ## 

tax <- assignTaxonomy (seqtab.nochim.all, "/isibhv/projects/FRAMdata/MolObs/SedTrap_AtlantECO/scripts/silva_nr99_v138.2_toSpecies_trainset.fa.gz", 
  tryRC = T, multithread = T)

tax2 <- assignTaxonomy (seqtab.nochim.all.foralph, "/isibhv/projects/FRAMdata/MolObs/SedTrap_AtlantECO/scripts/silva_nr99_v138.2_toSpecies_trainset.fa.gz", 
                       tryRC = T, multithread = T)


table(tax[, 1]) #[Domain]: Archaea 443 | Bacteria 14532 | Eukaryota 38
table(tax2[, 1]) #[Domain]: Archaea 464 | Bacteria 15828 | Eukaryota 42 

#------------------ Remove NA on phylum level 
#ACHTUNG:: [Avoid this step for alpha-div metrics]

sum(is.na(tax[, 2])) #how many rows has NA values at phylum level: 217

tax.good <- tax[!is.na(tax[, 2]),] #remove

seqtab.nochim.allgood <- seqtab.nochim.all[,rownames(tax.good)]

summary(rowSums(seqtab.nochim.allgood))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#6   64354  100923  124529  161827  643774 

#------------------ Format tables
seqtab.nochim.allprint <- seqtab.nochim.allgood
seqtab.print.foralph <- seqtab.nochim.all.foralph

tax.print <- tax.good
tax2.print <- tax2

all.equal(colnames(seqtab.nochim.allprint), rownames(tax.print)) #TRUE
all.equal(colnames(seqtab.print.foralph), rownames(tax2.print)) #TRUE 

#------------------ Changing the sequences to ASV+Number
colnames(seqtab.nochim.allprint) <- paste(
  "asv", 1:ncol(seqtab.nochim.allgood), sep = "")

rownames(tax.print) <- paste("asv", 1:ncol(seqtab.nochim.allgood), sep = "")

#For the alpha-div file
colnames(seqtab.print.foralph) <- paste(
  "asv", 1:ncol(seqtab.nochim.all.foralph), sep = "")

rownames(tax2.print) <- paste("asv", 1:ncol(seqtab.nochim.all.foralph), sep = "")

#------------------ Export

write.csv(seqtab.nochim.allprint,"seqtab_PROK_All.csv", sep="/t", quote=F) #ASV table
write.csv(tax.print,"tax_PROK_All.csv", sep="/t", quote=F) #Taxonomy Table

uniquesToFasta(seqtab.nochim.allgood,"uniques_PROK_All.fasta")

save.image("PROK_All.Rdata")


#The following files were not removed singletons or filtered (removed NA at phylum level)
#For Alpha-diversity metrics
write.csv(seqtab.print.foralph,"seqtab_PROK_FORALPHA.csv", sep="/t", quote=F) #ASV table for Alpha-div 
write.csv(tax2.print,"tax_PROK_FORALPHA.csv", sep="/t", quote=F) #Taxonomy Table for Alpha-div

uniquesToFasta(seqtab.nochim.all.foralph,"uniques_PROK_FORALPHA.fasta") 

save.image("PROK_FORALPHA.Rdata")


###############################################################################

