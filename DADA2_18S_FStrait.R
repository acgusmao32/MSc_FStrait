#DADA2 for Microbial Eukaryotes
#MSc Fram Strait - Project
#Modified from Stefan Neuhaus

#OBS: With previously clipped files using cutadapt

#[1] Sediment Traps
#[2] RAS

# load needed packages
library(dada2) #Version: 1.32.0
library(stringr) #Version: 1.5.1
library(ggplot2) #Version: 3.5.1
library(ShortRead)#Version 1.62.0
library(Biostrings)#Version: 2.72.1

#---------------------  PRIMERS

#Euk 18S V4 (Fadeev et al., 2018) -> 528iF/964iR
FWD_PRIMER = "GCGGTAATTCCAGCTCCAA"
REV_PRIMER = "ACTTTCGTTCTTGATYRR"

#---------------------ORGANIZING THE FILES

# set data directory

seq_data_dir1 <- "/isibhv/projects/FRAMdata/MolObs/SedTrap_AtlantECO/Clipped/Traps_Amplicon/Clipped_18S"
seq_data_dir2 <- "/isibhv/projects/FRAMdata/MolObs/SedTrap_AtlantECO/Clipped/RAS_Amplicon/Clipped_18S"

#Folder name of sequencing run to denoise 
# e.g. "M03457_0101_000000000-K3N5M" (AWI convention)
seq_run1 <- "MSC_Ana_Traps_FStrait_METAB_18S" #[1]
seq_run2 <- "MSC_Ana_RAS_FStrait_METAB_18S" #[2]

#Specify segment positions in file name of each set to be sample names
# "_" and "-" are delimiters
#example: RV-AWI5_S32_L001_R1_001.fastq.gz => sgmt=c(1) => "RV-AWI5"
#example2: HE542_St501_11mu_L001_R1_001.fastq.gz => sgmt=c(2,3) => St501_11mu

sgmt1 = c(4,5) #[1]
sgmt2 <- c(2, 3, 6) #[2]        

# define file name suffixes to remove them from the sample name

#[1]
R1_suffix1 <- "_R1_clean_clipped.fastq.gz"
R2_suffix1 <- "_R2_clean_clipped.fastq.gz"

#[2]
R1_suffix2 <- "_R1_001_trimmed.fastq.gz"
R2_suffix2 <- "_R2_001_trimmed.fastq.gz"

#PATHES FOR DIRECTORIES

#define working directory as project directory
project_dir <- file.path(getwd())

#path to raw data folder to work with
raw_dir1 <- file.path(seq_data_dir1) #[1]
raw_dir2 <- file.path(seq_data_dir2) #[2]

#temporary data dir (delete after analysis run) for results of intermidate steps 

temp_dir1 <- file.path(project_dir,paste(seq_run1,"temp",sep="_"))
qualFiltTrim_dir1 <- file.path(temp_dir1,"qualFiltTrim")

temp_dir2 <- file.path(project_dir,paste(seq_run2,"temp",sep="_"))
qualFiltTrim_dir2 <- file.path(temp_dir2,"qualFiltTrim")

#create directories for intermediate steps
if(!dir.exists(temp_dir1)) dir.create(temp_dir1)
if(!dir.exists(temp_dir2)) dir.create(temp_dir2)

if(!dir.exists(qualFiltTrim_dir1)) dir.create(qualFiltTrim_dir1)
if(!dir.exists(qualFiltTrim_dir2)) dir.create(qualFiltTrim_dir2)

#Construct needed file lists

#[1]
#basenames of all files
fnFs.files1 <- basename(sort(list.files(raw_dir1, pattern =R1_suffix1, full.names = TRUE)))
fnRs.files1 <- basename(sort(list.files(raw_dir1, pattern =R2_suffix1, full.names = TRUE)))
#pathes of raw files to analyse
fnFs.raw1 <- file.path(raw_dir1,fnFs.files1)
fnRs.raw1 <- file.path(raw_dir1,fnRs.files1)

# file pathes of intermediate result files
fnFs.qualFiltTrim1 <- file.path(qualFiltTrim_dir1,basename(fnFs.raw1))
fnRs.qualFiltTrim1 <- file.path(qualFiltTrim_dir1,basename(fnRs.raw1))
#create sample names
sample.names1 <- unname(sapply(fnFs.files1,function(filenames) 
  paste(strsplit(str_remove(filenames,R1_suffix1),"[_-]")[[1]][sgmt1], collapse="_")))
sample.names1

#create file to sample name mapping
file.list1 <- cbind(sample.names1,fnFs.files1, fnRs.files1)
file.list1

#[2]
#basenames of all files
fnFs.files2 <- basename(sort(list.files(raw_dir2, pattern =R1_suffix2, full.names = TRUE)))
fnRs.files2 <- basename(sort(list.files(raw_dir2, pattern =R2_suffix2, full.names = TRUE)))
#pathes of raw files to analyse
fnFs.raw2 <- file.path(raw_dir2,fnFs.files2)
fnRs.raw2 <- file.path(raw_dir2,fnRs.files2)

# file pathes of intermediate result files
fnFs.qualFiltTrim2 <- file.path(qualFiltTrim_dir2,basename(fnFs.raw2))
fnRs.qualFiltTrim2 <- file.path(qualFiltTrim_dir2,basename(fnRs.raw2))
#create sample names
sample.names2 <- unname(sapply(fnFs.files2,function(filenames) 
  paste(strsplit(str_remove(filenames,R1_suffix2),"[_-]")[[1]][sgmt2], collapse="_")))
sample.names2

#create file to sample name mapping
file.list2 <- cbind(sample.names2,fnFs.files2, fnRs.files2)
file.list2

#---------------QUALITY PLOTS OF "RAW" FILES

#Plot Phred quality score profiles; always subset to reduce computational costs
#Data are aggregate = instead of showing a plot for each read, it does an average of all them. 

#[1]
QualityProfileFs1 <- list()

for(i in 1:length(fnFs.files1)) { #for each file on Fs_16S
  QualityProfileFs1[[i]] <- list() #create a new empty list
  QualityProfileFs1[[i]][[1]] <- plotQualityProfile(file.path(raw_dir1,fnFs.files1[i]))} #and plot a quality profile
pdf("QualityProfileForward_TRAPS.pdf") #saving into .pdf

for(i in 1:length(fnFs.files1)) #organize the plots
{do.call("grid.arrange", QualityProfileFs1[[i]])}
dev.off()
rm(QualityProfileFs1)

#Reverse
QualityProfileRs1 <- list()

for(i in 1:length(fnRs.files1)) {
  QualityProfileRs1[[i]] <- list() #create a new empty list
  QualityProfileRs1[[i]][[1]] <- plotQualityProfile(file.path(raw_dir1,fnRs.files1[i]))}
pdf("QualityProfileReverse_TRAPS.pdf") #saving into .pdf

for(i in 1:length(fnRs.files1)) #organize the plots
{do.call("grid.arrange", QualityProfileRs1[[i]])}
dev.off()
rm(QualityProfileRs1)

#[2]
QualityProfileFs2 <- list()

for(i in 1:length(fnFs.files2)) { #for each file on Fs_16S
  QualityProfileFs2[[i]] <- list() #create a new empty list
  QualityProfileFs2[[i]][[1]] <- plotQualityProfile(file.path(raw_dir2,fnFs.files2[i]))} #and plot a quality profile
pdf("QualityProfileForward_RAS.pdf") #saving into .pdf

for(i in 1:length(fnFs.files2)) #organize the plots
{do.call("grid.arrange", QualityProfileFs2[[i]])}
dev.off()
rm(QualityProfileFs2)

#Reverse
QualityProfileRs2 <- list()

for(i in 1:length(fnRs.files2)) {
  QualityProfileRs2[[i]] <- list() #create a new empty list
  QualityProfileRs2[[i]][[1]] <- plotQualityProfile(file.path(raw_dir2,fnRs.files2[i]))}
pdf("QualityProfileReverse_RAS.pdf") #saving into .pdf

for(i in 1:length(fnRs.files2)) #organize the plots
{do.call("grid.arrange", QualityProfileRs2[[i]])}
dev.off()
rm(QualityProfileRs2)

#--------------- PRIMER CHECKS
#If necessary, to check cutadapt

#Create all possible primer versions (fwd, fwd_complement, rev, rev_complement)

allOrients <- function(primer) {
  #Create all orientations of the input sequence
  #require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  
}

#define and display all primer versions
FWD_PRIMER.orients <- allOrients(FWD_PRIMER)
REV_PRIMER.orients <- allOrients(REV_PRIMER)
FWD_PRIMER.orients 
REV_PRIMER.orients 

#check primer situation in sample subset of sequence files
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
for(i in 1:nrow(file.list)){
  colnames(file.list) <- NULL
  print(file.list[i,1])
  print(rbind(FWD_PRIMER.ForwardReads = sapply(FWD_PRIMER.orients, primerHits, fn = file.path(raw_dir,file.list[i,2])), 
              FWD_PRIMER.ReverseReads = sapply(FWD_PRIMER.orients, primerHits, fn = file.path(raw_dir,file.list[i,3])), 
              REV_PRIMER.ForwardReads = sapply(REV_PRIMER.orients, primerHits, fn = file.path(raw_dir,file.list[i,2])), 
              REV_PRIMER.ReverseReads = sapply(REV_PRIMER.orients, primerHits, fn = file.path(raw_dir,file.list[i,3]))))
}

# ---------------- FILTER AND TRIM SEQUENCES
#Total Amplicon Length = 964 (Rev) - 528 (fwd) = 436 bp

#Overlap = 260 + 200 (trunc Len values) - 436 = 24 bp of overlap #[1]
#Overlap = 230 + 216 (trunc Len values) - 436 = 10 bp of overlap #[2]

#truncLen parameters estimated by inspecting quality profiles
#roughly when green line (phred score mean) drops below 30
#maxEE is considered to be a function of truncLen; it is easiest to set a factor/proportion between both values, e.g. 100

#[1]
filterOut1 <- filterAndTrim(fnFs.raw1,fnFs.qualFiltTrim1,
                           fnRs.raw1,fnRs.qualFiltTrim1, 
                           maxN=0,maxEE=c(2.3,2.2),truncLen=c(260,200),
                           verbose = TRUE, rm.phix = TRUE, 
                           compress = TRUE, multithread = TRUE)

filterOut1 # display table of sequences lost by quality filtering
summary(filterOut1[, 2]/filterOut1[, 1]) #check results, should retain > 70% -- OK

#[2]
filterOut2 <- filterAndTrim(fnFs.raw2,fnFs.qualFiltTrim2,
                           fnRs.raw2,fnRs.qualFiltTrim2, 
                           maxN=0,maxEE=c(3,3),truncLen=c(230,216),
                           verbose = TRUE, rm.phix = TRUE, 
                           compress = TRUE, multithread = TRUE)

filterOut2 # display table of sequences lost by quality filtering
summary(filterOut[, 2]/filterOut[, 1]) #check results, should retain > 70% -- OK
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.3100  0.7709  0.7884  0.7632  0.8035  0.8332 

# ---------------- NEW QUALITY PLOTS - AFTER TRIMMED/FILTERED

#[1]
#Forward
QualityProfileTFs1 <- list()

for(i in 1:length(fnFs.files1)) { #for each file on Fs_16S
  QualityProfileTFs1[[i]] <- list() #create a new empty list
  QualityProfileTFs1[[i]][[1]] <- plotQualityProfile(file.path(qualFiltTrim_dir1,fnFs.files1[i]))} #and plot a quality profile
pdf("Quality Profile Forward After Trim_TRAPS.pdf") #saving into .pdf

for(i in 1:length(fnFs.files1)) #organize the plots
{do.call("grid.arrange", QualityProfileTFs1[[i]])}
dev.off()
rm(QualityProfileTFs1)

#Reverse
QualityProfileTRs1 <- list()

for(i in 1:length(fnRs.files1)) {
  QualityProfileTRs1[[i]] <- list() #create a new empty list
  QualityProfileTRs1[[i]][[1]] <- plotQualityProfile(file.path(qualFiltTrim_dir1,fnRs.files1[i]))}
pdf("Quality Profile Reverse After Trim_TRAPS.pdf") #saving into .pdf

for(i in 1:length(fnRs.files1)) #organize the plots
{do.call("grid.arrange", QualityProfileTRs1[[i]])}
dev.off()
rm(QualityProfileTRs1)


#[2]
#Forward
QualityProfileTFs2 <- list()

for(i in 1:length(fnFs.files2)) { #for each file on Fs_16S
  QualityProfileTFs2[[i]] <- list() #create a new empty list
  QualityProfileTFs2[[i]][[1]] <- plotQualityProfile(file.path(qualFiltTrim_dir2,fnFs.files2[i]))} #and plot a quality profile
pdf("Quality Profile Forward After Trim_RAS.pdf") #saving into .pdf

for(i in 1:length(fnFs.files2)) #organize the plots
{do.call("grid.arrange", QualityProfileTFs2[[i]])}
dev.off()
rm(QualityProfileTFs2)

#Reverse
QualityProfileTRs2 <- list()

for(i in 1:length(fnRs.files2)) {
  QualityProfileTRs2[[i]] <- list() #create a new empty list
  QualityProfileTRs2[[i]][[1]] <- plotQualityProfile(file.path(qualFiltTrim_dir2,fnRs.files2[i]))}
pdf("Quality Profile Reverse After Trim_RAS.pdf") #saving into .pdf

for(i in 1:length(fnRs.files2)) #organize the plots
{do.call("grid.arrange", QualityProfileTRs2[[i]])}
dev.off()
rm(QualityProfileTRs2)

# ---------------- DEREPLICATION
# Dereplicate identical sequences to reduce computational costs
# Input: Files; Output: R-object

#[1]
fnFs.deRep1 <- derepFastq(fnFs.qualFiltTrim1)
fnRs.deRep1 <- derepFastq(fnRs.qualFiltTrim1)

#don't use file names but corresponding sample names as the names here will be used to create the final ASV table later
names(fnFs.deRep1) <- sample.names1
names(fnRs.deRep1) <- sample.names1

#[2]
fnFs.deRep2 <- derepFastq(fnFs.qualFiltTrim2)
fnRs.deRep2 <- derepFastq(fnRs.qualFiltTrim2)

#don't use file names but corresponding sample names as the names here will be used to create the final ASV table later
names(fnFs.deRep2) <- sample.names2
names(fnRs.deRep2) <- sample.names2

# ---------------- ERROR RATES
#LEARN DADA2 ERROR RATES

#[1]
errFWDs1 <- learnErrors(fnFs.deRep1, multithread=TRUE,randomize=TRUE, nbases = 1e8)
errREVs1 <- learnErrors(fnRs.deRep1, multithread=TRUE,randomize=TRUE, nbases = 1e8)

#[2]
errFWDs2 <- learnErrors(fnFs.deRep2, multithread=TRUE,randomize=TRUE, nbases = 1e8)
errREVs2 <- learnErrors(fnRs.deRep2, multithread=TRUE,randomize=TRUE, nbases = 1e8)

#PLOT ERROR PROFILE
#black curve should fit black dots and monotonically decreasing

#[1]
pdf("ErrorProfiles_TRAPS.pdf")
plotErrors(errFWDs1, nominalQ=TRUE)
plotErrors(errREVs1, nominalQ=TRUE)
dev.off()

#[2]
pdf("ErrorProfiles_RAS.pdf")
plotErrors(errFWDs2, nominalQ=TRUE)
plotErrors(errREVs2, nominalQ=TRUE)
dev.off()

#DADA SAMPLE INFERENCE
#pool="FALSE" if computation too costly; "pseudo" as compromise between computational costs and ASV detection sensitivity

#[1]
dadaFWDs1 <- dada(fnFs.deRep1, err=errFWDs1, multithread=TRUE, pool="pseudo")
dadaREVs1 <- dada(fnRs.deRep1, err=errREVs1, multithread=TRUE, pool="pseudo")

#inspect dada inference objects
dadaFWDs1
dadaREVs1


#[2]
dadaFWDs2 <- dada(fnFs.deRep2, err=errFWDs2, multithread=TRUE, pool="pseudo")
dadaREVs2 <- dada(fnRs.deRep2, err=errREVs2, multithread=TRUE, pool="pseudo")

#inspect dada inference objects
dadaFWDs2
dadaREVs2

# ---------------- MERGE PAIRED ENDS

#[1]
mergers1 <- mergePairs(dadaFWDs1, fnFs.deRep1, dadaREVs1, fnRs.deRep1, minOverlap=25, maxMismatch=0, verbose=TRUE)
# display mergers object; check nmatch (length of actual overlap) to check if remaining sequence lengths and selected overlap were sufficient
mergers1

#[2]
mergers2 <- mergePairs(dadaFWDs2, fnFs.deRep2, dadaREVs2, fnRs.deRep2, minOverlap=25, maxMismatch=0, verbose=TRUE)
# display mergers object; check nmatch (length of actual overlap) to check if remaining sequence lengths and selected overlap were sufficient
mergers2

# ---------------- CONSTRUCT, INSPECT AND SAVE SEQUENCE TABLE

#create seqtab
seqtab1 <- makeSequenceTable(mergers1)
seqtab2 <- makeSequenceTable(mergers2)

#number samples
dim(seqtab1)[1] 
dim(seqtab2)[1] #104

#number of ASVs
dim(seqtab1)[2] 
dim(seqtab2)[2] #19282

#---------------------- EXPORT

#[1] save sequence tabs 
saveRDS(seqtab1, file = paste("TRAP_EUK_seqtab",".rda",sep=""))
write.csv(t(seqtab1), paste("TRAP_EUK_seqtab",".csv",sep=""), quote=FALSE )

#[2] save sequence tabs
saveRDS(seqtab2, file = paste("RAS_EUK_seqtab",".rda",sep=""))
write.csv(t(seqtab2), paste("RAS_EUK_seqtab",".csv",sep=""), quote=FALSE )

#---------------------- TRACK READS
#Obtaining some info about the samples

#[1]
getN1 <- function(x) sum(getUniques(x))
track1 <- cbind(filterOut1, sapply(dadaFWDs1, getN1), sapply(dadaREVs1, getN1), sapply(mergers1, getN1))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track1) <- c("raw", "qualFiltered", "denoisedF", "denoisedR", "merged")
rownames(track1) <- sample.names1
track1

saveRDS(track1, file = paste("TRAP_EUK_trackedReads",".rda",sep=""))
write.csv(track1, paste("TRAP_EUK_trackedReads",".csv",sep=""), quote=FALSE )

#[2]
getN2 <- function(x) sum(getUniques(x))
track2 <- cbind(filterOut2, sapply(dadaFWDs2, getN2), sapply(dadaREVs2, getN2), sapply(mergers2, getN2))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track2) <- c("raw", "qualFiltered", "denoisedF", "denoisedR", "merged")
rownames(track2) <- sample.names2
track2

saveRDS(track2, file = paste("RAS_EUK_trackedReads",".rda",sep=""))
write.csv(track2, paste("RAS_EUK_trackedReads",".csv",sep=""), quote=FALSE )

#---------------------- 2nd script - MERGE THE DATASETS and ANNOTATE
# To create one ASV table from multiple MiSeq runs


#PREPARE ENVIRONMENT
#define database directory
db_dir <- file.path("/isibhv/projects/p_bioinf2/dbs/dada2_dbs")

# select database and corresponding taxonomic level pattern
taxDB <- "pr2_version_5.0.0_SSU_dada2.fasta"
taxLvs <- c("Domain", "Kingdom", "Supergroup","Division", "Class", "Order", "Family", "Genus", "Species")

# Load seqtabs from both datasets
sq1 <- readRDS("TRAP_EUK_seqtab.rda")
sq2 <- readRDS("RAS_EUK_seqtab.rda")

## Merge
seqtab <- mergeSequenceTables(
  sq1,sq2, repeats="error")


#---------------------- REMOVE CHIMERAS
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)

summary(rowSums(seqtab.nochim)/rowSums(seqtab))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9477  0.9835  0.9919  0.9882  0.9965  1.0000

#---------------------- REMOVE SINGLETONS AND JUNK SEQUENCES
#Do not remove singletons for the file taht will be used for alpha-div analyses

# Determine amplicon length/size range 

#Merged dataset
size <- table(rep(nchar(colnames(seqtab.nochim)), colSums(seqtab.nochim))) #[All]
size

# "c" adjusted to size range of amplicons, have a look in the output of the last command

seqtab.nochim.all <- seqtab.nochim[, nchar(
  colnames(seqtab.nochim)) %in% c(318:403) & 
    colSums(seqtab.nochim) > 1]

dim(seqtab.nochim.all) #147 samples, 12911 ASVs

summary(rowSums(seqtab.nochim.all))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#457   84923  112698  131926  168567  616290 

summary(rowSums(seqtab.nochim.all)/rowSums(seqtab.nochim))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.5493  0.9955  0.9993  0.9901  0.9998  1.0000

#[For alpha metrics]
seqtab.nochim.all.foralph <- seqtab.nochim[, nchar(
  colnames(seqtab.nochim)) %in% c(318:403)]

dim(seqtab.nochim.all.foralph) #147 samples, 13037 ASVs

summary(rowSums(seqtab.nochim.all.foralph)/rowSums(seqtab.nochim))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.5493  0.9955  0.9993  0.9902  0.9998  1.0000 

#----------- SAVE CLEAN ASVs TABLES

write.csv(seqtab.nochim.all, paste("EUK_seqtab_nochim_all", ".csv",sep=""), quote=FALSE )
saveRDS(seqtab.nochim.all, file = paste("EUK_seqtab_nochim_all",".rda",sep=""))

#---------------------- UPDATE TO THE TRACKED READS FILE

#add sequence numbers after chimera removal to trackedReads dataframe
tracked.reads.nochim <- cbind(track, rowSums(seqtab.nochim.all))

#add column headers
colnames(tracked.reads.nochim) <- c("raw", "qualfiltered", "denoisedF", "denoisedR", "merged", "nonchim")

#display read track table
(tracked.reads.nochim)

#save read track table
saveRDS(tracked.reads.nochim, file = paste("EUK_trackedReads_nochim",".rda",sep=""))
write.csv(tracked.reads.nochim, paste("EUK_trackedReads_nochim",".csv",sep=""), quote=FALSE )


#---------------------- TAXONOMY - PR2 V 5.0.0

#Assign taxonomy and calculate boothsrap values
taxa1 <- assignTaxonomy(seqtab.nochim.all, file.path(db_dir,taxDB), outputBootstraps=TRUE, minBoot=70, taxLevels = taxLvs, tryRC=TRUE, multithread=TRUE)
taxa2 <- assignTaxonomy(seqtab.nochim.all.foralph, file.path(db_dir,taxDB), outputBootstraps=TRUE, minBoot=70, taxLevels = taxLvs, tryRC=TRUE, multithread=TRUE)

# Extract taxonomic assignments
taxonomy_table <- taxa1$tax
taxonomy_table2 <- taxa2$tax


# Extract bootstrap values
bootstrap_values <- taxa1$boot
bootstrap_values2 <- taxa2$boot

#How many Domain
table(taxonomy_table[, 1]) #Eukaryota 12596 |EukNucl 14 |Eukplas 1
table(taxonomy_table2[, 1]) #Eukaryota 12721 |EukNucl 14 |Eukplas 1


#---------------------- Check NA on Kingdom level

sum(is.na(taxonomy_table[, 2])) #how many rows has NA values at Kingdom level: 1588

seqtab.nochim.allgood <- seqtab.nochim.all[,rownames(taxonomy_table)] #to be sure that taxa and asv table match
seqtab.alph.good <- seqtab.nochim.all.foralph[,rownames(taxonomy_table2)]


#---------------------- Format tables
seqtab.nochim.allprint <- seqtab.nochim.allgood
seqtab.alph.print <- seqtab.alph.good

all.equal(colnames(seqtab.nochim.allprint), rownames(taxonomy_table)) #TRUE
all.equal(colnames(seqtab.alph.print), rownames(taxonomy_table2)) #TRUE


#Changing the sequences to ASV+Number
#[1]
colnames(seqtab.nochim.allprint) <- paste(
  "asv", 1:ncol(seqtab.nochim.allgood), sep = "")

taxa_new <- taxonomy_table

rownames(taxa_new) <- paste("asv", 1:ncol(seqtab.nochim.allgood), sep = "")

#[2]
colnames(seqtab.alph.print) <- paste(
  "asv", 1:ncol(seqtab.alph.good), sep = "")

taxa_new2 <- taxonomy_table2

rownames(taxa_new2) <- paste("asv", 1:ncol(seqtab.alph.good), sep = "")

#--------------------------- EXPORT

#save taxonomy table
write.csv(taxa_new, paste("All_EUK_taxa", ".csv",sep=""), quote=FALSE )
saveRDS(taxa_new, file = paste("All_EUK_taxa",".rda",sep=""))

write.csv(taxa_new2, paste("TAX_EUK_FALPHA", ".csv",sep=""), quote=FALSE )
saveRDS(taxa_new2, file = paste("TAX_EUK_FALPHA",".rda",sep=""))

#save ASV table
write.csv(cbind(t(seqtab.nochim.allprint)), paste("All_EUK_asv", ".csv",sep=""), quote=FALSE )
saveRDS(cbind(t(seqtab.nochim.allprint)), file = paste("All_EUK_asv",".rda",sep=""))

write.csv(cbind(t(seqtab.alph.print)), paste("ASV_EUK_FALPHA", ".csv",sep=""), quote=FALSE )
saveRDS(cbind(t(seqtab.alph.print)), file = paste("ASV_EUK_FALPHA",".rda",sep=""))
