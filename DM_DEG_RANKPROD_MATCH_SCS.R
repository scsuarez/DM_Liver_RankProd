# load packages needed for probe conversion to genes
library("rat2302.db")
library(annotate)

# read output significant gene output from RankProd Back into R
UPSIG_CarTet_7d_1175 <- read.delim("RP_Sig_Results/CarbonTet_1175mg_7d_UPSig.txt", header=TRUE, sep="\t", row.names = NULL)
DOWNSIG_CarTet_7d_1175 <- read.delim("RP_Sig_Results/CarbonTet_1175mg_7d_DOWNSig.txt", header=TRUE, sep="\t", row.names = NULL)
UPSIG_CarTet_7d_400 <- read.delim("RP_Sig_Results/CarbonTet_400mg_7d_UPSig.txt", header=TRUE, sep="\t", row.names = NULL)
DOWNSIG_CarTet_7d_400 <- read.delim("RP_Sig_Results/CarbonTet_400mg_7d_DOWNSig.txt", header=TRUE, sep="\t", row.names = NULL)
UPSIG_Clotri_3d_89 <- read.delim("RP_Sig_Results/Clotrimazole_89mg_3d_UPSig.txt", header=TRUE, sep="\t", row.names = NULL)
DOWNSIG_Clotri_3d_89 <- read.delim("RP_Sig_Results/Clotrimazole_89mg_3d_DOWNSig.txt", header=TRUE, sep="\t", row.names = NULL)
UPSIG_Clotri_5d_178 <- read.delim("RP_Sig_Results/Clotrimazole_178mg_5d_UPSig.txt", header=TRUE, sep="\t", row.names = NULL)
DOWNSIG_Clotri_5d_178 <- read.delim("RP_Sig_Results/Clotrimazole_178mg_5d_DOWNSig.txt", header=TRUE, sep="\t", row.names = NULL)
UPSIG_Erlot_1d_58 <- read.delim("RP_Sig_Results/Erlotinib_58mg_1d_UPSig.txt", header=TRUE, sep="\t", row.names = NULL)
DOWNSIG_Erlot_1d_58 <- read.delim("RP_Sig_Results/Erlotinib_58mg_1d_DOWNSig.txt", header=TRUE, sep="\t", row.names = NULL)
UPSIG_Hydraz_5d_45 <- read.delim("RP_Sig_Results/Hydrazine_45mg_5d_UPSig.txt", header=TRUE, sep="\t", row.names = NULL)
DOWNSIG_Hydraz_5d_45 <- read.delim("RP_Sig_Results/Hydrazine_45mg_5d_DOWNSig.txt", header=TRUE, sep="\t", row.names = NULL)
UPSIG_Phenobarb_1d_54 <- read.delim("RP_Sig_Results/Phenobarbital_54mg_1d_UPSig.txt", header=TRUE, sep="\t", row.names = NULL)
DOWNSIG_Phenobarb_1d_54 <- read.delim("RP_Sig_Results/Phenobarbital_54mg_1d_DOWNSig.txt", header=TRUE, sep="\t", row.names = NULL)

# for an alternative strategy we did not end up employing
# n <- max(length(UPSIG_CarTet_7d_1175[,1]), length(UPSIG_CarTet_7d_400[,1]), length(UPSIG_Clotri_3d_89[,1]), length(UPSIG_Clotri_5d_178[,1]), length(UPSIG_Erlot_1d_58[,1]), length(UPSIG_Hydraz_5d_45[,1]), length(UPSIG_Phenobarb_1d_54[,1]))
# length(UPSIG_CarTet_7d_1175_Probes) <- n
# length(UPSIG_CarTet_7d_400_Probes) <- n
# length(UPSIG_Clotri_3d_89_Probes) <- n
# length(UPSIG_Clotri_5d_178_Probes) <- n
# length(UPSIG_Erlot_1d_58_Probes) <- n
# length(UPSIG_Hydraz_5d_45_Probes) <- n
# length(UPSIG_Phenobarb_1d_54_Probes) <- n
# UPSig_Probes <- cbind(UPSIG_CarTet_7d_1175_Probes, UPSIG_CarTet_7d_400_Probes, UPSIG_Clotri_3d_89_Probes, UPSIG_Clotri_5d_178_Probes, UPSIG_Erlot_1d_58_Probes, UPSIG_Hydraz_5d_45_Probes, UPSIG_Phenobarb_1d_54_Probes)

# get lists of probes for all significant up conditions
UPSIGProbes_CarTet_7d_1175 <- as.vector(UPSIG_CarTet_7d_1175[,1])
UPSIGProbes_CarTet_7d_400 <- as.vector(UPSIG_CarTet_7d_400[,1])
UPSIGProbes_Clotri_3d_89 <- as.vector(UPSIG_Clotri_3d_89[,1])
UPSIGProbes_Clotri_5d_178 <- as.vector(UPSIG_Clotri_5d_178[,1])
UPSIGProbes_Erlot_1d_58 <- as.vector(UPSIG_Erlot_1d_58[,1])
UPSIGProbes_Hydraz_5d_45 <- as.vector(UPSIG_Hydraz_5d_45[,1])
UPSIGProbes_Phenobarb_1d_54 <- as.vector(UPSIG_Phenobarb_1d_54[,1])

# list all significant up probe lists
UpSigProbes <- ls(pattern = "UPSIGProbes")
# combine all up probe lists into a single vector
all.UpSig.Probes <- do.call("c", mget(UpSigProbes))
# creates a contingency table that counts the occurances of each probe, sorts the probes by frequency, and outputs as a data frame
count.UpSig.Probes <- as.data.frame(sort(table(all.UpSig.Probes), decreasing = TRUE))

# convert probes to Rat Genes
PROBES<- as.character(row.names(count.UpSig.Probes))
# Use select and chip.db to extract ids (instead of get/mget).Note: one to many mapping occurs
AOUT <- select(rat2302.db, PROBES, c("SYMBOL","ENTREZID", "GENENAME"))
DUP_AOUT <- AOUT[duplicated(AOUT$PROBEID)|duplicated(AOUT$PROBEID,fromLast=TRUE),]
DUP_AOUT_PROBE <- DUP_AOUT$PROBEID
AOUT_UNIQ <- AOUT[!(AOUT$PROBEID %in% DUP_AOUT_PROBE),]
# in duplicated, remove LOC probes + duplicate copy of probes
# Note: LOC probes are probes that map to gene symb starting with LOC... (duplicate).Same probe also map to another gene symb
DUP_AOUT1<-DUP_AOUT[!grepl("LOC",DUP_AOUT$SYMBOL)&!duplicated(DUP_AOUT$PROBEID),]
# rbind DUP_AOUT1 to AOUT_UNIQ
PROB_RID_OUT <- rbind(AOUT_UNIQ,DUP_AOUT1)
# bind to final data frame with probes, frequency, gene symbol, entrezid, and gene name
UP_sig_counts <- cbind(PROB_RID_OUT[,1], count.UpSig.Probes, PROB_RID_OUT[,2:4])
rownames(UP_sig_counts) <- NULL
names(UP_sig_counts)[1] <- "UP_SIG_PROBES"
names(UP_sig_counts)[2] <- "FREQUENCY"

# write significant down genes with counts to delimited text file
write.table(UP_sig_counts ,file="RP_Sig_Results/DM_RP_UPsig_counts.txt",sep="\t",row.names=FALSE)

# get lists of probes for all significant down conditions
DOWNSIGProbes_CarTet_7d_1175 <- as.vector(DOWNSIG_CarTet_7d_1175[,1])
DOWNSIGProbes_CarTet_7d_400 <- as.vector(DOWNSIG_CarTet_7d_400[,1])
DOWNSIGProbes_Clotri_3d_89 <- as.vector(DOWNSIG_Clotri_3d_89[,1])
DOWNSIGProbes_Clotri_5d_178 <- as.vector(DOWNSIG_Clotri_5d_178[,1])
DOWNSIGProbes_Erlot_1d_58 <- as.vector(DOWNSIG_Erlot_1d_58[,1])
DOWNSIGProbes_Hydraz_5d_45 <- as.vector(DOWNSIG_Hydraz_5d_45[,1])
DOWNSIGProbes_Phenobarb_1d_54 <- as.vector(DOWNSIG_Phenobarb_1d_54[,1])

# list all down probe lists
DownSigProbes <- ls(pattern = "DOWNSIGProbes")
# combine all down probe lists into a single vector
all.DownSig.Probes <- do.call("c", mget(DownSigProbes)) 
# creates a contingency table that counts the occurances of each probe, sorts the probes by frequency, and outputs as a data frame
count.DownSig.Probes <- as.data.frame(sort(table(all.DownSig.Probes), decreasing = TRUE))

# convert probes to Rat Genes
PROBES<- as.character(row.names(count.DownSig.Probes))
# Use select and chip.db to extract ids (instead of get/mget).Note: one to many mapping occurs
AOUT <- select(rat2302.db, PROBES, c("SYMBOL","ENTREZID", "GENENAME"))
DUP_AOUT <- AOUT[duplicated(AOUT$PROBEID)|duplicated(AOUT$PROBEID,fromLast=TRUE),]
DUP_AOUT_PROBE <- DUP_AOUT$PROBEID
AOUT_UNIQ <- AOUT[!(AOUT$PROBEID %in% DUP_AOUT_PROBE),]
# in duplicated, remove LOC probes + duplicate copy of probes
# Note: LOC probes are probes that map to gene symb starting with LOC... (duplicate).Same probe also map to another gene symb
DUP_AOUT1<-DUP_AOUT[!grepl("LOC",DUP_AOUT$SYMBOL)&!duplicated(DUP_AOUT$PROBEID),]
# rbind DUP_AOUT1 to AOUT_UNIQ
PROB_RID_OUT <- rbind(AOUT_UNIQ,DUP_AOUT1)
# bind to final data frame with probes, frequency, gene symbol, entrezid, and gene name
DOWN_sig_counts <- cbind(PROB_RID_OUT[,1], count.DownSig.Probes, PROB_RID_OUT[,2:4])
rownames(DOWN_sig_counts) <- NULL
names(DOWN_sig_counts)[1] <- "DOWN_SIG_PROBES"
names(DOWN_sig_counts)[2] <- "FREQUENCY"

# write significant down genes with counts to delimited text file
write.table(DOWN_sig_counts ,file="RP_Sig_Results/DM_RP_DOWNSig_counts.txt",sep="\t",row.names=FALSE)

DM_UPSig_counts_liver <- read.delim("RP_Sig_Results/DM_RP_UPsig_counts.txt", header=TRUE, sep="\t", row.names = NULL)
DM_DOWNSig_counts_liver <- read.delim("RP_Sig_Results/DM_RP_DOWNSig_counts.txt", header=TRUE, sep="\t", row.names = NULL)