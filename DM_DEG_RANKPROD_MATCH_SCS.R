# load packages needed for probe conversion to genes
library("rat2302.db")
library(annotate)
# file for conversion of Rat entrez gene id into human entrez gene id
# utilizes Bioconductor package biomaRt
library("biomaRt")
# utilizes complicated call structure, but this is the order you do it in
# generic - can be used for any of the list functions of biomaRt
# ensembl_us_east <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="useast.ensembl.org")
# if you want to use any of the functions of biomaRt you need to know, 
## the mart (database) you want to query
### obtained by listMarts(host="useast.ensembl.org")
## the dataset 
### available using listDatasets()
## the host
### default host ="www.ensembl.org" is no longer functional
### use a mirror like host="useast.ensembl.org"
# for joining and df maniulation
library(plyr)
# for human gene annotation
library("org.Hs.eg.db")
# making database simplification easier
library(sqldf)

# use the ensembl database for human genes for one part of our query
human <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="useast.ensembl.org")
# use the ensembl database for rat genes for the first part of our query
rat<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="rnorvegicus_gene_ensembl", host="useast.ensembl.org")

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
# Use select and chip.db to extract ids (instead of get/mget).Note: one to many mapping occurs if we had not filtered for entrezid duplicates
AOUT <- select(rat2302.db, PROBES, "ENTREZID")
count.UPSIG.entrez <- cbind(AOUT[1], count.UpSig.Probes, AOUT[2])
colnames(count.UPSIG.entrez) <- c("PROBEID", "FREQ", "ENTREZ_GENE_ID_RAT")
# get information on 2 linked datasets - homology mapping
# attributes available by listAttributes()
# filters available by listFilters()
RHIDs_UPSIG <- getLDS(attributes=c("entrezgene"), filters="entrezgene", values= AOUT$ENTREZID, mart=rat, attributesL=c("entrezgene"), martL=human, verbose = TRUE, uniqueRows = FALSE)
# set names of input and output columns
colnames(RHIDs_UPSIG) <- c("ENTREZ_GENE_ID_RAT", "ENTREZ_GENE_ID_HUMAN")
# join the 2 data frames including all values in count.UPSIG.entrez and matching values in RHIDs_UPSIG
count.UPSIG.RH <- join(count.UPSIG.entrez, RHIDs_UPSIG, by = "ENTREZ_GENE_ID_RAT", type = "left")
count.UPSIG.RH[ ,4] <- as.character(count.UPSIG.RH[ ,4])
# get gene annotation information
h.UPSIG.entrezid <- as.character(count.UPSIG.RH$ENTREZ_GENE_ID_HUMAN[!is.na(count.UPSIG.RH$ENTREZ_GENE_ID_HUMAN)])
h.UPSIG.entrezsymbol  <- as.character(unlist(mget(h.UPSIG.entrezid, envir=org.Hs.egSYMBOL, ifnotfound=NA)))
h.UPSIG.entrezgene  <- as.character(unlist(mget(h.UPSIG.entrezid, envir=org.Hs.egGENENAME, ifnotfound=NA)))

# don't use this now, but could be useful later for KEGG Pathway data
# h.UPSIG.entrezPATH  <- as.vector(unlist(mget(h.UPSIG.entrezid, envir=org.Hs.egPATH, ifnotfound=NA)))  

# bind entrez id, symbol, and gene information into 1 df
h.UPSIG <- as.data.frame(cbind(h.UPSIG.entrezid, h.UPSIG.entrezsymbol, h.UPSIG.entrezgene))
# set colnames of df so it can be joined
colnames(h.UPSIG) <- c("ENTREZ_GENE_ID_HUMAN", "SYMBOL", "GENENAME")
# set column classes to character for join
for(i in 1:dim(h.UPSIG)[2]){
  h.UPSIG[ ,i] <- as.character(h.UPSIG[ ,i])
}
# join df with human entrez information to our table with probes, frequency, rat entrezid
count.UPSIG.RHcomplete <- join(count.UPSIG.RH, h.UPSIG, by = "ENTREZ_GENE_ID_HUMAN", type = "left")
# use SQL and sqldf to subset for distinct cases (no repeated information, which happens for some reason with join)
count.UPSIG.RHcomplete <- sqldf('SELECT DISTINCT * FROM [count.UPSIG.RHcomplete]')

# write significant up genes with counts to delimited text file
write.table(count.UPSIG.RHcomplete ,file = "RP_Sig_Analysis/DM_RP_UPsig_entrez_counts.txt",sep="\t",row.names=FALSE)

# write probe and human entrez id to file for DAVID web analysis
write(count.UPSIG.RHcomplete$PROBEID, file = "RP_Sig_Analysis/DM_RP_UPsig_probes.txt")
write(h.UPSIG.entrezid, file = "RP_Sig_Analysis/DM_RP_UPsig_hentrezid.txt")

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
# Use select and chip.db to extract ids (instead of get/mget).Note: one to many mapping occurs if we had not filtered for entrezid duplicates
BOUT <- select(rat2302.db, PROBES, "ENTREZID")
count.DOWNSIG.entrez <- cbind(BOUT[1], count.DownSig.Probes, BOUT[2])
colnames(count.DOWNSIG.entrez) <- c("PROBEID", "FREQ", "ENTREZ_GENE_ID_RAT")
# get information on 2 linked datasets - homology mapping
# attributes available by listAttributes()
# filters available by listFilters()
RHIDs_DOWNSIG <- getLDS(attributes=c("entrezgene"), filters="entrezgene", values= BOUT$ENTREZID, mart=rat, attributesL=c("entrezgene"), martL=human, verbose = TRUE, uniqueRows = FALSE)
# set names of input and output columns
colnames(RHIDs_DOWNSIG) <- c("ENTREZ_GENE_ID_RAT", "ENTREZ_GENE_ID_HUMAN")
# join the 2 data frames including all values in count.UPSIG.entrez and matching values in RHIDs_UPSIG
count.DOWNSIG.RH <- join(count.DOWNSIG.entrez, RHIDs_DOWNSIG, by = "ENTREZ_GENE_ID_RAT", type = "left")
count.DOWNSIG.RH[ ,4] <- as.character(count.DOWNSIG.RH[ ,4])
# get gene annotation information
h.DOWNSIG.entrezid <- as.character(count.DOWNSIG.RH$ENTREZ_GENE_ID_HUMAN[!is.na(count.DOWNSIG.RH$ENTREZ_GENE_ID_HUMAN)])
h.DOWNSIG.entrezsymbol  <- as.character(unlist(mget(h.DOWNSIG.entrezid, envir=org.Hs.egSYMBOL, ifnotfound=NA)))
h.DOWNSIG.entrezgene  <- as.character(unlist(mget(h.DOWNSIG.entrezid, envir=org.Hs.egGENENAME, ifnotfound=NA)))

# bind entrez id, symbol, and gene information into 1 df
h.DOWNSIG <- as.data.frame(cbind(h.DOWNSIG.entrezid, h.DOWNSIG.entrezsymbol, h.DOWNSIG.entrezgene))
# set colnames of df so it can be joined
colnames(h.DOWNSIG) <- c("ENTREZ_GENE_ID_HUMAN", "SYMBOL", "GENENAME")
# set column classes to character for join
for(i in 1:dim(h.DOWNSIG)[2]){
  h.DOWNSIG[ ,i] <- as.character(h.DOWNSIG[ ,i])
}
# join df with human entrez information to our table with probes, frequency, rat entrezid
count.DOWNSIG.RHcomplete <- join(count.DOWNSIG.RH, h.DOWNSIG, by = "ENTREZ_GENE_ID_HUMAN", type = "left")
# use SQL and sqldf to subset for distinct cases (no repeated information, which happens for some reason with join)
count.DOWNSIG.RHcomplete <- sqldf('SELECT DISTINCT * FROM [count.DOWNSIG.RHcomplete]')

# write significant up genes with counts to delimited text file
write.table(count.DOWNSIG.RHcomplete ,file = "RP_Sig_Analysis/DM_RP_DOWNsig_entrez_counts.txt",sep="\t",row.names=FALSE)

# write probe and human entrez id to file for DAVID web analysis
write(count.DOWNSIG.RHcomplete$PROBEID, file = "RP_Sig_Analysis/DM_RP_DOWNsig_probes.txt")
write(h.DOWNSIG.entrezid, file = "RP_Sig_Analysis/DM_RP_DOWNsig_hentrezid.txt")