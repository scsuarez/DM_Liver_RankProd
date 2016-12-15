# library relevant packages for differentially expressed gene analysis
library(Biobase)
library(genefilter)
library(RankProd)
library(affy)
library("rat2302.db")
library(annotate)

# Read in the tables that data depends on

#read in fold change matrix
# 645 Chemical Exposure Conditions,  Probes(genes)
# read in treatment to control mappings
mappings <- read.delim("Liver_Treated_to_Control_Mapping.txt", header=TRUE, sep="\t")
# load eset for rankprod analysis
load("eset_liver_pc2_nsFilterNV.RData")


# from MDAH, file - DM_LIVER_DEG_RANKPROD_2.R
# get DEGs (DoD BHSAI)
 
RP_genes <- function(eset_name, TREATMENT_CHEMICAL, TREATMENT_DOSE, TREATMENT_DURATION, TREATMENT_VEHICLE, TREATMENT_ROUTE, filename_up, filename_down){
    eset <- eset_name
    pdata <- pData(eset)
    pdataSorted <- pdata[order(pdata$CHEMICAL, pdata$DOSE, pdata$DURATION, pdata$VEHICLE, pdata$ROUTE), ]
    sampleNamesSorted <- row.names(pdataSorted)
    esetSorted <- eset[,sampleNamesSorted]
    eset <- esetSorted
    pdata <- pdataSorted
    #^^^
    num_replicates_treatments <- numeric(0)
    num_replicates_controls <- numeric(0)
    #^^^
    pdata_c <- pdata[pdata$CHEMICAL == "", ] # controls
    pdata_t <- pdata[pdata$CHEMICAL != "", ] # treatments
    #^^
    # mappings <- read.delim("2_liver_treated_to_control_detailed_mapping.txt", header=TRUE, sep="\t", row.names = NULL)
    # this is MDAH Mappings file
    # Use ours as above
    
    treatments <- unique(pdata_t[, c("CHEMICAL", "DOSE", "DURATION", "VEHICLE", "ROUTE")])
    # treatments[treatments$CHEMICAL=="CARBON TETRACHLORIDE",]
    #^^
    pdata_t_current <- pdata_t[pdata_t$CHEMICAL== TREATMENT_CHEMICAL & pdata_t$DOSE == TREATMENT_DOSE & pdata_t$DURATION == TREATMENT_DURATION & pdata_t$VEHICLE== TREATMENT_VEHICLE & pdata_t$ROUTE== TREATMENT_ROUTE,]
    # enter data for TREATMENT_i, ie data for treatment of interest, in this case from PathResultp05
    treatmentNames <- row.names(pdata_t_current) 
    c.mappings <- mappings[mappings$TREATMENT_ARRAY_ID %in% treatmentNames,] # Map controls
    controlNames <- as.character(unique(c.mappings$CONTROL_ARRAY_ID))
    controlNames <- controlNames[controlNames %in% row.names(pdata_c)] # make sure these are in eset
    eset_i <- eset[,c(controlNames, treatmentNames)] # eset with one treatment and its controls
    # assign treatment (1) and control (0) classes
    m.c1 <- rep(c(0,1), c(length(controlNames), length(treatmentNames))) 
    m.RP.out <- RP(exprs(eset_i), m.c1, logged=TRUE) # compute DEG
    # plotRP(m.RP.out, cutoff=0.05) # cutoff for pfp which is same as FDR
    m.topgenes <- topGene(m.RP.out,cutoff=1.0,method="pval",logged=TRUE,logbase=2,gene.names=row.names(exprs(eset_i)))
    RP_P2M_CLO_UP <- m.topgenes$Table1
    RP_P2M_CLO_DOWN <- m.topgenes$Table2
    
    # UP
    RP_UP1<-as.data.frame(RP_P2M_CLO_UP[,"pfp",drop=FALSE])
    colnames(RP_UP1)<-"pfp"
    RP_UP1sig<-subset(RP_UP1,pfp<0.05)
    # DOWN
    RP_DOWN1<-as.data.frame(RP_P2M_CLO_DOWN[,"pfp",drop=FALSE])
    colnames(RP_DOWN1)<-"pfp"
    RP_DOWN1sig<-subset(RP_DOWN1,pfp<0.05)
    
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    PROBES<- as.character(rownames(RP_UP1sig))
    # Use select and chip.db to extract ids (instead of get/mget).Note: one to many mapping occurs
    AOUT <- select(rat2302.db, PROBES, c("SYMBOL","ENTREZID", "GENENAME"))
    DUP_AOUT <- AOUT[duplicated(AOUT$PROBEID)|duplicated(AOUT$PROBEID,fromLast=TRUE),]
    DUP_AOUT_PROBE <- DUP_AOUT$PROBEID
    AOUT_UNIQ <- AOUT[!(AOUT$PROBEID %in% DUP_AOUT_PROBE),]
    # in duplicated, remove LOC probes + duplicate copy of probes
    # Note: LOC probes are probes that map to gene symb starting with LOC... (duplicate).Same probe also map to another gene symb
    DUP_AOUT1<-DUP_AOUT[!grepl("LOC",DUP_AOUT$SYMBOL)&!duplicated(DUP_AOUT$PROBEID),]
    #rbind DUP_AOUT1 to AOUT_UNIQ
    PROB_RID_OUT <- rbind(AOUT_UNIQ,DUP_AOUT1)
    UP1_sig <- cbind(PROB_RID_OUT, RP_UP1sig) 
    
    
    PROBES<- as.character(rownames(RP_DOWN1sig))
    # Use select and chip.db to extract ids (instead of get/mget).Note: one to many mapping occurs
    AOUT <- select(rat2302.db, PROBES, c("SYMBOL","ENTREZID", "GENENAME"))
    DUP_AOUT <- AOUT[duplicated(AOUT$PROBEID)|duplicated(AOUT$PROBEID,fromLast=TRUE),]
    DUP_AOUT_PROBE <- DUP_AOUT$PROBEID
    AOUT_UNIQ <- AOUT[!(AOUT$PROBEID %in% DUP_AOUT_PROBE),]
    # in duplicated, remove LOC probes + duplicate copy of probes
    # Note: LOC probes are probes that map to gene symb starting with LOC... (duplicate).Same probe also map to another gene symb
    DUP_AOUT1<-DUP_AOUT[!grepl("LOC",DUP_AOUT$SYMBOL)&!duplicated(DUP_AOUT$PROBEID),]
    #rbind DUP_AOUT1 to AOUT_UNIQ
    PROB_RID_OUT <- rbind(AOUT_UNIQ,DUP_AOUT1)
    DOWN1_sig <- cbind(PROB_RID_OUT, RP_DOWN1sig)
    
    # write sig up and down as delimted text
    write.table(UP1_sig,file=filename_up, sep="\t",row.names=FALSE)
    write.table(DOWN1_sig,file=filename_down, sep="\t",row.names=FALSE)
}

# for all 9 conditions of interest
# CARBON TETRACHLORIDE  1175 mg/kg  CORN OIL 100 %  ORAL GAVAGE	7 d
TREATMENT_CHEMICAL <- "CARBON TETRACHLORIDE"
TREATMENT_DOSE <- "1175 mg/kg"
TREATMENT_DURATION <- "7 d"
TREATMENT_VEHICLE <- "CORN OIL 100 %"
TREATMENT_ROUTE <- "ORAL GAVAGE"
filename_up <- "RP_Sig_Results/CarbonTet_1175mg_7d_UPSig.txt"
filename_down <- "RP_Sig_Results/CarbonTet_1175mg_7d_DOWNSig.txt"

# CARBON TETRACHLORIDE  400 mg/kg	CORN OIL 100 %	ORAL GAVAGE	7 d
TREATMENT_CHEMICAL <- "CARBON TETRACHLORIDE"
TREATMENT_DOSE <- "400 mg/kg"
TREATMENT_DURATION <- "7 d"
TREATMENT_VEHICLE <- "CORN OIL 100 %"
TREATMENT_ROUTE <- "ORAL GAVAGE"
filename_up <- "RP_Sig_Results/CarbonTet_400mg_7d_UPSig.txt"
filename_down <- "RP_Sig_Results/CarbonTet_400mg_7d_DOWNSig.txt"

# CLOTRIMAZOLE  178 mg/kg	CORN OIL 100 %	ORAL GAVAGE	3 d
# all arrays for this chemical exposure were judged outliers by ArrayQualityMetrics
## TREATMENT_CHEMICAL <- "CLOTRIMAZOLE"
## TREATMENT_DOSE <- "178 mg/kg"
## TREATMENT_DURATION <- "3 d"
## TREATMENT_VEHICLE <- "CORN OIL 100 %"
## TREATMENT_ROUTE <- "ORAL GAVAGE"
## filename_up <- "RP_Sig_Results/Clotrimazole_178mg_3d_UPSig.txt"
## filename_down <- "RP_Sig_Results/Clotrimazole_178mg_3d_UPSig.txt"

# CLOTRIMAZOLE  178 mg/kg	CORN OIL 100 %	ORAL GAVAGE	5 d
TREATMENT_CHEMICAL <- "CLOTRIMAZOLE"
TREATMENT_DOSE <- "178 mg/kg"
TREATMENT_DURATION <- "5 d"
TREATMENT_VEHICLE <- "CORN OIL 100 %"
TREATMENT_ROUTE <- "ORAL GAVAGE"
filename_up <- "RP_Sig_Results/Clotrimazole_178mg_5d_UPSig.txt"
filename_down <- "RP_Sig_Results/Clotrimazole_178mg_5d_DOWNSig.txt"

# CLOTRIMAZOLE  89 mg/kg	CORN OIL 100 %	ORAL GAVAGE	3 d
TREATMENT_CHEMICAL <- "CLOTRIMAZOLE"
TREATMENT_DOSE <- "89 mg/kg"
TREATMENT_DURATION <- "3 d"
TREATMENT_VEHICLE <- "CORN OIL 100 %"
TREATMENT_ROUTE <- "ORAL GAVAGE"
filename_up <- "RP_Sig_Results/Clotrimazole_89mg_3d_UPSig.txt"
filename_down <- "RP_Sig_Results/Clotrimazole_89mg_3d_DOWNSig.txt"

# ERLOTINIB  58 mg/kg	WATER 100 %	ORAL GAVAGE	1 d
TREATMENT_CHEMICAL <- "ERLOTINIB"
TREATMENT_DOSE <- "58 mg/kg"
TREATMENT_DURATION <- "1 d"
TREATMENT_VEHICLE <- "WATER 100 %"
TREATMENT_ROUTE <- "ORAL GAVAGE"
filename_up <- "RP_Sig_Results/Erlotinib_58mg_1d_UPSig.txt"
filename_down <- "RP_Sig_Results/Erlotinib_58mg_1d_DOWNSig.txt"

# HYDRAZINE  45 mg/kg	WATER 100 %	ORAL GAVAGE	5 d
TREATMENT_CHEMICAL <- "HYDRAZINE"
TREATMENT_DOSE <- "45 mg/kg"
TREATMENT_DURATION <- "5 d"
TREATMENT_VEHICLE <- "WATER 100 %"
TREATMENT_ROUTE <- "ORAL GAVAGE"
filename_up <- "RP_Sig_Results/Hydrazine_45mg_5d_UPSig.txt"
filename_down <- "RP_Sig_Results/Hydrazine_45mg_5d_DOWNSig.txt"

# MICONAZOLE  920 mg/kg	CORN OIL 100 %	ORAL GAVAGE	5 d
# all arrays for this chemical exposure were judged outliers by ArrayQualityMetrics
## TREATMENT_CHEMICAL <- "MICONAZOLE"
## TREATMENT_DOSE <- "920 mg/kg"
## TREATMENT_DURATION <- "5 d"
## TREATMENT_VEHICLE <- "CORN OIL 100 %"
## TREATMENT_ROUTE <- "ORAL GAVAGE"
## filename_up <- "RP_Sig_Results/Miconozole_920mg_5d_UPSig.txt"
## filename_down <- "RP_Sig_Results/Miconozole_920mg_5d_DOWNSig.txt"

# PHENOBARBITAL  54 mg/kg	WATER 100 %	ORAL GAVAGE	1 d
TREATMENT_CHEMICAL <- "PHENOBARBITAL"
TREATMENT_DOSE <- "54 mg/kg"
TREATMENT_DURATION <- "1 d"
TREATMENT_VEHICLE <- "WATER 100 %"
TREATMENT_ROUTE <- "ORAL GAVAGE"
filename_up <- "RP_Sig_Results/Phenobarbital_54mg_1d_UPSig.txt"
filename_down <- "RP_Sig_Results/Phenobarbital_54mg_1d_DOWNSig.txt"

RP_genes(eset_liver_pc.2_nsFilterNV_rma_qc, TREATMENT_CHEMICAL, TREATMENT_DOSE, TREATMENT_DURATION, TREATMENT_VEHICLE, TREATMENT_ROUTE, filename_up, filename_down)
