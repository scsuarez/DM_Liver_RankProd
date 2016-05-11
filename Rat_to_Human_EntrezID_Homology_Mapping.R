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

# use the ensembl database for human genes for one part of our query
human <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="useast.ensembl.org")
# use the ensembl database for rat genes for the first part of our query
rat<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="rnorvegicus_gene_ensembl", host="useast.ensembl.org")
# get information on 2 linked datasets - homology mapping
# attributes available by listAttributes()
# filters available by listFilters()
RHIDs_UPSIG <- getLDS(attributes=c("entrezgene"), filters="entrezgene", values= DM_UPSig_counts_liver$ENTREZID, mart=rat, attributesL=c("entrezgene"), martL=human, verbose = TRUE, uniqueRows = FALSE)
# set names of input and output columns
colnames(RHIDs_UPSIG) <- c("ENTREZ_GENE_ID_RAT", "ENTREZ_GENE_ID_HUMAN")
countsRH_UPSIG <- DM_UPSig_counts_liver[DM_UPSig_counts_liver$ENTREZID %in% RHIDs_UPSIG$ENTREZ_GENE_ID_RAT, ]
countsRH_UPSIG <- countsRH_UPSIG[order(countsRH_UPSIG$ENTREZID), ]
RHIDs_UPSIG <- RHIDs_UPSIG[order(RHIDs_UPSIG$ENTREZ_GENE_ID_RAT), ]
summary(countsRH_UPSIG$ENTREZID %in% RHIDs_UPSIG$ENTREZ_GENE_ID_RAT)
