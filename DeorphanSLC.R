# 1 Before you start ####

# --------------------------- #
#     Set up environment      #
# --------------------------- #

# Create an empty workspace
rm(list = ls(all = TRUE))

# Create a new folder under current directory
dir.create("~/DeorphanSLC_manuscript")
setwd("~/DeorphanSLC_manuscript")

# Create subfolders
dir.create("Raw Data")
dir.create("Results")

# --------------------------- #
#      Download raw data      #
# --------------------------- #

## 1.1 Acquire CCLE 2019 omics ####
# 1: Go to https://depmap.org/portal/data_page/?tab=allData
# 2: In the drop-down list under "Select a file set to view", select "CCLE 2019"
# 3: Under Primary files, download "CCLE_metabolomics_20190502.csv" (CCLE2019 transcriptomics);
#    Under Supplemental Files, download "CCLE_RNAseq_genes_counts_20180929.gct.gz" (CCLE2019 metabolomics);
#    Under Supplemental Files, download "Cell_lines_annotations_20181226.txt"
# 4: Single click the .gct.gz file to release .gct
# 5: Drag the files downloaded under "~/DeorphanSLC_manuscript/Raw Data/" directory

## 1.2 Acquire NCI60 omics ####
# 1: To prepare NCI60 transcriptomics, lease refer to Methods provided in Perez & Sarkies (2023) (https://doi.org/10.1371/journal.pbio.3002354)
#    Specifically: Download raw RNA sequencing reads using “SRAtoolkit”, SRA identifier numbers available in table S5 of the above paper
# 2: Calculate read counts per gene by counting reads overlapping exons, here I cite:

biocmanager_packages = c("DESeq2",
                         "GenomicRanges",
                         "GenomicAlignments",
                         "MetaboSignal",
                         "KEGGREST",
                         "TxDb.Hsapiens.UCSC.hg38.knownGene",
                         "org.Hs.es.db")
lapply(biocmanager_packages, function(x) {
    if (!require(x, character.only = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)){
        install.packages("BiocManager")
      }
      BiocManager::install(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

exbygene = exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")
# point to directory with aligned sequencing bam files
thebams = list.files("insert_bam_directory", full = TRUE)

olap = summarizeOverlaps(exbygene, thebams)

# 3: Normalise the read counts with MRN normalisation, and save as .rds file, here I cite:
deseq = DESeqDataSet(olap, design= ~ 1)

deseq = estimateSizeFactors(deseq)
NCI60_normalised_counts = counts(deseq, normalized = TRUE)

NCI60RNA_metadata <- read.table("input/NCI60_RNAseq_SRA.csv", sep = ",", header = TRUE)
colnames(NCI60_normalised_counts) = NCI60RNA_metadata[match(str_remove(colnames(NCI60_normalised_counts), pattern = "\\.bam"), NCI60RNA_metadata$Run), "LibraryName"]

names_changed_to_match_format = str_remove(str_remove(str_replace(str_replace(str_remove(colnames(NCI60_normalised_counts), "^[A-Z]{2,3}:"), "/", "_"), " ", "_"), "\\(TB\\)"), "_ATCC")
names_changed_to_match_format = str_replace(names_changed_to_match_format, pattern = "MDA-N", replacement = "MDA-MB-468")

colnames(NCI60_normalised_counts) = names_changed_to_match_format

saveRDS(NCI60_normalised_counts, "Raw Data/NCI60_RNAseq_normalised_counts.rds")

# 4: To prepare NCI60 metabolomics, go to https://dtp.cancer.gov/, in the drop-down list under "Databases & Tools", select "DTP data search tools"
# 5: Under "Search Data of Molecular Characterization in NCI60 Tumor Cell Lines", select "Download NCI60 molecular target characterization data"
# 6: Under "Metabolomic Data From Metabolon- data averaged from triplicate experiments" subtitle, download the file "WEB_DATA_METABOLON.ZIP"
# 7: Unzip the "WEB_DATA_METABOLON.ZIP" file, drag the .csv file under "~/DeorphanSLC_manuscript/Raw Data/" directory

## 1.3 Acquire CCL180 omics ####
# 1: Go to https://doi.org/10.3929/ethz-b-000511784, download the file "primary analysis (metabolomics) (ZIP, 27.81Mb)"
# 2: Unzip the file and drag all the file contained under "~/DeorphanSLC_manuscript/Raw Data/" directory

## 1.4 Acquire CRISPR-Cas9 Gene dependency data ####
# 1: Go to https://depmap.org/portal/data_page/?tab=allData
# 2: In the drop-down list under "Select a file set to view" AND after "Version:", select "DepMap Public 23Q2 omics"
# 3: Under Primary files, download "CRISPRGeneEffect.csv" and "Model.csv"
# 4: Drag the "CRISPRGeneEffect.csv" and "Model.csv" under "~/DeorphanSLC_manuscript/Raw Data/" directory

## 1.5 Acquire PRISM drug repurposing screen data ####
# 1: Go to https://depmap.org/repurposing/
# 2: Under "Secondary screen" subtitle, download "secondary-screen-cell-line-info.csv";
#                                                "secondary-screen-replicate-collapsed-treatment-info.csv";
#                                                "secondary-screen-replicate-collapsed-logfold-change.csv"
# 3: Drag files downloaded under "~/DeorphanSLC_manuscript/Raw Data/" directory

## 1.6 Download supplementary tables  ####
# 1: Download Table S1 to your local computer, save the sheets with name "SLC transcripts" as a separate .csv file ("SLC transcripts.csv");
#                                              save the sheets with name "SLC annotation 2023" as a separate .csv file ("SLC annotation 2023.csv");
#             Table S5 to your local computer, save as a .csv file with the name "Known interaction" ("Known interaction.csv");
#             Table S6 to your local computer, save as a .csv file with the name "Translation book" ("Translation book.csv");
#             Table S8 to your local computer, save as a .csv file with the name "Adjacency matrix" ("Adjacency matrix.csv");
#             Table S9 to your local computer, save as a .csv file with the name "Pathway" ("Pathway.csv");
#             Table S11 to your local computer, save as a .csv file with the name "Known drug target" ("Known drug target.csv");
#             Table S12 to your local computer, save as a .csv file with the name "Filtering thresholds" ("Filtering thresholds.csv");
#             Table S13 to your local computer, save as a .csv file with the name "SLC-drug interaction predictions" ("SLC-drug interaction predictions.csv")
# 2: Drag these tables under "~/DeorphanSLC_manuscript/Raw Data/" directory

# 2 Normalisation, correlation and transformation ####

packages <- c("data.table",
              "tidyverse",
              "biomaRt", 
              "CePa",
              "progress",
              "imager",
              "foreach",
              "doParallel"
)
lapply(packages, function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Read SLC transcript annotation
SLC_annotation_2024 = fread("Raw Data/SLC transcripts.csv", select = c(1:6))

## 2.1 CCLE2019 omics ####
### §1 read files ####
CCLE2019_transcriptomics = read.gct(file="Raw Data/CCLE_RNAseq_genes_counts_20180929.gct")
CCLE2019_metabolomics = fread("Raw Data/CCLE_metabolomics_20190502.csv")

### §2 normalise transcript counts with MRN #####
# Create dds object in DESeq2
colData = data.frame(condition= colnames(CCLE2019_transcriptomics), row.names = colnames(CCLE2019_transcriptomics))
dds = DESeqDataSetFromMatrix(countData = CCLE2019_transcriptomics,
                              colData = colData,
                              design = ~ condition)

# Quality control: filtering rows with less than 10 reads
dds = dds[rowSums(counts(dds)) >= 10,]

# Normalisation with MRN
dds = estimateSizeFactors(dds)
CCLE2019_transcriptomics_normalised = counts(dds, normalized = TRUE)

# cleaning
rm(colData)
rm(dds)

### §3 format normalised transcript counts and metabolomics #####
# Edit cell line name with wrong format
CCLE2019_transcriptomics_normalised = as.data.frame(CCLE2019_transcriptomics_normalised)
colnames(CCLE2019_transcriptomics_normalised)[1:15] = substr(
  x = colnames(CCLE2019_transcriptomics_normalised)[1:15],
  start = 2, stop = nchar(colnames(CCLE2019_transcriptomics_normalised)[1:15])
)
# check cell lines specific to data frames (transcript counts and metabolomics)
print(paste0(
  length(colnames(CCLE2019_transcriptomics_normalised) [!(colnames(CCLE2019_transcriptomics_normalised) %in% CCLE2019_metabolomics$CCLE_ID)]),
  " cell lines are specific to CCLE2019 transcriptomics."
))

print(paste0(
  length(CCLE2019_metabolomics$CCLE_ID [!(CCLE2019_metabolomics$CCLE_ID %in% colnames(CCLE2019_transcriptomics_normalised))]),
  " cell lines are specific to CCLE2019 metabolomics."
))

# remove cell lines specific to data frames, 913 cell lines overlapping between transcript counts and metabolomics
CCLE2019_transcriptomics_normalised = CCLE2019_transcriptomics_normalised[,!colnames(CCLE2019_transcriptomics_normalised) %in%
                                                 colnames(CCLE2019_transcriptomics_normalised) [!(colnames(CCLE2019_transcriptomics_normalised) %in% CCLE2019_metabolomics$CCLE_ID)]]
CCLE2019_metabolomics = CCLE2019_metabolomics[!CCLE2019_metabolomics$CCLE_ID %in% 
                                                 CCLE2019_metabolomics$CCLE_ID [!(CCLE2019_metabolomics$CCLE_ID %in% colnames(CCLE2019_transcriptomics_normalised))],]
# match the order of cell lines in transcript counts and metabolomics
CCLE2019_metabolomics = t(CCLE2019_metabolomics)
colnames(CCLE2019_metabolomics) = CCLE2019_metabolomics[1,] # cell line name into column name
CCLE2019_metabolomics = CCLE2019_metabolomics[-c(1,2),]
reorder_idx = match(colnames(CCLE2019_transcriptomics_normalised), colnames(CCLE2019_metabolomics))
CCLE2019_metabolomics = CCLE2019_metabolomics[,reorder_idx]
all(colnames(CCLE2019_transcriptomics_normalised) == colnames(CCLE2019_metabolomics)) # the cell line order is aligned between transcript counts and metabolomics
CCLE2019_metabolomics <- data.table(Metabolite = rownames(CCLE2019_metabolomics), CCLE2019_metabolomics)

# cleaning
rm(reorder_idx)

### §4 translate Ensembl name to gene symbol #####
# translation using biomaRt
transcript_name = sub('\\.[0-9]*$', '', rownames(CCLE2019_transcriptomics_normalised)) # remove the version info from Ensembl name
rownames(CCLE2019_transcriptomics_normalised) = transcript_name
mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gene_name = getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                  filters = "ensembl_gene_id", values = transcript_name,
                  mart = mart)
# remove transcripts failed to translate
CCLE2019_transcriptomics_normalised = CCLE2019_transcriptomics_normalised[!rownames(CCLE2019_transcriptomics_normalised) %in% transcript_name [!(transcript_name %in% gene_name$ensembl_gene_id)],]

# returned 48431 rows in transcript counts, but gene name has 48432 rows
# the naming of one Ensembl name has changed across version and not updated: ENSG00000230417 can be translated to "LINC00595"(.12) and "LINC00856"(.6).
# Here "LINC00595" will be preserved.
gene_name = gene_name[-which(gene_name$hgnc_symbol == "LINC00856"),]

# attach translated gene symbols to the transcript counts
SYMBOL = c()
for (i in 1:nrow(CCLE2019_transcriptomics_normalised)){
  ensembl_name = rownames(CCLE2019_transcriptomics_normalised)[i]
  symbol_name = gene_name$hgnc_symbol[which(gene_name$ensembl_gene_id == ensembl_name)]
  SYMBOL = c(SYMBOL, symbol_name)
  print(i)
}
CCLE2019_transcriptomics_normalised$SYMBOL = SYMBOL
CCLE2019_transcriptomics_normalised = CCLE2019_transcriptomics_normalised %>% dplyr::select(SYMBOL, everything()) # move symbol column to the first column
CCLE2019_transcriptomics_normalised = CCLE2019_transcriptomics_normalised[-which(CCLE2019_transcriptomics_normalised$SYMBOL==""),] # remove blank symbol

# cleaning
rm(mart)
rm(gene_name)
rm(ensembl_name)
rm(SYMBOL)
rm(i)
rm(symbol_name)
rm(transcript_name)

### §5 separate SLC genes from translated gene symbol #####
CCLE2019_SLC_ntc = CCLE2019_transcriptomics_normalised[SYMBOL %in% SLC_annotation_2024$HGNC_symbol] # 516 SLC recovered automatically
SLC_annotation_2024[!HGNC_symbol %in% CCLE2019_SLC_ntc$SYMBOL] # 3 SLC is left behind: SLC25A53, SLC5A3, SLC6A14

### §6 CCLE2019 SLC TvsM #####
# create empty matrix to store statistical output
SLC_OutputCor = matrix(0, nrow=nrow(CCLE2019_SLC_ntc),ncol=nrow(CCLE2019_metabolomics))
SLC_OutputP = matrix(1, nrow=nrow(CCLE2019_SLC_ntc),ncol=nrow(CCLE2019_metabolomics))

for(i in 1:nrow(CCLE2019_SLC_ntc)){
  for(j in 1:nrow(CCLE2019_metabolomics)){
    CorIn = cor.test(as.numeric(CCLE2019_SLC_ntc[i,2:ncol(CCLE2019_SLC_ntc)]),
                     as.numeric(CCLE2019_metabolomics[j,2:ncol(CCLE2019_metabolomics)]), method="spearman", exact = FALSE)
    SLC_OutputCor[i,j] = CorIn$estimate
    SLC_OutputP[i,j] = CorIn$p.value
  }
  print(i)
}

# FDR for multiple correlation analysis
SLC_OutputPadj = matrix(1, nrow=nrow(CCLE2019_SLC_ntc),ncol=nrow(CCLE2019_metabolomics))
for (k in 1:nrow(SLC_OutputP)) {
  SLC_OutputPadj[k,] = p.adjust(SLC_OutputP[k,], method = "BH", n = length(SLC_OutputP[k,]))
}

SLC_OutputPadj = as.data.frame(SLC_OutputPadj)
SLC_OutputCor = as.data.frame(SLC_OutputCor)
colnames(SLC_OutputPadj) = CCLE2019_metabolomics$Metabolite #set up rowname and colname of adjusted p value
colnames(SLC_OutputCor) = CCLE2019_metabolomics$Metabolite #set up rowname and colname of correlation value
SLC_OutputPadj$SYMBOL = CCLE2019_SLC_ntc$SYMBOL
SLC_OutputCor$SYMBOL = CCLE2019_SLC_ntc$SYMBOL
SLC_OutputPadj = SLC_OutputPadj |> dplyr::select(SYMBOL, everything())
SLC_OutputCor = SLC_OutputCor |> dplyr::select(SYMBOL, everything())

# melt correlation matrices
mSLC_OutputPadj = melt(SLC_OutputPadj, id.vars = "SYMBOL")
mSLC_OutputCor = melt(SLC_OutputCor, id.vars = "SYMBOL")

# check if melted data frame has the same arrangement
all(mSLC_OutputPadj$SYMBOL == mSLC_OutputCor$SYMBOL) # TRUE
all(mSLC_OutputPadj$variable == mSLC_OutputCor$variable) # TRUE

SLC_interaction = mSLC_OutputPadj
colnames(SLC_interaction) = c("SYMBOL", "METABOLITE", "padj")
SLC_interaction$rho = mSLC_OutputCor$value

fwrite(SLC_interaction, "Results/CCLE2019 TvsM.csv")

## 2.2 NCI60 omics ####
### §1 read files #####
NCI60_transcriptomics = readRDS("Raw Data/NCI60_RNAseq_normalised_counts.rds")
NCI60_metabolomics = fread("01 Raw Data/01 Correlation Analysis/02 NCI60a/WEB_DATA_METABOLON.csv") # 352 metabolite over 58 cell lines (https://wiki.nci.nih.gov/display/ncidtpdata/molecular+target+data)

### §2 format transcriptomics and metabolomics #####
# dcast to form metabolite x cell line matrix
NCI60_metabolomics = NCI60_metabolomics[,c(2,7,8)]
NCI60_metabolomics = dcast(data = NCI60_metabolomics,formula = TITLE~cellname,fun.aggregate = sum,value.var = "VALUE")
# correct name format in metabolomics according to transcript data
overlap = intersect(colnames(NCI60_transcriptomics), colnames(NCI60_metabolomics)[2:ncol(NCI60_metabolomics)]) # 50 overlapping cellname
colnames(NCI60_transcriptomics) [!colnames(NCI60_transcriptomics) %in% overlap] # 10 transcript cell name not in overlap
colnames(NCI60_metabolomics)[2:ncol(NCI60_metabolomics)] [!colnames(NCI60_metabolomics)[2:ncol(NCI60_metabolomics)] %in% overlap] # 8 metabolomics cell name not in overlap
colnames(NCI60_metabolomics) = c(   "Metabolite",      "786-0",           "A498",            "A549",            "ACHN",            "BT-549",          "CAKI-1",          "CCRF-CEM",       
                                    "COLO_205",        "DU-145",          "EKVX",            "HCC-2998",        "HCT-116",         "HCT-15",          "HL-60",           "HOP-62",         
                                    "HOP-92",          "HS_578T",         "HT29",            "IGROV1",          "K-562",           "KM12",            "LOX_IMVI",        "M14",            
                                    "MALME-3M",        "MCF7",            "MDA-MB-231",      "MDA-MB-435",      "MOLT-4",          "NCI-H226",        "NCI-H23",         "NCI-H322M",      
                                    "NCI-H460",        "NCI-H522",        "NCI_ADR-RES",     "OVCAR-3",         "OVCAR-4",         "OVCAR-5",         "OVCAR-8",         "PC-3",           
                                    "RPMI-8226",       "RXF_393",         "SF-268",          "SF-295",          "SF-539",          "SK-MEL-28",       "SK-MEL-5",        "SK-OV-3",        
                                    "SN12C",           "SNB-19",          "SNB-75",          "SR",              "SW-620",          "T-47D",           "TK-10",           "U251",           
                                    "UACC-257",        "UACC-62",         "UO-31"   )

# record of change: A549/ATCC --> A549; COLO 205 --> COLO_205; HL-60(TB) --> HL-60; HS 578T --> HS_578T; LOX IMVI --> LOX_IMVI;
#                   MDA-MB-231/ATCC --> MDA-MB-231; NCI/ADR-RES --> NCI_ADR-RES; RXF-393 --> RXF_393
overlap = intersect(colnames(NCI60_transcriptomics), colnames(NCI60_metabolomics)[2:ncol(NCI60_metabolomics)]) # second verify, 58 cell name matched
# remove unfound cell name in MRN transcript counts (SK-MEL-2, MDA-MB-468)
NCI60_transcriptomics = as.data.frame(NCI60_transcriptomics[, -c(which(colnames(NCI60_transcriptomics) == "SK-MEL-2" |
                                                                       colnames(NCI60_transcriptomics) == "MDA-MB-468"))])
# translate entrez id to gene symbol
NCI60_transcriptomics$SYMBOL = mapIds(org.Hs.eg.db, keys = rownames(NCI60_transcriptomics), keytype = "ENTREZID", column = "SYMBOL")
NCI60_transcriptomics = NCI60_transcriptomics |> dplyr::select(SYMBOL, everything())
# reorder cell name to match transcript data
reorder_idx = match(colnames(NCI60_transcriptomics)[2:ncol(NCI60_transcriptomics)], colnames(NCI60_metabolomics)[2:ncol(NCI60_metabolomics)])
reorder_idx = c(1, reorder_idx+1)
NCI60_metabolomics = NCI60_metabolomics[, reorder_idx]
all(colnames(NCI60_transcriptomics)[2:ncol(NCI60_transcriptomics)] == colnames(NCI60_metabolomics)[2:ncol(NCI60_metabolomics)]) # the cell line order is aligned between transcript counts and metabolomics

# cleaning
rm(overlap)
rm(reorder_idx)

### §3 separate SLC genes from translated gene symbol #####
NCI60_SLC_ntc = NCI60_transcriptomics[which(NCI60_transcriptomics$SYMBOL %in% SLC_annotation_2024$HGNC_symbol),] # 518 SLC transcripts recovered
SLC_annotation_2024[!HGNC_symbol %in% NCI60_SLC_ntc$SYMBOL] # 1 SLC is left behind: SLCO1B7

### §4 NCI60 SLC TvsM #####

# create empty matrix to store statistical output
SLC_OutputCor = matrix(0, nrow=nrow(NCI60_SLC_ntc),ncol=nrow(NCI60_metabolomics))
SLC_OutputP = matrix(1, nrow=nrow(NCI60_SLC_ntc),ncol=nrow(NCI60_metabolomics))

for(i in 1:nrow(NCI60_SLC_ntc)){
  for(j in 1:nrow(NCI60_metabolomics)){
    CorIn = cor.test(as.numeric(NCI60_SLC_ntc[i,2:ncol(NCI60_SLC_ntc)]),
                     as.numeric(NCI60_metabolomics[j,2:ncol(NCI60_metabolomics)]), method="spearman", exact = FALSE)
    SLC_OutputCor[i,j] = CorIn$estimate
    SLC_OutputP[i,j] = CorIn$p.value
  }
  print(i)
}

# FDR for multiple correlation analysis
SLC_OutputPadj = matrix(1, nrow=nrow(NCI60_SLC_ntc),ncol=nrow(NCI60_metabolomics))
for (k in 1:nrow(SLC_OutputP)) {
  SLC_OutputPadj[k,] = p.adjust(SLC_OutputP[k,], method = "BH", n = length(SLC_OutputP[k,]))
}

SLC_OutputPadj = as.data.frame(SLC_OutputPadj)
SLC_OutputCor = as.data.frame(SLC_OutputCor)
colnames(SLC_OutputPadj) = NCI60_metabolomics$Metabolite #set up rowname and colname of adjusted p value
colnames(SLC_OutputCor) = NCI60_metabolomics$Metabolite #set up rowname and colname of correlation value
SLC_OutputPadj$SYMBOL = NCI60_SLC_ntc$SYMBOL
SLC_OutputCor$SYMBOL = NCI60_SLC_ntc$SYMBOL
SLC_OutputPadj = SLC_OutputPadj |> dplyr::select(SYMBOL, everything())
SLC_OutputCor = SLC_OutputCor |> dplyr::select(SYMBOL, everything())

# melt correlation matrices
mSLC_OutputPadj = melt(SLC_OutputPadj, id.vars = "SYMBOL")
mSLC_OutputCor = melt(SLC_OutputCor, id.vars = "SYMBOL")

# check if melted data frame has the same arrangement
all(mSLC_OutputPadj$SYMBOL == mSLC_OutputCor$SYMBOL) # TRUE
all(mSLC_OutputPadj$variable == mSLC_OutputCor$variable) # TRUE

SLC_interaction = mSLC_OutputPadj
colnames(SLC_interaction) = c("SYMBOL", "METABOLITE", "padj")
SLC_interaction$rho = mSLC_OutputCor$value

fwrite(SLC_interaction, "Results/NCI60 TvsM.csv")

## 2.3 CCL180 omics ####
### §1 read files #####
CCL180_ion_intensity = fread("Raw Data/meta_Intensities.csv")
CCL180_ion_annotations = fread("Raw Data/TopIds.csv")

### §2 format transcript counts and metabolomics #####
# construct metabolomics with ion intensity and ion annotation
CCL180_ion_intensity[,V3 := CCL180_ion_annotations$`Top annotation name`]
colnames(CCL180_ion_intensity) = gsub(".*\\/ ", "", colnames(CCL180_ion_intensity))
CCL180_metabolomics = data.table(matrix(ncol = 180, nrow = 1809))
colnames(CCL180_metabolomics) = unique(colnames(CCL180_ion_intensity)[4:ncol(CCL180_ion_intensity)])
CCL180_metabolomics = data.table(Metabolite = CCL180_ion_intensity$V3, CCL180_metabolomics)

for (idx in 2:ncol(CCL180_metabolomics)){
  cell_line_name = colnames(CCL180_metabolomics)[idx]
  sub_ion_intensity = CCL180_ion_intensity[, which(colnames(CCL180_ion_intensity) == cell_line_name), with = FALSE]
  CCL180_metabolomics[,(idx) := rowMeans(sub_ion_intensity)]
  print(idx)
}

# match cell line in metabolomics and transcript counts
CCL180_SLC_ntc = CCLE2019_SLC_ntc
colnames(CCL180_SLC_ntc) = gsub("_.*", "", colnames(CCL180_SLC_ntc))
# trim
CCL180_cell_lines = colnames(CCL180_metabolomics)[2:ncol(CCL180_metabolomics)]
CCLE2019_cell_lines = colnames(CCL180_SLC_ntc)[2:ncol(CCL180_SLC_ntc)]
overlap = intersect(CCL180_cell_lines, CCLE2019_cell_lines) # 151 cell name overlap
CCL180_cell_lines [!CCL180_cell_lines %in% overlap] # 29 CCL180 cell lines not sharing with CCLE2019
CCLE2019_cell_lines [!CCLE2019_cell_lines %in% overlap] # 762 CCLE2019 cell lines not sharing with CCL180

# the reason why only 151, not 153 cell lines claimed in original paper, is that the CCLE2019 transcript counts removed cell lines not in CCLE2019 metabolomics
# and the 2 cell lines removed are included in CCL180. For the sake of consistency, they will not be included in this study

# match
CCL180_metabolomics = CCL180_metabolomics[,c(TRUE, CCL180_cell_lines %in% overlap), with = FALSE]
CCL180_SLC_ntc = CCL180_SLC_ntc[,c(TRUE, CCLE2019_cell_lines %in% overlap), with = FALSE]
CCL180_cell_lines = colnames(CCL180_metabolomics)[2:ncol(CCL180_metabolomics)] # updated CCL180 cell lines
CCLE2019_cell_lines = colnames(CCL180_SLC_ntc)[2:ncol(CCL180_SLC_ntc)] # updated CCLE2019 cell lines
reorder_idx = match(CCLE2019_cell_lines, CCL180_cell_lines)
CCL180_metabolomics = CCL180_metabolomics[,c(1, reorder_idx+1), with = FALSE]
all(CCLE2019_cell_lines == colnames(CCL180_metabolomics)[2:ncol(CCL180_metabolomics)]) # the cell line order is aligned between transcript counts and metabolomics

# adjust duplicated terms
index_of_change = c(650, 890, 1137, 1334, 1402, 1468, 1554, 1670)
metabolite_of_change = c("5a-Androst-3-en-17-one.2", "Acetaminophen glucuronide.2", "Dihydrofukinolide.2", "Acuminoside.2",
                         "alpha-Tocopherol acetate.2", "Theasapogenol A.2", "Musabalbisiane C.2", "Opiorphin.2")
CCL180_metabolomics[index_of_change, Metabolite := metabolite_of_change]

# change ADP-ribose 1""-2"" cyclic phosphate
CCL180_metabolomics[1609, Metabolite := "ADP-ribose 1-2 cyclic phosphate"]

# cleaning
rm(CCL180_cell_lines)
rm(CCLE2019_cell_lines)
rm(cell_line_name)
rm(idx)
rm(overlap)
rm(reorder_idx)
rm(index_of_change)
rm(metabolite_of_change)

### §3 CCL180 SLC TvsM #####

# create empty matrix to store statistical output
SLC_OutputCor = matrix(0, nrow=nrow(CCL180_SLC_ntc),ncol=nrow(CCL180_metabolomics))
SLC_OutputP = matrix(1, nrow=nrow(CCL180_SLC_ntc),ncol=nrow(CCL180_metabolomics))

for(i in 1:nrow(CCL180_SLC_ntc)){
  for(j in 1:nrow(CCL180_metabolomics)){
    CorIn = cor.test(as.numeric(CCL180_SLC_ntc[i,2:ncol(CCL180_SLC_ntc)]),
                      as.numeric(CCL180_metabolomics[j,2:ncol(CCL180_metabolomics)]), method="spearman", exact = FALSE)
    SLC_OutputCor[i,j] = CorIn$estimate
    SLC_OutputP[i,j] = CorIn$p.value
  }
  print(i)
}

# FDR for multiple correlation analysis
SLC_OutputPadj = matrix(1, nrow=nrow(CCL180_SLC_ntc),ncol=nrow(CCL180_metabolomics))
for (k in 1:nrow(SLC_OutputP)) {
  SLC_OutputPadj[k,] = p.adjust(SLC_OutputP[k,], method = "BH", n = length(SLC_OutputP[k,]))
}

SLC_OutputPadj = as.data.frame(SLC_OutputPadj)
SLC_OutputCor = as.data.frame(SLC_OutputCor)
colnames(SLC_OutputPadj) = CCL180_metabolomics$Metabolite #set up rowname and colname of adjusted p value
colnames(SLC_OutputCor) = CCL180_metabolomics$Metabolite #set up rowname and colname of correlation value
SLC_OutputPadj$SYMBOL = CCL180_SLC_ntc$SYMBOL
SLC_OutputCor$SYMBOL = CCL180_SLC_ntc$SYMBOL
SLC_OutputPadj = SLC_OutputPadj |> dplyr::select(SYMBOL, everything())
SLC_OutputCor = SLC_OutputCor |> dplyr::select(SYMBOL, everything())

# melt correlation matrices
mSLC_OutputPadj = melt(SLC_OutputPadj, id.vars = "SYMBOL")
mSLC_OutputCor = melt(SLC_OutputCor, id.vars = "SYMBOL")

# check if melted data frame has the same arrangement
all(mSLC_OutputPadj$SYMBOL == mSLC_OutputCor$SYMBOL) # TRUE
all(mSLC_OutputPadj$variable == mSLC_OutputCor$variable) # TRUE

SLC_interaction = mSLC_OutputPadj
colnames(SLC_interaction) = c("SYMBOL", "METABOLITE", "padj")
SLC_interaction$rho = mSLC_OutputCor$value

fwrite(SLC_interaction, "Results/CCL180 TvsM.csv")

## 2.4 Transformation by Z-score based on metabolites ####
### §1 read files #### 
CCLE2019.raw = fread("Results/CCLE2019 TvsM.csv")
NCI60.raw    = fread("Results/NCI60 TvsM.csv")
CCL180.raw   = fread("Results/CCL180 TvsM.csv")

### §2 transform rho with pseudo Z-score by metabolite #####
Rho_transform_metabolite = function(dataset){
  dataset_metabolite = dataset[, unique(METABOLITE)]
  reference = data.table(Metabolite = dataset_metabolite,
                         Mean = unlist(lapply(dataset_metabolite, function(x){dataset[METABOLITE == x, mean(abs(rho))]})),
                         SD = unlist(lapply(dataset_metabolite, function(x){dataset[METABOLITE == x, sd(abs(rho))]})))
  new_dataset = dataset
  colnames(new_dataset)[4] = "rho_transformed"
  for (i in 1:nrow(reference)){
    rho_transformed = (abs(dataset[METABOLITE == reference[i, Metabolite], rho]) - reference[i, Mean])/reference[i, SD]
    new_dataset[METABOLITE == reference[i, Metabolite], rho_normalised := rho_transformed]
  }
  return(new_dataset)
}

CCLE2019.tra = Rho_transform_metabolite(CCLE2019.raw)
NCI60.tra    = Rho_transform_metabolite(NCI60.raw)
CCL180.tra   = Rho_transform_metabolite(CCL180.raw)

fwrite(CCLE2019.tra, "Results/CCLE2019 TvsM (transformed by metabolite).csv")
fwrite(NCI60.tra,    "Results/NCI60 TvsM (transformed by metabolite).csv")
fwrite(CCL180.tra,   "Results/CCL180 TvsM (transformed by metabolite).csv")

# --------------------------- #
#       Figure 1B data        #
# --------------------------- #

SLC6A6 = CCLE2019.tra[SYMBOL == "SLC6A6"]
ggplot() +
  geom_point(SLC6A6[order(rho_transformed, decreasing=TRUE)][3:nrow(SLC6A6)], mapping = aes(x = rho_transformed, y = -log10(padj)), color = "grey") +
  geom_point(SLC6A6[order(rho_transformed, decreasing=TRUE)][1:2], mapping = aes(x = rho_transformed, y = -log10(padj)), color = "#990000") +
  expand_limits(x = 10) +
  scale_x_continuous(breaks = seq(0,10,5)) +
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,20)) +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))
SLC6A8 = CCLE2019.tra[SYMBOL == "SLC6A8"]
ggplot() +
  geom_point(SLC6A8[order(rho_transformed, decreasing=TRUE)][3:nrow(SLC6A8)], mapping = aes(x = rho_transformed, y = -log10(padj)), color = "grey") +
  geom_point(SLC6A8[order(rho_transformed, decreasing=TRUE)][1:2], mapping = aes(x = rho_transformed, y = -log10(padj)), color = "#990000") +
  expand_limits(x = 10) +
  scale_x_continuous(breaks = seq(0,10,5)) +
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,20)) +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

# 3 CRISPR-Cas9 gene dependency ####
### 3.1 read files #####
CRISPR_GeneEffect = fread("Raw Data/CRISPRGeneEffect.csv")
CCLE2019_metabolomics = fread("Raw Data/CCLE_metabolomics_20190502.csv") # Re-load the original file into work space
Cell_line_metadata = fread("Raw Data/02 CRISPR dependency/Model.csv")

### 3.2 pull SLC_GeneEffect out from CRISPR_GeneEffect ####
# re-format gene column and extract SLC
colnames(CRISPR_GeneEffect) = gsub("\\s*\\([^\\)]+\\)","",colnames(CRISPR_GeneEffect))
SLC_GeneEffect = CRISPR_GeneEffect[, c(TRUE, colnames(CRISPR_GeneEffect)[2:ncol(CRISPR_GeneEffect)] %in% SLC_annotation_2024$HGNC_symbol), with = FALSE] # 500 SLC found
SLC_annotation_2024$HGNC_symbol[!SLC_annotation_2024$HGNC_symbol %in% colnames(SLC_GeneEffect)]

# re-format cell line rows in CERES file and create SLC_GeneEffect and nonSLC_GeneEffect
ach_match = c()
for (Ach in SLC_GeneEffect$ModelID){
  if (Ach %in% CCLE2019_metabolomics$DepMap_ID){
    ach_match = c(ach_match, which(SLC_GeneEffect$ModelID == Ach))
  }
}
SLC_GeneEffect = SLC_GeneEffect[ach_match, ] # 625 cell lines overlapped with CCLE 2019 metabolomics
SLC_GeneEffect = dcast(melt(SLC_GeneEffect, id.vars = "ModelID"), variable ~ ModelID)
colnames(SLC_GeneEffect)[1] = "SYMBOL"

# cleaning
rm(Ach)
rm(ach_match)

### 3.3 transform CCLE2019 metabolomics by z-score #####
# trim unfound cell lines in metadata
CCLE2019_metabolomics = CCLE2019_metabolomics[!is.na(DepMap_ID)]
all(CCLE2019_metabolomics$DepMap_ID %in% Cell_line_metadata$ModelID) # TRUE
Cell_line_metadata = Cell_line_metadata[Cell_line_metadata$ModelID %in% CCLE2019_metabolomics$DepMap_ID,]
# check tissue with only 1 cell lines in metadata
tissue_table = Cell_line_metadata[, .N, by = SampleCollectionSite]
bad_tissue = tissue_table[N == 1, SampleCollectionSite]
# 4 tissue, cervix, prostate, salivary_gland, unknown, only have 1 cell line linked with it
# remove these 4 tissue from both metadata and metabolomics
bad_cell_line = Cell_line_metadata$ModelID[Cell_line_metadata$SampleCollectionSite %in% bad_tissue]
CCLE2019_metabolomics = CCLE2019_metabolomics[-which(CCLE2019_metabolomics$DepMap_ID %in% bad_cell_line), ]
Cell_line_metadata = Cell_line_metadata[-which(Cell_line_metadata$ModelID %in% bad_cell_line),]
tissue_table = Cell_line_metadata[, .N, by = SampleCollectionSite]
# construct tissue information
tissue = tissue_table[, SampleCollectionSite]
tissue_metabolomics_mean = as.data.frame(matrix(nrow = length(tissue), ncol = ncol(CCLE2019_metabolomics)-1))
tissue_metabolomics_sd = as.data.frame(matrix(nrow = length(tissue), ncol = ncol(CCLE2019_metabolomics)-1))

colnames(tissue_metabolomics_mean) = c("Tissue", colnames(CCLE2019_metabolomics)[3:ncol(CCLE2019_metabolomics)])
colnames(tissue_metabolomics_sd) = c("Tissue", colnames(CCLE2019_metabolomics)[3:ncol(CCLE2019_metabolomics)])
tissue_metabolomics_mean$Tissue = tissue
tissue_metabolomics_sd$Tissue = tissue

for (i in 1:length(tissue)){
  cell_line_same_tissue = Cell_line_metadata$ModelID[Cell_line_metadata$SampleCollectionSite == tissue[i]]
  metabolomics_subset = CCLE2019_metabolomics[CCLE2019_metabolomics$DepMap_ID %in% cell_line_same_tissue, ]
  tissue_metabolomics_mean[i,2:ncol(tissue_metabolomics_mean)] = colMeans(metabolomics_subset[,3:ncol(metabolomics_subset)])
  tissue_metabolomics_sd[i,2:ncol(tissue_metabolomics_sd)] = apply(metabolomics_subset[,3:ncol(metabolomics_subset)],2,sd)
}

CCLE2019_metabolomics_z_score = as.data.frame(matrix(nrow = nrow(CCLE2019_metabolomics), ncol = ncol(CCLE2019_metabolomics)-1))
colnames(CCLE2019_metabolomics_z_score) = colnames(CCLE2019_metabolomics)[2:ncol(CCLE2019_metabolomics)]
CCLE2019_metabolomics_z_score$DepMap_ID = CCLE2019_metabolomics$DepMap_ID

CCLE2019_metabolomics = as.data.frame(CCLE2019_metabolomics)
for (i in 1:nrow(CCLE2019_metabolomics)){
  tissue_category = Cell_line_metadata$SampleCollectionSite[Cell_line_metadata$ModelID == CCLE2019_metabolomics$DepMap_ID[i]]
  this_tissue_mean = tissue_metabolomics_mean[which(tissue_metabolomics_mean$Tissue == tissue_category), ]
  this_tissue_sd = tissue_metabolomics_sd[which(tissue_metabolomics_sd$Tissue == tissue_category), ]
  z_score = c()
  for (j in 3:ncol(CCLE2019_metabolomics)){
    z_score = c(z_score, (CCLE2019_metabolomics[i, j] - this_tissue_mean[, j-1])/this_tissue_sd[, j-1])
  }
  CCLE2019_metabolomics_z_score[i, 2:ncol(CCLE2019_metabolomics_z_score)] = z_score
  print(i)
}

### 3.4 re-format CCLE2019 metabolomics to match with SLC_GeneEffect #####

# re-format cell line rows in CCLE2019 metabolomics
ach_match = c()
for (Ach in CCLE2019_metabolomics_z_score$DepMap_ID){
  if (Ach %in% colnames(SLC_GeneEffect)){
    ach_match = c(ach_match, which(CCLE2019_metabolomics_z_score$DepMap_ID == Ach))
  }
}

CCLE2019_metabolomics_z_score = CCLE2019_metabolomics_z_score[ach_match, ]
CCLE2019_metabolomics_z_score = as.data.table(CCLE2019_metabolomics_z_score)
CCLE2019_metabolomics_z_score = dcast(melt(CCLE2019_metabolomics_z_score, id.vars = "DepMap_ID"), variable ~ DepMap_ID)
colnames(CCLE2019_metabolomics_z_score)[1] = "Metabolite"

ach_match = c()
for (Ach in colnames(CCLE2019_metabolomics_z_score)){
  if (Ach %in% colnames(SLC_GeneEffect)){
    ach_match = c(ach_match, which(colnames(SLC_GeneEffect) == Ach))
  }
}

SLC_GeneEffect = SLC_GeneEffect[, c(1, ach_match), with = FALSE]
all(colnames(SLC_GeneEffect)[2:ncol(SLC_GeneEffect)] == colnames(CCLE2019_metabolomics_z_score)[2:ncol(CCLE2019_metabolomics_z_score)]) #TRUE

### 3.5 wilcoxon test between top 20% (more negative) and bottom 20% (less negative) #####
# remove any SLC that have less 20 cell lines in top 20% and negative 20% after removing positive cell lines (i.e. less than 100 negative cell lines)
cell_line_checker = data.frame(SLC_name = SLC_GeneEffect$SYMBOL,
                               Negative_cell_line = unlist(lapply(1:nrow(SLC_GeneEffect), function(x){length(which(as.numeric(SLC_GeneEffect[x, 2:ncol(SLC_GeneEffect)])<0))})))
SLC_GeneEffect = SLC_GeneEffect[!(SLC_GeneEffect$SYMBOL %in% cell_line_checker$SLC_name[which(cell_line_checker$Negative_cell_line < 100)]),]

# Wilcoxon test
SLC_metab_wcx_table = as.data.frame(matrix(nrow = nrow(SLC_GeneEffect), ncol = nrow(CCLE2019_metabolomics_z_score)))
rownames(SLC_metab_wcx_table) = SLC_GeneEffect$SYMBOL
colnames(SLC_metab_wcx_table) = CCLE2019_metabolomics_z_score$Metabolite
SLC_metab_fc_table = SLC_metab_wcx_table

for (i in 1:nrow(SLC_GeneEffect)){
  # select out negative cell line
  SLC_GeneEffect_negative = SLC_GeneEffect[i, (which(as.numeric(SLC_GeneEffect[i,2:ncol(SLC_GeneEffect)])<0)+1), with = FALSE]
  # select out top 20% (more negative) and bottom 20% (less negative) cell lines
  Top20Per = colnames(SLC_GeneEffect_negative[, which(as.numeric(SLC_GeneEffect_negative) < quantile(as.numeric(SLC_GeneEffect_negative), probs = c(0.2))), with = FALSE])
  Bot20Per = colnames(SLC_GeneEffect_negative[, which(as.numeric(SLC_GeneEffect_negative) > quantile(as.numeric(SLC_GeneEffect_negative), probs = c(0.8))), with = FALSE])
  # match to metabolomics
  Top20Per_meta = CCLE2019_metabolomics_z_score[,colnames(CCLE2019_metabolomics_z_score) %in% Top20Per, with = FALSE]
  Bot20Per_meta = CCLE2019_metabolomics_z_score[,colnames(CCLE2019_metabolomics_z_score) %in% Bot20Per, with = FALSE]
  for (j in 1:nrow(CCLE2019_metabolomics_z_score)){
    SLC_metab_wcx_table[i,j] = wilcox.test(as.numeric(Top20Per_meta[j,]), as.numeric(Bot20Per_meta[j,]), exact = FALSE)$p.value
    SLC_metab_fc_table[i,j] = mean(as.numeric(Top20Per_meta[j,])) - mean(as.numeric(Bot20Per_meta[j,]))
  }
  print(i)
}

# adjust p-value
SLC_metab_wcx_table_padj = SLC_metab_wcx_table
for (i in 1:nrow(SLC_metab_wcx_table)){
  SLC_metab_wcx_table_padj[i,] = p.adjust(SLC_metab_wcx_table[i,], method = "BH", n = length(SLC_metab_wcx_table[i,]))
}

# melt data table
SLC_metab_wcx_table_padj_for_melt = SLC_metab_wcx_table_padj
SLC_metab_wcx_table_padj_for_melt$SYMBOL = rownames(SLC_metab_wcx_table_padj_for_melt)
SLC_metab_wcx_table_padj_melted = melt(SLC_metab_wcx_table_padj_for_melt, id.vars = "SYMBOL")

SLC_metab_fc_table_for_melt = SLC_metab_fc_table
SLC_metab_fc_table_for_melt$SYMBOL = rownames(SLC_metab_fc_table_for_melt)
SLC_metab_fc_table_melted = melt(SLC_metab_fc_table_for_melt, id.vars = "SYMBOL")

all(SLC_metab_fc_table_melted$SYMBOL == SLC_metab_wcx_table_padj_melted$SYMBOL) # TRUE
all(SLC_metab_fc_table_melted$variable == SLC_metab_wcx_table_padj_melted$variable) # TRUE

colnames(SLC_metab_wcx_table_padj_melted) = c("SYMBOL", "METABOLITE", "padj")
colnames(SLC_metab_fc_table_melted) = c("SYMBOL", "METABOLITE", "FC")
SLC_dependency_wcx_table = data.frame(SYMBOL = SLC_metab_wcx_table_padj_melted$SYMBOL, METABOLITE = SLC_metab_wcx_table_padj_melted$METABOLITE, 
                                      padj = SLC_metab_wcx_table_padj_melted$padj, FC = SLC_metab_fc_table_melted$FC)

fwrite(SLC_dependency_wcx_table, "Results/CRISPR-Cas9 gene dependency (metabolite in Z-score).csv")

# 4 Construct Evidence Table (et) for every pair ####
## 4.1 read files ####
tb = fread("Raw Data/Translation book.csv")

## 4.2 define functions ####
# ConvertNames #
ConvertNames = function(name, source){
  # Process source input #
  if (is.null(source)){
    stop("Source of metabolite name required!")
  }
  else if (source %in% colnames(tb)){
    # Process name input #
    if (is.null(name)){
      stop("Name of metabolite required!")
    }else if (name %in% tb[,get(source)]){
      ans = tb[get(source) == name]
    }else{
      stop("Name entered not recognised!")
    }
  }else{
    stop("Source entered not recognised!")
  }
  
  # return ans as result #
  return(ans)
}
# BuildEvidenceTable #
BuildEvidenceTable = function(gene, formula){
  # build table #
  ans = data.table(
    SYMBOL = gene, Formula = formula, 
    NCI60.raw = NA, CCLE2019.raw = NA, CCL180.raw = NA,
    NCI60.tra = NA, CCLE2019.tra = NA, CCL180.tra = NA
  )
  # process formula #
  line = try(ConvertNames(formula, "Formula"), silent = TRUE)
  if (grepl("Error", line[[1]]) == TRUE){
    stop("Error in converting formula!")
  }else{
    compound_info = ConvertNames(formula, "Formula")
  }
  info_source = c(
    "NCI60.raw", "CCLE2019.raw", "CCL180.raw", 
    "NCI60.tra", "CCLE2019.tra", "CCL180.tra"
  )
  for (i in 1:length(info_source)){
    if (is.na(compound_info[,get(info_source[i])]) | compound_info[,get(info_source[i])] == ""){
      ans[, info_source[i]:=NA]
    }else{
      type = substr(info_source[i], nchar(info_source[i]) - 2, nchar(info_source[i]))
      if (type == "raw") {
        ans[, info_source[i]:=get(info_source[i])[SYMBOL == gene & METABOLITE == compound_info[,get(info_source[i])]][,rho]]
      }else if (type == "tra"){
        ans[, info_source[i]:=get(info_source[i])[SYMBOL == gene & METABOLITE == compound_info[,get(info_source[i])]][,rho_transformed]]
      }
    }
  }
  return(ans)
}
# Outersect #
outersect = function(x, y) {
  sort(c(x[!x%in%y],
         y[!y%in%x]))
}

## 4.3 Produce Evidence Table (et) ####
specific_SLC_pool = outersect(unique(CCLE2019.raw$SYMBOL), unique(NCI60.raw$SYMBOL))
overlapped_SLC_pool = intersect(unique(CCLE2019.raw$SYMBOL), unique(NCI60.raw$SYMBOL))

# overlapping SLC #
et = data.table(
  SYMBOL = character(), Formula = character(), 
  NCI60.raw = numeric(), CCLE2019.raw = numeric(), CCL180.raw = numeric(),
  NCI60.tra = numeric(), CCLE2019.tra = numeric(), CCL180.tra = numeric()
)
for (i in 1:length(overlapped_SLC_pool)){
  for (j in 1:length(tb$Formula)){
    temp = BuildEvidenceTable(overlapped_SLC_pool[i], tb$Formula[j])
    et = rbind(et, temp)
  }
  print(i)
}

# specific SLC #
"SLC25A53" %in% unique(NCI60.raw$SYMBOL)
"SLC5A3" %in% unique(NCI60.raw$SYMBOL)
"SLC6A14" %in% unique(NCI60.raw$SYMBOL)
"SLCO1B7" %in% unique(CCLE2019.raw$SYMBOL)

et_specific = data.table(
  SYMBOL = character(), Formula = character(), 
  NCI60.raw = numeric(), CCLE2019.raw = numeric(), CCL180.raw = numeric(),
  NCI60.tra = numeric(), CCLE2019.tra = numeric(), CCL180.tra = numeric()
)

# specific SLC - specific to NCI60 #
for (i in 1:(length(specific_SLC_pool)-1)){
  for (j in 1:length(tb$Formula)){
    temp = data.table(
      SYMBOL = specific_SLC_pool[i], Formula = tb$Formula[j], 
      NCI60.raw = NA, CCLE2019.raw = NA, CCL180.raw = NA,
      NCI60.tra = NA, CCLE2019.tra = NA, CCL180.tra = NA
    )
    compound_info = ConvertNames(tb$Formula[j], "Formula")
    if (is.na(compound_info[,NCI60]) | compound_info[,NCI60] == ""){
      temp[, NCI60.raw:=NA]
    }else{
      temp[, NCI60.raw:=NCI60.raw[SYMBOL == specific_SLC_pool[i] & METABOLITE == compound_info[,NCI60]][,rho]]
      temp[, NCI60.tra:=NCI60.tra[SYMBOL == specific_SLC_pool[i] & METABOLITE == compound_info[,NCI60]][,rho_normalised]]
    }
    et_specific = rbind(et_specific, temp)
    print(paste(i, ":", j/length(tb$Formula)))
  }
}
# specific SLC - specific to CCLE2019 #
for (i in 4){
  for (j in 1:length(tb$Formula)){
    temp = data.table(
      SYMBOL = specific_SLC_pool[i], Formula = tb$Formula[j], 
      NCI60.raw = NA, CCLE2019.raw = NA, CCL180.raw = NA,
      NCI60.tra = NA, CCLE2019.tra = NA, CCL180.tra = NA
    )
    compound_info = ConvertNames(tb$Formula[j], "Formula")
    if (is.na(compound_info[,CCLE2019]) | compound_info[,CCLE2019] == ""){
      temp[, NCI60.raw:=NA]
    }else{
      temp[, CCLE2019.raw:=CCLE2019.raw[SYMBOL == specific_SLC_pool[i] & METABOLITE == compound_info[,CCLE2019]][,rho]]
      temp[, CCLE2019.tra:=CCLE2019.tra[SYMBOL == specific_SLC_pool[i] & METABOLITE == compound_info[,CCLE2019]][,rho_normalised]]
    }
    et_specific = rbind(et_specific, temp)
    print(paste(i, ":", j/length(tb$Formula)))
  }
}

# merge et_specific to et #
et_specific = et_specific[rowSums(is.na(et_specific[,3:8]))<6] # remove all interaction with NA in all three matrices
et = rbind(et, et_specific)

# merge with CRISPR dependency wilcoxon test padj
et$CRISPRDEP = vector('numeric',nrow(et))
for (i in 1:nrow(et)){
  compound_info = ConvertNames(et$Formula[i], "Formula")
  if (is.na(compound_info[,CCLE2019]) | compound_info[,CCLE2019] == ""){
    et[i, CRISPRDEP:= NA]
    print(i/nrow(et))
  }else if (!et$SYMBOL[i] %in% unique(SLC_dependency_wcx_table$SYMBOL)){
    et[i, CRISPRDEP:= NA]
    print(i/nrow(et))
  }else{
    et[i, CRISPRDEP:= SLC_dependency_wcx_table[SYMBOL == et$SYMBOL[i] & METABOLITE == compound_info[,CCLE2019]][,padj]]
  }
  print(i/nrow(et))
}

fwrite(et, "Results/Evidence table.csv")

# 5 Concordance and benchmark ####

# --------------------------- #
#       Figure 1C data        #
# --------------------------- #

# Evidence Table (rho).csv
et.raw = et |>
  dplyr::select(SYMBOL, Formula, NCI60.raw, CCLE2019.raw, CCL180.raw) |>
  rename_with(~ c("SYMBOL", "Formula", "NCI60", "CCLE2019", "CCL180"))

# Evidence Table Known Iteractions (rho).csv
ki = fread("Raw Data/Known interaction.csv")[1:837]
ki = ki[!is.na(Formula)][Status != "Derivative"]
etki.raw = ki
colnames(etki.raw)[1] = "SYMBOL"
etki.raw = etki.raw |>
  mutate(
    NCI60 = vector("numeric", nrow(etki.raw)),
    CCLE2019 = vector("numeric", nrow(etki.raw)),
    CCL180 = vector("numeric", nrow(etki.raw)),
    CRISPRDEP = vector("numeric", nrow(etki.raw))
  )
for (i in 1:nrow(etki.raw)){
  if (etki.raw$SYMBOL[i] %in% unique(et.raw$SYMBOL)){
    etki.raw[i, c("NCI60", "CCLE2019", "CCL180", "CRISPRDEP") :=
           as.list(et.raw[SYMBOL == etki.raw$SYMBOL[i] & Formula == etki.raw$Formula[i]][,NCI60:CRISPRDEP])]
  }else{
    etki.raw[i, c("NCI60", "CCLE2019", "CCL180", "CRISPRDEP") := list(NA, NA, NA, NA)]
  }
}
etki.raw = etki.raw[rowSums(is.na(etki.raw[,8:11]))<4]
  
## 5.1 Concordance ####
### §1 All pairs ####
#### CCLE2019 and NCI60 ####
et_1 <- et.raw[rowSums(is.na(et.raw[,3:4])) == 0][, !5:6]
length(unique(et_1$Formula))
length(unique(et_1$SYMBOL))
cor.test(et_1$CCLE2019, et_1$NCI60, method = "spearman", exact = FALSE)$p.value
cor.test(et_1$CCLE2019, et_1$NCI60, method = "spearman", exact = FALSE)$estimate

#### CCLE2019 and CCL180 ####
et_2 = et.raw[rowSums(is.na(et.raw[,4:5])) == 0][, !c(3,6)]
length(unique(et_2$Formula))
length(unique(et_2$SYMBOL))
cor.test(et_2$CCLE2019, et_2$CCL180, method = "spearman", exact = FALSE)$p.value
cor.test(et_2$CCLE2019, et_2$CCL180, method = "spearman", exact = FALSE)$estimate

#### NCI60 and CCL180 ####
et_3 = et.raw[rowSums(is.na(et.raw[,c(3,5)])) == 0][, !c(4,6)]
length(unique(et_3$Formula))
length(unique(et_3$SYMBOL))
cor.test(et_3$NCI60, et_3$CCL180, method = "spearman", exact = FALSE)$p.value
cor.test(et_3$NCI60, et_3$CCL180, method = "spearman", exact = FALSE)$estimate

### §2 Known pairs ####
#### CCLE2019 and NCI60 ####
etki_1 <- etki.raw[rowSums(is.na(etki.raw[,3:4])) == 0]
length(unique(etki_1$Formula))
length(unique(etki_1$SYMBOL))
cor.test(etki_1$CCLE2019, etki_1$NCI60, method = "spearman", exact = FALSE)$p.value
cor.test(etki_1$CCLE2019, etki_1$NCI60, method = "spearman", exact = FALSE)$estimate

#### CCLE2019 and CCL180 ####
etki_2 = etki.raw[rowSums(is.na(etki.raw[,4:5])) == 0]
length(unique(etki_2$Formula))
length(unique(etki_2$SYMBOL))
cor.test(etki_2$CCLE2019, etki_2$CCL180, method = "spearman", exact = FALSE)$p.value
cor.test(etki_2$CCLE2019, etki_2$CCL180, method = "spearman", exact = FALSE)$estimate

#### NCI60 and CCL180 ####
etki_3 = etki.raw[rowSums(is.na(etki.raw[,c(3,5)])) == 0]
length(unique(etki_3$Formula))
length(unique(etki_3$SYMBOL))
cor.test(etki_3$NCI60, etki_3$CCL180, method = "spearman", exact = FALSE)$p.value
cor.test(etki_3$NCI60, etki_3$CCL180, method = "spearman", exact = FALSE)$estimate

## 5.2 Benchmark ####
et.tra = et |>
  dplyr::select(SYMBOL, Formula, NCI60.tra, CCLE2019.tra, CCL180.tra) |>
  rename_with(~ c("SYMBOL", "Formula", "NCI60", "CCLE2019", "CCL180"))
# Generate etki with transformed rho
etki.tra = ki
colnames(etki.tra)[1] = "SYMBOL"
etki.tra = etki.tra |>
  mutate(
    NCI60 = vector("numeric", nrow(etki.tra)),
    CCLE2019 = vector("numeric", nrow(etki.tra)),
    CCL180 = vector("numeric", nrow(etki.tra)),
    CRISPRDEP = vector("numeric", nrow(etki.tra))
  )

for (i in 1:nrow(etki.tra)){
  if (etki.tra$SYMBOL[i] %in% unique(et.tra$SYMBOL)){
    etki.tra[i, c("NCI60", "CCLE2019", "CCL180", "CRISPRDEP") :=
           as.list(et.tra[SYMBOL == etki.tra$SYMBOL[i] & Formula == etki.tra$Formula[i]][,NCI60:CRISPRDEP])]
  }else{
    etki.tra[i, c("NCI60", "CCLE2019", "CCL180", "CRISPRDEP") := list(NA, NA, NA, NA)]
  }
}

etki.tra = etki.tra[rowSums(is.na(etki.tra[,8:11]))<4]

# Define a function for benchmark and distribution plotter
benchmark = function(dataset, source, n_sample = 100){
  etki_gene_pool = unique(dataset$SYMBOL)
  etki_metabolite_pool = unique(dataset$Formula)
  
  fake_etki_collection = lapply(1:n_sample, function(sample_time){
    # initialise fake dataset #
    fake_etki = setNames(data.table(matrix(ncol = ncol(dataset), nrow = nrow(dataset))), colnames(dataset))
    fake_etki[, Substrate := NULL]
    char.col = colnames(fake_etki)[1:3]
    num.col = colnames(fake_etki)[4]
    fake_etki[, (char.col) := lapply(.SD, as.character), .SDcols = char.col]
    fake_etki[, (num.col) := lapply(.SD, as.character), .SDcols = num.col]
    fake_etki[, c("SYMBOL", "Gene_type") := as.list(dataset[, c("SYMBOL", "Gene_type")])]
    # for each SLC, shuffle fake substrate formula #
    for (j in 1:nrow(fake_etki)){
      true_substrate = dataset[SYMBOL == fake_etki$SYMBOL[j], Formula]
      fake_substrate_pool = etki_metabolite_pool[!etki_metabolite_pool %in% true_substrate]
      fake_substrate = sample(fake_substrate_pool, 1)
      while (nrow(et[SYMBOL == fake_etki$SYMBOL[j] & Formula == fake_substrate][,NCI60:CRISPRDEP]) == 0){
        fake_substrate = sample(fake_substrate_pool, 1)
      }
      fake_etki[j, Formula := fake_substrate]
      fake_etki[j, (source) := 
                  et[SYMBOL == fake_etki$SYMBOL[j] & Formula == fake_etki$Formula[j], source, with = FALSE]]
    }
    print(paste(sample_time, "done"))
    fake_etki = fake_etki |> 
      mutate(
        Iteration = sample_time,
        .before = everything()
      )
    return(fake_etki)
  }) |> rbindlist()
  
  return(fake_etki_collection)
}

if(!require("xtable")){
  devtools::install_github("psyteachr/introdataviz")
}
analyseDist = function(data.simu, data.true, dataset){
  colnames(data.true)[5] = "value"
  colnames(data.simu)[5] = "value"
  
  if (!dataset %in% c("CCLE2019", "NCI60", "CCL180")){
    stop("Please enter valid dataset: CCLE2019 | NCI60 | CCL180.")
  }else if (dataset == "CCLE2019") {
    colorset = c("known" = "#8ABFBD", "simulated" = "grey")
  }else if (dataset == "NCI60") {
    colorset = c("known" = "#FDF18A", "simulated" = "grey")
  }else if (dataset == "CCL180") {
    colorset = c("known" = "#987B9F", "simulated" = "grey")
  }
  
  # Bootstrap
  boot = sapply(1:100, function(s){
    resample = sample(data.true[, as.numeric(value)], size = nrow(data.true), replace = TRUE)
    return(mean(resample))
  })
  
  Dist = data.table(
    type = rep(c("simulated", "known"), each = 100),
    value = c(data.simu[, mean(value), by = Iteration][, V1], boot)
  )
  
  rain_height = .1
  
  g = ggplot(Dist, aes(x = "", y = value, fill = type)) +
    # mountain
    introdataviz::geom_flat_violin(trim=FALSE, alpha = 0.4, position = position_nudge(x = -0.6)) +
    # rain
    geom_point(aes(colour = type), size = 1, alpha = 0.4, show.legend = FALSE, position = position_jitter(width = rain_height, height = 0)) +
    geom_boxplot(width = rain_height, alpha = 0.4, show.legend = FALSE, outlier.shape = NA) +
    scale_y_continuous(name = "Average Transformed Spearman's Rho",
                       breaks = seq(-0.2, 0.7, 0.1), 
                       limits = c(-0.3, 0.7)) +
    # Manual color settings
    scale_fill_manual(values = colorset) + # For violin plots
    scale_color_manual(values = colorset) + # For points
    coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
  
  p = t.test(
    x = Dist[type == "simulated", value],
    y = Dist[type == "known", value],
    alternative = "less"
  )$p.value
  
  return(list(dist.data = Dist, dist.plot = g, dist.pvalue = p))
}

# separate into individual etki
NCI60_etki = etki.tra[!is.na(NCI60)][, c(1,2,4,5,8)]
CCLE2019_etki = etki.tra[!is.na(CCLE2019)][, c(1,2,4,5,9)]
CCL180_etki = etki.tra[!is.na(CCL180)][, c(1,2,4,5,10)]
CRISPR_etki = etki.tra[!is.na(CRISPRDEP)][, c(1,2,4,5,11)]

# generate individual simulation
NCI60_simulated_random = benchmark(NCI60_etki, "NCI60")
CCLE2019_simulated_random = benchmark(CCLE2019_etki, "CCLE2019")
CCL180_simulated_random = benchmark(CCL180_etki, "CCL180")
CRISPR_simulated_random = benchmark(CRISPR_etki, "CRISPRDEP")

### §1 Cell panel benchmarking ####
NCI60_dist = analyseDist(data.simu = NCI60_simulated_random, data.true = NCI60_etki, dataset = "NCI60")
NCI60_dist$dist.plot
NCI60_dist$dist.pvalue

CCLE2019_dist = analyseDist(data.simu = CCLE2019_simulated_random, data.true = NCI60_etki, dataset = "CCLE2019")
CCLE2019_dist$dist.plot
CCLE2019_dist$dist.pvalue

CCL180_dist = analyseDist(data.simu = CCL180_simulated_random, data.true = NCI60_etki, dataset = "CCL180")
CCL180_dist$dist.plot
CCL180_dist$dist.pvalue

# --------------------------- #
#       Figure 1D data        #
# --------------------------- #

# Fig 1D (cell panel).csv
benchmark_dist = rbind(
  NCI60_dist$dist.data |> mutate(dataset = "NCI60", .before = everything()),
  CCLE2019_dist$dist.data |> mutate(dataset = "CCLE2019", .before = everything()),
  CCL180_dist$dist.data |> mutate(dataset = "CCL180", .before = everything())
)

### §2 Gene dependency benchmarking ####
CRISPR_report = data.table(
  Type = c("Known", rep("Simulated", 100)),
  Mean = c(CRISPR_etki[, mean(CRISPRDEP)], 
           unlist(lapply(CRISPR_simulated_random, function(x){x[, mean(as.numeric(CRISPRDEP))]}))),
  Fraction = c(CRISPR_etki[CRISPRDEP < 0.05, .N]/CRISPR_etki[,.N], 
               unlist(lapply(CRISPR_simulated_random, function(x){x[CRISPRDEP < 0.05, .N]/x[,.N]})))
)

# CRISPR Dependency recovery across alpha
known_recovery = data.table(Type = rep("Known", 100),
                            AlphaCutoff = seq(0.001, 0.1, 0.001),
                            Fraction_recovery = vector("numeric", 100)
)
for (i in 1:nrow(known_recovery)){
  nr = CRISPR_etki[as.numeric(CRISPRDEP) < known_recovery[i, AlphaCutoff], .N]/CRISPR_etki[,.N]
  known_recovery[i, Fraction_recovery := nr]
}

simulated_recovery = data.table(Type = rep("Simulated", 10000),
                                AlphaCutoff = rep(seq(0.001, 0.1, 0.001), each = 100),
                                Fraction_recovery = vector("numeric", 10000)
)
for (i in seq(0.001, 0.1, 0.001)){
  nr = unlist(lapply(CRISPR_simulated_random, function(x){x[CRISPRDEP < i, .N]/x[,.N]}))
  simulated_recovery[AlphaCutoff == i, Fraction_recovery := nr]
}

# --------------------------- #
#       Figure 2B data        #
# --------------------------- #

# CRISPRDEP_FractionRecoveryAcrossAlpha.csv
recovery = rbind(known_recovery, simulated_recovery)

# 6 Adjacency ####
library(MetaboSignal)
## 6.1 The generation of adjacency matrix ####
# Note: adjacency matrix is provided as a supplementary table (S8), please refer to section 1.5.
#       The code below demonstrates how adjacency matrix is produced.
# get human pathways
human_pathways = MS_getPathIds(organism_code = "hsa")
human_metabolic_pathways = human_pathways[human_pathways[, "Path_type"] == "metabolic",]
human_signaling_pathways = human_pathways[human_pathways[, "Path_type"] == "signaling",]

# build reaction network
hmn = MS_reactionNetwork(human_metabolic_pathways[, "Path_id"]) # human metabolic network

# build adjacency matrix
compound = unique(hmn[,"source"][grep("cpd", hmn[,"source"])])
AM = matrix(nrow = length(compound), ncol = length(compound))
colnames(AM) = compound
rownames(AM) = compound
pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                      total = length(compound),
                      complete = "=",   # Completion bar character
                      incomplete = "-", # Incomplete bar character
                      current = ">",    # Current bar character
                      clear = FALSE,    # If TRUE, clears the bar when finish
                      width = 100)      # Width of the progress bar
for (i in 1:length(compound)){
  pb$tick()
  for (j in 1:length(compound)){
    if (i == j){
      AM[i,j] = 0
    }else{
      distance = MS_shortestPaths(hmn, compound[i], compound[j], mode = "out", type = "first")
      if (is.null(distance)){
        AM[i,j] = Inf
      }else{
        AM[i,j] = (length(distance)-1)/2
      }
    }
  }
}

MS_exportCytoscape(hmn, "hsa")

AM = as.data.table(AM)
AM[, COMPOUND := compound]
setcolorder(AM, "COMPOUND")

fwrite(AM, "Results/Adjacency matrix.csv")

## 6.2 General behavior of metabolites that are proximal and distant to the metabolite of interest ####
colnames(AM) = gsub("cpd:", "", colnames(AM))
AM$COMPOUND = gsub("cpd:", "", AM$COMPOUND)

### §1 Define functions ####
# Remember to load ConvertNames() function under 6.2 if you jump start #

# getAM #

# this function finds the derivative compound of input formula, return a named list of kegg id, with the name indicating the steps
# required to reach the vertices

# step.inf = NULL returns NA (not passing qc) or finite derivatives;
# step.inf = TRUE returns NA (not passing qc) or infinite derivatives;

getAM = function(formula, step.inf = NULL, is.kegg = FALSE){
  # if the formula entered is already a kegg id (is.kegg = TRUE)
  if (isTRUE(is.kegg)){
    kegg_id = formula
  }else{
    kegg_id = as.character(ConvertNames(formula, "Formula")[,"KEGG ID"])
  }
  qc = TRUE # logical term of quality control
  generic.lipid = c("C00157", "C00165", "C00195", "C00350", "C00416", "C00422", "C00550", "C02530", "C04230", "C04438") # kegg id belonging to lipidomics
  # if kegg_id is NA in translation table #
  if(is.na(kegg_id)){
    qc = FALSE
    # if kegg_id is not NA, but not found in the adjacency matrix #
  }else if (nchar(kegg_id) == 6 & !kegg_id %in% AM$COMPOUND){
    qc = FALSE
    # if kegg_id is not NA, only has one associated kegg id, but is lipidomic #
  }else if (nchar(kegg_id) == 6 & kegg_id %in% generic.lipid){
    qc = FALSE
    # if kegg_id is not NA, only has one associated kegg id, and is found in the adjacency matrix #
  }else if (nchar(kegg_id) == 6 & kegg_id %in% AM$COMPOUND){
    kegg_id_qc = kegg_id
    # if kegg_id is not NA, but contain multiple kegg id linked with the input formula #
  }else if (nchar(kegg_id) > 6){
    kegg_id = gsub(" ", "", unlist(strsplit(kegg_id, ";"))) # break the kegg id into individual id and remove space
    kegg_id_qc = kegg_id[kegg_id %in% AM$COMPOUND & !kegg_id %in% generic.lipid] # only select the kegg id found in AM and not lipid
    # if all kegg ids are not found in AM
    if (rlang::is_empty(kegg_id_qc)){
      qc = FALSE
    }
  }else{
    qc = FALSE
  }
  
  # if compound passed qc and there's no requirement of finding infinite metabolites #
  if(qc & is.null(step.inf)){
    ans = c()
    # gather adjacent metabolite for the given kegg id #
    for (i in 1:length(kegg_id_qc)){
      AM_related = AM[COMPOUND == kegg_id_qc[i]]
      AM_related = AM_related[,.SD,.SDcols=which(AM_related[,.SD] != Inf & AM_related[,.SD] != 0)]
      # if there is no adjacent metabolite derivative for the given kegg id #
      if (ncol(AM_related) == 1){
        ans_temp = character(0)
        # if there is adjacent metabolite derivative for the given kegg id #
      }else{
        ans_temp = colnames(AM_related)[2:ncol(AM_related)]
        names(ans_temp) = as.numeric(AM_related[,2:ncol(AM_related)])
        ans = c(ans, ans_temp)
      }
    }
    ans = ans[!ans %in% generic.lipid] # remove lipid derivatives in the answer
    # if the accumulated number of adjacent metabolite derivative equals to 0 #
    if (is.null(ans)){
      return(NA)
      # if the accumulated number of adjacent metabolite derivative does not equals to 0 #
    }else{
      return(ans)
    }
    # if compound passed qc and there's requirement of finding infinite metabolites #
  }else if (qc & !is.null(step.inf)){
    ans = vector("list", length(kegg_id_qc))
    # gather infinite metabolite for the given kegg id #
    for (i in 1:length(kegg_id_qc)){
      AM_related = AM[COMPOUND == kegg_id_qc[i]]
      AM_related = AM_related[,.SD,.SDcols=which(AM_related[,.SD] == Inf & AM_related[,.SD] != 0)]
      ans[[i]] = colnames(AM_related)[2:ncol(AM_related)]
    }
    # if there is only 1 kegg id passing qc, return the info of infinite metabolites #
    if (length(ans) == 1){
      return(ans[[1]])
      # if there is multiple kegg ids passing qc, return the intersect of the infinite metabolites #
    }else{
      ans_intersect = intersect(ans[[1]], ans[[2]])
      if (length(ans) == 2){
        return(ans_intersect)
      }else{
        for (j in 3:length(ans)){
          ans_intersect = intersect(ans_intersect, ans[[j]])
        }
        return(ans_intersect)
      }
    }
    # if compound did not pass qc #
  }else{
    return(NA)
  }
}

# translate() capture and translate the kegg id from the output of getAM()
translate = function(KEGG.id){
  kegg_info = data.table(KEGG = as.character(KEGG.id),
                         Formula = vector("character", length(KEGG.id)))
  for (j in 1:nrow(kegg_info)){
    info = grep(kegg_info$KEGG[j], tb$`KEGG ID`)
    if (length(info) == 0){
      kegg_info[j, Formula := NA]
    }else if (length(info) == 1){
      kegg_info[j, Formula := tb$Formula[info]]
    }else if (length(info) > 1){
      kegg_info[j, Formula := paste(tb$Formula[info], collapse = "<>")]
    }
  }
  # Remove the derivative that's not found in the translation book #
  kegg_info = na.omit(kegg_info)
  # if no derivative is found in the translation book, stop translation #
  if (nrow(kegg_info) == 0){
    return(NA)
    # if there's only one formula matching the given kegg
  }else if (length(grep("<>", kegg_info$Formula)) == 0){
    return(kegg_info)
  }else{
    # Break the <> and insert a new row #
    kegg_info = kegg_info |> 
      mutate(Formula = strsplit(as.character(Formula), "<>")) |> 
      unnest(Formula) |>
      group_by(Formula) |>
      summarise(KEGG = toString(sort(unique(KEGG)), .groups='drop'))
    kegg_info = as.data.table(kegg_info)
    return(kegg_info)
  }
}

# getNR() gathers the information for formula after the formula is translated by translate()
getNR = function(gene, translated.formula){
  data = data.table(Formula = translated.formula[, Formula],
                    NCI60 = vector("numeric", nrow(translated.formula)),
                    CCLE2019 = vector("numeric", nrow(translated.formula)),
                    CCL180 = vector("numeric", nrow(translated.formula)))
  for (k in 1:nrow(data)){
    if(nrow(et[SYMBOL == gene & Formula == data[k, Formula]]) != 0){
      data[k, c("NCI60", "CCLE2019", "CCL180") :=
             as.list(et[SYMBOL == gene & Formula == data[k, Formula]][,NCI60:CCL180])]
    }else{
      data[k, c("NCI60", "CCLE2019", "CCL180") := as.list(NA, NA, NA)]
    }
  }
  return(data)
}

# getSource() acquires the metabolomics which contains the concentration info for this formula
getSource = function(formula){
  logical = as.logical((is.na(tb[Formula == formula, NCI60:CCL180])))
  ans = colnames(tb[,NCI60:CCL180])[!logical]
  return(ans)
}

# check_NA can determine if the result is an NA value
check_NA = function(x){
  if (length(x) == 1){
    if (is.na(x)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}

# channel calculates the PS difference between compounds and derivatives
channel = function(derivative_table, substrate_table, s_source){
  # collecting rho and p value #
  rho_diff_val = c()
  for (m in 1:nrow(derivative_table)){
    d_source = getSource(derivative_table$Formula[m]) # the source of derivatives
    channeling = length(s_source)*length(d_source) # calculate channeling value
    ### 1 substrate source and 1 derivative source ###
    if (channeling == 1){
      rho_diff_val = c(rho_diff_val, as.numeric(substrate_table[, s_source, with = FALSE]) - 
                         as.numeric(derivative_table[m, d_source, with = FALSE])
      )
    }
    ### 1 comp source and 2 deri source | 2 comp source and 1 deri source ###
    else if(channeling == 2 | channeling == 3){
      same_source = intersect(s_source, d_source)
      ## if there is one omic source overlapped ##
      if (length(same_source) != 0){
        rho_diff_val = c(rho_diff_val, as.numeric(substrate_table[, same_source, with = FALSE]) - 
                           as.numeric(derivative_table[m, same_source, with = FALSE])
        )
        ## if there is no omic source overlapped, calculate the average of correlation statistics ##
      }else{
        diff = c()
        if (length(s_source) > length(d_source)){
          for (n in 1:length(s_source)){
            diff = c(diff, as.numeric(substrate_table[, s_source[n], with = FALSE]) -
                       as.numeric(derivative_table[m, d_source, with = FALSE]))
          }
          rho_diff_val = c(rho_diff_val, mean(diff, na.rm = TRUE))
        }else{
          for (n in 1:length(d_source)){
            diff = c(diff, as.numeric(substrate_table[, s_source, with = FALSE]) -
                       as.numeric(derivative_table[m, d_source[n], with = FALSE]))
          }
          rho_diff_val = c(rho_diff_val, mean(diff, na.rm = TRUE))
        }
      }
    }
    ### 2 comp source and 2 deri source, 3 comp source and 3 deri source ###
    else if (channeling == 4 | channeling == 9){
      same_source = intersect(s_source, d_source)
      ## if there is 1 omics source shared ##
      if (length(same_source) == 1){
        rho_diff_val = c(rho_diff_val, as.numeric(substrate_table[, same_source, with = FALSE]) - 
                           as.numeric(derivative_table[m, same_source, with = FALSE])
        )
        ## if there are 2 or 3 omics source shared, calculate the average of correlation statistics ##
      }else if (length(same_source) > 1){
        diff = c()
        for (n in 1:length(same_source)){
          diff = c(diff, as.numeric(substrate_table[, same_source[n], with = FALSE]) - 
                     as.numeric(derivative_table[m, same_source[n], with = FALSE]))
        }
        rho_diff_val = c(rho_diff_val, mean(diff, na.rm = TRUE))
      }
    }
    ### 2 comp source and 3 deri source, 3 comp source and 2 deri source ###
    else if (channeling == 6){
      ## however change, there will be 2 overlapping omics ##
      same_source = intersect(s_source, d_source)
      diff = c()
      for (n in 1:length(same_source)){
        diff = c(diff, as.numeric(substrate_table[, same_source[n], with = FALSE]) - 
                   as.numeric(derivative_table[m, same_source[n], with = FALSE]))
      }
      rho_diff_val = c(rho_diff_val, mean(diff, na.rm = TRUE))
    }
  }
  return(rho_diff_val)
}

### §2 Investigate behavior ####
# Note: Both etki.raw and etki.tra are investigated, remove comment mark below to produce data of interest:

# etki = etki.raw
# etki = etki.tra

#### Quality control of the etki ####
# QC: filter out compounds with derivatives of limited quality #
QC_category = c()
for (i in 1:nrow(etki)){
  # acquire derivative for compound input
  temp = getAM(etki$Formula[i])
  
  # *QC1*: filter out blindly all compound with a) compound input not found in the given matrix;                #
  #                                             b) derivative not found in the given matrix;                    #
  #                                             c) lesser amount of derivatives are found in the given matrix.  #
  if (length(temp) < 3){
    QC_category = c(QC_category, i)
  }
  print(i)
}

print(
  paste0(
    round(100*(1-length(QC_category)/nrow(etki)), 2),
    "% interaction passess QC, covering ",
    round(100*length(etki[-QC_category][, unique(Formula)])/length(etki[, unique(Formula)]), 2),
    "% of all metabolites."
  )
)

etki_qc = etki[-QC_category]

#### Calculating WT.pvalue for etki_qc across steps ####
step_table = data.table(Steps = seq(1, 10, 1),
                        WT.pvalue = vector("numeric", length(seq(1, 10, 1))))
Inf_across_step = vector("list", nrow(step_table))
for (i in 1:nrow(step_table)){
  step = step_table[i, Steps]
  # make dNR_step #
  dNR_step = data.table(SYMBOL = etki_qc[,SYMBOL],
                        Formula = etki_qc[,Formula],
                        dNR = vector("numeric", nrow(etki_qc))
  )
  for (j in 1:nrow(etki_qc)){
    temp = getAM(etki_qc[j, Formula])
    proximal = temp[which(names(temp) == step)]
    if (length(proximal) == 0){             
      dNR_step[j, dNR := NaN]                
    }else{                                   
      proximal_temp = translate(proximal)   
      if (length(proximal_temp) == 1){      
        if (is.na(proximal_temp)){           
          dNR_step[j, dNR := NaN]           
        }else{
          proximal_temp = getNR(etki_qc[j, SYMBOL], proximal_temp)
          proximal_dNR = channel(proximal_temp, etki_qc[j], getSource(etki_qc[j, Formula]))
          dNR_step[j, dNR := mean(abs(proximal_dNR))]
        }
      }else{
        proximal_temp = getNR(etki_qc[j, SYMBOL], proximal_temp)
        proximal_dNR = channel(proximal_temp, etki_qc[j], getSource(etki_qc[j, Formula]))
        dNR_step[j, dNR := mean(abs(proximal_dNR))]
      }
    }
    print(paste("step", i ,":", j))
  }
  dNR_step = na.omit(dNR_step)
  
  # take the dNR_step and make dNR_inf #
  dNR_inf = data.table(SYMBOL = dNR_step[,SYMBOL],
                       Formula = dNR_step[,Formula],
                       dNR = vector("numeric", nrow(dNR_step))
  )
  for (j in 1:nrow(dNR_inf)){
    infinite = getAM(dNR_inf[j, Formula], step.inf = T)
    if (length(infinite) == 0){            
      dNR_inf[j, dNR := NaN]               
    }else{                                  
      infinite_temp = translate(infinite)   
      if (length(infinite_temp) == 1){       
        if (is.na(infinite_temp)){         
          dNR_inf[j, dNR := NaN]            
        }else{
          infinite_temp = infinite_temp[sample(1:nrow(infinite_temp), 20)]
          infinite_temp = getNR(dNR_inf[j, SYMBOL], infinite_temp)
          etki_qc_reference = etki_qc[SYMBOL == dNR_inf[j, SYMBOL] & Formula == dNR_inf[j, Formula]]
          if(nrow(etki_qc_reference) > 1){
            etki_qc_reference = etki_qc_reference[1]
          }
          infinite_dNR = channel(infinite_temp, etki_qc_reference, getSource(dNR_inf[j, Formula]))
          dNR_inf[j, dNR := mean(abs(infinite_dNR))]
        }
      }else{
        infinite_temp = infinite_temp[sample(1:nrow(infinite_temp), 20)]
        infinite_temp = getNR(dNR_inf[j, SYMBOL], infinite_temp)
        etki_qc_reference = etki_qc[SYMBOL == dNR_inf[j, SYMBOL] & Formula == dNR_inf[j, Formula]]
        if(nrow(etki_qc_reference) > 1){
          etki_qc_reference = etki_qc_reference[1]
        }
        infinite_dNR = channel(infinite_temp, etki_qc_reference, getSource(dNR_inf[j, Formula]))
        dNR_inf[j, dNR := mean(abs(infinite_dNR))]
      }
    }
    print(paste("step inf:", j))
  }
  dNR_inf = na.omit(dNR_inf)
  WT = wilcox.test(dNR_step$dNR, dNR_inf$dNR, alternative = "less")$p.value
  # store data #
  Inf_across_step[[i]] = dNR_inf
  step_table[i, WT.pvalue := WT]
  print(paste("step", i, "completed"))
}

#### Simulated etki_qc compare to random #####
step_table_benchmarking = data.table(Steps = rep(seq(1,10,1), each = 100),
                                     WT.pvalue = vector("numeric", length(rep(seq(1,10,1), each = 100)))
)

for (i in 1:nrow(step_table_benchmarking)){
  step = step_table_benchmarking[i, Steps]
  # make dNR_step #
  dNR_step = data.table(SYMBOL = etki_qc[,SYMBOL],
                        Formula = etki_qc[,Formula],
                        dNR = vector("numeric", nrow(etki_qc))
  )
  for (j in 1:nrow(etki_qc)){
    temp = getAM(etki_qc[j, Formula])
    proximal = temp[which(names(temp) == step)]
    if (length(proximal) == 0){             
      dNR_step[j, dNR := NaN]               
    }else{                                  
      proximal_temp = translate(proximal)    
      if (length(proximal_temp) == 1){      
        if (is.na(proximal_temp)){          
          dNR_step[j, dNR := NaN]           
        }else{
          temp_fake = getAM(etki_qc[j, Formula], step.inf = T)
          temp_fake = translate(temp_fake)
          temp_fake = temp_fake[sample(1:nrow(temp_fake), length(proximal_temp))]
          temp_fake = getNR(etki_qc[j, SYMBOL], temp_fake)
          fake_dNR = channel(temp_fake, etki_qc[j], getSource(etki_qc[j, Formula]))
          dNR_step[j, dNR := mean(abs(fake_dNR))]
        }
      }else{
        temp_fake = getAM(etki_qc[j, Formula], step.inf = T)
        temp_fake = translate(temp_fake)
        temp_fake = temp_fake[sample(1:nrow(temp_fake), length(proximal_temp))]
        temp_fake = getNR(etki_qc[j, SYMBOL], temp_fake)
        fake_dNR = channel(temp_fake, etki_qc[j], getSource(etki_qc[j, Formula]))
        dNR_step[j, dNR := mean(abs(fake_dNR))]
      }
    }
    print(paste("row", i ,":", j))
  }
  dNR_step = na.omit(dNR_step)
  
  # compare with stored dNR_inf #
  dNR_inf = Inf_across_step[[(step)]]
  WT = wilcox.test(dNR_step$dNR, dNR_inf$dNR)$p.value
  # store data #
  step_table_benchmarking[i, WT.pvalue := WT]
  print(paste("row", i, "completed"))
}

# adjust p-value for every step due to multiple comparison #
step_table_benchmarking[, class := rep("simulated", nrow(step_table_benchmarking))]
step_table[, class := rep("known", nrow(step_table))]

# If you chose etki.raw for etki under 6.2$2
# --------------------------- #
#       Figure 3B data        #
# --------------------------- #

# If you chose etki.tra for etki under 6.2$2
# --------------------------- #
#      Figure S2A data        #
# --------------------------- #

step_assembly = rbind(step_table_benchmarking, step_table)

for (step in seq(1, 10, 1)){
  p = step_assembly[Steps == step, WT.pvalue]
  step_assembly[Steps == step, WT.pvalue := p.adjust(p, method = "BH", n = length(p))]
  print(step)
}

step_simulated = step_assembly[class == "simulated"]
step_known = step_assembly[class == "known"]

ggplot() +
  # box plot of fecs
  geom_boxplot(data = step_simulated, 
               aes(x = Steps, y= -log10(WT.pvalue), group = Steps), fill = "grey") +
  # points of ecs
  geom_line(data = step_known,
            aes(x= Steps, y= -log10(WT.pvalue)),
            color = "#990000",
            linetype="dashed") + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks=seq(0,10,1))

# 7 Confidence scores ####
ki = fread("01 Raw Data/00 SLC annotation 2024/known interaction/Known interaction formula 20240321.csv")[1:837]
ki = ki[!is.na(Status)][!is.na(Formula)]

etki = ki
colnames(etki)[1] = "SYMBOL"
etki = etki |>
  mutate(
    NCI60 = vector("numeric", nrow(etki)),
    CCLE2019 = vector("numeric", nrow(etki)),
    CCL180 = vector("numeric", nrow(etki)),
    CRISPRDEP = vector("numeric", nrow(etki))
  )

for (i in 1:nrow(etki)){
  if (etki$SYMBOL[i] %in% unique(et.tra$SYMBOL)){
    etki[i, c("NCI60", "CCLE2019", "CCL180", "CRISPRDEP") :=
           as.list(et.tra[SYMBOL == etki$SYMBOL[i] & Formula == etki$Formula[i]][,NCI60:CRISPRDEP])]
  }else{
    etki[i, c("NCI60", "CCLE2019", "CCL180", "CRISPRDEP") := list(NA, NA, NA, NA)]
  }
}

etki = etki[rowSums(is.na(etki[,8:11]))<4]
etki_Tr = etki

## 7.1 Simulation generation ####
fake_etki_Tr = lapply(1:100, function(sample_time){
  etki_gene_pool = unique(etki_Tr$SYMBOL)
  etki_metabolite_pool = unique(etki_Tr$Formula)
  # initialise fake etki_Tr data table #
  fake_etki = setNames(data.table(matrix(ncol = ncol(etki_Tr), nrow = nrow(etki_Tr))), colnames(etki_Tr))
  fake_etki[, Substrate := NULL]
  char.col = colnames(fake_etki)[1:6]
  num.col = colnames(fake_etki)[7:10]
  fake_etki[, (char.col) := lapply(.SD, as.character), .SDcols = char.col]
  fake_etki[, (num.col) := lapply(.SD, as.character), .SDcols = num.col]
  fake_etki[, c("SYMBOL", "Gene_type","Localisation") := as.list(etki_Tr[, c("SYMBOL", "Gene_type","Localisation")])]
  fake_etki[, Status := rep("Fake", nrow(fake_etki))]
  # for each SLC, shuffle fake substrate formula #
  for (j in 1:nrow(fake_etki)){
    true_substrate = etki_Tr[SYMBOL == fake_etki$SYMBOL[j]][, Formula]
    fake_substrate_pool = etki_metabolite_pool[!etki_metabolite_pool %in% true_substrate]
    fake_substrate = sample(fake_substrate_pool, 1)
    while (nrow(et.tra[SYMBOL == fake_etki$SYMBOL[j] & Formula == fake_substrate][,NCI60:CRISPRDEP]) == 0){
      fake_substrate = sample(fake_substrate_pool, 1)
    }
    fake_etki[j, Formula := fake_substrate]
    fake_etki[j, c("NCI60", "CCLE2019", "CCL180", "CRISPRDEP") := 
                as.list(et.tra[SYMBOL == fake_etki$SYMBOL[j] & Formula == fake_etki$Formula[j]][,NCI60:CRISPRDEP])]
  }
  print(paste(sample_time, "done"))
  fake_etki = fake_etki |>
    mutate(
      Iteration = sample_time,
      .before = everything()
    )
  return(fake_etki)
}) |> rbindlist()

## 7.2 Cell panel parameter (ToD) optimisation ####
ConfidenceScore_A_by_set = function(data, ki, a, b, c, score.only = FALSE){
  # Factor A: Predictive Scores
  ## quantile ##
  ns_dist1 = na.omit(ki$NCI60) # NCI60 #
  ns_dist2 = na.omit(ki$CCLE2019) # CCLE2019 #
  ns_dist3 = na.omit(ki$CCL180) # CCL180 #
  q1 = quantile(ns_dist1[ns_dist1 > as.numeric(a)], seq(0, 1, by = 0.1))
  q2 = quantile(ns_dist2[ns_dist2 > as.numeric(b)], seq(0, 1, by = 0.1))
  q3 = quantile(ns_dist3[ns_dist3 > as.numeric(c)], seq(0, 1, by = 0.1))
  
  Factor_A_info = data.table(Factor_A_name = c("NCI60", "CCLE2019", "CCL180"),
                             Factor_A_dist = c("q1", "q2", "q3"),
                             Factor_A_threshold = c(as.numeric(a), as.numeric(b), as.numeric(c)))
  Score = sapply(1:nrow(data), function(i){
    # score prediction for 3 matrices #
    Factor_A = vector("integer", nrow(Factor_A_info))
    for (j in 1:nrow(Factor_A_info)){
      n = as.numeric(data[i, get(Factor_A_info[j, Factor_A_name])])
      if (is.na(n) | n < Factor_A_info[j, Factor_A_threshold]){
        Factor_A[j] = 0
      } else if (n > Factor_A_info[j, Factor_A_threshold]){
        Factor_A[j] = findInterval(n, get(Factor_A_info[j, Factor_A_dist]))*3
      }
    }
    Factor_A = sum(Factor_A)
    return(as.numeric(Factor_A))
  })
  
  if (score.only){
    return(Score)
  }else{
    score_table = data |>
      mutate(
        Score = Score
      )
    return(score_table)
  }
}
Best_threshold_detect_A_by_set = function(ki.data, si.data, from.range){
  # Setting up parallel
  numCores = 8  # Use all available cores except one for other processes
  cl = makeCluster(numCores)
  
  # a: range for NCI60
  registerDoParallel(cl)
  range_a = foreach(a = from.range, .combine = rbind, .export = ls(globalenv())) %dopar% {
    library(data.table)
    library(tidyverse)
    # calculate scores for etki and fake etki under threshold a
    # i.e. if we only consider NCI60 transformed-rho > a, and assign confidence score based on a
    #      what would be the confidence score for each SLC-metabolite pair?
    ecs = ConfidenceScore_A_by_set(data = etki_Tr, ki = ki.data, a = a, b = Inf, c = Inf, score.only = TRUE)
    fecs = lapply(1:100, function(f){
      simuScore = ConfidenceScore_A_by_set(data = si.data[Iteration == f], ki = ki.data, a = a, b = Inf, c = Inf, score.only = TRUE)
      return(simuScore)
    })
    
    # calculate fraction recovery (TP - FP) for each given score cutoff
    # i.e. how many SLC-metabolite pair can be discovered, if we set score cutoff from 0 (all) - 34 (none)
    fr = lapply(0:34, function(cutoff){
      TP = length(ecs[ecs >= cutoff])/length(ecs)
      FP = sapply(fecs, function(f){length(f[f >= cutoff])}/length(f))
      
      if (TP == 0 & mean(FP) == 0){
        pvalue = 1
      }else if (TP == 1 & mean(FP) == 1) {
        pvalue = 1
      }else{
        pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
      }
      
      return(
        data.table(
          Cutoff = cutoff,
          TP = TP,
          FP.mean = mean(FP),
          FP.sd = sd(FP),
          TTP.literal = TP - mean(FP),
          TTP.in.sd = (TP - mean(FP))/sd(FP),
          pvalue = pvalue
        )
      )
    }) |> rbindlist() |>
      mutate(
        Threshold = a,
        .before = everything()
      )
    
    return(fr)
    
  }
  
  range_a_detected = range_a[, mean(TTP.literal), by = Threshold][order(V1, decreasing = TRUE)][1, Threshold]
  
  message("range detection a finished at ", Sys.time())
  
  # b: range for CCLE2019
  range_b = foreach(b = from.range, .combine = rbind, .export = ls(globalenv())) %dopar% {
    library(data.table)
    library(tidyverse)
    # calculate scores for etki and fake etki under threshold b
    # i.e. if we only consider CCLE2019 transformed-rho > b, and assign confidence score based on b
    #      what would be the confidence score for each SLC-metabolite pair?
    ecs = ConfidenceScore_A_by_set(data = etki_Tr, ki = ki.data, a = Inf, b = b, c = Inf, score.only = TRUE)
    fecs = lapply(1:100, function(f){
      simuScore = ConfidenceScore_A_by_set(data = si.data[Iteration == f], ki = ki.data, a = Inf, b = b, c = Inf, score.only = TRUE)
      return(simuScore)
    })
    
    # calculate fraction recovery (TP - FP) for each given score cutoff
    # i.e. how many SLC-metabolite pair can be discovered, if we set score cutoff from 0 (all) - 34 (none)
    fr = lapply(0:34, function(cutoff){
      TP = length(ecs[ecs >= cutoff])/length(ecs)
      FP = sapply(fecs, function(f){length(f[f >= cutoff])}/length(f))
      
      if (TP == 0 & mean(FP) == 0){
        pvalue = 1
      }else if (TP == 1 & mean(FP) == 1) {
        pvalue = 1
      }else{
        pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
      }
      
      return(
        data.table(
          Cutoff = cutoff,
          TP = TP,
          FP.mean = mean(FP),
          FP.sd = sd(FP),
          TTP.literal = TP - mean(FP),
          TTP.in.sd = (TP - mean(FP))/sd(FP),
          pvalue = pvalue
        )
      )
    }) |> rbindlist() |>
      mutate(
        Threshold = b,
        .before = everything()
      )
    
    return(fr)
    
  }
  
  range_b_detected = range_b[, mean(TTP.literal), by = Threshold][order(V1, decreasing = TRUE)][1, Threshold]
  
  message("range detection b finished at ", Sys.time())
  
  # c: range for CCL180
  
  range_c = foreach(c = from.range, .combine = rbind, .export = ls(globalenv())) %dopar% {
    library(data.table)
    library(tidyverse)
    # calculate scores for etki and fake etki under threshold c
    # i.e. if we only consider CCLE2019 transformed-rho > c, and assign confidence score based on c
    #      what would be the confidence score for each SLC-metabolite pair?
    ecs = ConfidenceScore_A_by_set(data = etki_Tr, ki = ki.data, a = Inf, b = Inf, c = c, score.only = TRUE)
    fecs = lapply(1:100, function(f){
      simuScore = ConfidenceScore_A_by_set(data = si.data[Iteration == f], ki = ki.data, a = Inf, b = Inf, c = c, score.only = TRUE)
      return(simuScore)
    })
    
    # calculate fraction recovery (TP - FP) for each given score cutoff
    # i.e. how many SLC-metabolite pair can be discovered, if we set score cutoff from 0 (all) - 34 (none)
    fr = lapply(0:34, function(cutoff){
      TP = length(ecs[ecs >= cutoff])/length(ecs)
      FP = sapply(fecs, function(f){length(f[f >= cutoff])}/length(f))
      
      if (TP == 0 & mean(FP) == 0){
        pvalue = 1
      }else if (TP == 1 & mean(FP) == 1) {
        pvalue = 1
      }else{
        pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
      }
      
      return(
        data.table(
          Cutoff = cutoff,
          TP = TP,
          FP.mean = mean(FP),
          FP.sd = sd(FP),
          TTP.literal = TP - mean(FP),
          TTP.in.sd = (TP - mean(FP))/sd(FP),
          pvalue = pvalue
        )
      )
    }) |> rbindlist() |>
      mutate(
        Threshold = c,
        .before = everything()
      )
    
    return(fr)
    
  }
  
  range_c_detected = range_c[, mean(TTP.literal), by = Threshold][order(V1, decreasing = TRUE)][1, Threshold]
  
  message("range detection c finished at ", Sys.time())
  
  stopCluster(cl)
  
  return(list(a = range_a_detected, b = range_b_detected, c = range_c_detected, a.record = range_a, b.record = range_b, c.record = range_c))
}

ConfA.threshold.detection = Best_threshold_detect_A_by_set(ki.data = etki_Tr, si.data = fake_etki_Tr, from.range = seq(0,1,0.1))

# --------------------------- #
#      Figure S1A data        #
# --------------------------- #

# Fig S1A.csv
ToD = ConfA.threshold.detection$a.record
ggplot() +
  geom_line(data = ToD, mapping = aes(x = Cutoff, y = TTP.literal, group = Threshold), col = "grey") +
  geom_line(data = ToD[Threshold == 0.0], mapping = aes(x = Cutoff, y = TTP.literal, group = Threshold), col = "#BF9000") +
  scale_y_continuous(limits = c(0,0.11), breaks = seq(0,0.1,0.02)) +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  )

# --------------------------- #
#      Figure S1B data        #
# --------------------------- #
# Fig S1B.csv
ToD = ConfA.threshold.detection$b.record
ggplot() +
  geom_line(data = ToD, mapping = aes(x = Cutoff, y = TTP.literal, group = Threshold), col = "grey") +
  geom_line(data = ToD[Threshold == 0.0], mapping = aes(x = Cutoff, y = TTP.literal, group = Threshold), col = "#BF9000") +
  scale_y_continuous(limits = c(0,0.11), breaks = seq(0,0.1,0.02)) +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  )

# --------------------------- #
#      Figure S1C data        #
# --------------------------- #
# Fig S1C.csv
ToD = ConfA.threshold.detection$c.record
ggplot() +
  geom_line(data = ToD, mapping = aes(x = Cutoff, y = TTP.literal, group = Threshold), col = "grey") +
  geom_line(data = ToD[Threshold == 0.0], mapping = aes(x = Cutoff, y = TTP.literal, group = Threshold), col = "#BF9000") +
  scale_y_continuous(limits = c(0,0.11), breaks = seq(0,0.1,0.02)) +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  )

## 7.3 Gene dependency parameter (ToD) optimisation ####
Best_threshold_detect_B_by_set = function(ki.data, si.data, from.range){
  fr = lapply(from.range, function(cutoff){
    TP = nrow(ki.data[CRISPRDEP <= cutoff])/nrow(ki.data)
    FP = sapply(1:100, function(f){nrow(si.data[Iteration == f & CRISPRDEP <= cutoff])/nrow(si.data[Iteration == f])})
    
    if (TP == 0 & mean(FP) == 0){
      pvalue = 1
    }else{
      pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
    }
    
    return(
      data.table(
        Cutoff = cutoff,
        TP = TP,
        FP.mean = mean(FP),
        FP.sd = sd(FP),
        TTP.literal = TP - mean(FP),
        TTP.in.sd = (TP - mean(FP))/sd(FP),
        pvalue = pvalue
      )
    )
  }) |> rbindlist()
  
  message("range detection d finished at ", Sys.time())
  
  # Make decisions based on range report
  range_d_detected = fr[, mean(TTP.literal), by = Cutoff][order(V1, decreasing = TRUE)][1, Cutoff]
  return(list(d = range_d_detected, d.record = fr))
}

ConfB.threshold.detection = Best_threshold_detect_B_by_set(ki.data = etki_Tr, si.data = fake_etki_Tr, from.range = seq(0,0.2,0.01))
# --------------------------- #
#      Figure S1D data        #
# --------------------------- #
# Fig S1D.csv
ToD = ConfB.threshold.detection$d.record
ggplot() +
  geom_point(data = ToD, mapping = aes(x = Cutoff, y = TTP.literal), size = 3, color = "grey") +  
  geom_point(data = ToD[Cutoff == 0.16], mapping = aes(x = Cutoff, y = TTP.literal), size = 3, color = "#990000") +
  scale_y_continuous(limits = c(0,0.04), breaks = seq(0,0.04,0.01)) +
  scale_x_continuous(limits = c(0,0.2), breaks = seq(0,0.2,0.02)) +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  )

## 7.4 Gene dependency parameter (weight b) optimisation ####
ConfidenceScore_B_by_set = function(data, ki, d, e, score.only = FALSE){
  # Factor B quantile
  ns_dist4 = -log10(na.omit(ki$CRISPRDEP))
  q4 = quantile(ns_dist4[ns_dist4 > -log10(d)], seq(0, 1, by = 0.1))
  
  # Calculate scores
  Score = sapply(1:nrow(data), function(i){
    # Factor B: Dependency Scores #
    if (is.na(data[i, CRISPRDEP]) | data[i, -log10(as.numeric(CRISPRDEP))] < -log10(d)){
      Factor_B = 0
    } else if (data[i, -log10(as.numeric(CRISPRDEP))] > -log10(d)){
      Factor_B = findInterval(data[i, -log10(as.numeric(CRISPRDEP))], q4)*as.numeric(e)
    }
    return(as.numeric(Factor_B))
  })
  
  if (score.only){
    return(Score)
  }else{
    score_table = data |>
      mutate(
        Score = Score
      )
    return(score_table)
  }
}

Best_scalar_detect_B_by_set = function(ki.data, si.data, d, from.range){
  # e: range of scalar reward for CRISPR dependency
  range_e = lapply(from.range, function(e){
    # calculate scores for etki and fake etki given scalar e
    # i.e. if we only consider CRISPR dependency pvalue < d, and assign confidence score scaled by e
    #      what would be the confidence score for each SLC-metabolite pair?
    ecs = ConfidenceScore_B_by_set(data = etki_Tr, ki = ki.data, d = d, e = e, score.only = TRUE)
    fecs = lapply(1:100, function(f){
      simuScore = ConfidenceScore_B_by_set(data = si.data[Iteration == f], ki = ki.data, d = d, e = e, score.only = TRUE)
      return(simuScore)
    })
    
    # calculate fraction recovery (TP - FP) for each given score cutoff
    # i.e. how many SLC-metabolite pair can be discovered, if we set score cutoff from 0 (all) - 12 (none)
    fr = lapply(0:12, function(cutoff){
      TP = length(ecs[ecs >= cutoff])/length(ecs)
      FP = sapply(fecs, function(f){length(f[f >= cutoff])}/length(f))
      
      if (TP == 0 & mean(FP) == 0){
        pvalue = 1
      }else if (TP == 1 & mean(FP) == 1) {
        pvalue = 1
      }else{
        pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
      }
      
      return(
        data.table(
          Cutoff = cutoff,
          TP = TP,
          FP.mean = mean(FP),
          FP.sd = sd(FP),
          TTP.literal = TP - mean(FP),
          TTP.in.sd = (TP - mean(FP))/sd(FP),
          pvalue = pvalue
        )
      )
    }) |> rbindlist() |>
      mutate(
        Scalar = e,
        .before = everything()
      )
    
    return(fr)
    
  }) |> rbindlist()
  
  message("range detection e finished at ", Sys.time())
  
  # Make decisions based on e range report
  e.report = range_e[, mean(TTP.literal), by = Scalar][order(V1, decreasing = TRUE)][1, Scalar]
  
  return(list(e = e.report, e.record = range_e))
}

ConfB.scalar.detection = Best_scalar_detect_B_by_set(ki.data = etki_Tr, si.data = fake_etki_Tr, d = as.numeric(ConfB.threshold.detection$d), from.range = seq(0.1, 1, 0.1))

# --------------------------- #
#      Figure S1E data        #
# --------------------------- #
parameter_b = ConfB.scalar.detection$e.record
ggplot() +
  geom_line(data = Scalar_b, mapping = aes(x = Cutoff, y = TTP.literal, group = Scalar), col = "grey") +
  geom_line(data = Scalar_e[Scalar == 1.0], mapping = aes(x = Cutoff, y = TTP.literal, group = Scalar), col = "#990000") +
  scale_y_continuous(limits = c(0,0.11), breaks = seq(0,0.1,0.02)) +
  scale_x_continuous(limits = c(0,12), breaks = seq(0,12,2)) +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  )

## 7.5 Unit of adjacency optimisation ####
ConfidenceScore_ABC_adj1 = function(data, ki = etki_Tr, score.only = FALSE){
  ## constant ##
  ## loading final constant ##
  constant = c(
    a = 0.00, b = 0.00, c = 1.0, 
    d = 0.16, e = 1.00,
    o = 0.00, p = 0.00, q = 0.00, m = 1.00
  )
  
  ## quantile ##
  ns_dist1 = na.omit(ki$NCI60) # NCI60 #
  ns_dist2 = na.omit(ki$CCLE2019) # CCLE2019 #
  ns_dist3 = na.omit(ki$CCL180) # CCL180 #
  ns_dist4 = -log10(na.omit(ki$CRISPRDEP)) # CRISPR dependency
  q1 = quantile(ns_dist1[ns_dist1 > as.numeric(constant["a"])], seq(0, 1, by = 0.1))
  q2 = quantile(ns_dist2[ns_dist2 > as.numeric(constant["b"])], seq(0, 1, by = 0.1))
  q3 = quantile(ns_dist3[ns_dist3 > as.numeric(constant["c"])], seq(0, 1, by = 0.1))
  q4 = quantile(ns_dist4[ns_dist4 > -log10(constant["d"])], seq(0, 1, by = 0.1))
  q5 = quantile(ns_dist1[ns_dist1 > as.numeric(constant["o"])], seq(0, 1, by = 0.1))
  q6 = quantile(ns_dist2[ns_dist1 > as.numeric(constant["p"])], seq(0, 1, by = 0.1))
  q7 = quantile(ns_dist3[ns_dist1 > as.numeric(constant["q"])], seq(0, 1, by = 0.1))
  
  data_info = data.table(
    data_name = c("NCI60", "CCLE2019", "CCL180"),
    FactorA_dist = c("q1", "q2", "q3"),
    FactorC_dist = c("q5", "q6", "q7"),
    FactorA_threshold = c(as.numeric(constant["a"]), as.numeric(constant["b"]), as.numeric(constant["c"])),
    FactorC_threshold = c(as.numeric(constant["o"]), as.numeric(constant["p"]), as.numeric(constant["q"]))
  )
  
  cl = makeCluster(8)
  registerDoParallel(cl)
  
  Score = foreach(i = 1:nrow(data), .combine = c, .export = unique(c(ls(globalenv()), ls(environment())))) %dopar% {
    library(data.table)
    library(tidyverse)
    
    # Factor A #
    Factor_A = vector("integer", nrow(data_info))
    for (j in 1:nrow(data_info)){
      n = as.numeric(data[i, get(data_info[j, data_name])])
      if (is.na(n) | n < data_info[j, FactorA_threshold]){
        Factor_A[j] = 0
      } else if (n > data_info[j, FactorA_threshold]){
        Factor_A[j] = findInterval(n, get(data_info[j, FactorA_dist]))*3
      }
    }
    Factor_A = as.numeric(sum(Factor_A))
    
    # Factor B #
    if (is.na(data[i, CRISPRDEP]) | data[i, -log10(as.numeric(CRISPRDEP))] < -log10(constant["d"])){
      Factor_B = 0
    } else if (data[i, -log10(as.numeric(CRISPRDEP))] > -log10(constant["d"])){
      Factor_B = findInterval(data[i, -log10(as.numeric(CRISPRDEP))], q4)*as.numeric(constant["e"])
    }
    
    # Factor C #
    if (is.na(data$Formula[i])){
      temp = getAM(data$KEGG_ID[i], is.kegg = TRUE)
    }else{
      temp = getAM(data$Formula[i])
    }
    if (check_NA(temp)){
      Factor_C = 0
    }else if (!1 %in% names(temp)){
      Factor_C = 0
    }else{
      proximal = translate(temp[which(names(temp) %in% 1)])
      if (check_NA(proximal)){
        Factor_C = 0
      }else{
        proximal_ns = getNR(data[i, SYMBOL], proximal)
        if (nrow(proximal_ns) != 0) {
          Factor_C = 0
          for (l in 1:nrow(data_info)){
            N = unique(proximal_ns[, get(data_info[l, data_name])])
            for (n in N){
              if (is.na(n) | n < data_info[l, FactorC_threshold]){
                Factor_C = Factor_C + 0
              }else if (n > data_info[l, FactorC_threshold] & n < data_info[l, FactorA_threshold]){
                Factor_C = Factor_C + 0.5*constant["m"]
              }else if (n > data_info[l, FactorA_threshold]){
                Factor_C = Factor_C + findInterval(n, get(data_info[l, FactorC_dist]))*constant["m"]
              }
            }
            print(Factor_C)
          }
        }else{
          Factor_C = 0
        }
      }
    }
    
    # Final score #
    #finalscore = Factor_A + Factor_B + Factor_C
    finalscore = Factor_C
    return(finalscore)
  }
  
  stopCluster(cl)
  
  if (score.only){
    return(as.numeric(Score))
  }else{
    score_table = data |>
      mutate(
        Score = as.numeric(Score)
      )
    return(score_table)
  }
}
ConfidenceScore_ABC_adj2 = function(data, ki = etki_Tr, score.only = FALSE){
  ## constant ##
  ## loading final constant ##
  constant = c(
    a = 0.00, b = 0.00, c = 1.0, 
    d = 0.16, e = 1.00,
    o = 0.00, p = 0.00, q = 0.00, m = 1.00
  )
  
  ## quantile ##
  ns_dist1 = na.omit(ki$NCI60) # NCI60 #
  ns_dist2 = na.omit(ki$CCLE2019) # CCLE2019 #
  ns_dist3 = na.omit(ki$CCL180) # CCL180 #
  ns_dist4 = -log10(na.omit(ki$CRISPRDEP)) # CRISPR dependency
  q1 = quantile(ns_dist1[ns_dist1 > as.numeric(constant["a"])], seq(0, 1, by = 0.1))
  q2 = quantile(ns_dist2[ns_dist2 > as.numeric(constant["b"])], seq(0, 1, by = 0.1))
  q3 = quantile(ns_dist3[ns_dist3 > as.numeric(constant["c"])], seq(0, 1, by = 0.1))
  q4 = quantile(ns_dist4[ns_dist4 > -log10(constant["d"])], seq(0, 1, by = 0.1))
  q5 = quantile(ns_dist1[ns_dist1 > as.numeric(constant["o"])], seq(0, 1, by = 0.1))
  q6 = quantile(ns_dist2[ns_dist1 > as.numeric(constant["p"])], seq(0, 1, by = 0.1))
  q7 = quantile(ns_dist3[ns_dist1 > as.numeric(constant["q"])], seq(0, 1, by = 0.1))
  
  data_info = data.table(
    data_name = c("NCI60", "CCLE2019", "CCL180"),
    FactorA_dist = c("q1", "q2", "q3"),
    FactorC_dist = c("q5", "q6", "q7"),
    FactorA_threshold = c(as.numeric(constant["a"]), as.numeric(constant["b"]), as.numeric(constant["c"])),
    FactorC_threshold = c(as.numeric(constant["o"]), as.numeric(constant["p"]), as.numeric(constant["q"]))
  )
  
  cl = makeCluster(8)
  registerDoParallel(cl)
  
  Score = foreach(i = 1:nrow(data), .combine = c, .export = unique(c(ls(globalenv()), ls(environment())))) %dopar% {
    library(data.table)
    library(tidyverse)
    
    # Factor A #
    Factor_A = vector("integer", nrow(data_info))
    for (j in 1:nrow(data_info)){
      n = as.numeric(data[i, get(data_info[j, data_name])])
      if (is.na(n) | n < data_info[j, FactorA_threshold]){
        Factor_A[j] = 0
      } else if (n > data_info[j, FactorA_threshold]){
        Factor_A[j] = findInterval(n, get(data_info[j, FactorA_dist]))*3
      }
    }
    Factor_A = as.numeric(sum(Factor_A))
    
    # Factor B #
    if (is.na(data[i, CRISPRDEP]) | data[i, -log10(as.numeric(CRISPRDEP))] < -log10(constant["d"])){
      Factor_B = 0
    } else if (data[i, -log10(as.numeric(CRISPRDEP))] > -log10(constant["d"])){
      Factor_B = findInterval(data[i, -log10(as.numeric(CRISPRDEP))], q4)*as.numeric(constant["e"])
    }
    
    # Factor C #
    if (is.na(data$Formula[i])){
      temp = getAM(data$KEGG_ID[i], is.kegg = TRUE)
    }else{
      temp = getAM(data$Formula[i])
    }
    if (check_NA(temp)){
      Factor_C = 0
    }else if (!1 %in% names(temp) & !2 %in% names(temp)){
      Factor_C = 0
    }else{
      proximal = translate(temp[which(names(temp) %in% c(1,2))])
      if (check_NA(proximal)){
        Factor_C = 0
      }else{
        proximal_ns = getNR(data[i, SYMBOL], proximal)
        if (nrow(proximal_ns) != 0) {
          Factor_C = 0
          for (l in 1:nrow(data_info)){
            N = unique(proximal_ns[, get(data_info[l, data_name])])
            for (n in N){
              if (is.na(n) | n < data_info[l, FactorC_threshold]){
                Factor_C = Factor_C + 0
              }else if (n > data_info[l, FactorC_threshold] & n < data_info[l, FactorA_threshold]){
                Factor_C = Factor_C + 0.5*constant["m"]
              }else if (n > data_info[l, FactorA_threshold]){
                Factor_C = Factor_C + findInterval(n, get(data_info[l, FactorC_dist]))*constant["m"]
              }
            }
            print(Factor_C)
          }
        }else{
          Factor_C = 0
        }
      }
    }
    
    # Final score #
    #finalscore = Factor_A + Factor_B + Factor_C
    finalscore = Factor_C
    return(finalscore)
  }
  
  stopCluster(cl)
  
  if (score.only){
    return(as.numeric(Score))
  }else{
    score_table = data |>
      mutate(
        Score = as.numeric(Score)
      )
    return(score_table)
  }
}
ConfidenceScore_ABC_adj3 = function(data, ki = etki_Tr, score.only = FALSE){
  ## constant ##
  ## loading final constant ##
  constant = c(
    a = 0.00, b = 0.00, c = 1.0, 
    d = 0.16, e = 1.00,
    o = 0.00, p = 0.00, q = 0.00, m = 1.00
  )
  
  ## quantile ##
  ns_dist1 = na.omit(ki$NCI60) # NCI60 #
  ns_dist2 = na.omit(ki$CCLE2019) # CCLE2019 #
  ns_dist3 = na.omit(ki$CCL180) # CCL180 #
  ns_dist4 = -log10(na.omit(ki$CRISPRDEP)) # CRISPR dependency
  q1 = quantile(ns_dist1[ns_dist1 > as.numeric(constant["a"])], seq(0, 1, by = 0.1))
  q2 = quantile(ns_dist2[ns_dist2 > as.numeric(constant["b"])], seq(0, 1, by = 0.1))
  q3 = quantile(ns_dist3[ns_dist3 > as.numeric(constant["c"])], seq(0, 1, by = 0.1))
  q4 = quantile(ns_dist4[ns_dist4 > -log10(constant["d"])], seq(0, 1, by = 0.1))
  q5 = quantile(ns_dist1[ns_dist1 > as.numeric(constant["o"])], seq(0, 1, by = 0.1))
  q6 = quantile(ns_dist2[ns_dist1 > as.numeric(constant["p"])], seq(0, 1, by = 0.1))
  q7 = quantile(ns_dist3[ns_dist1 > as.numeric(constant["q"])], seq(0, 1, by = 0.1))
  
  data_info = data.table(
    data_name = c("NCI60", "CCLE2019", "CCL180"),
    FactorA_dist = c("q1", "q2", "q3"),
    FactorC_dist = c("q5", "q6", "q7"),
    FactorA_threshold = c(as.numeric(constant["a"]), as.numeric(constant["b"]), as.numeric(constant["c"])),
    FactorC_threshold = c(as.numeric(constant["o"]), as.numeric(constant["p"]), as.numeric(constant["q"]))
  )
  
  cl = makeCluster(8)
  registerDoParallel(cl)
  
  Score = foreach(i = 1:nrow(data), .combine = c, .export = unique(c(ls(globalenv()), ls(environment())))) %dopar% {
    library(data.table)
    library(tidyverse)
    
    # Factor A #
    Factor_A = vector("integer", nrow(data_info))
    for (j in 1:nrow(data_info)){
      n = as.numeric(data[i, get(data_info[j, data_name])])
      if (is.na(n) | n < data_info[j, FactorA_threshold]){
        Factor_A[j] = 0
      } else if (n > data_info[j, FactorA_threshold]){
        Factor_A[j] = findInterval(n, get(data_info[j, FactorA_dist]))*3
      }
    }
    Factor_A = as.numeric(sum(Factor_A))
    
    # Factor B #
    if (is.na(data[i, CRISPRDEP]) | data[i, -log10(as.numeric(CRISPRDEP))] < -log10(constant["d"])){
      Factor_B = 0
    } else if (data[i, -log10(as.numeric(CRISPRDEP))] > -log10(constant["d"])){
      Factor_B = findInterval(data[i, -log10(as.numeric(CRISPRDEP))], q4)*as.numeric(constant["e"])
    }
    
    # Factor C #
    if (is.na(data$Formula[i])){
      temp = getAM(data$KEGG_ID[i], is.kegg = TRUE)
    }else{
      temp = getAM(data$Formula[i])
    }
    if (check_NA(temp)){
      Factor_C = 0
    }else if (!1 %in% names(temp) & !2 %in% names(temp) & !3 %in% names(temp)){
      Factor_C = 0
    }else{
      proximal = translate(temp[which(names(temp) %in% c(1,2,3))])
      if (check_NA(proximal)){
        Factor_C = 0
      }else{
        proximal_ns = getNR(data[i, SYMBOL], proximal)
        if (nrow(proximal_ns) != 0) {
          Factor_C = 0
          for (l in 1:nrow(data_info)){
            N = unique(proximal_ns[, get(data_info[l, data_name])])
            for (n in N){
              if (is.na(n) | n < data_info[l, FactorC_threshold]){
                Factor_C = Factor_C + 0
              }else if (n > data_info[l, FactorC_threshold] & n < data_info[l, FactorA_threshold]){
                Factor_C = Factor_C + 0.5*constant["m"]
              }else if (n > data_info[l, FactorA_threshold]){
                Factor_C = Factor_C + findInterval(n, get(data_info[l, FactorC_dist]))*constant["m"]
              }
            }
            print(Factor_C)
          }
        }else{
          Factor_C = 0
        }
      }
    }
    
    # Final score #
    # finalscore = Factor_A + Factor_B + Factor_C
    finalscore = Factor_C
    return(finalscore)
  }
  
  stopCluster(cl)
  
  if (score.only){
    return(as.numeric(Score))
  }else{
    score_table = data |>
      mutate(
        Score = as.numeric(Score)
      )
    return(score_table)
  }
}

ConfABC_adj1 = ConfidenceScore_ABC_adj1(etki_Tr)
ConfABC_adj2 = ConfidenceScore_ABC_adj2(etki_Tr)
ConfABC_adj3 = ConfidenceScore_ABC_adj3(etki_Tr)

simulated_ConfABC_adj1 = lapply(fake_etki_Tr, ConfidenceScore_ABC_adj1)
simulated_ConfABC_adj2 = lapply(fake_etki_Tr, ConfidenceScore_ABC_adj2)
simulated_ConfABC_adj3 = lapply(fake_etki_Tr, ConfidenceScore_ABC_adj3)

Fr = lapply(seq(0,1,0.001), function(Th){
  Th1 = ConfABC_adj1[, max(as.numeric(Score))] * Th
  Th2 = ConfABC_adj2[, max(as.numeric(Score))] * Th
  Th3 = ConfABC_adj3[, max(as.numeric(Score))] * Th
  Adj1 = nrow(ConfABC_adj1[as.numeric(Score) > Th1])/nrow(ConfABC_adj1) - sapply(simulated_ConfABC_adj1, function(f){nrow(f[as.numeric(Score) > Th1])/nrow(f)})
  Adj2 = nrow(ConfABC_adj2[as.numeric(Score) > Th2])/nrow(ConfABC_adj2) - sapply(simulated_ConfABC_adj2, function(f){nrow(f[as.numeric(Score) > Th2])/nrow(f)})
  Adj3 = nrow(ConfABC_adj3[as.numeric(Score) > Th3])/nrow(ConfABC_adj3) - sapply(simulated_ConfABC_adj3, function(f){nrow(f[as.numeric(Score) > Th3])/nrow(f)})
  message(Th)
  return(
    data.table(
      Th = Th,
      Class = rep(c("Adj1", "Adj2", "Adj3"), each = 100),
      Fr = c(Adj1, Adj2, Adj3)
    )
  )
}) |> rbindlist()

# --------------------------- #
#      Figure S2B data        #
# --------------------------- #
Fr_mean = Fr[, mean(Fr), by = .(Th, Class)][Th <= 0.30]
ggplot(Fr_mean, aes(x = Th, y = V1, color = Class)) + 
  geom_line() +
  scale_color_manual(values=c("lightgreen", "#990000", "lightblue")) +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  )

## 7.6 Metabolite adjacency parameter (ToD) optimisation ####
ConfidenceScore_ABC_by_set = function(data, ki, o, p, q, m, score.only = FALSE){
  
  constant = c(
    # Threshold for ConfA, scaled by 3
    a = ConfA.threshold.detection$a,
    b = ConfA.threshold.detection$b,
    c = ConfA.threshold.detection$c,
    # Threshold for ConfB
    d = ConfB.threshold.detection$d,
    # Scalar for ConfB
    e = ConfB.scalar.detection$e,
    # Threshold for ConfC
    o = o,
    p = p,
    q = q,
    # Scalar for ConfC
    m = m
  )
  
  ## quantile ##
  ns_dist1 = na.omit(ki$NCI60) # NCI60 #
  ns_dist2 = na.omit(ki$CCLE2019) # CCLE2019 #
  ns_dist3 = na.omit(ki$CCL180) # CCL180 #
  ns_dist4 = -log10(na.omit(ki$CRISPRDEP)) # CRISPR dependency
  q1 = quantile(ns_dist1[ns_dist1 > as.numeric(constant["a"])], seq(0, 1, by = 0.1))
  q2 = quantile(ns_dist2[ns_dist2 > as.numeric(constant["b"])], seq(0, 1, by = 0.1))
  q3 = quantile(ns_dist3[ns_dist3 > as.numeric(constant["c"])], seq(0, 1, by = 0.1))
  q4 = quantile(ns_dist4[ns_dist4 > -log10(constant["d"])], seq(0, 1, by = 0.1))
  q5 = quantile(ns_dist1[ns_dist1 > as.numeric(constant["o"])], seq(0, 1, by = 0.1))
  q6 = quantile(ns_dist2[ns_dist1 > as.numeric(constant["p"])], seq(0, 1, by = 0.1))
  q7 = quantile(ns_dist3[ns_dist1 > as.numeric(constant["q"])], seq(0, 1, by = 0.1))
  
  data_info = data.table(
    data_name = c("NCI60", "CCLE2019", "CCL180"),
    FactorA_dist = c("q1", "q2", "q3"),
    FactorC_dist = c("q5", "q6", "q7"),
    FactorA_threshold = c(as.numeric(constant["a"]), as.numeric(constant["b"]), as.numeric(constant["c"])),
    FactorC_threshold = c(as.numeric(constant["o"]), as.numeric(constant["p"]), as.numeric(constant["q"]))
  )
  
  Score = lapply(1:nrow(data), function(i){
    library(data.table)
    library(tidyverse)
    
    # Factor A #
    Factor_A = vector("integer", nrow(data_info))
    for (j in 1:nrow(data_info)){
      n = as.numeric(data[i, get(data_info[j, data_name])])
      if (is.na(n) | n < data_info[j, FactorA_threshold]){
        Factor_A[j] = 0
      } else if (n > data_info[j, FactorA_threshold]){
        Factor_A[j] = findInterval(n, get(data_info[j, FactorA_dist]))*3
      }
    }
    Factor_A = as.numeric(sum(Factor_A))
    
    # Factor B #
    if (is.na(data[i, CRISPRDEP]) | data[i, -log10(as.numeric(CRISPRDEP))] < -log10(constant["d"])){
      Factor_B = 0
    } else if (data[i, -log10(as.numeric(CRISPRDEP))] > -log10(constant["d"])){
      Factor_B = findInterval(data[i, -log10(as.numeric(CRISPRDEP))], q4)*as.numeric(constant["e"])
    }
    
    # Factor C #
    if (is.na(data$Formula[i])){
      temp = getAM(data$KEGG_ID[i], is.kegg = TRUE)
    }else{
      temp = getAM(data$Formula[i])
    }
    if (check_NA(temp)){
      Factor_C = 0
    }else if (!1 %in% names(temp) & !2 %in% names(temp)){
      Factor_C = 0
    }else{
      proximal = translate(temp[which(names(temp) %in% c(1,2))])
      if (check_NA(proximal)){
        Factor_C = 0
      }else{
        proximal_ns = getNR(data[i, SYMBOL], proximal)
        if (nrow(proximal_ns) != 0) {
          Factor_C = 0
          for (l in 1:nrow(data_info)){
            N = unique(proximal_ns[, get(data_info[l, data_name])])
            for (n in N){
              if (is.na(n) | n < data_info[l, FactorC_threshold]){
                Factor_C = Factor_C + 0
              }else if (n > data_info[l, FactorC_threshold] & n < data_info[l, FactorA_threshold]){
                Factor_C = Factor_C + constant["m"]/2
              }else if (n > data_info[l, FactorA_threshold]){
                Factor_C = Factor_C + findInterval(n, get(data_info[l, FactorA_dist]))*constant["m"]
              }
            }
          }
        }else{
          Factor_C = 0
        }
      }
    }
    
    # Final score #
    finalscore = Factor_A + Factor_B + Factor_C
    return(finalscore)
  }) |> unlist()
  
  if (score.only){
    return(as.numeric(Score))
  }else{
    score_table = data |>
      mutate(
        Score = as.numeric(Score)
      )
    return(score_table)
  }
}

Best_threshold_detect_C_by_set = function(ki.data, si.data, from.range){
  # Setting up parallel
  numCores = 8  # Use all available cores except one for other processes
  cl = makeCluster(numCores)
  
  # o: derivative range for NCI60
  registerDoParallel(cl)
  range_o = foreach(o = from.range, .combine = rbind, .export = ls(globalenv())) %dopar% {
    library(data.table)
    library(tidyverse)
    # calculate scores for etki and fake etki under threshold o
    # i.e. if we only consider NCI60 derivative transformed-rho > o, and assign confidence score based on o
    #      what would be the confidence score for each SLC-metabolite pair?
    ecs = ConfidenceScore_ABC_by_set(data = etki_Tr, ki = ki.data, o = o, p = Inf, q = Inf, m = 1, score.only = TRUE)
    
    fecs = lapply(1:100, function(f){
      simuScore = ConfidenceScore_ABC_by_set(data = si.data[Iteration == f], ki = ki.data, o = o, p = Inf, q = Inf, m = 1, score.only = TRUE)
      message(f)
      return(simuScore)
    })
    
    # calculate fraction recovery (TP - FP) for each given score cutoff
    # i.e. how many SLC-metabolite pair can be discovered, if we set score cutoff from 0 (all) - 60 (none)
    fr = lapply(0:150, function(cutoff){
      TP = length(ecs[ecs >= cutoff])/length(ecs)
      FP = sapply(fecs, function(f){length(f[f >= cutoff])}/length(f))
      
      if (TP == 0 & mean(FP) == 0){
        pvalue = 1
      }else if (TP == 1 & mean(FP) == 1) {
        pvalue = 1
      }else{
        pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
      }
      
      return(
        data.table(
          Cutoff = cutoff,
          TP = TP,
          FP.mean = mean(FP),
          FP.sd = sd(FP),
          TTP.literal = TP - mean(FP),
          TTP.in.sd = (TP - mean(FP))/sd(FP),
          pvalue = pvalue
        )
      )
    }) |> rbindlist() |>
      mutate(
        Threshold = o,
        .before = everything()
      )
    
    return(fr)
    
  }
  message("Range o finished at ", Sys.time())
  
  # p: derivative range for CCLE2019
  range_p = foreach(p = from.range, .combine = rbind, .export = ls(globalenv())) %dopar% {
    library(data.table)
    library(tidyverse)
    # calculate scores for etki and fake etki under threshold p
    # i.e. if we only consider CCLE2019 transformed-rho > p, and assign confidence score based on p
    #      what would be the confidence score for each SLC-metabolite pair?
    ecs = ConfidenceScore_ABC_by_set(data = etki_Tr, ki = ki.data, o = Inf, p = p, q = Inf, m = 1, score.only = TRUE)
    fecs = lapply(1:100, function(f){
      simuScore = ConfidenceScore_ABC_by_set(data = si.data[Iteration == f], ki = ki.data, o = Inf, p = p, q = Inf, m = 1, score.only = TRUE)
      return(simuScore)
    })
    
    # calculate fraction recovery (TP - FP) for each given score cutoff
    # i.e. how many SLC-metabolite pair can be discovered, if we set score cutoff from 0 (all) - 34 (none)
    fr = lapply(0:150, function(cutoff){
      TP = length(ecs[ecs >= cutoff])/length(ecs)
      FP = sapply(fecs, function(f){length(f[f >= cutoff])}/length(f))
      
      if (TP == 0 & mean(FP) == 0){
        pvalue = 1
      }else if (TP == 1 & mean(FP) == 1) {
        pvalue = 1
      }else{
        pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
      }
      
      return(
        data.table(
          Cutoff = cutoff,
          TP = TP,
          FP.mean = mean(FP),
          FP.sd = sd(FP),
          TTP.literal = TP - mean(FP),
          TTP.in.sd = (TP - mean(FP))/sd(FP),
          pvalue = pvalue
        )
      )
    }) |> rbindlist() |>
      mutate(
        Threshold = p,
        .before = everything()
      )
    
    return(fr)
    
  }
  message("Range p finished at ", Sys.time())
  
  # q: derivative range for CCL180
  range_q = foreach(q = from.range, .combine = rbind, .export = ls(globalenv())) %dopar% {
    library(data.table)
    library(tidyverse)
    # calculate scores for etki and fake etki under threshold c
    # i.e. if we only consider CCLE2019 transformed-rho > c, and assign confidence score based on c
    #      what would be the confidence score for each SLC-metabolite pair?
    ecs = ConfidenceScore_ABC_by_set(data = etki_Tr, ki = ki.data, o = Inf, p = Inf, q = q, m = 1, score.only = TRUE)
    fecs = lapply(1:100, function(f){
      simuScore = ConfidenceScore_ABC_by_set(data = si.data[Iteration == f], ki = ki.data, o = Inf, p = Inf, q = q, m = 1, score.only = TRUE)
      return(simuScore)
    })
    
    # calculate fraction recovery (TP - FP) for each given score cutoff
    # i.e. how many SLC-metabolite pair can be discovered, if we set score cutoff from 0 (all) - 34 (none)
    fr = lapply(0:150, function(cutoff){
      TP = length(ecs[ecs >= cutoff])/length(ecs)
      FP = sapply(fecs, function(f){length(f[f >= cutoff])}/length(f))
      
      if (TP == 0 & mean(FP) == 0){
        pvalue = 1
      }else if (TP == 1 & mean(FP) == 1) {
        pvalue = 1
      }else{
        pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
      }
      
      return(
        data.table(
          Cutoff = cutoff,
          TP = TP,
          FP.mean = mean(FP),
          FP.sd = sd(FP),
          TTP.literal = TP - mean(FP),
          TTP.in.sd = (TP - mean(FP))/sd(FP),
          pvalue = pvalue
        )
      )
    }) |> rbindlist() |>
      mutate(
        Threshold = q,
        .before = everything()
      )
    
    return(fr)
    
  }
  
  message("Range q finished at ", Sys.time())
  
  stopCluster(cl)
  
  # Make decisions on o, p, and q based on range report
  o.report = range_o[, mean(TTP.literal), by = Threshold][order(V1, decreasing = TRUE)][1, Threshold]
  p.report = range_p[, mean(TTP.literal), by = Threshold][order(V1, decreasing = TRUE)][1, Threshold]
  q.report = range_q[, mean(TTP.literal), by = Threshold][order(V1, decreasing = TRUE)][1, Threshold]
  
  return(list(o = o.report, p = p.report, q = q.report, o.record = range_o, p.record = range_p, q.record = range_q))
}

ConfC.threshold.detection = Best_threshold_detect_C_by_set(ki.data = etki_Tr, si.data = fake_etki_Tr, from.range = seq(0,1,0.1))

# --------------------------- #
#      Figure S1F data        #
# --------------------------- #
ToD = ConfC.threshold.detection$o.record
ggplot() +
  geom_line(data = ToD, mapping = aes(x = Cutoff, y = TTP.literal, group = Threshold), col = "grey") +
  geom_line(data = ToD[Threshold == 0.1], mapping = aes(x = Cutoff, y = TTP.literal, group = Threshold), col = "#BF9000") +
  scale_y_continuous(limits = c(-0.01,0.15), breaks = seq(0,0.2,0.05)) +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  )

# --------------------------- #
#      Figure S1G data        #
# --------------------------- #
ToD = ConfC.threshold.detection$p.record
ggplot() +
  geom_line(data = ToD, mapping = aes(x = Cutoff, y = TTP.literal, group = Threshold), col = "grey") +
  geom_line(data = ToD[Threshold == 0.1], mapping = aes(x = Cutoff, y = TTP.literal, group = Threshold), col = "#8ABFBD") +
  scale_y_continuous(limits = c(-0.01,0.15), breaks = seq(0,0.2,0.05)) +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  )

# --------------------------- #
#      Figure S1H data        #
# --------------------------- #
ToD = ConfC.threshold.detection$q.record
ggplot() +
  geom_line(data = ToD, mapping = aes(x = Cutoff, y = TTP.literal, group = Threshold), col = "grey") +
  geom_line(data = ToD[Threshold == 0.2], mapping = aes(x = Cutoff, y = TTP.literal, group = Threshold), col = "#987B9F") +
  scale_y_continuous(limits = c(-0.01,0.15), breaks = seq(0,0.2,0.05)) +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  )

## 7.7 Metabolite adjacency parameter (weight c) optimisation ####
Best_scalar_detect_C_by_set = function(ki.data, si.data, o, p, q, from.range){
  # Setting up parallel
  numCores = 8  # Use all available cores except one for other processes
  cl = makeCluster(numCores)
  
  # m: range of scalar reward for Adjacency
  registerDoParallel(cl)
  range_m = foreach(m = from.range, .combine = rbind, .export = ls(globalenv())) %dopar% {
    library(data.table)
    library(tidyverse)
    # calculate scores for etki and fake etki given scalar m
    # i.e. if we only consider transformed-rho > o, p, q, and assign confidence score scaled by m
    #      what would be the confidence score for each SLC-metabolite pair?
    ecs = ConfidenceScore_ABC_by_set(data = etki_Tr, ki = ki.data, o = o, p = p, q = q, m = m, score.only = TRUE)
    fecs = lapply(1:100, function(f){
      simuScore = ConfidenceScore_ABC_by_set(data = si.data[Iteration == f], ki = ki.data, o = o, p = p, q = q, m = m, score.only = TRUE)
      return(simuScore)
    })
    
    # calculate fraction recovery (TP - FP) for each given score cutoff
    # i.e. how many SLC-metabolite pair can be discovered, if we set score cutoff from 0 (all) - 12 (none)
    fr = lapply(0:150, function(cutoff){
      TP = length(ecs[ecs >= cutoff])/length(ecs)
      FP = sapply(fecs, function(f){length(f[f >= cutoff])}/length(f))
      
      if (TP == 0 & mean(FP) == 0){
        pvalue = 1
      }else if (TP == 1 & mean(FP) == 1) {
        pvalue = 1
      }else{
        pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
      }
      
      return(
        data.table(
          Cutoff = cutoff,
          TP = TP,
          FP.mean = mean(FP),
          FP.sd = sd(FP),
          TTP.literal = TP - mean(FP),
          TTP.in.sd = (TP - mean(FP))/sd(FP),
          pvalue = pvalue
        )
      )
    }) |> rbindlist() |>
      mutate(
        Scalar = m,
        .before = everything()
      )
    
    return(fr)
    
  }
  
  message("Range m finished at ", Sys.time())
  
  stopCluster(cl)
  
  # Make decisions based on m range report
  range_m_detected = range_m[, mean(TTP.literal), by = Scalar][order(V1, decreasing = TRUE)][1, Scalar]
  
  return(list(m = range_m_detected, m.record = range_m))
}

ConfC.scalar.detection = Best_scalar_detect_C_by_set(ki.data = etki_Tr, si.data = fake_etki_Tr,
                                                     o = ConfC.threshold.detection$o,
                                                     p = ConfC.threshold.detection$p,
                                                     q = ConfC.threshold.detection$q,
                                                     from.range = seq(0,1,0.1))
# --------------------------- #
#      Figure S1I data        #
# --------------------------- #
ToD = ConfC.threshold.detection$m.record
ggplot() +
  geom_line(data = Scalar_c, mapping = aes(x = Cutoff, y = TTP.literal, group = Scalar), col = "grey") +
  geom_line(data = Scalar_c[Scalar == 1.0], mapping = aes(x = Cutoff, y = TTP.literal, group = Scalar), col = "#990000") +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  )

# 8 Compare confidence score outcomes ####
## 8.1 Define ConfidenceScore ####
ConfidenceScore_A = function(data, ki = etki, score.only = FALSE){
  ## loading final constant ##
  constant = c(
    a = 0.00, b = 0.00, c = 1.0, 
    d = 0.16, e = 1.00,
    o = 0.10, p = 0.10, q = 0.20, m = 1.00
  )
  
  ## quantile ##
  ns_dist1 = na.omit(ki$NCI60) # NCI60 #
  ns_dist2 = na.omit(ki$CCLE2019) # CCLE2019 #
  ns_dist3 = na.omit(ki$CCL180) # CCL180 #
  # ns_dist4 = -log10(na.omit(ki$CRISPRDEP)) # CRISPR dependency
  q1 = quantile(ns_dist1[ns_dist1 > as.numeric(constant["a"])], seq(0, 1, by = 0.1))
  q2 = quantile(ns_dist2[ns_dist2 > as.numeric(constant["b"])], seq(0, 1, by = 0.1))
  q3 = quantile(ns_dist3[ns_dist3 > as.numeric(constant["c"])], seq(0, 1, by = 0.1))
  # q4 = quantile(ns_dist4[ns_dist4 > -log10(constant["d"])], seq(0, 1, by = 0.1))
  # q5 = quantile(ns_dist1[ns_dist1 > as.numeric(constant["o"])], seq(0, 1, by = 0.1))
  # q6 = quantile(ns_dist2[ns_dist1 > as.numeric(constant["p"])], seq(0, 1, by = 0.1))
  # q7 = quantile(ns_dist3[ns_dist1 > as.numeric(constant["q"])], seq(0, 1, by = 0.1))
  
  data_info = data.table(
    data_name = c("NCI60", "CCLE2019", "CCL180"),
    FactorA_dist = c("q1", "q2", "q3"),
    FactorC_dist = c("q5", "q6", "q7"),
    FactorA_threshold = c(as.numeric(constant["a"]), as.numeric(constant["b"]), as.numeric(constant["c"])),
    FactorC_threshold = c(as.numeric(constant["o"]), as.numeric(constant["p"]), as.numeric(constant["q"]))
  )
  
  # Factor A #
  Score = sapply(1:nrow(data), function(i){
    # score prediction for 3 matrices #
    Factor_A = vector("integer", nrow(data_info))
    for (j in 1:nrow(data_info)){
      n = as.numeric(data[i, get(data_info[j, data_name])])
      if (is.na(n) | n < data_info[j, FactorA_threshold]){
        Factor_A[j] = 0
      } else if (n > data_info[j, FactorA_threshold]){
        Factor_A[j] = findInterval(n, get(data_info[j, FactorA_dist]))*3
      }
    }
    Factor_A = as.numeric(sum(Factor_A))
    return(as.numeric(Factor_A))
  })
  
  if (score.only){
    return(Score)
  }else{
    score_table = data |>
      mutate(
        Score = Score
      )
    return(score_table)
  }
}
ConfidenceScore_AB = function(data, ki = etki, score.only = FALSE){
  ## loading final constant ##
  constant = c(
    a = 0.00, b = 0.00, c = 1.0, 
    d = 0.16, e = 1.00,
    o = 0.10, p = 0.10, q = 0.20, m = 1.00
  )
  
  ## quantile ##
  ns_dist1 = na.omit(ki$NCI60) # NCI60 #
  ns_dist2 = na.omit(ki$CCLE2019) # CCLE2019 #
  ns_dist3 = na.omit(ki$CCL180) # CCL180 #
  ns_dist4 = -log10(na.omit(ki$CRISPRDEP)) # CRISPR dependency
  q1 = quantile(ns_dist1[ns_dist1 > as.numeric(constant["a"])], seq(0, 1, by = 0.1))
  q2 = quantile(ns_dist2[ns_dist2 > as.numeric(constant["b"])], seq(0, 1, by = 0.1))
  q3 = quantile(ns_dist3[ns_dist3 > as.numeric(constant["c"])], seq(0, 1, by = 0.1))
  q4 = quantile(ns_dist4[ns_dist4 > -log10(constant["d"])], seq(0, 1, by = 0.1))
  # q5 = quantile(ns_dist1[ns_dist1 > as.numeric(constant["o"])], seq(0, 1, by = 0.1))
  # q6 = quantile(ns_dist2[ns_dist1 > as.numeric(constant["p"])], seq(0, 1, by = 0.1))
  # q7 = quantile(ns_dist3[ns_dist1 > as.numeric(constant["q"])], seq(0, 1, by = 0.1))
  
  data_info = data.table(
    data_name = c("NCI60", "CCLE2019", "CCL180"),
    FactorA_dist = c("q1", "q2", "q3"),
    FactorC_dist = c("q5", "q6", "q7"),
    FactorA_threshold = c(as.numeric(constant["a"]), as.numeric(constant["b"]), as.numeric(constant["c"])),
    FactorC_threshold = c(as.numeric(constant["o"]), as.numeric(constant["p"]), as.numeric(constant["q"]))
  )
  
  Score = sapply(1:nrow(data), function(i){
    # Factor A #
    Factor_A = vector("integer", nrow(data_info))
    for (j in 1:nrow(data_info)){
      n = as.numeric(data[i, get(data_info[j, data_name])])
      if (is.na(n) | n < data_info[j, FactorA_threshold]){
        Factor_A[j] = 0
      } else if (n > data_info[j, FactorA_threshold]){
        Factor_A[j] = findInterval(n, get(data_info[j, FactorA_dist]))*3
      }
    }
    Factor_A = as.numeric(sum(Factor_A))
    
    # Factor B #
    if (is.na(data[i, CRISPRDEP]) | data[i, -log10(as.numeric(CRISPRDEP))] < -log10(constant["d"])){
      Factor_B = 0
    } else if (data[i, -log10(as.numeric(CRISPRDEP))] > -log10(constant["d"])){
      Factor_B = findInterval(data[i, -log10(as.numeric(CRISPRDEP))], q4)*as.numeric(constant["e"])
    }
    
    finalscore = as.numeric(Factor_A) + as.numeric(Factor_B)
  })
  
  if (score.only){
    return(Score)
  }else{
    score_table = data |>
      mutate(
        Score = Score
      )
    return(score_table)
  }
}
ConfidenceScore_ABC = function(data, ki = etki, score.only = FALSE){
  ## constant ##
  ## loading final constant ##
  constant = c(
    a = 0.00, b = 0.00, c = 1.0, 
    d = 0.16, e = 1.00,
    o = 0.10, p = 0.10, q = 0.20, m = 1.00
  )
  
  ## quantile ##
  ns_dist1 = na.omit(ki$NCI60) # NCI60 #
  ns_dist2 = na.omit(ki$CCLE2019) # CCLE2019 #
  ns_dist3 = na.omit(ki$CCL180) # CCL180 #
  ns_dist4 = -log10(na.omit(ki$CRISPRDEP)) # CRISPR dependency
  q1 = quantile(ns_dist1[ns_dist1 > as.numeric(constant["a"])], seq(0, 1, by = 0.1))
  q2 = quantile(ns_dist2[ns_dist2 > as.numeric(constant["b"])], seq(0, 1, by = 0.1))
  q3 = quantile(ns_dist3[ns_dist3 > as.numeric(constant["c"])], seq(0, 1, by = 0.1))
  q4 = quantile(ns_dist4[ns_dist4 > -log10(constant["d"])], seq(0, 1, by = 0.1))
  q5 = quantile(ns_dist1[ns_dist1 > as.numeric(constant["o"])], seq(0, 1, by = 0.1))
  q6 = quantile(ns_dist2[ns_dist1 > as.numeric(constant["p"])], seq(0, 1, by = 0.1))
  q7 = quantile(ns_dist3[ns_dist1 > as.numeric(constant["q"])], seq(0, 1, by = 0.1))
  
  data_info = data.table(
    data_name = c("NCI60", "CCLE2019", "CCL180"),
    FactorA_dist = c("q1", "q2", "q3"),
    FactorC_dist = c("q5", "q6", "q7"),
    FactorA_threshold = c(as.numeric(constant["a"]), as.numeric(constant["b"]), as.numeric(constant["c"])),
    FactorC_threshold = c(as.numeric(constant["o"]), as.numeric(constant["p"]), as.numeric(constant["q"]))
  )
  
  cl = makeCluster(8)
  registerDoParallel(cl)
  
  Score = foreach(i = 1:nrow(data), .combine = c, .export = unique(c(ls(globalenv()), ls(environment())))) %dopar% {
    library(data.table)
    library(tidyverse)
    
    # Factor A #
    Factor_A = vector("integer", nrow(data_info))
    for (j in 1:nrow(data_info)){
      n = as.numeric(data[i, get(data_info[j, data_name])])
      if (is.na(n) | n < data_info[j, FactorA_threshold]){
        Factor_A[j] = 0
      } else if (n > data_info[j, FactorA_threshold]){
        Factor_A[j] = findInterval(n, get(data_info[j, FactorA_dist]))*3
      }
    }
    Factor_A = as.numeric(sum(Factor_A))
    
    # Factor B #
    if (is.na(data[i, CRISPRDEP]) | data[i, -log10(as.numeric(CRISPRDEP))] < -log10(constant["d"])){
      Factor_B = 0
    } else if (data[i, -log10(as.numeric(CRISPRDEP))] > -log10(constant["d"])){
      Factor_B = findInterval(data[i, -log10(as.numeric(CRISPRDEP))], q4)*as.numeric(constant["e"])
    }
    
    # Factor C #
    if (is.na(data$Formula[i])){
      temp = getAM(data$KEGG_ID[i], is.kegg = TRUE)
    }else{
      temp = getAM(data$Formula[i])
    }
    if (check_NA(temp)){
      Factor_C = 0
    }else if (!1 %in% names(temp) & !2 %in% names(temp)){
      Factor_C = 0
    }else{
      proximal = translate(temp[which(names(temp) %in% c(1,2))])
      if (check_NA(proximal)){
        Factor_C = 0
      }else{
        proximal_ns = getNR(data[i, SYMBOL], proximal)
        if (nrow(proximal_ns) != 0) {
          Factor_C = 0
          for (l in 1:nrow(data_info)){
            N = unique(proximal_ns[, get(data_info[l, data_name])])
            for (n in N){
              if (is.na(n) | n < data_info[l, FactorC_threshold]){
                Factor_C = Factor_C + 0
              }else if (n > data_info[l, FactorC_threshold] & n < data_info[l, FactorA_threshold]){
                Factor_C = Factor_C + 0.5*constant["m"]
              }else if (n > data_info[l, FactorA_threshold]){
                Factor_C = Factor_C + findInterval(n, get(data_info[l, FactorC_dist]))*constant["m"]
              }
            }
            print(Factor_C)
          }
        }else{
          Factor_C = 0
        }
      }
    }
    
    # Final score #
    finalscore = Factor_A + Factor_B + Factor_C
    return(finalscore)
  }
  
  stopCluster(cl)
  
  if (score.only){
    return(as.numeric(Score))
  }else{
    score_table = data |>
      mutate(
        Score = as.numeric(Score)
      )
    return(score_table)
  }
}

## 8.2 Compare ConfidenceScore across incorporation ####
# Known data #
ConfA = ConfidenceScore_A(etki_Tr)
ConfAB = ConfidenceScore_AB(etki_Tr)
ConfABC = ConfidenceScore_ABC(etki_Tr)

# Simulated data #
simulated_ConfA = lapply(fake_etki_Tr, ConfidenceScore_A)
simulated_ConfAB = lapply(fake_etki_Tr, ConfidenceScore_AB)
simulated_ConfABC = lapply(fake_etki_Tr, ConfidenceScore_ABC)

# --------------------------- #
#       Figure 1D data        #
# --------------------------- #
analyseDist = function(data.simu, data.true){
  # Bootstrap
  boot = sapply(1:100, function(s){
    resample = sample(data.true[, as.numeric(Score)], size = nrow(data.true), replace = TRUE)
    return(mean(resample))
  })
  
  Dist = data.table(
    type = rep(c("simulated", "known"), each = 100),
    Score = c(data.simu[, mean(as.numeric(Score)), by = Iteration][, V1], boot)
  )
  
  rain_height = .1
  colorset = c("known" = "#990000", "simulated" = "grey")
  g = ggplot(Dist, aes(x = "", y = Score, fill = type)) +
    # mountain
    introdataviz::geom_flat_violin(trim=FALSE, alpha = 0.4, position = position_nudge(x = -0.6)) +
    # rain
    geom_point(aes(colour = type), size = 1, alpha = 0.4, show.legend = FALSE, position = position_jitter(width = rain_height, height = 0)) +
    geom_boxplot(width = rain_height, alpha = 0.4, show.legend = FALSE, outlier.shape = NA) +
    scale_y_continuous(name = "Average Confidence Score",
                       breaks = seq(6, 20, 2), 
                       limits = c(6, 20)) +
    # Manual color settings
    scale_fill_manual(values = colorset) + # For violin plots
    scale_color_manual(values = colorset) + # For points
    coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
  
  p = t.test(
    x = Dist[type == "simulated", Score],
    y = Dist[type == "known", Score],
    alternative = "less"
  )$p.value
  
  return(list(dist.data = Dist, dist.plot = g, dist.pvalue = p))
}
simulated_ConfA_4dist = lapply(1:100, function(f){
  data = simulated_ConfA[[f]] |>
    mutate(
      Iteration = f,
      .before = everything()
    ) |>
    mutate(
      Score = as.numeric(Score)
    )
  return(data)
}) |> rbindlist()

conf.dist = analyseDist(data.simu = simulated_ConfA_4dist, data.true = ConfA)
conf.dist$dist.plot
conf.dist$dist.pvalue
# Fig 1D (conf score).csv
conf.dist$dist.data

# --------------------------- #
#       Figure 1E data        #
# --------------------------- #

Fr = lapply(0:100, function(th){
  true_fr = nrow(ConfA[as.numeric(Score) > th])/nrow(ConfA)
  fake_fr = sapply(1:100, function(f){
    this_iteration = simulated_ConfA[[f]]
    return(nrow(this_iteration[as.numeric(Score) > th])/nrow(this_iteration))
  })
  message(th)
  return(
    data.table(
      Th = th,
      Class = c("Known", rep("Simulated random", 100)),
      Fraction = c(true_fr, fake_fr)
    )
  )
}) |> rbindlist()

wilcox.test(Fr[Class == "Simulated random", Fraction], Fr[Class == "Known", Fraction], alternative = "less")

ggplot() +
  geom_violin(data = Fr[Class == "Simulated random"], 
              mapping = aes(x = Th, y = Fraction, group = Th),
              color = "grey") +
  geom_line(data = Fr[Class == "Known"],
            mapping = aes(x = Th, y = Fraction),
            color = "#990000",
            linetype="dashed") + 
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks=seq(0,80,20)) +
  scale_y_continuous(limits = c(0, 0.65), breaks = seq(0,0.6,0.1))

# --------------------------- #
#     Figure 2C & 3C data     #
# --------------------------- #
ConfA_dcs = ConfA[, mean(as.numeric(Score))] - unlist(lapply(simulated_ConfA, function(x){mean(as.numeric(x[,Score]))}))
ConfAB_dcs = ConfAB[, mean(as.numeric(Score))] - unlist(lapply(simulated_ConfAB, function(x){mean(as.numeric(x[,Score]))}))
ConfABC_dcs = ConfABC[, mean(as.numeric(Score))] - unlist(lapply(simulated_ConfABC, function(x){mean(as.numeric(x[,Score]))}))

dCS = data.table(dCS = c(ConfA_dcs, ConfAB_dcs, ConfABC_dcs),
                 Factor = rep(c("A", "AB", "ABC"), each = 100))

shapiro.test(dCS[Factor == "A", dCS])
shapiro.test(dCS[Factor == "AB", dCS])
shapiro.test(dCS[Factor == "ABC", dCS])

# Figure 2C
ggplot(dCS[Factor %in% c("A", "AB")], aes(x = Factor, y = dCS, fill = Factor)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#3D5C6F", "#F9AE78")) +
  expand_limits(y = 6) +
  scale_y_continuous(breaks = seq(4,6,1)) +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

t.test(dCS[Factor == "A", dCS], dCS[Factor == "AB", dCS], alternative = "l")$p.value

# Figure 3C
mean((ConfABC_dcs - ConfAB_dcs)/ConfAB_dcs)
ggplot(dCS, aes(x = Factor, y = dCS, fill = Factor)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#3D5C6F", "#F9AE78","#E47159")) +
  expand_limits(y = 14) +
  scale_y_continuous(breaks = seq(4,14,2)) +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

t.test(dCS[Factor == "AB", dCS], dCS[Factor == "ABC", dCS], alternative = "l")$p.value

# --------------------------- #
#      Figure 3D data         #
# --------------------------- #

Fr = lapply(0:190, function(Th){
  FrA = nrow(ConfA[as.numeric(Score) > Th])/nrow(ConfA) - sapply(simulated_ConfA, function(f){nrow(f[as.numeric(Score) > Th])/nrow(f)})
  FrAB = nrow(ConfAB[as.numeric(Score) > Th])/nrow(ConfAB) - sapply(simulated_ConfAB, function(f){nrow(f[as.numeric(Score) > Th])/nrow(f)})
  FrABC = nrow(ConfABC[as.numeric(Score) > Th])/nrow(ConfABC) - sapply(simulated_ConfABC, function(f){nrow(f[as.numeric(Score) > Th])/nrow(f)})
  message(Th)
  return(
    data.table(
      Th = Th,
      Class = rep(c("A", "AB", "ABC"), each = 100),
      Fr = c(FrA, FrAB, FrABC)
    )
  )
}) |> rbindlist()

ggplot(Fr, aes(x = Th, y = Fr, group = interaction(Th, Class), fill = Class, color = Class)) +
  geom_violin(alpha = 0.2) +
  scale_color_manual(values=c("#3D5C6F", "#F9AE78","#E47159")) +
  scale_fill_manual(values=c("#3D5C6F", "#F9AE78","#E47159")) +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  # scale_x_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  scale_y_continuous(limits = c(0, 0.2), breaks = seq(0, 0.2, 0.05))

wilcox.test(Fr[Class == "AB", Fr], Fr[Class == "ABC", Fr], alternative = "l")$p.value

# --------------------------- #
#     Figure S3A data         #
# --------------------------- #

roc = function(data.simu, data.true){
  boundary = c(
    true.upper = data.true[, max(as.numeric(Score))],
    fake.upper = data.simu[, max(as.numeric(Score))]
  )
  boundary.upper = ceiling(max(boundary))
  boundary.lower = 0
  threshold_ticking = seq(boundary.upper, boundary.lower, -0.1)
  ROCcurve = lapply(threshold_ticking, function(t){
    true_pos  = nrow(data.true[Score >= t])/nrow(data.true)
    false_pos = sapply(1:100, function(i){
      this_iteration = data.simu[Iteration == i]
      return(nrow(this_iteration[Score >= t])/nrow(this_iteration))
    })
    progress = round(which(threshold_ticking == t)/length(threshold_ticking), digits = 3) * 100
    if (progress %% 10 == 0) {
      message("Generating ROC ", progress, " %")
    }
    return(
      data.table(
        threshold = t,
        true.positive.rate = true_pos,
        false.positive.rate = mean(false_pos)
      )
    )
  }) |> rbindlist()
  # Generates plot
  ROCcurve.plot = ggplot(ROCcurve) +
    geom_point(aes(x = false.positive.rate, y = true.positive.rate)) +
    geom_abline(col = "grey", linetype = "dashed") +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  # Calculate AUC with the Trapezoidal rule
  sensitivity = ROCcurve$true.positive.rate
  specificity = ROCcurve$false.positive.rate
  height = (sensitivity[-1]+sensitivity[-length(sensitivity)])/2
  width = diff(specificity)
  AUC = sum(height*width)
  
  return(
    list(
      ROC.data = ROCcurve,
      ROC.plot = ROCcurve.plot,
      AUC = AUC
    )
  )
}
simulated_ConfA.roc = lapply(1:100, function(s){
  data = simulated_ConfA[[s]] |>
    mutate(
      Iteration = s,
      .before = everything()
    ) |>
    mutate(
      Score = as.numeric(Score)
    )
  return(data)
}) |> rbindlist()
ConfA.roc = ConfA |>
  mutate(
    Score = as.numeric(Score)
  )
# Fig S3A (left).csv
simulated_ConfA.roc.object = roc(data.simu = simulated_ConfA.roc, data.true = ConfA.roc)
simulated_ConfA.roc.object$ROC.plot
simulated_ConfA.roc.object$ROC.data
simulated_ConfA.roc.object$AUC

simulated_ConfAB.roc = lapply(1:100, function(s){
  data = simulated_ConfAB[[s]] |>
    mutate(
      Iteration = s,
      .before = everything()
    )  |>
    mutate(
      Score = as.numeric(Score)
    )
  return(data)
}) |> rbindlist()
ConfAB.roc = ConfAB |>
  mutate(
    Score = as.numeric(Score)
  )
# Fig S3A (middle).csv
simulated_ConfAB.roc.object = roc(data.simu = simulated_ConfAB.roc, data.true = ConfAB.roc)
simulated_ConfAB.roc.object$ROC.plot
simulated_ConfAB.roc.object$ROC.data
simulated_ConfAB.roc.object$AUC

simulated_ConfABC.roc = lapply(1:100, function(s){
  data = simulated_ConfABC[[s]] |>
    mutate(
      Iteration = s,
      .before = everything()
    )  |>
    mutate(
      Score = as.numeric(Score)
    )
  return(data)
}) |> rbindlist()
ConfABC.roc = ConfABC |>
  mutate(
    Score = as.numeric(Score)
  )
# Fig S3A (right).csv
simulated_ConfABC.roc.object = roc(data.simu = simulated_ConfABC.roc, data.true = ConfABC.roc)
simulated_ConfABC.roc.object$ROC.plot
simulated_ConfABC.roc.object$ROC.data
simulated_ConfABC.roc.object$AUC

# --------------------------- #
#     Figure S3B data         #
# --------------------------- #

prc = function(data.simu, data.true){
  boundary = c(
    true.upper = data.true[, max(as.numeric(Score))],
    fake.upper = data.simu[, max(as.numeric(Score))]
  )
  boundary.upper = ceiling(max(boundary))
  boundary.lower = 0
  threshold_ticking = seq(boundary.lower, boundary.upper, 0.1)
  PRCcurve = lapply(threshold_ticking, function(t){
    true_pos  = nrow(data.true[Score >= t])
    false_pos = sapply(1:100, function(i){
      this_iteration = data.simu[Iteration == i]
      return(nrow(this_iteration[Score >= t]))
    })
    false_neg = nrow(data.true[Score < t])
    if (true_pos == 0 & mean(false_pos) == 0){
      precision = 1
    }else{
      precision = true_pos/(true_pos + mean(false_pos))
    }
    recall = true_pos/(true_pos + false_neg)
    progress = round(which(threshold_ticking == t)/length(threshold_ticking), digits = 3) * 100
    if (progress %% 10 == 0) {
      message("Generating PRC ", progress, " %")
    }
    return(
      data.table(
        threshold = t,
        true.positive = true_pos,
        false.positive = mean(false_pos),
        false.negative = false_neg,
        Precision = precision,
        Recall = recall
      )
    )
  }) |> rbindlist()
  # Generates plot
  PRCcurve.plot = ggplot(PRCcurve) +
    geom_point(aes(x = Recall, y = Precision)) +
    geom_line(aes(x = Recall, y = Precision)) +
    geom_hline(col = "grey", linetype = "dashed", yintercept = 0.5) +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  # Calculate AUC with the Trapezoidal rule
  sensitivity = PRCcurve$Precision
  specificity = PRCcurve$Recall
  height = (sensitivity[-1]+sensitivity[-length(sensitivity)])/2
  width = -diff(specificity)
  AUC = sum(height*width)
  
  return(
    list(
      PRC.data = PRCcurve,
      PRC.plot = PRCcurve.plot,
      AUC = AUC
    )
  )
}

simulated_ConfA.prc = lapply(1:100, function(s){
  data = simulated_ConfA[[s]] |>
    mutate(
      Iteration = s,
      .before = everything()
    ) |>
    mutate(
      Score = as.numeric(Score)
    )
  return(data)
}) |> rbindlist()
ConfA.prc = ConfA |>
  mutate(
    Score = as.numeric(Score)
  )
# Fig S3B (left).csv
simulated_ConfA.prc.object = prc(data.simu = simulated_ConfA.prc, data.true = ConfA.prc)
simulated_ConfA.prc.object$PRC.plot
simulated_ConfA.prc.object$PRC.data
simulated_ConfA.prc.object$AUC

simulated_ConfAB.prc = lapply(1:100, function(s){
  data = simulated_ConfAB[[s]] |>
    mutate(
      Iteration = s,
      .before = everything()
    )  |>
    mutate(
      Score = as.numeric(Score)
    )
  return(data)
}) |> rbindlist()
ConfAB.prc = ConfAB |>
  mutate(
    Score = as.numeric(Score)
  )
# Fig S3B (middle).csv
simulated_ConfAB.prc.object = prc(data.simu = simulated_ConfAB.prc, data.true = ConfAB.prc)
simulated_ConfAB.prc.object$PRC.plot
simulated_ConfAB.prc.object$PRC.data
simulated_ConfAB.prc.object$AUC

simulated_ConfABC.prc = lapply(1:100, function(s){
  data = simulated_ConfABC[[s]] |>
    mutate(
      Iteration = s,
      .before = everything()
    )  |>
    mutate(
      Score = as.numeric(Score)
    )
  return(data)
}) |> rbindlist()
ConfABC.prc = ConfABC |>
  mutate(
    Score = as.numeric(Score)
  )
# Fig S3B (right).csv
simulated_ConfABC.prc.object = prc(data.simu = simulated_ConfABC.prc, data.true = ConfABC.prc)
simulated_ConfABC.prc.object$PRC.plot
simulated_ConfABC.prc.object$PRC.data
simulated_ConfABC.prc.object$AUC

# 9 Partition, training and testing ####
dir.create("Results/AutoOpt")
## 9.1 Partition and training ####
n_Opt = 10
for (Opt_time in 1:n_Opt){
  ### §1 Known interaction partition ####
  # 70% for training
  training_set = sample(1:nrow(etki), 0.7*nrow(etki))
  etki_Tr = etki[training_set]
  # 30% for testing
  entire_set = 1:nrow(etki)
  testing_set = entire_set[!entire_set %in% training_set]
  etki_Te = etki[testing_set]
  
  # Simulation generation
  fake_etki_Tr = lapply(1:100, function(sample_time){
    etki_gene_pool = unique(etki_Tr$SYMBOL)
    etki_metabolite_pool = unique(etki_Tr$Formula)
    # initialise fake etki_Tr data table #
    fake_etki = setNames(data.table(matrix(ncol = ncol(etki_Tr), nrow = nrow(etki_Tr))), colnames(etki_Tr))
    fake_etki[, Substrate := NULL]
    char.col = colnames(fake_etki)[1:6]
    num.col = colnames(fake_etki)[7:10]
    fake_etki[, (char.col) := lapply(.SD, as.character), .SDcols = char.col]
    fake_etki[, (num.col) := lapply(.SD, as.character), .SDcols = num.col]
    fake_etki[, c("SYMBOL", "Gene_type","Localisation") := as.list(etki_Tr[, c("SYMBOL", "Gene_type","Localisation")])]
    fake_etki[, Status := rep("Fake", nrow(fake_etki))]
    # for each SLC, shuffle fake substrate formula #
    for (j in 1:nrow(fake_etki)){
      true_substrate = etki_Tr[SYMBOL == fake_etki$SYMBOL[j]][, Formula]
      fake_substrate_pool = etki_metabolite_pool[!etki_metabolite_pool %in% true_substrate]
      fake_substrate = sample(fake_substrate_pool, 1)
      while (nrow(et[SYMBOL == fake_etki$SYMBOL[j] & Formula == fake_substrate][,NCI60:CRISPRDEP]) == 0){
        fake_substrate = sample(fake_substrate_pool, 1)
      }
      fake_etki[j, Formula := fake_substrate]
      fake_etki[j, c("NCI60", "CCLE2019", "CCL180", "CRISPRDEP") := 
                  as.list(et[SYMBOL == fake_etki$SYMBOL[j] & Formula == fake_etki$Formula[j]][,NCI60:CRISPRDEP])]
    }
    print(paste(sample_time, "done"))
    fake_etki = fake_etki |>
      mutate(
        Iteration = sample_time,
        .before = everything()
      )
    return(fake_etki)
  }) |> rbindlist()
  
  ## §2 Confidence Score A Optimisation ####
  ConfidenceScore_A_by_set = function(data, ki, a, b, c, score.only = FALSE){
    # Factor A: Predictive Scores
    ## quantile ##
    ns_dist1 = na.omit(ki$NCI60) # NCI60 #
    ns_dist2 = na.omit(ki$CCLE2019) # CCLE2019 #
    ns_dist3 = na.omit(ki$CCL180) # CCL180 #
    q1 = quantile(ns_dist1[ns_dist1 > as.numeric(a)], seq(0, 1, by = 0.1))
    q2 = quantile(ns_dist2[ns_dist2 > as.numeric(b)], seq(0, 1, by = 0.1))
    q3 = quantile(ns_dist3[ns_dist3 > as.numeric(c)], seq(0, 1, by = 0.1))
    
    Factor_A_info = data.table(Factor_A_name = c("NCI60", "CCLE2019", "CCL180"),
                               Factor_A_dist = c("q1", "q2", "q3"),
                               Factor_A_threshold = c(as.numeric(a), as.numeric(b), as.numeric(c)))
    Score = sapply(1:nrow(data), function(i){
      # score prediction for 3 matrices #
      Factor_A = vector("integer", nrow(Factor_A_info))
      for (j in 1:nrow(Factor_A_info)){
        n = as.numeric(data[i, get(Factor_A_info[j, Factor_A_name])])
        if (is.na(n) | n < Factor_A_info[j, Factor_A_threshold]){
          Factor_A[j] = 0
        } else if (n > Factor_A_info[j, Factor_A_threshold]){
          Factor_A[j] = findInterval(n, get(Factor_A_info[j, Factor_A_dist]))*3
        }
      }
      Factor_A = sum(Factor_A)
      return(as.numeric(Factor_A))
    })
    
    if (score.only){
      return(Score)
    }else{
      score_table = data |>
        mutate(
          Score = Score
        )
      return(score_table)
    }
  }
  Best_threshold_detect_A_by_set = function(ki.data, si.data, from.range){
    # Setting up parallel
    numCores = 8  # Use all available cores except one for other processes
    cl = makeCluster(numCores)
    
    # a: range for NCI60
    registerDoParallel(cl)
    range_a = foreach(a = from.range, .combine = rbind, .export = ls(globalenv())) %dopar% {
      library(data.table)
      library(tidyverse)
      # calculate scores for etki and fake etki under threshold a
      # i.e. if we only consider NCI60 transformed-rho > a, and assign confidence score based on a
      #      what would be the confidence score for each SLC-metabolite pair?
      ecs = ConfidenceScore_A_by_set(data = etki_Tr, ki = ki.data, a = a, b = Inf, c = Inf, score.only = TRUE)
      fecs = lapply(1:100, function(f){
        simuScore = ConfidenceScore_A_by_set(data = si.data[Iteration == f], ki = ki.data, a = a, b = Inf, c = Inf, score.only = TRUE)
        return(simuScore)
      })
      
      # calculate fraction recovery (TP - FP) for each given score cutoff
      # i.e. how many SLC-metabolite pair can be discovered, if we set score cutoff from 0 (all) - 34 (none)
      fr = lapply(0:34, function(cutoff){
        TP = length(ecs[ecs >= cutoff])/length(ecs)
        FP = sapply(fecs, function(f){length(f[f >= cutoff])}/length(f))
        
        if (TP == 0 & mean(FP) == 0){
          pvalue = 1
        }else if (TP == 1 & mean(FP) == 1) {
          pvalue = 1
        }else{
          pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
        }
        
        return(
          data.table(
            Cutoff = cutoff,
            TP = TP,
            FP.mean = mean(FP),
            FP.sd = sd(FP),
            TTP.literal = TP - mean(FP),
            TTP.in.sd = (TP - mean(FP))/sd(FP),
            pvalue = pvalue
          )
        )
      }) |> rbindlist() |>
        mutate(
          Threshold = a,
          .before = everything()
        )
      
      return(fr)
      
    }
    
    message("range detection a finished at ", Sys.time())
    
    # b: range for CCLE2019
    range_b = foreach(b = from.range, .combine = rbind, .export = ls(globalenv())) %dopar% {
      library(data.table)
      library(tidyverse)
      # calculate scores for etki and fake etki under threshold b
      # i.e. if we only consider CCLE2019 transformed-rho > b, and assign confidence score based on b
      #      what would be the confidence score for each SLC-metabolite pair?
      ecs = ConfidenceScore_A_by_set(data = etki_Tr, ki = ki.data, a = Inf, b = b, c = Inf, score.only = TRUE)
      fecs = lapply(1:100, function(f){
        simuScore = ConfidenceScore_A_by_set(data = si.data[Iteration == f], ki = ki.data, a = Inf, b = b, c = Inf, score.only = TRUE)
        return(simuScore)
      })
      
      # calculate fraction recovery (TP - FP) for each given score cutoff
      # i.e. how many SLC-metabolite pair can be discovered, if we set score cutoff from 0 (all) - 34 (none)
      fr = lapply(0:34, function(cutoff){
        TP = length(ecs[ecs >= cutoff])/length(ecs)
        FP = sapply(fecs, function(f){length(f[f >= cutoff])}/length(f))
        
        if (TP == 0 & mean(FP) == 0){
          pvalue = 1
        }else if (TP == 1 & mean(FP) == 1) {
          pvalue = 1
        }else{
          pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
        }
        
        return(
          data.table(
            Cutoff = cutoff,
            TP = TP,
            FP.mean = mean(FP),
            FP.sd = sd(FP),
            TTP.literal = TP - mean(FP),
            TTP.in.sd = (TP - mean(FP))/sd(FP),
            pvalue = pvalue
          )
        )
      }) |> rbindlist() |>
        mutate(
          Threshold = b,
          .before = everything()
        )
      
      return(fr)
      
    }
    
    message("range detection b finished at ", Sys.time())
    
    # c: range for CCL180
    
    range_c = foreach(c = from.range, .combine = rbind, .export = ls(globalenv())) %dopar% {
      library(data.table)
      library(tidyverse)
      # calculate scores for etki and fake etki under threshold c
      # i.e. if we only consider CCLE2019 transformed-rho > c, and assign confidence score based on c
      #      what would be the confidence score for each SLC-metabolite pair?
      ecs = ConfidenceScore_A_by_set(data = etki_Tr, ki = ki.data, a = Inf, b = Inf, c = c, score.only = TRUE)
      fecs = lapply(1:100, function(f){
        simuScore = ConfidenceScore_A_by_set(data = si.data[Iteration == f], ki = ki.data, a = Inf, b = Inf, c = c, score.only = TRUE)
        return(simuScore)
      })
      
      # calculate fraction recovery (TP - FP) for each given score cutoff
      # i.e. how many SLC-metabolite pair can be discovered, if we set score cutoff from 0 (all) - 34 (none)
      fr = lapply(0:34, function(cutoff){
        TP = length(ecs[ecs >= cutoff])/length(ecs)
        FP = sapply(fecs, function(f){length(f[f >= cutoff])}/length(f))
        
        if (TP == 0 & mean(FP) == 0){
          pvalue = 1
        }else if (TP == 1 & mean(FP) == 1) {
          pvalue = 1
        }else{
          pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
        }
        
        return(
          data.table(
            Cutoff = cutoff,
            TP = TP,
            FP.mean = mean(FP),
            FP.sd = sd(FP),
            TTP.literal = TP - mean(FP),
            TTP.in.sd = (TP - mean(FP))/sd(FP),
            pvalue = pvalue
          )
        )
      }) |> rbindlist() |>
        mutate(
          Threshold = c,
          .before = everything()
        )
      
      return(fr)
      
    }
    
    message("range detection c finished at ", Sys.time())
    
    stopCluster(cl)
    
    # Make decisions on a, b, and c based on range report
    a.report = range_a[order(TTP.literal, decreasing = TRUE)][1, Threshold]
    b.report = range_b[order(TTP.literal, decreasing = TRUE)][1, Threshold]
    c.report = range_c[order(TTP.literal, decreasing = TRUE)][1, Threshold]
    
    return(list(a = a.report, b = b.report, c = c.report, a.record = range_a, b.record = range_b, c.record = range_c))
  }
  
  ConfA.threshold.detection = Best_threshold_detect_A_by_set(ki.data = etki_Tr, si.data = fake_etki_Tr, from.range = seq(0,1,0.1))
  
  ## §3 Confidence Score B Optimisation ####
  #### 3.1 Determine a p-value threshold for CRISPR dependency screen ####
  # i.e. Only p-value that's below the threshold will be considered for score reward
  Best_threshold_detect_B_by_set = function(ki.data, si.data, from.range){
    fr = lapply(from.range, function(cutoff){
      TP = nrow(ki.data[CRISPRDEP <= cutoff])/nrow(ki.data)
      FP = sapply(1:100, function(f){nrow(si.data[Iteration == f & CRISPRDEP <= cutoff])/nrow(si.data[Iteration == f])})
      
      if (TP == 0 & mean(FP) == 0){
        pvalue = 1
      }else{
        pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
      }
      
      return(
        data.table(
          Cutoff = cutoff,
          TP = TP,
          FP.mean = mean(FP),
          FP.sd = sd(FP),
          TTP.literal = TP - mean(FP),
          TTP.in.sd = (TP - mean(FP))/sd(FP),
          pvalue = pvalue
        )
      )
    }) |> rbindlist()
    
    message("range detection d finished at ", Sys.time())
    
    # Make decisions based on range report
    d.report = fr[order(TTP.literal, decreasing = TRUE)][1, Cutoff]
    return(list(d = d.report, d.record = fr))
  }
  
  ConfB.threshold.detection = Best_threshold_detect_B_by_set(ki.data = etki_Tr, si.data = fake_etki_Tr, from.range = seq(0,0.2,0.01))
  
  #### 3.2 Determine a score reward for CRISPR dependency screen ####
  ConfidenceScore_B_by_set = function(data, ki, d, e, score.only = FALSE){
    # Factor B quantile
    ns_dist4 = -log10(na.omit(ki$CRISPRDEP))
    q4 = quantile(ns_dist4[ns_dist4 > -log10(d)], seq(0, 1, by = 0.1))
    
    # Calculate scores
    Score = sapply(1:nrow(data), function(i){
      # Factor B: Dependency Scores #
      if (is.na(data[i, CRISPRDEP]) | data[i, -log10(as.numeric(CRISPRDEP))] < -log10(d)){
        Factor_B = 0
      } else if (data[i, -log10(as.numeric(CRISPRDEP))] > -log10(d)){
        Factor_B = findInterval(data[i, -log10(as.numeric(CRISPRDEP))], q4)*as.numeric(e)
      }
      return(as.numeric(Factor_B))
    })
    
    if (score.only){
      return(Score)
    }else{
      score_table = data |>
        mutate(
          Score = Score
        )
      return(score_table)
    }
  }
  
  Best_scalar_detect_B_by_set = function(ki.data, si.data, d, from.range){
    # e: range of scalar reward for CRISPR dependency
    range_e = lapply(from.range, function(e){
      # calculate scores for etki and fake etki given scalar e
      # i.e. if we only consider CRISPR dependency pvalue < d, and assign confidence score scaled by e
      #      what would be the confidence score for each SLC-metabolite pair?
      ecs = ConfidenceScore_B_by_set(data = etki_Tr, ki = ki.data, d = d, e = e, score.only = TRUE)
      fecs = lapply(1:100, function(f){
        simuScore = ConfidenceScore_B_by_set(data = si.data[Iteration == f], ki = ki.data, d = d, e = e, score.only = TRUE)
        return(simuScore)
      })
      
      # calculate fraction recovery (TP - FP) for each given score cutoff
      # i.e. how many SLC-metabolite pair can be discovered, if we set score cutoff from 0 (all) - 12 (none)
      fr = lapply(0:12, function(cutoff){
        TP = length(ecs[ecs >= cutoff])/length(ecs)
        FP = sapply(fecs, function(f){length(f[f >= cutoff])}/length(f))
        
        if (TP == 0 & mean(FP) == 0){
          pvalue = 1
        }else if (TP == 1 & mean(FP) == 1) {
          pvalue = 1
        }else{
          pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
        }
        
        return(
          data.table(
            Cutoff = cutoff,
            TP = TP,
            FP.mean = mean(FP),
            FP.sd = sd(FP),
            TTP.literal = TP - mean(FP),
            TTP.in.sd = (TP - mean(FP))/sd(FP),
            pvalue = pvalue
          )
        )
      }) |> rbindlist() |>
        mutate(
          Scalar = e,
          .before = everything()
        )
      
      return(fr)
      
    }) |> rbindlist()
    
    message("range detection e finished at ", Sys.time())
    
    # Make decisions based on e range report
    e.report = range_e[order(TTP.literal, decreasing = TRUE)][1, Scalar]
    
    return(list(e = e.report, e.record = range_e))
  }
  
  ConfB.scalar.detection = Best_scalar_detect_B_by_set(ki.data = etki_Tr, si.data = fake_etki_Tr, d = as.numeric(ConfB.threshold.detection$d), from.range = seq(0.1, 3, 0.1))
  
  ### §4 Confidence Score C Optimisation ####
  #### 4.1 Determine thresholds for adjacency matrix ####
  ConfidenceScore_C_by_set = function(data, ki, o, p, q, m, score.only = FALSE){
    ## quantile ##
    ns_dist1 = na.omit(ki$NCI60) # NCI60 #
    ns_dist2 = na.omit(ki$CCLE2019) # CCLE2019 #
    ns_dist3 = na.omit(ki$CCL180) # CCL180 #
    q1 = quantile(ns_dist1[ns_dist1 > as.numeric(o)], seq(0, 1, by = 0.1))
    q2 = quantile(ns_dist2[ns_dist2 > as.numeric(p)], seq(0, 1, by = 0.1))
    q3 = quantile(ns_dist3[ns_dist3 > as.numeric(q)], seq(0, 1, by = 0.1))
    
    data_info = data.table(data_name = c("NCI60", "CCLE2019", "CCL180"),
                           data_dist = c("q1", "q2", "q3"),
                           Factor_C_threshold = c(as.numeric(o), as.numeric(p), as.numeric(q)))
    Score = sapply(1:nrow(data), function(i){
      
      if (is.na(data$Formula[i])){
        temp = getAM(data$KEGG_ID[i], is.kegg = TRUE)
      }else{
        temp = getAM(data$Formula[i])
      }
      if (check_NA(temp)){
        Factor_C = 0
      }else if (!1 %in% names(temp) & !2 %in% names(temp)){
        Factor_C = 0
      }else{
        proximal = translate(temp[which(names(temp) %in% c(1,2))])
        if (check_NA(proximal)){
          Factor_C = 0
        }else{
          proximal_ns = getNR(data[i, SYMBOL], proximal)
          if (nrow(proximal_ns) != 0) {
            Factor_C = 0
            for (l in 1:nrow(data_info)){
              N = unique(proximal_ns[, get(data_info[l, data_name])])
              for (n in N){
                if (is.na(n) | n < data_info[l, Factor_C_threshold]){
                  Factor_C = Factor_C + 0
                }else if (n > data_info[l, Factor_C_threshold]){
                  Factor_C = Factor_C + findInterval(n, get(data_info[l, data_dist]))*m
                }
              }
            }
          }else{
            Factor_C = 0
          }
        }
      }
      
      return(as.numeric(Factor_C))
    })
    
    if (score.only){
      return(Score)
    }else{
      score_table = data |>
        mutate(
          Score = Score
        )
      return(score_table)
    }
  }
  
  Best_threshold_detect_C_by_set = function(ki.data, si.data, from.range){
    # Setting up parallel
    numCores = 8  # Use all available cores except one for other processes
    cl = makeCluster(numCores)
    
    # o: derivative range for NCI60
    registerDoParallel(cl)
    range_o = foreach(o = from.range, .combine = rbind, .export = ls(globalenv())) %dopar% {
      library(data.table)
      library(tidyverse)
      # calculate scores for etki and fake etki under threshold o
      # i.e. if we only consider NCI60 derivative transformed-rho > o, and assign confidence score based on o
      #      what would be the confidence score for each SLC-metabolite pair?
      ecs = ConfidenceScore_C_by_set(data = etki_Tr, ki = ki.data, o = o, p = Inf, q = Inf, m = 1, score.only = TRUE)
      
      fecs = lapply(1:100, function(f){
        simuScore = ConfidenceScore_C_by_set(data = si.data[Iteration == f], ki = ki.data, o = o, p = Inf, q = Inf, m = 1, score.only = TRUE)
        message(f)
        return(simuScore)
      })
      
      # calculate fraction recovery (TP - FP) for each given score cutoff
      # i.e. how many SLC-metabolite pair can be discovered, if we set score cutoff from 0 (all) - 60 (none)
      fr = lapply(0:60, function(cutoff){
        TP = length(ecs[ecs >= cutoff])/length(ecs)
        FP = sapply(fecs, function(f){length(f[f >= cutoff])}/length(f))
        
        if (TP == 0 & mean(FP) == 0){
          pvalue = 1
        }else if (TP == 1 & mean(FP) == 1) {
          pvalue = 1
        }else{
          pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
        }
        
        return(
          data.table(
            Cutoff = cutoff,
            TP = TP,
            FP.mean = mean(FP),
            FP.sd = sd(FP),
            TTP.literal = TP - mean(FP),
            TTP.in.sd = (TP - mean(FP))/sd(FP),
            pvalue = pvalue
          )
        )
      }) |> rbindlist() |>
        mutate(
          Threshold = o,
          .before = everything()
        )
      
      return(fr)
      
    }
    message("Range o finished at ", Sys.time())
    
    # p: derivative range for CCLE2019
    range_p = foreach(p = from.range, .combine = rbind, .export = ls(globalenv())) %dopar% {
      library(data.table)
      library(tidyverse)
      # calculate scores for etki and fake etki under threshold p
      # i.e. if we only consider CCLE2019 transformed-rho > p, and assign confidence score based on p
      #      what would be the confidence score for each SLC-metabolite pair?
      ecs = ConfidenceScore_C_by_set(data = etki_Tr, ki = ki.data, o = Inf, p = p, q = Inf, m = 1, score.only = TRUE)
      fecs = lapply(1:100, function(f){
        simuScore = ConfidenceScore_C_by_set(data = si.data[Iteration == f], ki = ki.data, o = Inf, p = p, q = Inf, m = 1, score.only = TRUE)
        return(simuScore)
      })
      
      # calculate fraction recovery (TP - FP) for each given score cutoff
      # i.e. how many SLC-metabolite pair can be discovered, if we set score cutoff from 0 (all) - 34 (none)
      fr = lapply(0:60, function(cutoff){
        TP = length(ecs[ecs >= cutoff])/length(ecs)
        FP = sapply(fecs, function(f){length(f[f >= cutoff])}/length(f))
        
        if (TP == 0 & mean(FP) == 0){
          pvalue = 1
        }else if (TP == 1 & mean(FP) == 1) {
          pvalue = 1
        }else{
          pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
        }
        
        return(
          data.table(
            Cutoff = cutoff,
            TP = TP,
            FP.mean = mean(FP),
            FP.sd = sd(FP),
            TTP.literal = TP - mean(FP),
            TTP.in.sd = (TP - mean(FP))/sd(FP),
            pvalue = pvalue
          )
        )
      }) |> rbindlist() |>
        mutate(
          Threshold = p,
          .before = everything()
        )
      
      return(fr)
      
    }
    message("Range p finished at ", Sys.time())
    
    # q: derivative range for CCL180
    range_q = foreach(q = from.range, .combine = rbind, .export = ls(globalenv())) %dopar% {
      library(data.table)
      library(tidyverse)
      # calculate scores for etki and fake etki under threshold c
      # i.e. if we only consider CCLE2019 transformed-rho > c, and assign confidence score based on c
      #      what would be the confidence score for each SLC-metabolite pair?
      ecs = ConfidenceScore_C_by_set(data = etki_Tr, ki = ki.data, o = Inf, p = Inf, q = q, m = 1, score.only = TRUE)
      fecs = lapply(1:100, function(f){
        simuScore = ConfidenceScore_C_by_set(data = si.data[Iteration == f], ki = ki.data, o = Inf, p = Inf, q = q, m = 1, score.only = TRUE)
        return(simuScore)
      })
      
      # calculate fraction recovery (TP - FP) for each given score cutoff
      # i.e. how many SLC-metabolite pair can be discovered, if we set score cutoff from 0 (all) - 34 (none)
      fr = lapply(0:34, function(cutoff){
        TP = length(ecs[ecs >= cutoff])/length(ecs)
        FP = sapply(fecs, function(f){length(f[f >= cutoff])}/length(f))
        
        if (TP == 0 & mean(FP) == 0){
          pvalue = 1
        }else if (TP == 1 & mean(FP) == 1) {
          pvalue = 1
        }else{
          pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
        }
        
        return(
          data.table(
            Cutoff = cutoff,
            TP = TP,
            FP.mean = mean(FP),
            FP.sd = sd(FP),
            TTP.literal = TP - mean(FP),
            TTP.in.sd = (TP - mean(FP))/sd(FP),
            pvalue = pvalue
          )
        )
      }) |> rbindlist() |>
        mutate(
          Threshold = q,
          .before = everything()
        )
      
      return(fr)
      
    }
    
    message("Range q finished at ", Sys.time())
    
    stopCluster(cl)
    
    # Make decisions on o, p, and q based on range report
    o.report = range_o[order(TTP.literal, decreasing = TRUE)][1, Threshold]
    p.report = range_p[order(TTP.literal, decreasing = TRUE)][1, Threshold]
    q.report = range_q[order(TTP.literal, decreasing = TRUE)][1, Threshold]
    
    return(list(o = o.report, p = p.report, q = q.report, o.record = range_o, p.record = range_p, q.record = range_q))
  }
  
  ConfC.threshold.detection = Best_threshold_detect_C_by_set(ki.data = etki_Tr, si.data = fake_etki_Tr, from.range = seq(0,1,0.1))
  
  #### 4.2 Determine score reward for adjacency matrix ####
  Best_scalar_detect_C_by_set = function(ki.data, si.data, o, p, q, from.range){
    # Setting up parallel
    numCores = 8  # Use all available cores except one for other processes
    cl = makeCluster(numCores)
    
    # m: range of scalar reward for Adjacency
    registerDoParallel(cl)
    range_m = foreach(m = from.range, .combine = rbind, .export = ls(globalenv())) %dopar% {
      library(data.table)
      library(tidyverse)
      # calculate scores for etki and fake etki given scalar m
      # i.e. if we only consider transformed-rho > o, p, q, and assign confidence score scaled by m
      #      what would be the confidence score for each SLC-metabolite pair?
      ecs = ConfidenceScore_C_by_set(data = etki_Tr, ki = ki.data, o = o, p = p, q = q, m = m, score.only = TRUE)
      fecs = lapply(1:100, function(f){
        simuScore = ConfidenceScore_C_by_set(data = si.data[Iteration == f], ki = ki.data, o = o, p = p, q = q, m = m, score.only = TRUE)
        return(simuScore)
      })
      
      # calculate fraction recovery (TP - FP) for each given score cutoff
      # i.e. how many SLC-metabolite pair can be discovered, if we set score cutoff from 0 (all) - 12 (none)
      fr = lapply(0:60, function(cutoff){
        TP = length(ecs[ecs >= cutoff])/length(ecs)
        FP = sapply(fecs, function(f){length(f[f >= cutoff])}/length(f))
        
        if (TP == 0 & mean(FP) == 0){
          pvalue = 1
        }else if (TP == 1 & mean(FP) == 1) {
          pvalue = 1
        }else{
          pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
        }
        
        return(
          data.table(
            Cutoff = cutoff,
            TP = TP,
            FP.mean = mean(FP),
            FP.sd = sd(FP),
            TTP.literal = TP - mean(FP),
            TTP.in.sd = (TP - mean(FP))/sd(FP),
            pvalue = pvalue
          )
        )
      }) |> rbindlist() |>
        mutate(
          Scalar = m,
          .before = everything()
        )
      
      return(fr)
      
    }
    
    message("Range m finished at ", Sys.time())
    
    stopCluster(cl)
    
    # Make decisions based on m range report
    m.report = range_m[order(TTP.literal, decreasing = TRUE)][1, Scalar]
    
    return(list(m = m.report, m.record = range_m))
  }
  
  ConfC.scalar.detection = Best_scalar_detect_C_by_set(ki.data = etki_Tr, si.data = fake_etki_Tr,
                                                       o = ConfC.threshold.detection$o,
                                                       p = ConfC.threshold.detection$p,
                                                       q = ConfC.threshold.detection$q,
                                                       from.range = seq(0,1,0.1))
  
  ### §5 Finalising the ConfABC computation ####
  ConstantIndex = c(
    # Threshold for ConfA, scaled by 3
    a = ConfA.threshold.detection$a,
    b = ConfA.threshold.detection$b,
    c = ConfA.threshold.detection$c,
    # Threshold for ConfB
    d = ConfB.threshold.detection$d,
    # Scalar for ConfB
    e = ConfB.scalar.detection$e,
    # Threshold for ConfC
    o = ConfC.threshold.detection$o,
    p = ConfC.threshold.detection$p,
    q = ConfC.threshold.detection$q,
    # Scalar for ConfC
    m = ConfC.scalar.detection$m
  )
  
  output = list(etki_Tr = etki_Tr, etki_Te = etki_Te, constant = ConstantIndex)
  saveRDS(output, paste0("AutoOpt/output", Opt_time, ".rds"))
  message("Optimisation ", Opt_time, "finished at ", Sys.time())
}

## 9.2 Testing ####
ConfidenceScore_ABC_by_set = function(data, ki, constant, score.only = FALSE){
  ## quantile ##
  ns_dist1 = na.omit(ki$NCI60) # NCI60 #
  ns_dist2 = na.omit(ki$CCLE2019) # CCLE2019 #
  ns_dist3 = na.omit(ki$CCL180) # CCL180 #
  ns_dist4 = -log10(na.omit(ki$CRISPRDEP)) # CRISPR dependency
  q1 = quantile(ns_dist1[ns_dist1 > as.numeric(constant["a"])], seq(0, 1, by = 0.1))
  q2 = quantile(ns_dist2[ns_dist2 > as.numeric(constant["b"])], seq(0, 1, by = 0.1))
  q3 = quantile(ns_dist3[ns_dist3 > as.numeric(constant["c"])], seq(0, 1, by = 0.1))
  q4 = quantile(ns_dist4[ns_dist4 > -log10(constant["d"])], seq(0, 1, by = 0.1))
  q5 = quantile(ns_dist1[ns_dist1 > as.numeric(constant["o"])], seq(0, 1, by = 0.1))
  q6 = quantile(ns_dist2[ns_dist1 > as.numeric(constant["p"])], seq(0, 1, by = 0.1))
  q7 = quantile(ns_dist3[ns_dist1 > as.numeric(constant["q"])], seq(0, 1, by = 0.1))
  
  data_info = data.table(
    data_name = c("NCI60", "CCLE2019", "CCL180"),
    FactorA_dist = c("q1", "q2", "q3"),
    FactorC_dist = c("q5", "q6", "q7"),
    FactorA_threshold = c(as.numeric(constant["a"]), as.numeric(constant["b"]), as.numeric(constant["c"])),
    FactorC_threshold = c(as.numeric(constant["o"]), as.numeric(constant["p"]), as.numeric(constant["q"]))
  )
  
  cl = makeCluster(8)
  registerDoParallel(cl)
  
  Score = foreach(i = 1:nrow(data), .combine = c, .export = unique(c(ls(globalenv()), ls(environment())))) %dopar% {
    library(data.table)
    library(tidyverse)
    
    # Factor A #
    Factor_A = vector("integer", nrow(data_info))
    for (j in 1:nrow(data_info)){
      n = as.numeric(data[i, get(data_info[j, data_name])])
      if (is.na(n) | n < data_info[j, FactorA_threshold]){
        Factor_A[j] = 0
      } else if (n > data_info[j, FactorA_threshold]){
        Factor_A[j] = findInterval(n, get(data_info[j, FactorA_dist]))*3
      }
    }
    Factor_A = as.numeric(sum(Factor_A))
    
    # Factor B #
    if (is.na(data[i, CRISPRDEP]) | data[i, -log10(as.numeric(CRISPRDEP))] < -log10(constant["d"])){
      Factor_B = 0
    } else if (data[i, -log10(as.numeric(CRISPRDEP))] > -log10(constant["d"])){
      Factor_B = findInterval(data[i, -log10(as.numeric(CRISPRDEP))], q4)*as.numeric(constant["e"])
    }
    
    # Factor C #
    if (is.na(data$Formula[i])){
      temp = getAM(data$KEGG_ID[i], is.kegg = TRUE)
    }else{
      temp = getAM(data$Formula[i])
    }
    if (check_NA(temp)){
      Factor_C = 0
    }else if (!1 %in% names(temp) & !2 %in% names(temp)){
      Factor_C = 0
    }else{
      proximal = translate(temp[which(names(temp) %in% c(1,2))])
      if (check_NA(proximal)){
        Factor_C = 0
      }else{
        proximal_ns = getNR(data[i, SYMBOL], proximal)
        if (nrow(proximal_ns) != 0) {
          Factor_C = 0
          for (l in 1:nrow(data_info)){
            N = unique(proximal_ns[, get(data_info[l, data_name])])
            for (n in N){
              if (is.na(n) | n < data_info[l, FactorC_threshold]){
                Factor_C = Factor_C + 0
              }else if (n > data_info[l, FactorC_threshold]){
                Factor_C = Factor_C + findInterval(n, get(data_info[l, FactorC_dist]))*constant["m"]
              }
            }
          }
        }else{
          Factor_C = 0
        }
      }
    }
    
    # Final score #
    finalscore = Factor_A + Factor_B + Factor_C
    return(finalscore)
  }
  
  stopCluster(cl)
  
  if (score.only){
    return(as.numeric(Score))
  }else{
    score_table = data |>
      mutate(
        Score = as.numeric(Score)
      )
    return(score_table)
  }
}
analyseDist = function(data.true, data.simu){
  data.true = sapply(1:100, function(s){
    resample = sample(data.true, size = length(data.true), replace = TRUE)
    return(mean(resample))
  })
  data.simu = sapply(data.simu, mean)
  
  Dist = data.table(
    type = rep(c("simulated", "known"), each = 100),
    Score = c(data.simu, data.true)
  )
  
  rain_height = .1
  colorset = c("known" = "#990000", "simulated" = "grey")
  g = ggplot(Dist, aes(x = "", y = Score, fill = type)) +
    # mountain
    introdataviz::geom_flat_violin(trim=FALSE, alpha = 0.4, position = position_nudge(x = -0.6)) +
    # rain
    geom_point(aes(colour = type), size = 1, alpha = 0.4, show.legend = FALSE, position = position_jitter(width = rain_height, height = 0)) +
    geom_boxplot(width = rain_height, alpha = 0.4, show.legend = FALSE, outlier.shape = NA) +
    scale_y_continuous(name = "Average Confidence Score",
                       breaks = seq(15, 40, 5), 
                       limits = c(15, 40)) +
    # Manual color settings
    scale_fill_manual(values = colorset) + # For violin plots
    scale_color_manual(values = colorset) + # For points
    coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
  
  p = t.test(
    x = Dist[type == "simulated", Score],
    y = Dist[type == "known", Score],
    alternative = "less"
  )$p.value
  
  return(list(Dist.data = Dist, Dist.plot = g, pvalue = p))
}


### §1 Read saved Tr/Te data ####
TrTe_list = list.files(path = "Results/AutoOpt/")
TrTe = lapply(
  TrTe_list, 
  function(path){readRDS(paste0("Results/AutoOpt/", path))}
)

### §2 Validation ####
Name_collector = c()
Dist_collector = vector("list", 10)
Discoverable_interactions_number_collector = c()
Discoverable_interactions_within_top_collector = c()

for (t in 1:10){
  # Retrieve the etki data #
  TrTe_name = str_remove(TrTe_list[t], ".rds")
  etki_Tr = TrTe[[t]]$etki_Tr
  etki_Te = TrTe[[t]]$etki_Te
  fake_etki_Tr = TrTe[[t]]$fake_etki_Tr
  ConstantIndex = TrTe[[t]]$constant
  Name_collector = c(Name_collector, TrTe_name)
  
  # Calculate the scores #
  optimised_real = ConfidenceScore_ABC_by_set(data = etki_Tr, ki = etki_Tr, constant = ConstantIndex, score.only = TRUE)
  optimised_fake = lapply(1:100, function(f){
    simuScore = ConfidenceScore_ABC_by_set(data = fake_etki_Tr[Iteration == f], ki = etki_Tr, constant = ConstantIndex, score.only = TRUE)
    message(f)
    return(simuScore)
  })
  
  # Reporter 1: Distributions #
  Dist.object = analyseDist(data.true = optimised_real, data.simu = optimised_fake)
  Dist_collector[[t]] = Dist.object
  
  # Decide the best rank cutoff to use #
  fr = lapply(0:150, function(cutoff){
    TP = length(optimised_real[optimised_real >= cutoff])/length(optimised_real)
    FP = sapply(optimised_fake, function(f){length(f[f >= cutoff])/length(f)})
    
    if (TP == 0 & mean(FP) == 0){
      pvalue = 1
    }else if (TP == 1 & mean(FP) == 1) {
      pvalue = 1
    }else{
      pvalue = t.test(x = FP, mu = TP, alternative = "less")$p.value
    }
    
    return(
      data.table(
        Cutoff = cutoff,
        TP = TP,
        FP.mean = mean(FP),
        FP.sd = sd(FP),
        TTP.literal = TP - mean(FP),
        TTP.in.sd = (TP - mean(FP))/sd(FP),
        pvalue = pvalue
      )
    )
  }) |> rbindlist()
  
  optimised_cutoff = fr[order(TTP.literal, decreasing = TRUE)][1, Cutoff]
  
  # Testing #
  optimised_testing = ConfidenceScore_ABC_by_set(data = etki_Te, ki = etki_Tr, constant = ConstantIndex, score.only = FALSE)
  
  # Reporter 2: How many fractions are above the optimised confidence score cutoff? #
  Discoverable_interactions_number = optimised_testing[Score >= optimised_cutoff, .N]/optimised_testing[, .N]
  Discoverable_interactions_number_collector = c(Discoverable_interactions_number_collector, Discoverable_interactions_number)
  
  # Reporter 3: How many predictions can be found in the top 4% of each SLC against 1980 metabolites? #
  SLC_being_tested = unique(optimised_testing$SYMBOL)
  RankPer_record = lapply(SLC_being_tested, function(slc){
    et_atomic = et[SYMBOL == slc]
    et_atomic_score = ConfidenceScore_ABC_by_set(data = et_atomic, ki = etki_Tr, constant = ConstantIndex, score.only = FALSE)
    message(slc)
    return(et_atomic_score)
  }) |> rbindlist()
  RankPer_testing_set = sapply(1:nrow(optimised_testing), function(i){
    RankPer_record_atomic = RankPer_record[SYMBOL == optimised_testing[i, SYMBOL]][order(Score, decreasing = TRUE)]
    first_zero = min(which(RankPer_record_atomic$Score == 0))
    RankPer_record_atomic = RankPer_record_atomic |>
      mutate(
        RankPer = c((1:first_zero)/first_zero, rep(1, times = nrow(RankPer_record_atomic) - first_zero))
      )
    return(RankPer_record_atomic[Formula == optimised_testing[i, Formula], RankPer])
  })
  Discoverable_interactions_within_top = length(RankPer_testing_set[RankPer_testing_set <= 0.03])/length(RankPer_testing_set)
  Discoverable_interactions_within_top_collector = c(Discoverable_interactions_within_top_collector, Discoverable_interactions_within_top)
  
  # Save image #
  save.image(paste0("03 Results/05 Building ConfidenceScore/TrTe/Validated ", TrTe_name, ".RData"))
}

names(Dist_collector) = Name_collector
names(Discoverable_interactions_number_collector) = Name_collector
names(Discoverable_interactions_within_top_collector) = Name_collector

# Summarise discovery in testing sets
Discoverable_summary = data.table(
  Iteration = as.numeric(str_remove(Name_collector, "output")),
  Discoverable.total = Discoverable_interactions_number_collector,
  Discoverable.in.top = Discoverable_interactions_within_top_collector
) |>
  mutate(
    Discoverable.out.top = Discoverable.total - Discoverable.in.top,
    Not.Discoverable = 1-Discoverable.total
  ) |>
  melt(id.vars = "Iteration")

# --------------------------- #
#       Figure 4D data        #
# --------------------------- #

ggplot(Discoverable_summary[variable != "Discoverable.total"],
       aes(fill = factor(variable, levels = c("Not.Discoverable", "Discoverable.out.top", "Discoverable.in.top")), 
           x = factor(Iteration, levels = seq(10,1,-1)), 
           y = value)) + 
  scale_fill_manual(values = c("lightgrey", "pink", "#990000")) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = "Iteration", y = "Fraction") +
  scale_y_continuous(limits = c(0,1.0), breaks = seq(0,1.0,0.2)) +
  theme(
    plot.title = element_text(size=11),legend.position="none",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  ) + 
  geom_hline(yintercept = Discoverable_summary[variable == "Discoverable.in.top", mean(value)], col = "black", linetype="dashed") +
  coord_flip()

# --------------------------- #
#       Figure S4 data        #
# --------------------------- #

for (i in 1:10){
  print(Dist_collector[[i]]$Dist.plot)
  print(Dist_collector[[i]]$Dist.data)
}

# 10 Calculating Confidence score for all pairs ####
constant = c(
  a = 0.00, b = 0.00, c = 1.00, 
  d = 0.16, e = 1.00,
  o = 0.10, p = 0.10, q = 0.20, m = 1.00
)

ns_dist1 = na.omit(etki$NCI60) # NCI60 #
ns_dist2 = na.omit(etki$CCLE2019) # CCLE2019 #
ns_dist3 = na.omit(etki$CCL180) # CCL180 #
ns_dist4 = -log10(na.omit(etki$CRISPRDEP)) # CRISPR dependency
q1 = quantile(ns_dist1[ns_dist1 > as.numeric(constant["a"])], seq(0, 1, by = 0.1))
q2 = quantile(ns_dist2[ns_dist2 > as.numeric(constant["b"])], seq(0, 1, by = 0.1))
q3 = quantile(ns_dist3[ns_dist3 > as.numeric(constant["c"])], seq(0, 1, by = 0.1))
q4 = quantile(ns_dist4[ns_dist4 > -log10(constant["d"])], seq(0, 1, by = 0.1))
q5 = quantile(ns_dist1[ns_dist1 > as.numeric(constant["o"])], seq(0, 1, by = 0.1))
q6 = quantile(ns_dist2[ns_dist1 > as.numeric(constant["p"])], seq(0, 1, by = 0.1))
q7 = quantile(ns_dist3[ns_dist1 > as.numeric(constant["q"])], seq(0, 1, by = 0.1))

data_info = data.table(
  data_name = c("NCI60", "CCLE2019", "CCL180"),
  FactorA_dist = c("q1", "q2", "q3"),
  FactorC_dist = c("q5", "q6", "q7"),
  FactorA_threshold = c(as.numeric(constant["a"]), as.numeric(constant["b"]), as.numeric(constant["c"])),
  FactorC_threshold = c(as.numeric(constant["o"]), as.numeric(constant["p"]), as.numeric(constant["q"]))
)

ConfidenceScore_ABC = function(row, ki, score.only = FALSE){
  data_info = data.table(
    data_name = c("NCI60", "CCLE2019", "CCL180"),
    FactorA_dist = c("q1", "q2", "q3"),
    FactorC_dist = c("q5", "q6", "q7"),
    FactorA_threshold = c(as.numeric(constant["a"]), as.numeric(constant["b"]), as.numeric(constant["c"])),
    FactorC_threshold = c(as.numeric(constant["o"]), as.numeric(constant["p"]), as.numeric(constant["q"]))
  )
  
  # Factor A #
  Factor_A = vector("integer", nrow(data_info))
  for (j in 1:nrow(data_info)){
    n = as.numeric(row[, get(data_info[j, data_name])])
    if (is.na(n) | n < data_info[j, FactorA_threshold]){
      Factor_A[j] = 0
    } else if (n > data_info[j, FactorA_threshold]){
      Factor_A[j] = findInterval(n, get(data_info[j, FactorA_dist]))*3
    }
  }
  Factor_A = as.numeric(sum(Factor_A))
  
  # Factor B #
  if (is.na(row[, CRISPRDEP]) | row[, -log10(as.numeric(CRISPRDEP))] < -log10(constant["d"])){
    Factor_B = 0
  } else if (row[, -log10(as.numeric(CRISPRDEP))] > -log10(constant["d"])){
    Factor_B = findInterval(row[, -log10(as.numeric(CRISPRDEP))], q4)*as.numeric(constant["e"])
  }
  
  # Factor C #
  if (is.na(row$Formula)){
    temp = getAM(row$KEGG_ID, is.kegg = TRUE)
  }else{
    temp = getAM(row$Formula)
  }
  if (check_NA(temp)){
    Factor_C = 0
  }else if (!1 %in% names(temp) & !2 %in% names(temp)){
    Factor_C = 0
  }else{
    proximal = translate(temp[which(names(temp) %in% c(1,2))])
    if (check_NA(proximal)){
      Factor_C = 0
    }else{
      proximal_ns = getNR(row[, SYMBOL], proximal)
      if (nrow(proximal_ns) != 0) {
        Factor_C = 0
        for (l in 1:nrow(data_info)){
          N = unique(proximal_ns[, get(data_info[l, data_name])])
          for (n in N){
            if (is.na(n) | n < data_info[l, FactorC_threshold]){
              Factor_C = Factor_C + 0
            }else if (n > data_info[l, FactorC_threshold] & n < data_info[l, FactorA_threshold]){
              Factor_C = Factor_C + constant["m"]/2
            }else if (n > data_info[l, FactorA_threshold]){
              Factor_C = Factor_C + findInterval(n, get(data_info[l, FactorC_dist]))*constant["m"]
            }
          }
        }
      }else{
        Factor_C = 0
      }
    }
  }
  
  # Final score #
  finalscore = Factor_A + Factor_B + Factor_C
  return(list(A = Factor_A, B = Factor_B, C = as.numeric(Factor_C), Score = as.numeric(finalscore)))
}

et_ConfidenceScore = data.table(et.tra, 
                                A = vector("numeric", nrow(et.tra)),
                                B = vector("numeric", nrow(et.tra)),
                                C = vector("numeric", nrow(et.tra)),
                                Score = vector("numeric", nrow(et.tra)))
pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                      total = nrow(et_ConfidenceScore),
                      complete = "=",   # Completion bar character
                      incomplete = "-", # Incomplete bar character
                      current = ">",    # Current bar character
                      clear = TRUE,    # Clears the bar when finish
                      width = 100)      # Width of the progress bar
for (i in 1:nrow(et_ConfidenceScore)){
  et_ConfidenceScore[i, c("A", "B", "C", "Score") := ConfidenceScore_ABC(et_ConfidenceScore[i], ki = etki)]
  pb$tick()
}

fwrite(et_ConfidenceScore, "Results/Evidence table (Confidence Score).csv")

# 11 Survey Confidence score ####
etcs = fread("Results/Evidence table (Confidence Score new).csv")
known_list = fread("Raw Data/Known interaction.csv")[1:837]

etcs_ki = known_list[!is.na(Formula)]
colnames(etcs_ki)[1] = "SYMBOL"
etcs_ki = data.table(etcs_ki,
                     RankPer = vector("numeric", nrow(etcs_ki)),
                     A = vector("numeric", nrow(etcs_ki)),
                     B = vector("numeric", nrow(etcs_ki)),
                     C = vector("numeric", nrow(etcs_ki)),
                     Score = vector("numeric", nrow(etcs_ki)))
for (i in 1:nrow(etcs_ki)){
  if (etcs_ki[i, SYMBOL] %in% unique(etcs$SYMBOL)){
    #etcs_atomic = etcs[SYMBOL == etcs_ki$SYMBOL[i]][order(Score, A, decreasing = TRUE)]
    etcs_atomic = etcs[SYMBOL == etcs_ki$SYMBOL[i]][order(Score, decreasing = TRUE)]
    etcs_ki[i, c("A", "B", "C", "Score") := etcs_atomic[Formula == etcs_ki$Formula[i], c("A", "B", "C", "Score")]]
    # Find when score starts to equal to 0
    etcs_atomic_first_zero = min(which(etcs_atomic$Score == 0))
    # Calculate a rank percentile
    etcs_atomic = etcs_atomic |>
      mutate(
        RankPer = case_when(
          Score != 0 ~ which(etcs_atomic$Formula == etcs_ki$Formula[i])/(etcs_atomic_first_zero - 1),
          Score == 0 ~ 1
        )
      )
    etcs_ki[i, RankPer := etcs_atomic[Formula == etcs_ki$Formula[i], RankPer]]
  }else{
    etcs_ki[i, RankPer := NaN]
  }
}

etcs_ki = etcs_ki[!is.na(RankPer)]

## Generate simulated etcs_ki ##
fake_collection = vector("list", 100)
etcs_ki_gene_pool = unique(etcs_ki$SYMBOL)
etcs_ki_metabolite_pool = unique(etcs_ki$Formula)

for (sample_time in 1:100){
  # initialise fake etcs_ki data table #
  fake_etcs_ki = setNames(data.table(matrix(ncol = ncol(etcs_ki), nrow = nrow(etcs_ki))), colnames(etcs_ki))
  fake_etcs_ki[, c("SYMBOL", "Gene_type","Localisation") := as.list(etcs_ki[, c("SYMBOL", "Gene_type","Localisation")])]
  fake_etcs_ki[, Status := rep("Fake", nrow(fake_etcs_ki))]
  fake_etcs_ki = fake_etcs_ki |>
    select(-Substrate) |>
    mutate(
      Formula = as.character(Formula),
      KEGG_ID = as.character(KEGG_ID),
      Score = as.numeric(Score),
      RankPer = as.numeric(RankPer)
    )
  
  # for each SLC, shuffle fake substrate formula #
  for (j in 1:nrow(fake_etcs_ki)){
    if (!fake_etcs_ki[j, SYMBOL] %in% unique(etcs[, SYMBOL])){
      fake_etcs_ki[j, Rank := NaN]
    }else{
      true_substrate = etcs_ki[SYMBOL == fake_etcs_ki$SYMBOL[j]][, Formula]
      fake_substrate_pool = etcs_ki_metabolite_pool[!etcs_ki_metabolite_pool %in% true_substrate]
      fake_substrate = sample(fake_substrate_pool, 1)
      while (nrow(etcs[SYMBOL == fake_etcs_ki$SYMBOL[j] & Formula == fake_substrate]) == 0){
        fake_substrate = sample(fake_substrate_pool, 1)
      }
      fake_etcs_ki[j, Formula := fake_substrate]
      etcs_atomic = etcs[SYMBOL == fake_etcs_ki$SYMBOL[j]][order(Score, decreasing = TRUE)]
      # Find when score starts to equal to 0
      etcs_atomic_first_zero = min(which(etcs_atomic$Score == 0))
      # Calculate a rank percentile
      etcs_atomic = etcs_atomic |>
        mutate(
          RankPer = case_when(
            Score != 0 ~ which(etcs_atomic$Formula == fake_etcs_ki$Formula[j])/(etcs_atomic_first_zero - 1),
            Score == 0 ~ 1
          )
        )
      fake_etcs_ki[j, RankPer := etcs_atomic[Formula == fake_etcs_ki$Formula[j], RankPer]]
    }
  }
  fake_collection[[sample_time]] = fake_etcs_ki
  print(paste(sample_time, "done"))
}

# Distributions
analyseDist = function(data.simu, data.true){
  # Bootstrap
  boot = sapply(1:100, function(s){
    resample = sample(data.true[, as.numeric(RankPer)], size = nrow(data.true), replace = TRUE)
    return(median(resample))
  })
  
  Dist = data.table(
    type = rep(c("simulated", "known"), each = 100),
    RankPer = c(data.simu[, median(as.numeric(RankPer)), by = Iteration][, V1], boot)
  )
  
  rain_height = .1
  colorset = c("known" = "#990000", "simulated" = "grey")
  g = ggplot(Dist, aes(x = "", y = RankPer, fill = type)) +
    # mountain
    introdataviz::geom_flat_violin(trim=FALSE, alpha = 0.4, position = position_nudge(x = -0.6)) +
    # rain
    geom_point(aes(colour = type), size = 1, alpha = 0.4, show.legend = FALSE, position = position_jitter(width = rain_height, height = 0)) +
    geom_boxplot(width = rain_height, alpha = 0.4, show.legend = FALSE, outlier.shape = NA) +
    scale_y_continuous(name = "Median Rank Percentile",
                       breaks = seq(0.2, 1, 0.1), 
                       limits = c(0.2, 1)) +
    # Manual color settings
    scale_fill_manual(values = colorset) + # For violin plots
    scale_color_manual(values = colorset) + # For points
    coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
  
  p = t.test(
    x = Dist[type == "simulated", RankPer],
    y = Dist[type == "known", RankPer],
    alternative = "greater"
  )$p.value
  
  return(list(dist.data = Dist, dist.plot = g, dist.pvalue = p))
}

fake_collection = lapply(1:100, function(f){
  data = fake_collection[[f]] |>
    mutate(
      Iteration = f,
      .before = everything()
    ) |>
    mutate(
      RankPer = as.numeric(RankPer)
    )
  return(data)
}) |> rbindlist()

# --------------------------- #
#       Figure 4A data        #
# --------------------------- #

rank.object = analyseDist(data.simu = fake_collection[!is.na(RankPer)], data.true = etcs_ki[!is.na(RankPer)])
rank.object$dist.plot
rank.object$dist.pvalue
rank.object$dist.data

# --------------------------- #
#       Figure 4B data        #
# --------------------------- #
fr = lapply(seq(0.001,1,0.001), function(x){
  fr_true = nrow(etcs_ki[RankPer <= x])/nrow(etcs_ki)
  fr_simu = sapply(1:100, function(y){nrow(fake_collection[Iteration == y & RankPer <= x])/nrow(fake_collection[Iteration == y])})
  message(x)
  return(
    data.table(
      Th = x,
      fr = fr_true - fr_simu
    )
  )
}) |> rbindlist()

ggplot() +
  geom_violin(fr, mapping = aes(x = Th, y = fr, group = Th, alpha = 0.2), color = "grey", fill = "grey") +
  geom_line(data.table(Th = seq(0.001,1,0.001), fr = fr[, median(fr), by = Th][, V1]), mapping = aes(x = Th, y = fr)) +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))


fr_median = data.table(Th = seq(0.001,1,0.001), fr = fr[, median(fr), by = Th][, V1])

real = nrow(etcs_ki[RankPer <= 0.253])/nrow(etcs_ki)
fake = sapply(1:100, function(y){nrow(fake_collection[Iteration == y & RankPer <= 0.253])/nrow(fake_collection[Iteration == y])})
median(fake)
sd(fake)

## use rankper < 0.253 ##

# 12 Over-representation analysis ####
library(KEGGREST)

## 12.1 generate list of human metabolic pathways and linked compound to each pathways #####
# get human pathways (metabolic pathway definition from CCL180 paper)
human_metabolic_pathways = fread("01 Raw Data/05 Building ConfidenceScore/Pathways.csv")

# in data table format for re-visit
human_metabolic_pathways = data.table(human_metabolic_pathways,
                                      linked_compound = vector("character", length = nrow(human_metabolic_pathways)),
                                      compound_no = vector("numeric", length = nrow(human_metabolic_pathways))
)
for (i in 1:nrow(human_metabolic_pathways)){
  q = keggGet(human_metabolic_pathways[i, "Pathway id"])[[1]][["COMPOUND"]]
  id = paste(names(q), collapse = ";")
  no = length(q)
  human_metabolic_pathways[i, linked_compound := id]
  human_metabolic_pathways[i, compound_no := no]
  print(i)
}

# in list format for ORA (function category list for ORA)
human_metabolic_pathways_list = vector("list", nrow(human_metabolic_pathways))
names(human_metabolic_pathways_list) = unlist(human_metabolic_pathways[, "Pathway id"])
for (i in 1:length(human_metabolic_pathways_list)){
  human_metabolic_pathways_list[[i]] = names(keggGet(names(human_metabolic_pathways_list)[i])[[1]][["COMPOUND"]])
  print(i)
}

## 12.2 gather feature parameter for ORA #####
known_list = fread("Raw Data/Known interaction.csv")[1:837][!is.na(Status)]
orphan_list = fread("Raw Data/SLC transcript.csv", select = c(1:7))[Knowledge_type == "Orphan"]

## 12.3 ORA #####
# Lazy fisher's exact test
# test against the hypothesis that the feature does not positively linked to category
fisherExactTest = function(features, category, universe){
  contingencyTable = matrix(nrow = 2, ncol = 2)
  rownames(contingencyTable) = c("inCategory", "outCategory")
  colnames(contingencyTable) = c("inFeature", "outFeature")
  # prepare values
  x = sum(features %in% category)
  k = length(category)
  m = length(features)
  N = length(universe)
  n = N - m
  # fill contingency table
  contingencyTable["inCategory",  "inFeature"]  = x
  contingencyTable["inCategory",  "outFeature"] = k - x
  contingencyTable["outCategory", "inFeature"]  = m - x
  contingencyTable["outCategory", "outFeature"] = n - (k - x)
  # apply fisher exact test
  p = stats::fisher.test(x = contingencyTable, alternative = "greater")$p.value 
  # there are MORE number in the overlap of feature and category (top left cells)
  return(p)
}
# Over-representation analysis
overRepresentationAnalysis = function(features, funCatList, universe,
                                      minSize = 1, maxSize = Inf,
                                      removeNoOverlap = TRUE, pAdjustMethod = "BH") {
  
  # if missing universe, use all features in the functional category list
  if(missing(universe)) {
    universe = unique(unlist(funCatList))
  }
  # remove categories out of size
  size = lapply(funCatList, length) %>%
    unlist()
  funCatList = funCatList[size >= minSize & size <= maxSize]
  # calculate overlap and p value for each funCat, returning a one-row data frame with results
  resDf = lapply(1:length(funCatList), function(x) {
    funCat = funCatList[[x]]
    funCat = funCat[funCat %in% universe]
    overlappingFeatures = features[features %in% funCat]
    nOverlap = length(overlappingFeatures)
    if(nOverlap == 0) {
      overlappingFeatures = "No overlap"
    } else {
      overlappingFeatures = paste(unique(overlappingFeatures), collapse = ", ")
    }
    # get enrichment p values
    fisherP = fisherExactTest(features = features, category = funCat, universe = universe)
    # return out df
    outDT = data.table(
      functionalCategory = names(funCatList)[x],
      overlap = overlappingFeatures, 
      nOverlap = nOverlap, 
      pValue = fisherP
    )
    return(outDT)
  }) |> rbindlist()
  # adjust p values
  if(removeNoOverlap) {
    resDf = subset(resDf, resDf$nOverlap != 0)
  }
  resDf[ , "pAdj"] = stats::p.adjust(resDf$pValue, method = pAdjustMethod)
  # return data frame
  return(resDf)
}
# toolkit
ConvertNames = function(name, source){
  # Process source input #
  if (is.null(source)){
    stop("Source of metabolite name required!")
  }
  else if (source %in% colnames(tb)){
    # Process name input #
    if (is.null(name)){
      stop("Name of metabolite required!")
    }else if (name %in% tb[,get(source)]){
      ans = tb[get(source) == name]
    }else{
      stop("Name entered not recognised!")
    }
  }else{
    stop("Source entered not recognised!")
  }
  
  # return ans as result #
  return(ans)
}
check_NA = function(x){
  if (length(x) == 1){
    if (is.na(x)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}

# Universe
linked_to_metabolomics = unique(str_trim(unlist(str_split(unlist(tb[, "KEGG ID"]), ";"))))
linked_to_pathway = unique(str_trim(unlist(str_split(unlist(human_metabolic_pathways[, linked_compound]), ";"))))
univ = intersect(linked_to_metabolomics, linked_to_pathway)

## build known list for ORA ##
known_list_ORA = data.table(SYMBOL = known_list[, unique(HGNC_symbol)],
                            Substrate = vector("character", length(known_list[, unique(HGNC_symbol)])),
                            KEGG_ID = vector("character", length(known_list[, unique(HGNC_symbol)])))
for (i in 1:nrow(known_list_ORA)){
  s = paste(known_list[HGNC_symbol == known_list_ORA[i, SYMBOL], Substrate], collapse = "; ")
  f = known_list[HGNC_symbol == known_list_ORA[i, SYMBOL], Formula]
  if (length(which(is.na(f))) == 0){
    k = paste(na.omit(unlist(lapply(f, function(x){as.character(ConvertNames(x, "Formula")[,"KEGG ID"])}))), collapse = "; ")
  }else{
    kna_YES = known_list[HGNC_symbol == known_list_ORA[i, SYMBOL] & is.na(Formula), KEGG_ID]
    kna_NO = paste(na.omit(unlist(lapply(f[-which(is.na(f))], function(x){as.character(ConvertNames(x, "Formula")[,"KEGG ID"])}))), collapse = "; ")
    k = paste(c(kna_YES, kna_NO), collapse = "; ")
  }
  if (nchar(k) == 8 & str_sub(k, 7,8) == "; "){
    k = str_remove(k, "; ")
  }
  known_list_ORA[i, Substrate := s]
  known_list_ORA[i, KEGG_ID := k]
}

bulk.overRepresentationAnalysis.known = function(known_list_ORA, up.to.rank.per, set.p.value = 0.05){
  # Calculate ORA
  known_list_ORA_result = vector("list", nrow(known_list_ORA))
  names(known_list_ORA_result) = known_list_ORA[, SYMBOL]
  for (i in 1:nrow(known_list_ORA)){
    etcs_atomic = etcs[SYMBOL == known_list_ORA[i, SYMBOL]][order(Score, decreasing = TRUE)][Score != 0]
    etcs_atomic = etcs_atomic[1:(nrow(etcs_atomic)*up.to.rank.per)]
    etcs_atomic_feature = sapply(
      etcs_atomic[,na.omit(Formula)], function(x){as.character(ConvertNames(x, "Formula")[,"KEGG ID"])}
    ) |>
      str_split(";") |>
      unlist() |>
      str_trim()
    etcs_atomic_feature = etcs_atomic_feature[etcs_atomic_feature %in% univ]
    ora_reporting = overRepresentationAnalysis(features = unique(etcs_atomic_feature), 
                                               funCatList = human_metabolic_pathways_list, 
                                               universe = univ) |> as.data.table()
    ora_reporting = inner_join(x = ora_reporting[pAdj <= set.p.value][order(pAdj)], 
                               y = human_metabolic_pathways[, .(functionalCategory = `Pathway id`, `Pathway name`, Class, compound_no)], 
                               by = "functionalCategory")
    setcolorder(ora_reporting, c("functionalCategory", "Pathway name", "Class", "pAdj", "compound_no", "nOverlap", "overlap"))
    
    known_list_ORA_result[[i]] = ora_reporting
  }
  # Specify result
  known_list_ORA = data.table(known_list_ORA,
                              Status = vector("character", nrow(known_list_ORA)))
  for (i in 1:length(known_list_ORA_result)){
    pathway = inner_join(x = data.table("Pathway id" = known_list_ORA_result[[i]][,functionalCategory]),
                         y = human_metabolic_pathways[, .(`Pathway id`, linked_compound)],
                         by = "Pathway id")
    if (nrow(pathway) == 0){
      known_list_ORA[i, Status := "No Pathway Enriched"]
    }else{
      pathway_linked_compound = str_split(unlist(pathway[, linked_compound]), ";")[[1]]
      substrate_linked_compound = known_list_ORA[i, KEGG_ID]
      if (nchar(substrate_linked_compound) == 0){
        known_list_ORA[i, Status := "No Substrate Found"]
      }else{
        substrate_linked_compound = str_split(substrate_linked_compound, "; ")[[1]]
        known_list_ORA[i, Status := length(intersect(pathway_linked_compound, substrate_linked_compound))]
      }
    }
  }
  return(known_list_ORA)
}

# --------------------------- #
#       Figure 4C data        #
# --------------------------- #

plot_table = lapply(seq(0.01, 0.30, 0.01), function(k){
  a = bulk.overRepresentationAnalysis.known(known_list_ORA = known_list_ORA, up.to.rank.per = k)
  n = a[Status != "No Pathway Enriched" & Status != "No Substrate Found" & Status != "0", .N]
  m = mean(as.numeric(a[Status != "No Pathway Enriched" & Status != "No Substrate Found" & Status != "0", Status]))
  message(k)
  return(
    data.table(
      RankPer = k,
      Number = n,
      Mean = m
    )
  )
}) |> rbindlist()

plot(x = plot_table[, RankPer], y = plot_table[, Number])
ggplot(plot_table, aes(x = RankPer, y = Number)) +
  geom_point() +
  scale_y_continuous(limits = c(48, 66), breaks = seq(45, 65, 5)) +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

# 13 Generating prediction for orphans ####
# Note: the data is provided in Table S10, here's demonstration of how it was produced

orphan_list_ORA = data.table(SYMBOL = orphan_list[, unique(HGNC_symbol)],
                             Enriched_pathway = vector("character", length(orphan_list[, unique(HGNC_symbol)])))
orphan_list_ORA_result = vector("list", nrow(orphan_list_ORA))
names(orphan_list_ORA_result) = orphan_list_ORA[, SYMBOL]
for (i in 1:nrow(orphan_list_ORA)){
  etcs_atomic = etcs[SYMBOL == orphan_list_ORA[i, SYMBOL]][order(Score, decreasing = TRUE)][Score != 0]
  etcs_atomic = etcs_atomic[1:(nrow(etcs_atomic)*0.03)]
  etcs_atomic_feature = 
    lapply(etcs_atomic[,na.omit(Formula)], function(x){as.character(ConvertNames(x, "Formula")[,"KEGG ID"])}) |>
    unlist() |>
    str_split(";") |>
    unlist() |>
    str_trim()
  
  ora_reporting = overRepresentationAnalysis(features = unique(etcs_atomic_feature), 
                                             funCatList = human_metabolic_pathways_list, 
                                             universe = univ)
  ora_reporting = inner_join(x = data.table(ora_reporting)[pAdj <= 0.05][order(pAdj)][,pValue := NULL], 
                             y = human_metabolic_pathways[, .(functionalCategory = `Pathway id`, `Pathway name`, Class, compound_no)], 
                             by = "functionalCategory")
  setcolorder(ora_reporting, c("functionalCategory", "Pathway name", "Class", "pAdj", "compound_no", "nOverlap", "overlap"))
  
  orphan_list_ORA_result[[i]] = ora_reporting
  print(i)
}

for (i in 1:length(orphan_list_ORA_result)){
  pathway = inner_join(x = data.table("Pathway id" = orphan_list_ORA_result[[i]][,functionalCategory]),
                       y = human_metabolic_pathways[, .(`Pathway id`, linked_compound)],
                       by = "Pathway id")
  if (nrow(pathway) == 0){
    orphan_list_ORA[i, Enriched_pathway := "No Pathway Enriched"]
  }else{
    ep = paste(unlist(pathway[, `Pathway id`]), collapse = ";")
    orphan_list_ORA[i, Enriched_pathway := ep]
  }
}
orphan_list_ORA = as.data.table(separate_longer_delim(orphan_list_ORA, Enriched_pathway, delim = ";"))
orphan_list_ORA = data.table(orphan_list_ORA,
                             Linked_Compound = vector("character", nrow(orphan_list_ORA)))
for (i in 1:nrow(orphan_list_ORA)){
  if(orphan_list_ORA[i, Enriched_pathway] != 	"No Pathway Enriched"){
    ep = orphan_list_ORA[i, Enriched_pathway]
    lp = human_metabolic_pathways[`Pathway id` == ep, linked_compound]
    orphan_list_ORA[i, Linked_Compound := lp]
  }
}

orphan_prediction = lapply(1:nrow(orphan_list), function(o){
  orphan_gene = orphan_list[o, HGNC_symbol]
  etcs_atomic = etcs[SYMBOL == orphan_gene][order(Score, decreasing = TRUE)][Score != 0]
  etcs_atomic = etcs_atomic[1:(nrow(etcs_atomic)*0.04)]
  etcs_atomic_feature = lapply(
    etcs_atomic[,na.omit(Formula)], function(x){as.character(ConvertNames(x, "Formula")[,"KEGG ID"])}
  ) |>
    unlist() |>
    str_split(";") |>
    unlist() |>
    str_trim()
  
  ep = paste(orphan_list_ORA[SYMBOL == orphan_gene, Enriched_pathway], collapse = ";")
  if (ep != "No Pathway Enriched"){
    orphan_list_ORA_atomic = orphan_list_ORA[SYMBOL == orphan_gene]
    ep_status = vector("list", nrow(etcs_atomic))
    for (j in 1:nrow(etcs_atomic)){
      kegg_to_test = str_trim(unlist(str_split(ConvertNames(etcs_atomic[j, Formula], "Formula")[,"KEGG ID"], ";")))
      if (check_NA(kegg_to_test)){
        ep_status[[j]] = ""
      }else{
        path_loc = c()
        for (k in 1:length(kegg_to_test)){
          path_loc = c(grep(kegg_to_test[k], orphan_list_ORA_atomic[, Linked_Compound]), path_loc)
        }
        if (length(path_loc) == 0){
          ep_status[[j]] = ""
        }else{
          ep_status[[j]] = paste("Enriched in pathway", paste(unique(orphan_list_ORA_atomic[path_loc, Enriched_pathway]), collapse = ", "))
        }
      }
    }
    this_orphan_prediction = data.table(
      SYMBOL = orphan_gene,
      Prediction = etcs_atomic$Formula,
      A = etcs_atomic$A,
      B = etcs_atomic$B,
      C = etcs_atomic$C,
      Score = etcs_atomic$Score,
      Pathway_enriched = unlist(ep_status)
    )
  }else{
    this_orphan_prediction = data.table(
      SYMBOL = orphan_gene,
      Prediction = etcs_atomic$Formula,
      A = etcs_atomic$A,
      B = etcs_atomic$B,
      C = etcs_atomic$C,
      Score = etcs_atomic$Score,
      Pathway_enriched = ""
    )
  }
  return(this_orphan_prediction)
}) |> rbindlist()


# lag() shifts the row down by 1 row, and then the change is detected by when the lagged value does not equate the original value
change_indices = which(orphan_prediction$SYMBOL != dplyr::lag(orphan_prediction$SYMBOL, default = orphan_prediction$SYMBOL[1]))

# define blank row
blank_row = data.table(matrix("", ncol = ncol(orphan_prediction), nrow = 2))
colnames(blank_row) = colnames(orphan_prediction)

# add blank row
add_blank_rows = function(data, change_indices) {
  change_indices = change_indices - 1
  new_data = data.table()
  for (i in 1:nrow(data)) {
    new_data = rbind(new_data, data[i, ])
    if (i %in% change_indices) {  # Exclude first index as it is not a change
      new_data = rbind(new_data, blank_row)
    }
  }
  return(new_data)
}

orphan_prediction_edited = add_blank_rows(orphan_prediction, change_indices)
fwrite(orphan_prediction_edited, "Results/Orphan Prediction Result.csv")

# 14 Drug sensitivities ####
## 14.1 Calculate dufference in drug dose-responses ####
### §1 read files #####
CCLE_metadata = read.delim("Raw Data/Cell_lines_annotations_20181226.txt")
Cell_line_metadata = fread("Raw Data/secondary-screen-cell-line-info.csv")
Drug_treatment_metadata = fread("Raw Data/secondary-screen-replicate-collapsed-treatment-info.csv")
Drug_screen = fread("Raw Data/secondary-screen-replicate-collapsed-logfold-change.csv")

### §2 Primary treatment of file #####

# remove cell lines failed STR profiling #
Drug_screen = Drug_screen[-grep("FAILED_STR", V1),]

# translate the CCLE name to ACH code in CCLE2019_SLC_ntc #
CCLE_ach = unlist(lapply(colnames(CCLE2019_SLC_ntc)[2:ncol(CCLE2019_SLC_ntc)], function(x){CCLE_metadata$depMapID[which(CCLE_metadata$CCLE_ID == x)]}))
colnames(CCLE2019_SLC_ntc)[2:ncol(CCLE2019_SLC_ntc)] = CCLE_ach

# match cell lines in CCLE2019_SLC_ntc and drug_screen #
shared_cell_line = intersect(Drug_screen$V1, colnames(CCLE2019_SLC_ntc)[2:ncol(CCLE2019_SLC_ntc)])
CCLE2019_SLC_ntc = CCLE2019_SLC_ntc[, c(TRUE, colnames(CCLE2019_SLC_ntc)[2:ncol(CCLE2019_SLC_ntc)] %in% shared_cell_line), with = FALSE]
Drug_screen = Drug_screen[V1 %in% shared_cell_line, ]
setcolorder(CCLE2019_SLC_ntc, c("SYMBOL", Drug_screen$V1))

### §3 Fit curve to SLC-drug pairs with loess ####
Curve = function(S, D){
  # get the metadata for desired drug, produce the screen profile for that drug #
  D_metadata = Drug_treatment_metadata[name == D]
  D_screen = Drug_screen[,c(TRUE, colnames(Drug_screen)[2:ncol(Drug_screen)] %in% D_metadata$column_name), with = FALSE]
  
  # get the expression profile of desired SLC, acquire cell lines of high and low expression #
  SLC_expression = as.numeric(unlist(CCLE2019_SLC_ntc[SYMBOL == S, 3:ncol(CCLE2019_SLC_ntc)]))
  # if the SLC is not found in the transcriptomics, terminate the function and return NA #
  if (length(SLC_expression) == 0){
    output = list(NaN, NaN, NaN)
    names(output) = c("Best.log.Dose", "Minimum.Wilcox.p", "Absolute.Median.Difference")
    return(output)
  }
  q = quantile(SLC_expression, c(0.2, 0.8))
  # if the lower quantile is not 0 #
  if (q[1] != 0){
    low_expression_cell = colnames(CCLE2019_SLC_ntc)[which(SLC_expression < q[1]) + 2]
    # if the lower quantile is 0 #
  }else{
    low_expression_cell = colnames(CCLE2019_SLC_ntc)[which(SLC_expression == q[1]) + 2]
  }
  high_expression_cell = colnames(CCLE2019_SLC_ntc)[which(SLC_expression > q[2]) + 2]
  
  # pool dose-response curve data for both high and low expression #
  cell_viability_for_D = melt(D_screen, id.vars = "V1")
  colnames(cell_viability_for_D) = c("cell", "dose", "viability")
  cell_viability_for_D$dose = as.numeric(str_match(cell_viability_for_D$dose, "::\\s*(.*?)\\s*::")[,2])
  cell_viability_for_D$log.dose = round(log10(as.numeric(cell_viability_for_D$dose)), 2)
  cell_viability_for_D$expression_type = vector("character", nrow(cell_viability_for_D))
  
  for (i in 1:nrow(cell_viability_for_D)){
    if (cell_viability_for_D$cell[i] %in% low_expression_cell){
      cell_viability_for_D[i, expression_type := "low expression"] # low expression
    }else if (cell_viability_for_D$cell[i] %in% high_expression_cell){
      cell_viability_for_D[i, expression_type := "high expression"] # high expression
    }else{
      cell_viability_for_D[i, expression_type := "NE"] # not examined #
    }
  }
  
  cell_viability_for_D = cell_viability_for_D[expression_type != "NE"][order(expression_type)]
  fit_high = loess(viability ~ log.dose, data = cell_viability_for_D[expression_type == "high expression"])
  fit_low = loess(viability ~ log.dose, data = cell_viability_for_D[expression_type == "low expression"])
  new_dose = expand.grid(log.dose = seq(-3.21, 1, 0.01))
  preds_high = predict(fit_high, newdata = new_dose)
  preds_low = predict(fit_low, newdata = new_dose)
  t_statistic = t.test(preds_high, preds_low)
  output = list(p_value = t_statistic$p.value,
                md = as.numeric(t_statistic$estimate[1] - t_statistic$estimate[2]))
  return(output)
}
plotCurve = function(S, D, ylim = NULL){
  # get the metadata for desired drug, produce the screen profile for that drug #
  D_metadata = Drug_treatment_metadata[name == D]
  D_screen = Drug_screen[,c(TRUE, colnames(Drug_screen)[2:ncol(Drug_screen)] %in% D_metadata$column_name), with = FALSE]
  
  # get the expression profile of desired SLC, acquire cell lines of high and low expression #
  SLC_expression = as.numeric(unlist(CCLE2019_SLC_ntc[SYMBOL == S, 3:ncol(CCLE2019_SLC_ntc)]))
  # if the SLC is not found in the transcriptomics, terminate the function and return NA #
  if (length(SLC_expression) == 0){
    output = list(NaN, NaN, NaN)
    names(output) = c("Best.log.Dose", "Minimum.Wilcox.p", "Absolute.Median.Difference")
    return(output)
  }
  q = quantile(SLC_expression, c(0.2, 0.8))
  # if the lower quantile is not 0 #
  if (q[1] != 0){
    low_expression_cell = colnames(CCLE2019_SLC_ntc)[which(SLC_expression < q[1]) + 2]
    # if the lower quantile is 0 #
  }else{
    low_expression_cell = colnames(CCLE2019_SLC_ntc)[which(SLC_expression == q[1]) + 2]
  }
  high_expression_cell = colnames(CCLE2019_SLC_ntc)[which(SLC_expression > q[2]) + 2]
  
  # pool dose-response curve data for both high and low expression #
  cell_viability_for_D = melt(D_screen, id.vars = "V1")
  colnames(cell_viability_for_D) = c("cell", "dose", "viability")
  cell_viability_for_D$dose = as.numeric(str_match(cell_viability_for_D$dose, "::\\s*(.*?)\\s*::")[,2])
  cell_viability_for_D$log.dose = as.numeric(round(log10(as.numeric(cell_viability_for_D$dose)), 2))
  cell_viability_for_D$Expression_type = vector("character", nrow(cell_viability_for_D))
  
  for (i in 1:nrow(cell_viability_for_D)){
    if (cell_viability_for_D$cell[i] %in% low_expression_cell){
      cell_viability_for_D[i, Expression_type := "low expression"] # low expression
    }else if (cell_viability_for_D$cell[i] %in% high_expression_cell){
      cell_viability_for_D[i, Expression_type := "high expression"] # high expression
    }else{
      cell_viability_for_D[i, Expression_type := "NE"] # not examined #
    }
  }
  
  cell_viability_for_D = cell_viability_for_D[Expression_type != "NE"][order(Expression_type)]
  if (is.null(ylim)){
    p = ggplot(cell_viability_for_D, aes(x=log.dose, y=viability, fill = Expression_type, color = Expression_type)) + 
      geom_point(position = position_dodge(0.2)) +
      geom_smooth(method = 'loess') + 
      theme_minimal()
  }else{
    p = ggplot(cell_viability_for_D, aes(x=log.dose, y=viability, fill = Expression_type, color = Expression_type)) + 
      geom_point(position = position_dodge(0.2)) +
      geom_smooth(method = 'loess') + 
      scale_y_continuous(limits = ylim) +
      theme_minimal()
  }
  return(p)
}

# --------------------------- #
#        Figure 5B data       #
# --------------------------- #

plotCurve("SLC35F2", "YM155", ylim = c(-7.5, 2.5))
plotCurve("PTCHD4", "idasanutlin", ylim = c(-7.5, 2.5))
plotCurve("SLC19A2", "thiamine", ylim = c(-7.5, 2.5))

### §4 All interaction calculation ####
SLC_list = CCLE2019_SLC_ntc$SYMBOL
Compound_list = unique(Drug_treatment_metadata$name)

# SLC-Compound Interaction (SCI)
SCI = data.table(SLC = rep(SLC_list, each = length(Compound_list)),
                 Compound = rep(Compound_list, length(SLC_list)),
                 pvalue = vector("numeric", 747168),
                 mean.diff = vector("numeric", 747168))
for (i in 1:nrow(SCI)){
  ans = Curve(SCI[i, SLC], SCI[i, Compound])
  SCI[i, pvalue := ans$p_value]
  SCI[i, mean.diff := ans$md]
  print(i)
}

fwrite(SCI, "Results/SCI.csv")

## 14.2 Benchmark: Is SLC known to transport the cytotoxic drug affect viability more? ####
### §1 Known drug target ####
kdt = fread("Raw Data/Known drug target.csv")
kdt = kdt[Name %in% Drug_treatment_metadata$name]
# eliminate imatinib, which is a drug targeting specific mutation in cancer (BCR-ABL fusion gene & KIT)
kdt = kdt[Compound != "Imatinib"]
# build kdt Curve data table #
kdt_curve = kdt
kdt_curve$pvalue = vector("numeric", nrow(kdt_curve)) # p-value of paired t-test
kdt_curve$mean.diff = vector("numeric", nrow(kdt_curve)) # absolute mean difference by paired t-test

for (i in 1:nrow(kdt_curve)){
  gene = kdt_curve[i, SYMBOL]
  drug = kdt_curve[i, Name]
  kdt_curve[i, pvalue := SCI[SLC == gene & Compound == drug, pvalue]]
  kdt_curve[i, mean.diff := SCI[SLC == gene & Compound == drug, mean.diff]]
  print(i)
}

# remove those that are not apparent in the known drug target
kdt_curve = kdt_curve[pvalue < 0.05]

### §2 Simulated drug target ####
fake_kdt_curve_collection = vector("list", 100)
kdt_curve_gene_pool = unique(kdt_curve$SYMBOL)
kdt_curve_drug_pool = unique(kdt_curve$Name)

for (sample_time in 1:100){
  # initialise fake kdt_curve data table #
  fake_kdt_curve = setNames(data.table(matrix(ncol = ncol(kdt_curve), nrow = nrow(kdt_curve))), colnames(kdt_curve))
  char.col = colnames(fake_kdt_curve)[1:5]
  num.col = colnames(fake_kdt_curve)[6:7]
  fake_kdt_curve[, (char.col) := lapply(.SD, as.character), .SDcols = char.col]
  fake_kdt_curve[, (num.col) := lapply(.SD, as.character), .SDcols = num.col]
  fake_kdt_curve[, Action := rep("Fake", nrow(fake_kdt_curve))]
  fake_kdt_curve[, c("SYMBOL", "Localisation", "Compartment") := as.list(kdt_curve[, c("SYMBOL", "Localisation", "Compartment")])]
  # for each SLC, shuffle fake substrate formula #
  for (j in 1:nrow(fake_kdt_curve)){
    true_drug = kdt_curve[SYMBOL == fake_kdt_curve$SYMBOL[j]][, Name]
    fake_drug_pool = kdt_curve_drug_pool[!kdt_curve_drug_pool %in% true_drug]
    fake_drug = sample(fake_drug_pool, 1)
    fake_kdt_curve[j, Name := fake_drug]
    fake_kdt_curve[j, c("pvalue", "mean.diff") := 
                     SCI[SLC == fake_kdt_curve$SYMBOL[j] & Compound == fake_drug, c("pvalue", "mean.diff")]]
  }
  fake_kdt_curve_collection[[sample_time]] = fake_kdt_curve
  print(paste(sample_time, "done"))
}

### §3 Distributions ####

# --------------------------- #
#    Figure 5C (left) data    #
# --------------------------- #

analyseDist = function(data.simu, data.true, measurement){
  if (!measurement %in% c("logp", "abs.mean.diff")){
    stop("Measurement not valid. Choose between: logp; abs.mean.diff")
  }else if (measurement == "logp"){
    boot = sapply(1:100, function(s){
      resample = sample(data.true[, as.numeric(get(measurement))], size = nrow(data.true), replace = TRUE)
      return(median(resample))
    })
    Dist = data.table(
      type = rep(c("simulated", "known"), each = 100),
      value = c(data.simu[, median(as.numeric(get(measurement))), by = Iteration][, V1], boot)
    )
  }else if (measurement == "abs.mean.diff"){
    boot = sapply(1:100, function(s){
      resample = sample(data.true[, as.numeric(get(measurement))], size = nrow(data.true), replace = TRUE)
      return(mean(resample))
    })
    Dist = data.table(
      type = rep(c("simulated", "known"), each = 100),
      value = c(data.simu[, mean(as.numeric(get(measurement))), by = Iteration][, V1], boot)
    )
  }
  
  rain_height = .1
  colorset = c("known" = "#990000", "simulated" = "grey")
  if(measurement == "logp"){
    breaks = seq(0, 20, 2)
    name = "Median p-value"
  }else if (measurement == "abs.mean.diff"){
    breaks = seq(0, 0.4, 0.05)
    name = "Average Absolute Mean Difference"
  }
  
  g = ggplot(Dist, aes(x = "", y = value, fill = type)) +
    # mountain
    introdataviz::geom_flat_violin(trim=FALSE, alpha = 0.4, position = position_nudge(x = -0.6)) +
    # rain
    geom_point(aes(colour = type), size = 1, alpha = 0.4, show.legend = FALSE, position = position_jitter(width = rain_height, height = 0)) +
    geom_boxplot(width = rain_height, alpha = 0.4, show.legend = FALSE, outlier.shape = NA) +
    scale_y_continuous(name = name,
                       breaks = breaks, 
                       limits = c(min(breaks), max(breaks))) +
    # Manual color settings
    scale_fill_manual(values = colorset) + # For violin plots
    scale_color_manual(values = colorset) + # For points
    coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
  
  p = t.test(
    x = Dist[type == "simulated", value],
    y = Dist[type == "known", value],
    alternative = "less"
  )$p.value
  
  return(list(dist.data = Dist, dist.plot = g, dist.pvalue = p))
}

fake_kdt_curve_collection = lapply(1:100, function(f){
  data = fake_kdt_curve_collection[[f]] |>
    mutate(
      Iteration = f,
      .before = everything()
    )
  return(data)
}) |> rbindlist()

fake_kdt_curve_collection = fake_kdt_curve_collection |>
  mutate(
    logp = -log10(as.numeric(pvalue)),
    abs.mean.diff = abs(as.numeric(mean.diff))
  )
kdt_curve = kdt_curve |>
  mutate(
    logp = -log10(as.numeric(pvalue)),
    abs.mean.diff = abs(as.numeric(mean.diff))
  )

model.dist = analyseDist(data.simu = fake_kdt_curve_collection, data.true = kdt_curve, measurement = "abs.mean.diff")
model.dist$dist.plot
model.dist$dist.data
model.dist$dist.pvalue

# --------------------------- #
#   Figure S3C (left) data    #
# --------------------------- #

roc = function(data.simu, data.true){
  boundary = c(
    true.upper = data.true[, max(abs(as.numeric(abs.mean.diff)))],
    fake.upper = data.simu[, max(abs(as.numeric(abs.mean.diff)))],
    true.lower = data.true[, min(abs(as.numeric(abs.mean.diff)))],
    fake.lower = data.simu[, min(abs(as.numeric(abs.mean.diff)))]
  )
  boundary.upper = ceiling(max(boundary))
  boundary.lower = 0
  threshold_ticking = seq(boundary.upper, boundary.lower, -0.01)
  ROCcurve = lapply(threshold_ticking, function(t){
    true_pos  = nrow(data.true[abs.mean.diff > t])/nrow(data.true)
    false_pos = sapply(1:100, function(i){
      this_iteration = data.simu[Iteration == i]
      return(nrow(this_iteration[abs.mean.diff > t])/nrow(this_iteration))
    })
    progress = round(which(threshold_ticking == t)/length(threshold_ticking), digits = 3) * 100
    if (progress %% 10 == 0) {
      message("Generating ROC ", progress, " %")
    }
    return(
      data.table(
        threshold = t,
        true.positive.rate = true_pos,
        false.positive.rate = mean(false_pos)
      )
    )
  }) |> rbindlist()
  # Generates plot
  ROCcurve.plot = ggplot(ROCcurve) +
    geom_point(aes(x = false.positive.rate, y = true.positive.rate)) +
    geom_line(aes(x = false.positive.rate, y = true.positive.rate)) +
    geom_abline(col = "grey", linetype = "dashed") +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  # Calculate AUC with the Trapezoidal rule
  sensitivity = ROCcurve$true.positive.rate
  specificity = ROCcurve$false.positive.rate
  height = (sensitivity[-1]+sensitivity[-length(sensitivity)])/2
  width = diff(specificity)
  AUC = sum(height*width)
  
  return(
    list(
      ROC.data = ROCcurve,
      ROC.plot = ROCcurve.plot,
      AUC = AUC
    )
  )
}

roc.object.d = roc(data.simu = fake_kdt_curve_collection, data.true = kdt_curve)

roc.object.d$ROC.plot
roc.object.d$ROC.data
roc.object.d$AUC

# --------------------------- #
#  Figure S3C (right) data    #
# --------------------------- #

prc = function(data.simu, data.true){
  boundary = c(
    true.upper = data.true[, max(as.numeric(abs.mean.diff))],
    fake.upper = data.simu[, max(as.numeric(abs.mean.diff))]
  )
  boundary.upper = ceiling(max(boundary))
  boundary.lower = 0
  threshold_ticking = seq(boundary.lower, boundary.upper, 0.01)
  PRCcurve = lapply(threshold_ticking, function(t){
    true_pos  = nrow(data.true[abs.mean.diff >= t])
    false_pos = sapply(1:100, function(i){
      this_iteration = data.simu[Iteration == i]
      return(nrow(this_iteration[abs.mean.diff >= t]))
    })
    false_neg = nrow(data.true[abs.mean.diff < t])
    if (true_pos == 0 & mean(false_pos) == 0){
      precision = 1
    }else{
      precision = true_pos/(true_pos + mean(false_pos))
    }
    progress = round(which(threshold_ticking == t)/length(threshold_ticking), digits = 3) * 100
    if (progress %% 10 == 0) {
      message("Generating PRC ", progress, " %")
    }
    return(
      data.table(
        threshold = t,
        true.positive = true_pos,
        false.positive = mean(false_pos),
        false.negative = false_neg,
        Precision = precision,
        Recall = true_pos/(true_pos + false_neg)
      )
    )
  }) |> rbindlist()
  # Generates plot
  PRCcurve.plot = ggplot(PRCcurve) +
    geom_point(aes(x = Recall, y = Precision)) +
    geom_line(aes(x = Recall, y = Precision)) +
    geom_line(aes(x = Recall, y = Precision)) +
    geom_hline(col = "grey", linetype = "dashed", yintercept = 0.5) +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  # Calculate AUC with the Trapezoidal rule
  sensitivity = PRCcurve$Precision
  specificity = PRCcurve$Recall
  height = (sensitivity[-1]+sensitivity[-length(sensitivity)])/2
  width = -diff(specificity)
  AUC = sum(height*width)
  
  return(
    list(
      PRC.data = PRCcurve,
      PRC.plot = PRCcurve.plot,
      AUC = AUC
    )
  )
}

prc.object.d = prc(data.simu = fake_kdt_curve_collection, data.true = kdt_curve)

prc.object.d$PRC.plot
prc.object.d$PRC.data
prc.object.d$AUC


## 14.3 Benchmark: Is SLC known to transport the cytotoxic drug shows better correlation between expression and dose? ####
### §1 Known drug target ####
CheckDoseEffect = function(S, D, original.data = FALSE){
  # get the metadata for desired drug, produce the screen profile for that drug #
  D_metadata = Drug_treatment_metadata[name == D]
  D_screen = Drug_screen[,c(TRUE, colnames(Drug_screen)[2:ncol(Drug_screen)] %in% D_metadata$column_name), with = FALSE]
  
  # get the expression profile of desired SLC, acquire cell lines of high and low expression #
  SLC_expression = as.numeric(unlist(SLC_transcriptomics[SYMBOL == S, 3:ncol(SLC_transcriptomics)]))
  # if the SLC is not found in the transcriptomics, terminate the function and return NA #
  if (length(SLC_expression) == 0){
    output = list(NaN, NaN)
    names(output) = c("p_value", "md")
    return(output)
  }
  q = quantile(SLC_expression, c(0.2, 0.8))
  # if the lower quantile is not 0 #
  if (q[1] != 0){
    low_expression_cell = colnames(SLC_transcriptomics)[which(SLC_expression < q[1]) + 2]
    # if the lower quantile is 0 #
  }else{
    low_expression_cell = colnames(SLC_transcriptomics)[which(SLC_expression == q[1]) + 2]
  }
  high_expression_cell = colnames(SLC_transcriptomics)[which(SLC_expression > q[2]) + 2]
  
  # pool dose-response curve data for both high and low expression #
  cell_viability_for_D = melt(D_screen, id.vars = "V1")
  colnames(cell_viability_for_D) = c("cell", "dose", "viability")
  cell_viability_for_D$dose = as.numeric(str_match(cell_viability_for_D$dose, "::\\s*(.*?)\\s*::")[,2])
  cell_viability_for_D$log.dose = round(log10(as.numeric(cell_viability_for_D$dose)), 2)
  cell_viability_for_D$expression_type = vector("character", nrow(cell_viability_for_D))
  
  for (i in 1:nrow(cell_viability_for_D)){
    if (cell_viability_for_D$cell[i] %in% low_expression_cell){
      cell_viability_for_D[i, expression_type := "low expression"] # low expression
    }else if (cell_viability_for_D$cell[i] %in% high_expression_cell){
      cell_viability_for_D[i, expression_type := "high expression"] # high expression
    }else{
      cell_viability_for_D[i, expression_type := "NE"] # not examined #
    }
  }
  
  cell_viability_for_D = cell_viability_for_D[expression_type != "NE"][order(expression_type)]
  fit_high = loess(viability ~ log.dose, data = cell_viability_for_D[expression_type == "high expression"])
  fit_low = loess(viability ~ log.dose, data = cell_viability_for_D[expression_type == "low expression"])
  if (isTRUE(original.data)){
    original.dose = c(-3.214, -2.612, -2.010, -1.408, -0.806, -0.204, 0.398, 1.000)
    new_dose = expand.grid(log.dose = original.dose)
    preds_high = predict(fit_high, newdata = new_dose)
    preds_low = predict(fit_low, newdata = new_dose)
    CorIn = cor.test(original.dose, abs(as.numeric(preds_high - preds_low)), method = "spearman")
  }else{
    new_dose = expand.grid(log.dose = seq(-3.21, 1, 0.01))
    preds_high = predict(fit_high, newdata = new_dose)
    preds_low = predict(fit_low, newdata = new_dose)
    CorIn = cor.test(seq(-3.21, 1, 0.01), abs(as.numeric(preds_high - preds_low)), method = "spearman")
  }
  output = list(p_value = as.numeric(CorIn$p.value),
                rho = as.numeric(CorIn$estimate))
  return(output)
}
kdt_CDE = kdt_curve[,1:5]
kdt_CDE$pvalue = vector("numeric", nrow(kdt_CDE)) # p-value of spearman's correlation
kdt_CDE$rho = vector("numeric", nrow(kdt_CDE)) # rho of spearman's correlation
for (i in 1:nrow(kdt_CDE)){
  gene = kdt_CDE[i, SYMBOL]
  drug = kdt_CDE[i, Name]
  statistics = CheckDoseEffect(gene, drug)
  kdt_CDE[i, pvalue := statistics$p_value]
  kdt_CDE[i, rho := statistics$rho]
  print(i)
}
### §2 Simulated drug target ####
fake_kdt_CDE_collection = vector("list", 100)
kdt_CDE_gene_pool = unique(kdt_CDE$SYMBOL)
kdt_CDE_drug_pool = unique(kdt_CDE$Name)

for (sample_time in 1:100){
  # initialise fake kdt_CDE data table #
  fake_kdt_CDE = setNames(data.table(matrix(ncol = ncol(kdt_CDE), nrow = nrow(kdt_CDE))), colnames(kdt_CDE))
  char.col = colnames(fake_kdt_CDE)[1:5]
  num.col = colnames(fake_kdt_CDE)[6:7]
  fake_kdt_CDE[, (char.col) := lapply(.SD, as.character), .SDcols = char.col]
  fake_kdt_CDE[, (num.col) := lapply(.SD, as.character), .SDcols = num.col]
  fake_kdt_CDE[, c("SYMBOL", "Localisation", "Compartment") := as.list(kdt_CDE[, c("SYMBOL", "Localisation", "Compartment")])]
  # for each SLC, shuffle fake substrate formula #
  for (j in 1:nrow(fake_kdt_CDE)){
    true_drug = kdt_CDE[SYMBOL == fake_kdt_CDE$SYMBOL[j]][, Name]
    fake_drug_pool = kdt_CDE_drug_pool[!kdt_CDE_drug_pool %in% true_drug]
    fake_drug = sample(fake_drug_pool, 1)
    fake_kdt_CDE[j, Name := fake_drug]
    #statistics = CheckDoseEffect(fake_kdt_CDE[j, SYMBOL], fake_drug, original.data = TRUE)
    statistics = CheckDoseEffect(fake_kdt_CDE[j, SYMBOL], fake_drug)
    fake_kdt_CDE[j, c("pvalue", "rho") := list(statistics$p_value, statistics$rho)]
    print(paste(sample_time, ":", round(j/nrow(fake_kdt_CDE),2)))
  }
  fake_kdt_CDE_collection[[sample_time]] = fake_kdt_CDE
  print(paste(sample_time, "done"))
}

### §3 Distributions ####

# --------------------------- #
#   Figure 5C (right) data    #
# --------------------------- #

fake_kdt_CDE_collection = lapply(1:100, function(f){
  data = fake_kdt_CDE_collection[[f]] |>
    mutate(
      Iteration = f,
      .before = everything()
    )
  return(data)
}) |> rbindlist()

analyseDist = function(data.simu, data.true){
  boot = sapply(1:100, function(s){
    resample = sample(data.true[, abs(as.numeric(rho))], size = nrow(data.true), replace = TRUE)
    return(mean(resample))
  })
  Dist = data.table(
    type = rep(c("simulated", "known"), each = 100),
    value = c(data.simu[, mean(abs(as.numeric(rho))), by = Iteration][, V1], boot)
  )
  
  rain_height = .1
  colorset = c("known" = "#990000", "simulated" = "grey")
  breaks = seq(0, 1, 0.2)
  name = "Absolute Spearman's rho"
  
  g = ggplot(Dist, aes(x = "", y = value, fill = type)) +
    # mountain
    introdataviz::geom_flat_violin(trim=FALSE, alpha = 0.4, position = position_nudge(x = -0.6)) +
    # rain
    geom_point(aes(colour = type), size = 1, alpha = 0.4, show.legend = FALSE, position = position_jitter(width = rain_height, height = 0)) +
    geom_boxplot(width = rain_height, alpha = 0.4, show.legend = FALSE, outlier.shape = NA) +
    scale_y_continuous(name = name,
                       breaks = breaks, 
                       limits = c(min(breaks), max(breaks))) +
    # Manual color settings
    scale_fill_manual(values = colorset) + # For violin plots
    scale_color_manual(values = colorset) + # For points
    coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
  
  p = t.test(
    x = Dist[type == "simulated", value],
    y = Dist[type == "known", value],
    alternative = "less"
  )$p.value
  
  return(list(dist.data = Dist, dist.plot = g, dist.pvalue = p))
}

CDE.dist = analyseDist(data.simu = fake_kdt_CDE_collection, data.true = kdt_CDE)
CDE.dist$dist.plot
CDE.dist$dist.data
CDE.dist$dist.pvalue

# --------------------------- #
#   Figure S3D (left) data    #
# --------------------------- #

roc = function(data.simu, data.true){
  boundary = c(
    true.upper = data.true[, max(as.numeric(abs.rho))],
    fake.upper = data.simu[, max(as.numeric(abs.rho))],
    true.lower = data.true[, min(as.numeric(abs.rho))],
    fake.lower = data.simu[, min(as.numeric(abs.rho))]
  )
  boundary.upper = ceiling(max(boundary))
  boundary.lower = 0
  threshold_ticking = seq(boundary.upper, boundary.lower, -0.01)
  
  ROCcurve = lapply(threshold_ticking, function(t){
    true_pos  = nrow(data.true[abs.rho >= t])/nrow(data.true)
    false_pos = sapply(1:100, function(i){
      this_iteration = data.simu[Iteration == i]
      return(nrow(this_iteration[abs.rho >= t])/nrow(this_iteration))
    })
    progress = round(which(threshold_ticking == t)/length(threshold_ticking), digits = 3) * 100
    if (progress %% 10 == 0) {
      message("Generating ROC ", progress, " %")
    }
    return(
      data.table(
        threshold = t,
        true.positive.rate = true_pos,
        false.positive.rate = mean(false_pos)
      )
    )
  }) |> rbindlist()
  # Generates plot
  ROCcurve.plot = ggplot(ROCcurve) +
    geom_point(aes(x = false.positive.rate, y = true.positive.rate)) +
    geom_line(aes(x = false.positive.rate, y = true.positive.rate)) +
    geom_abline(col = "grey", linetype = "dashed") +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  # Calculate AUC with the Trapezoidal rule
  sensitivity = ROCcurve$true.positive.rate
  specificity = ROCcurve$false.positive.rate
  height = (sensitivity[-1]+sensitivity[-length(sensitivity)])/2
  width = diff(specificity)
  AUC = sum(height*width)
  
  return(
    list(
      ROC.data = ROCcurve,
      ROC.plot = ROCcurve.plot,
      AUC = AUC
    )
  )
}
roc.object.cdt = roc(data.simu = fake_kdt_CDE_collection, data.true = kdt_CDE)

roc.object.cdt$ROC.plot
roc.object.cdt$ROC.data
roc.object.cdt$AUC

# --------------------------- #
#  Figure S3D (right) data    #
# --------------------------- #

prc = function(data.simu, data.true){
  boundary = c(
    true.upper = data.true[, max(as.numeric(abs.rho))],
    fake.upper = data.simu[, max(as.numeric(abs.rho))]
  )
  boundary.upper = ceiling(max(boundary))
  boundary.lower = 0
  threshold_ticking = seq(boundary.lower, boundary.upper, 0.01)
  PRCcurve = lapply(threshold_ticking, function(t){
    true_pos  = nrow(data.true[abs.rho >= t])
    false_pos = sapply(1:100, function(i){
      this_iteration = data.simu[Iteration == i]
      return(nrow(this_iteration[abs.rho >= t]))
    })
    false_neg = nrow(data.true[abs.rho < t])
    if (true_pos == 0 & mean(false_pos) == 0){
      precision = 1
    }else{
      precision = true_pos/(true_pos + mean(false_pos))
    }
    progress = round(which(threshold_ticking == t)/length(threshold_ticking), digits = 3) * 100
    if (progress %% 10 == 0) {
      message("Generating PRC ", progress, " %")
    }
    return(
      data.table(
        threshold = t,
        true.positive = true_pos,
        false.positive = mean(false_pos),
        false.negative = false_neg,
        Precision = precision,
        Recall = true_pos/(true_pos + false_neg)
      )
    )
  }) |> rbindlist()
  # Generates plot
  PRCcurve.plot = ggplot(PRCcurve) +
    geom_point(aes(x = Recall, y = Precision)) +
    geom_line(aes(x = Recall, y = Precision)) +
    geom_line(aes(x = Recall, y = Precision)) +
    geom_hline(col = "grey", linetype = "dashed", yintercept = 0.5) +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  # Calculate AUC with the Trapezoidal rule
  sensitivity = PRCcurve$Precision
  specificity = PRCcurve$Recall
  height = (sensitivity[-1]+sensitivity[-length(sensitivity)])/2
  width = -diff(specificity)
  AUC = sum(height*width)
  
  return(
    list(
      PRC.data = PRCcurve,
      PRC.plot = PRCcurve.plot,
      AUC = AUC
    )
  )
}
prc.object.cdt = prc(data.simu = fake_kdt_CDE_collection, data.true = kdt_CDE)

prc.object.cdt$PRC.plot
prc.object.cdt$PRC.data
prc.object.cdt$AUC

### §4 calculate for every interaction ####
SCI$dose.pvalue = vector("numeric", nrow(SCI))
SCI$dose.rho = vector("numeric", nrow(SCI))

for (i in 1:nrow(SCI)){
  gene = SCI[i, SLC]
  drug = SCI[i, Compound]
  statistics = CheckDoseEffect(gene, drug)
  SCI[i, dose.pvalue := statistics$p_value]
  SCI[i, dose.rho := statistics$rho]
  print(i)
}

fwrite(SCI, "Results/SCI.csv")

## 14.3 Filtering thresholds for SLC-drug predictions ####
# Note: Filtering thresholds are provided as a supplementary table (S12).
#       The code below demonstrates how filtering thresholds are produced.
drug_list = unique(Drug_treatment_metadata[, name])
pseudoSLC_benchmark = vector("list", length(drug_list))
names(pseudoSLC_benchmark) = drug_list
for (drug_index in 1:length(drug_list)){
  drug = drug_list[drug_index]
  benchmarking = vector("list", 100)
  for (sample_time in 1:100){
    benchmarking[[sample_time]] = Curve(drug)
    print(paste(drug, ":", sample_time))
  }
  pseudoSLC_benchmark[[drug_index]] = benchmarking
}

pseudo_benchmark = data.table(Drug = rep(drug_list, each = 100),
                              logp = vector("numeric", length(drug_list)*100),
                              amd = vector("numeric", length(drug_list)*100))
for (drug_index in 1:length(drug_list)){
  p = unlist(lapply(1:100, function(x){pseudoSLC_benchmark[[drug_index]][[x]][[1]]}))
  md = unlist(lapply(1:100, function(x){pseudoSLC_benchmark[[drug_index]][[x]][[2]]}))
  pseudo_benchmark[Drug == drug_list[drug_index], logp := -log10(p)]
  pseudo_benchmark[Drug == drug_list[drug_index], amd := abs(md)]
}

boxplot(pseudo_benchmark[, mean(logp), by = Drug][, V1])
boxplot(pseudo_benchmark[, mean(amd), by = Drug][, V1])

pseudo_significance = data.table(Drug = drug_list,
                                 logp_threshold = vector("numeric", length(drug_list)),
                                 amd_threshold = vector("numeric", length(drug_list)))
for (drug_index in 1:length(drug_list)){
  pM = pseudo_benchmark[Drug == drug_list[drug_index], mean(logp)]
  pSD = 2*pseudo_benchmark[Drug == drug_list[drug_index], sd(logp)]
  pseudo_significance[Drug == drug_list[drug_index], logp_threshold := pM + pSD]
  
  amdM = pseudo_benchmark[Drug == drug_list[drug_index], mean(amd)]
  amdSD = 2*pseudo_benchmark[Drug == drug_list[drug_index], sd(amd)]
  pseudo_significance[Drug == drug_list[drug_index], amd_threshold := amdM + amdSD]
}

fwrite(pseudo_significance, "Results/Filtering thresholds.csv")

## 14.4 Generating SLC-drug interaction predictions ####
# Note: SLC-drug predictions are provided as a supplementary table (S13).
#       The code below demonstrates how predictions are produced.

# remove interactions that are likely to fall in pseudo significance #
pseudo_significance = fread("Raw Data/Filtering thresholds.csv")
SCI_filtered = setNames(data.table(matrix(ncol = ncol(SCI), nrow = 0)), colnames(SCI))
for (drug in pseudo_significance[, Drug]){
  logpTh = pseudo_significance[Drug == drug, logp_threshold]
  amdTh = pseudo_significance[Drug == drug, amd_threshold]
  SCI_atomic = SCI[Compound == drug][-log10(pvalue) > logpTh & abs(mean.diff) > amdTh]
  SCI_filtered = rbind(SCI_filtered, SCI_atomic)
}

# taking into account of dose-viability correlation
SCI_filtered = SCI_filtered[dose.pvalue < 0.05][order(-abs(mean.diff), -abs(dose.rho))]

# remove mutation specific drugs
SCI_filtered = SCI_filtered[Compound != "poziotinib"]

# generate prediction
drug_prediction = SCI_filtered[1:50]
fwrite(drug_prediction, "Results/SLC-drug interaction predictions.csv")

## 14.5 Interaction with RITA ####
dir.create("Results/RITA/")

drug_prediction = fread("Raw Data/SLC-drug interaction predictions.csv")
Curve = function(S, D){
  # get the metadata for desired drug, produce the screen profile for that drug #
  D_metadata = Drug_treatment_metadata[name == D]
  D_screen = Drug_screen[,c(TRUE, colnames(Drug_screen)[2:ncol(Drug_screen)] %in% D_metadata$column_name), with = FALSE]
  
  # get the expression profile of desired SLC, acquire cell lines of high and low expression #
  SLC_expression = as.numeric(unlist(CCLE2019_SLC_ntc[SYMBOL == S, 3:ncol(CCLE2019_SLC_ntc)]))
  # if the SLC is not found in the transcriptomics, terminate the function and return NA #
  if (length(SLC_expression) == 0){
    output = list(NaN, NaN, NaN)
    names(output) = c("Best.log.Dose", "Minimum.Wilcox.p", "Absolute.Median.Difference")
    return(output)
  }
  q = quantile(SLC_expression, c(0.2, 0.8))
  # if the lower quantile is not 0 #
  if (q[1] != 0){
    low_expression_cell = colnames(CCLE2019_SLC_ntc)[which(SLC_expression < q[1]) + 2]
    # if the lower quantile is 0 #
  }else{
    low_expression_cell = colnames(CCLE2019_SLC_ntc)[which(SLC_expression == q[1]) + 2]
  }
  high_expression_cell = colnames(CCLE2019_SLC_ntc)[which(SLC_expression > q[2]) + 2]
  
  # pool dose-response curve data for both high and low expression #
  cell_viability_for_D = melt(D_screen, id.vars = "V1")
  colnames(cell_viability_for_D) = c("cell", "dose", "viability")
  cell_viability_for_D$dose = as.numeric(str_match(cell_viability_for_D$dose, "::\\s*(.*?)\\s*::")[,2])
  cell_viability_for_D$log.dose = round(log10(as.numeric(cell_viability_for_D$dose)), 2)
  cell_viability_for_D$expression_type = vector("character", nrow(cell_viability_for_D))
  
  for (i in 1:nrow(cell_viability_for_D)){
    if (cell_viability_for_D$cell[i] %in% low_expression_cell){
      cell_viability_for_D[i, expression_type := "low expression"] # low expression
    }else if (cell_viability_for_D$cell[i] %in% high_expression_cell){
      cell_viability_for_D[i, expression_type := "high expression"] # high expression
    }else{
      cell_viability_for_D[i, expression_type := "NE"] # not examined #
    }
  }
  
  cell_viability_for_D = cell_viability_for_D[expression_type != "NE"][order(expression_type)]
  
  return(cell_viability_for_D)
}
RITA_prediction = durg_prediction[Compound == "RITA"]

for (i in 1:nrow(RITA_prediction)){
  temp = Curve(S = RITA_prediction[i, SLC], D = RITA_prediction[i, Compound])
  fwrite(temp, paste0("Results/RITA/", RITA_prediction[i, SLC], ".csv"))
}

temp = list.files(pattern="Results/RITA/")
sensitising = temp[!temp %in% "SLC3A1.csv"]
attenuating = "SLC3A1.csv"
SLCs = lapply(sensitising, function(x){fread(paste0("Results/RITA/", x))}) # read through all the files
SLCs_high_expression = lapply(SLCs, function(x){unlist(x[expression_type == "high expression", cell, by = cell][,1])})
SLCs_low_expression = lapply(SLCs, function(x){unlist(x[expression_type == "low expression", cell, by = cell][,1])})

draw = function(g){
  ggplot(g, aes(x=log10(dose), y=viability, fill =expression_type, color = expression_type)) + 
    scale_fill_manual(values=c("#63BBB0", "#AC667E")) +
    scale_color_manual(values=c("#63BBB0", "#AC667E")) +
    geom_point(position = position_dodge(0.2), alpha = 0.5) + 
    geom_smooth(method = 'loess') + 
    scale_y_continuous(limits = c(-7.5, 2.5), breaks = seq(-7.5, 2.5, 2.5)) +
    theme(
      legend.position="none", plot.title = element_text(size=11),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"))
}
glist = lapply(SLCs, draw)

# --------------------------- #
#       Figure S5A data       #
# --------------------------- #

for (i in 1:length(glist)){
  print(glist[[i]])
}

# --------------------------- #
#       Figure S5B data       #
# --------------------------- #

celllinecount_high = as.data.table(table(unlist(SLCs_high_expression)))
celllinecount_low = as.data.table(table(unlist(SLCs_low_expression)))

ggplot(a, aes(x = N)) +
  geom_histogram(fill = "#63BBB0", color="#e9ecef", binwidth = 1) +
  scale_x_continuous(breaks = seq(1,8,1)) +
  scale_y_continuous(limits = c(0, 150)) +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(b, aes(x = N)) +
  geom_histogram(fill = "#AC667E", color="#e9ecef", binwidth = 1) +
  scale_x_continuous(breaks = seq(1,8,1)) +
  scale_y_continuous(limits = c(0, 150)) +
  theme(
    legend.position="none", plot.title = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

# 15 Evaluating previous prediction ####
# Note: Prediction evaluations are provided as a supplementary table (S14).
#       The code below demonstrates how predictions are evaluated.

# The previous predictions are derived from from Meixner et al. (2020) (https://www.embopress.org/doi/full/10.15252/msb.20209652)
# 1: Download Figure EV4, save as "heatmap.png", drag under "~/DeorphanSLC_manuscript/Raw Data/" directory
# 2: Download Table EV1, save the sheet "slc_annotation" as a separate .csv file "SLC annotation 2020.csv",
#    drag under "~/DeorphanSLC_manuscript/Raw Data/" directory

heatmap = load.image("Raw Data/heatmap.png")
# divide channels
channel1 = heatmap[,,1]
channel2 = heatmap[,,2]
channel3 = heatmap[,,3]
channel4 = heatmap[,,4]

# Function to extract color of each cell
extract_cell_color = function(x, y) {
  # Define cell dimensions and grid line width
  cell_width = 78
  cell_height = 18
  grid_line_width = 3.37
  # Calculate coordinates of the cell
  cell_x = (x-1) * (cell_height + grid_line_width) + grid_line_width
  cell_y = (y-1) * (cell_width + grid_line_width) + grid_line_width
  # Extract color of the cell
  cell_color = rgb(red = median(channel1[cell_y:(cell_y + cell_width), cell_x:(cell_x + cell_height)]),
                   green = median(channel2[cell_y:(cell_y + cell_width),cell_x:(cell_x + cell_height)]),
                   blue = median(channel3[cell_y:(cell_y + cell_width),cell_x:(cell_x + cell_height)]),
                   alpha = median(channel4[cell_y:(cell_y + cell_width),cell_x:(cell_x + cell_height)]))
  return(cell_color)
}

# Convert colors to values
Spectrum_value = seq(1,0,-0.04)
Spectrum_color = c("#D40A00","#DB3826","#E24B32","#EC5D3B","#F57747",
                   "#F78C51","#F8A45C","#FCB768","#FBC879","#FCDC8C",
                   "#FCE69A","#FEF2AA","#FDF7BB","#F9F7CA","#F1F5D5",
                   "#E8F4EA","#DAF0F6","#C9E8EF","#B6DEEA","#A6D3E4",
                   "#92C3DB","#80B6D5","#6AA1C9","#5C91C3","#4C7DB8",
                   "#3670B3")

# Function to calculate color distance
color_distance <- function(color1, color2) {
  sum((as.integer(grDevices::col2rgb(color1)) - as.integer(grDevices::col2rgb(color2)))^2)
}

# Function to map cell color to spectrum value
map_color_to_spectrum <- function(cell_color) {
  distances <- sapply(Spectrum_color, function(spectrum_color) color_distance(cell_color, spectrum_color))
  closest_index <- which.min(distances)
  return(Spectrum_value[closest_index])
}

# For each cell of table, record a number
cells_matrix = matrix(NaN, nrow = 115, ncol = 18)

for (x in 1:115) {
  for (y in 1:18) {
    cell_color = extract_cell_color(x, y)
    cells_matrix[x, y] = map_color_to_spectrum(cell_color)
  }
}

y_name = c("MFSD9","SLC35A5","NIPAL3","SLC35F5","SLC35F1",
           "SLC35F2","MFSD11","SLC35C2","MFSD5","SLC35G2",
           "CLN3","SLC35G4","SLC39A11","SLC35F6","TMEM165",
           "SLC35E1","SLC10A4","SLC26A10","SLC10A3","LETMD1",
           "SLC25A51","MFSD14A","MFSD6L","SLC35E2A","SLC35G5",
           "SLC35G6","SLC25A52","SLC25A34","SLC25A35","SFXN5",
           "SLC25A40","SLC23A3","MFSD13A","UNC93B1","SLC35E4",
           "SLC50A1","MFSD3","SLC35A4","SLC25A46","SLC27A3",
           "SLC44A3","ANKH","LETM2","SLC38A10","NIPAL2",
           "SLC10A5","SLC10A7","SLC44A5","SLC22A23","SLC35F4",
           "MFSD1","SLC22A10","SLC36A3","SLC38A11","MFSD8",
           "SLC38A6","SLC9B1","SLC16A6","SPNS3","SLC45A3",
           "MFSD14B","SLC35G1","SLC16A4","SLC38A8","SVOP",
           "SPNS1","DIRC2","SLC16A14","SLC7A4","MFSD4B",
           "SLC22A15","SFXN4","SLC25A39","SLC25A48","SLC35D3",
           "MFSD12","SLC22A17","SLC22A31","SLC35G3","TMEM241",
           "SLC16A13","MFSD4A","SLC46A2","SLC2A6","SLC22A18",
           "SLC35E3","SLC16A5","SLC49A3","SLC2A8","SLC7A13",
           "SLC25A14","SLC25A43","UNC93A","SLC22A25","SLC25A45",
           "SLC35E2B","SLCO5A1","SLCO6A1","SLC6A16","SLC6A17",
           "MFSD6","SLC9C2","OCA2","SLC9A4","SLC45A4",
           "SLC12A9","SV2A","SLC12A8","SV2B","SV2C",
           "SLC7A14","SVOPL","TMEM104","SLC15A5","SLC22A14")
x_name = c("Metal cation","Sodium(+)","Alkali metal cation","L-alpha-amino acid","Elemental halogen",
           "Inorganic anion","Alkaline earth cation","Magnesium(2+)","Potassium(+)","Carbohydrate",
           "Calcium(2+)","Carboxylic acid","Lipid","Steroid","Nucleobase-containing",
           "Divalent metal cation","Transition element cation","Zinc(2+)")

prediction = data.table(Orphan_SLC = y_name,
                        cells_matrix)
colnames(prediction) = c("Orphan SLC", x_name)

review = fread("Raw Data/SLC annotation 2023.csv")
old_review = fread("Raw Data/SLC annotation 2020.csv")

orphan_set = c()
for (i in 1:nrow(prediction)){
  find = grep(prediction[i, `Orphan SLC`], review[, HGNC_symbol])
  if (length(find) == 1){
    orphan_set = c(orphan_set, find)
  }else if (prediction[i, `Orphan SLC`] == "MFSD1"){
    find = which(review[, HGNC_symbol] == "MFSD1")
    orphan_set = c(orphan_set, find)
  }else if (prediction[i, `Orphan SLC`] == "SVOP"){
    find = which(review[, HGNC_symbol] == "SVOP (SLC22B4)")
    orphan_set = c(orphan_set, find)
  }else if (prediction[i, `Orphan SLC`] == "MFSD6"){
    find = which(review[, HGNC_symbol] == "MFSD6")
    orphan_set = c(orphan_set, find)
  }
  print(paste(prediction[i, `Orphan SLC`], ":", find))
}

orphan_set = review[orphan_set]

orphan_set_old = c()
for (i in 1:nrow(prediction)){
  find = which(old_review[, HGNC_symbol] == prediction[i, `Orphan SLC`])
  orphan_set_old = c(orphan_set_old, find)
}

old_orphan_set = old_review[orphan_set_old]

term_predict = c()
for (i in 1:nrow(prediction)){
  atomic_data = paste(colnames(prediction[,2:ncol(prediction)])[prediction[i,2:ncol(prediction)]>0.5], collapse = ";")
  term_predict = c(term_predict, atomic_data)
}

prediction_table = data.table(`Orphan SLC` = prediction[, `Orphan SLC`],
                              Prediction = term_predict)
orphan_set[, Prediction := term_predict]
fwrite(orphan_set, "Results/Prediction evaluations.csv")




