suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(ggplot2))


samples_desc <- read_csv("../../../EXP22002069_auto/EXP22002069_samplesheet_auto.csv")
# samples_desc$batch <- factor(samples_desc$batch,unique(samples_desc$batch))

g.count <- readRDS("../../../EXP22002069_auto/Counts/salmon.merged.gene_counts.rds")
colnames(g.count) <- as.character(lapply(colnames(g.count),function(x) str_replace(x,"X","")) )
g.count.matrix <- assay(g.count,"counts")
head(g.count.matrix)

t.count <- readRDS("../../../EXP22002069_auto/Counts/salmon.merged.transcript_counts.rds")
colnames(t.count) <- as.character(lapply(colnames(t.count),function(x) str_replace(x,"X","")) )
t.count.matrix <- assay(t.count,"counts")
head(t.count.matrix)





data = read.table("../../../../HUVEC_RNA_SEQ/Ensemble_to_geneid_BIOMART.txt",sep="\t",check.names = F,header = T)
data <- data[data$transcript_biotype %in% c("protein_coding","lncRNA"),]

tx2gene <- read.csv("../../../EXP22002069_auto/Counts/salmon_tx2gene.tsv", sep="\t", header = FALSE)
colnames(tx2gene) = c("tx", "gene_id", "gene_name")

tx2gene <- merge(tx2gene,data,by.x="gene_id",by.y="ensembl_gene_id",all.x=T)  
# tx2gene <- tx2gene[!duplicated(tx2gene[,c("gene_id","gene_name")]),]






library(tidyverse)
library(ggpubr)

samplesheet <- read_csv("~/Google Drive/Shared drives/Raj Lab/Lab stuff/ProjectFluorescenceSurrogate/Data/I25.Line1ASO_HUVEC_Oct2022/StuffQi/SampleSheet_qy.csv")
samplesheet <- samplesheet[samplesheet$Basename %in% samples_desc$methylation_file,]

load("~/Google Drive/Shared drives/Raj Lab/Lab stuff/ProjectFluorescenceSurrogate/Data/I25.Line1ASO_HUVEC_Oct2022/StuffQi/dat0Noob.RData")
dat0 <- dat0 %>% tibble::column_to_rownames(var = "ID")
dat0 <- dat0[, names(dat0) %in% samplesheet$Basename]

epic <- data.table::fread("../../../../EPIC/EPICAnnotation_age_DNAmTL.txt") 
epic <- epic[,c("CpG","annotation","ENSEMBL","SYMBOL")]
# %>%
#   mutate(age_associated=(P.age.JHS < pvalue_cut)) %>%
#   mutate(methylated=ifelse((Z.age.JHS>0 & P.age.JHS < pvalue_cut),"HYPER",ifelse((Z.age.JHS<0 & P.age.JHS < pvalue_cut),"HYPO","NOCHANGE"))) 

# HUVEC ten age dataset

huvec_10_samples_desc <- read_csv("../../../../HUVEC_RNA_SEQ/nf-run-auto/huvec_metadata_samplesheet_auto.csv")
# samples_desc$batch <- factor(samples_desc$batch,unique(samples_desc$batch))

huvec_10.g.count <- readRDS("../../../../HUVEC_RNA_SEQ/nf-run-auto/Counts/salmon.merged.gene_counts.rds")
colnames(huvec_10.g.count) <- as.character(lapply(colnames(huvec_10.g.count),function(x) str_replace(x,"X","")) )
huvec_10.g.count.matrix <- assay(huvec_10.g.count,"counts")
head(huvec_10.g.count.matrix)


# HUVEC proteomics preliminary run dataset

library(tibble)
library(reshape2)
library(dplyr)
library(ggpubr)
p_meta <- readxl::read_xlsx("../../../PROTEOMICS/ExperimentalDesign_Template MS SK.xlsx",na = "N/A")
p_meta$VG_LABEL <- paste0(p_meta$`Cell Line`,"_",p_meta$`Generic Condition`,"_",p_meta$Compound)
proteo <- data.table::fread("../../../PROTEOMICS/20230224_KeRa1_simplified_PG_Matrix.txt",na.strings = "NaN")
proteo <- proteo[-1,]


save(list=c("samples_desc","g.count.matrix","t.count.matrix","dat0","tx2gene","epic","huvec_10_samples_desc","huvec_10.g.count.matrix","p_meta","proteo"),file = "../app_db/app_data.Rdata")
