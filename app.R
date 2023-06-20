# unloadNamespace("shiny")
# renv::init(project = "/Users/vipulgupta/Google Drive/Shared drives/Raj Lab/VipulGupta/HUVEC_RNA_SEQ/SHINYAPPS/HUVEC/")
# library(BiocManager)
# options(repos = BiocManager::repositories())
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(shinydashboard))# library(shinydashboard)
suppressPackageStartupMessages(library(shinythemes))# library(shinythemes)
suppressPackageStartupMessages(library(leaflet))
suppressPackageStartupMessages(library(shinyWidgets))
# suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(AnnotationHub))
suppressPackageStartupMessages(library(karyoploteR))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(rtracklayer))
# renv::snapshot()

########################################################### 
#Read the Meta data file
########################################################### 


load(file = "app_db/app_data.Rdata")
load(file = "app_db/browser_data.Rdata")

########################################################### 
#query gene list
########################################################### 
genes <- tx2gene[!duplicated(tx2gene[,c("gene_id","gene_name","transcript_biotype")]),]
duplicated_genes <- genes$gene_name[duplicated(genes$gene_name)]
genes <- genes[!genes$gene_name %in% duplicated_genes,]$gene_name

########################################################### 
#gene_level_expression v/s age plot for huvec young/middle/old
########################################################### 
genes1 <- tx2gene[!duplicated(tx2gene[,c("gene_id","gene_name","transcript_biotype")]),]

exp_line_plot <- function(gene_query){
  
  tx2gene <- tx2gene[tx2gene$transcript_biotype %in% c("protein_coding","lncRNA"),]
  counts <- g.count.matrix[unique(tx2gene[tx2gene$gene_name %in% gene_query,]$gene_id),]
  rownames(counts) <- plyr::mapvalues(rownames(counts),genes1$gene_id,genes1$gene_name,warn_missing = F)
  
  df.m <- reshape2::melt(as.matrix(log2(counts+1)),value.name = "counts")
  df.m$samples <- plyr::mapvalues(df.m$Var2,samples_desc$sample,samples_desc$Group)
  # df.m$samples <-factor(df.m$samples,levels=c("CPD19","CPD127","CPD275")) 
  df.m$samples <-factor(df.m$samples,levels=c("YOUNG","MEDIUM","OLD")) 
  df.m$age <- plyr::mapvalues(df.m$Var2,samples_desc$sample,as.numeric(samples_desc$mAge))
  df.m$age <- as.numeric(levels(df.m$age))[df.m$age]
  
  ggpubr::ggline(df.m, x = "samples", y = "counts", color = "Var1",
                 add = c("mean_se", "jitter"),palette="npg",numeric.x.axis = F,ylab="log2(Counts+1)", xlab="HUVEC cells")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits = c(min(log2(g.count.matrix+1)),max(log2(g.count.matrix+1)))) +
    theme_bw()+
    theme(text = element_text(size=16),axis.text.x = element_text())
  
}

########################################################### 
#gene_level_expression v/s age plot for huvec 10 time points
########################################################### 


exp_line_plot_huvec10 <- function(gene_query){
  
  tx2gene <- tx2gene[tx2gene$transcript_biotype %in% c("protein_coding","lncRNA"),]
  counts <- huvec_10.g.count.matrix[unique(tx2gene[tx2gene$gene_name %in% gene_query,]$gene_id),]
  rownames(counts) <- plyr::mapvalues(rownames(counts),genes1$gene_id,genes1$gene_name,warn_missing = F)
  
  df.m <- reshape2::melt(as.matrix(log2(counts+1)),value.name = "counts")
  df.m$samples <- plyr::mapvalues(df.m$Var2,huvec_10_samples_desc$sample,huvec_10_samples_desc$CPD)
  # df.m$samples <-factor(df.m$samples,levels=c("CPD19","CPD127","CPD275")) 
  # df.m$samples <-factor(df.m$samples,levels=c("Young","Middle","Old")) 
  df.m$age <- plyr::mapvalues(df.m$Var2,huvec_10_samples_desc$sample,as.numeric(huvec_10_samples_desc$mAge))
  df.m$age <- as.numeric(levels(df.m$age))[df.m$age]
  
  ggpubr::ggline(df.m, x = "age", y = "counts", color = "Var1",
                 add = c("mean_se", "jitter"),palette="npg",numeric.x.axis = T,ylab="log2(Counts+1)", xlab="HUVEC mAge")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits = c(min(log2(huvec_10.g.count.matrix+1)),max(log2(huvec_10.g.count.matrix+1)))) +
    theme_bw()+
    theme(text = element_text(size=16),axis.text.x = element_text())
  
}


########################################################### 
#EZH2-inhibition by GSK in HUVEC
########################################################### 


ezh2_inihibition_plot <- function(gene_query){
  
  choose_samples <- c("CTRL_GSK","GSK")
  meta <- meta_exo_ezh2 %>% filter(.,VG_label %in% choose_samples)
  tx2gene <- tx2gene[tx2gene$transcript_biotype %in% c("protein_coding","lncRNA"),]
  counts <- ezh2_exo_count[unique(tx2gene[tx2gene$gene_name %in% gene_query,]$gene_id),colnames(ezh2_exo_count) %in% meta$sample]
  rownames(counts) <- plyr::mapvalues(rownames(counts),genes1$gene_id,genes1$gene_name,warn_missing = F)
  
  df.m <- reshape2::melt(as.matrix(log2(counts+1)),value.name = "counts")
  df.m$CPD <- plyr::mapvalues(df.m$Var2,meta$sample,meta$CPD) %>% as.character() %>% as.numeric()
  df.m$VG_label <- plyr::mapvalues(df.m$Var2,meta$sample,meta$VG_label)
  df.m$gene_name <- plyr::mapvalues(df.m$Var1,tx2gene$gene_id,tx2gene$gene_name,warn_missing = F)
  
  ggpubr::ggline(df.m, x = "CPD", y = "counts", color = "VG_label",facet.by = "Var1",
                 palette="npg",numeric.x.axis = T,ylab="log2(Counts+1)", xlab="CPD",title = "EZH2-inhibitor treated HUVECs")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits = c(min(log2(ezh2_exo_count+1)),max(log2(ezh2_exo_count+1)))) +
    theme_bw()+
    theme(text = element_text(size=16),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
}

########################################################### 
#EZH2-inhibition by GSK in HUVEC
########################################################### 


exo_inihibition_plot <- function(gene_query){
  
  choose_samples <- c("CTRL_EXO","EXO")
  meta <- meta_exo_ezh2 %>% filter(.,VG_label %in% choose_samples)
  tx2gene <- tx2gene[tx2gene$transcript_biotype %in% c("protein_coding","lncRNA"),]
  counts <- ezh2_exo_count[unique(tx2gene[tx2gene$gene_name %in% gene_query,]$gene_id),colnames(ezh2_exo_count) %in% meta$sample]
  rownames(counts) <- plyr::mapvalues(rownames(counts),genes1$gene_id,genes1$gene_name,warn_missing = F)
  
  df.m <- reshape2::melt(as.matrix(log2(counts+1)),value.name = "counts")
  df.m$CPD <- plyr::mapvalues(df.m$Var2,meta$sample,meta$CPD) %>% as.character() %>% as.numeric()
  df.m$VG_label <- plyr::mapvalues(df.m$Var2,meta$sample,meta$VG_label)
  df.m$gene_name <- plyr::mapvalues(df.m$Var1,tx2gene$gene_id,tx2gene$gene_name,warn_missing = F)
  
  ggpubr::ggline(df.m, x = "CPD", y = "counts", color = "VG_label",facet.by = "Var1",
                 palette="npg",numeric.x.axis = T,ylab="log2(Counts+1)", xlab="CPD",title = "Exosome treated HUVECs")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits = c(min(log2(ezh2_exo_count+1)),max(log2(ezh2_exo_count+1)))) +
    theme_bw()+
    theme(text = element_text(size=16),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
}

########################################################### 
#REST-OE in HUVEC
########################################################### 


rest_oe_plot <- function(gene_query){
  
  choose_samples <- c("REST_CTRL","REST_OE")
  meta <- meta_exo_ezh2 %>% filter(.,VG_label %in% choose_samples)
  tx2gene <- tx2gene[tx2gene$transcript_biotype %in% c("protein_coding","lncRNA"),]
  counts <- ezh2_exo_count[unique(tx2gene[tx2gene$gene_name %in% gene_query,]$gene_id),colnames(ezh2_exo_count) %in% meta$sample]
  rownames(counts) <- plyr::mapvalues(rownames(counts),genes1$gene_id,genes1$gene_name,warn_missing = F)
  
  df.m <- reshape2::melt(as.matrix(log2(counts+1)),value.name = "counts")
  df.m$CPD <- plyr::mapvalues(df.m$Var2,meta$sample,meta$CPD) %>% as.character() %>% as.numeric()
  df.m$VG_label <- plyr::mapvalues(df.m$Var2,meta$sample,meta$VG_label)
  df.m$VG_label <- factor(df.m$VG_label,levels = choose_samples)
  df.m$gene_name <- plyr::mapvalues(df.m$Var1,tx2gene$gene_id,tx2gene$gene_name,warn_missing = F)
  
  ggpubr::ggboxplot(df.m, x = "VG_label", y = "counts", color = "VG_label",facet.by = "Var1",
                 palette="npg",numeric.x.axis = T,ylab="log2(Counts+1)", xlab="CPD",title = "REST-OE in HUVECs")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits = c(min(log2(ezh2_exo_count+1)),max(log2(ezh2_exo_count+1)))) +
    theme_bw()+
    theme(text = element_text(size=16),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
}

########################################################### 
#MVEC old middle and neonatal primary cells 2 biological replicates
########################################################### 


mvec_plot <- function(gene_query){
  
  choose_samples <- c("MVEC_NEONATAL","MVEC_MIDDLE","MVEC_OLD")
  meta <- meta_exo_ezh2 %>% filter(.,VG_label %in% choose_samples)
  tx2gene <- tx2gene[tx2gene$transcript_biotype %in% c("protein_coding","lncRNA"),]
  counts <- ezh2_exo_count[unique(tx2gene[tx2gene$gene_name %in% gene_query,]$gene_id),colnames(ezh2_exo_count) %in% meta$sample]
  rownames(counts) <- plyr::mapvalues(rownames(counts),genes1$gene_id,genes1$gene_name,warn_missing = F)
  
  df.m <- reshape2::melt(as.matrix(log2(counts+1)),value.name = "counts")
  df.m$CPD <- plyr::mapvalues(df.m$Var2,meta$sample,meta$CPD) %>% as.character() %>% as.numeric()
  df.m$VG_label <- plyr::mapvalues(df.m$Var2,meta$sample,meta$VG_label)
  df.m$VG_label <- factor(df.m$VG_label,levels = choose_samples)
  df.m$gene_name <- plyr::mapvalues(df.m$Var1,tx2gene$gene_id,tx2gene$gene_name,warn_missing = F)
  
  ggpubr::ggboxplot(df.m, x = "VG_label", y = "counts", color = "VG_label",facet.by = "Var1",
                    palette="npg",numeric.x.axis = T,ylab="log2(Counts+1)", xlab="CPD",title = "MVEC RNAseq")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits = c(min(log2(ezh2_exo_count+1)),max(log2(ezh2_exo_count+1)))) +
    theme_bw()+
    theme(text = element_text(size=16),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
}


########################################################### 
#Keratinocytes Rapamycin, Forskolin or Rapa+Forskolin RNAseq
########################################################### 


kera_rapa_fors_plot <- function(gene_query){
  
  choose_samples <- c("FSK_ctr","FSK_Rapa","FSK_Forsk","FSK_R_F")
  meta <- meta_exo_ezh2 %>% filter(.,VG_label %in% choose_samples)
  tx2gene <- tx2gene[tx2gene$transcript_biotype %in% c("protein_coding","lncRNA"),]
  counts <- ezh2_exo_count[unique(tx2gene[tx2gene$gene_name %in% gene_query,]$gene_id),colnames(ezh2_exo_count) %in% meta$sample]
  rownames(counts) <- plyr::mapvalues(rownames(counts),genes1$gene_id,genes1$gene_name,warn_missing = F)
  
  df.m <- reshape2::melt(as.matrix(log2(counts+1)),value.name = "counts")
  df.m$CPD <- plyr::mapvalues(df.m$Var2,meta$sample,meta$CPD) %>% as.character() %>% as.numeric()
  df.m$VG_label <- plyr::mapvalues(df.m$Var2,meta$sample,meta$VG_label)
  df.m$VG_label <- factor(df.m$VG_label,levels = choose_samples)
  df.m$gene_name <- plyr::mapvalues(df.m$Var1,tx2gene$gene_id,tx2gene$gene_name,warn_missing = F)
  
  ggpubr::ggboxplot(df.m, x = "VG_label", y = "counts", color = "VG_label",facet.by = "Var1",
                    palette="npg",numeric.x.axis = T,ylab="log2(Counts+1)", xlab="CPD",title = "Kerationocytes")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits = c(min(log2(ezh2_exo_count+1)),max(log2(ezh2_exo_count+1)))) +
    theme_bw()+
    theme(text = element_text(size=16),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
}

########################################################### 
#transcript_level_expression v/s age plot  for huvec young/middle/old
########################################################### 

transcript_age_plot <- function(gene_query){
  
  counts <- t.count.matrix[tx2gene[tx2gene$gene_name %in% gene_query,]$tx,] 
  df.m <- reshape2::melt(as.matrix(log2(counts+1)),value.name = "counts")
  # df.m$samples <- plyr::mapvalues(df.m$Var2,samples_desc$sample,paste0(samples_desc$CPD,"_",samples_desc$Repeat))
  df.m$samples <- plyr::mapvalues(df.m$Var2,samples_desc$sample,samples_desc$Group)
  # df.m$samples <-factor(df.m$samples,levels=c("CPD19","CPD127","CPD275")) 
  df.m$samples <-factor(df.m$samples,levels=c("YOUNG","MEDIUM","OLD")) 
  df.m$age <- plyr::mapvalues(df.m$Var2,samples_desc$sample,as.numeric(samples_desc$mAge))
  df.m$age <- as.numeric(levels(df.m$age))[df.m$age]
  ggpubr::ggline(df.m, x = "samples", y = "counts", facet.by = "Var1",
                 add = c("mean_se", "jitter"),palette="npg",numeric.x.axis = T,ylab="log2(Counts+1)", xlab="HUVECs Age")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits = c(min(log2(t.count.matrix+1)),max(log2(t.count.matrix+1)))) +
    theme_bw()+
    theme(text = element_text(size=16),axis.text.x = element_text())

} 


########################################################### 
#Methylation v/s age plot for huvec young/middle/old
########################################################### 

meth_age_plot <- function(gene_query){
  
  
  
  gene_cpg <- epic[epic$SYMBOL %in% gene_query,]
  counts <- dat0[rownames(dat0) %in% gene_cpg$CpG,]
  df.m <- reshape2::melt(as.matrix(counts),value.name = "counts")
  df.m$samples <- plyr::mapvalues(df.m$Var2,samples_desc$methylation_file,samples_desc$Group)
  df.m$samples <-factor(df.m$samples,levels=c("YOUNG","MEDIUM","OLD"))  
  df.m$age <- plyr::mapvalues(df.m$Var2,samples_desc$methylation_file,as.numeric(samples_desc$mAge))
  df.m$age <- as.numeric(levels(df.m$age))[df.m$age]
  df.m$annotation <- plyr::mapvalues(df.m$Var1,gene_cpg$CpG,paste0(gene_cpg$CpG,":",gene_cpg$annotation))
  ggpubr::ggline(df.m, x = "samples", y = "counts",facet.by = "annotation",
                 add = c("mean_se", "jitter"),palette="npg",numeric.x.axis = T,ylab="Beta values", xlab="HUVECs Age")+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits = c(0,1)) +
    theme_bw()+
    theme(text = element_text(size=7),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
 
}


########################################################### 
#Proteomics dataset plot for huvec and keratinocytes including HUVEC young/middle/old
########################################################### 

prot_plot_fun <- function(query_gene){
  
  # proteo_query <- proteo[proteo$Genes %in% query_gene,] %>% as.data.frame()
  # proteo_query <- reshape2::melt(t(proteo_query[,names(proteo_query)[-(61:65)]]),value.name = "proteomics_intensity")
  # 
  # proteo_query$VG_LABEL <- plyr::mapvalues(proteo_query$Var1,p_meta$`MS Hub ID`,p_meta$VG_LABEL)
  # proteo_query$proteomics_intensity <- as.numeric(proteo_query$proteomics_intensity )
  # proteo_query <- proteo_query %>% mutate(proteomics_intensity = if_else(is.na(proteomics_intensity), 0, proteomics_intensity))
  # proteo_query$`Tissue Type` <- plyr::mapvalues(proteo_query$Var1,p_meta$`MS Hub ID`,p_meta$`Tissue Type`)
  # 
  # ggpubr::ggboxplot(proteo_query %>% as.data.frame(),x="VG_LABEL",y="proteomics_intensity",color = "Tissue Type",add = "jitter",title = query_gene,
  #                   palette="npg",ylab="log2(Protein abundance)", xlab="Samples")+
  #   theme_bw()+
  #   theme(text = element_text(size=12),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  proteo_query <- proteo_tmt[proteo_tmt$`Gene Symbol` %in% query_gene,] %>% as.data.frame()
  rownames(proteo_query) <- proteo_query$`Gene Symbol`
  proteo_query <- proteo_query[,!names(proteo_query) %in% c("Gene Symbol")]
  proteo_query <- reshape2::melt(t(proteo_query),value.name = "proteomics_intensity",na.rm = F)
  
  
  meta_exo_ezh2_prot <- meta_exo_ezh2 %>% separate_rows(Proteomics_sample, sep=";", convert = TRUE)
  proteo_query$VG_LABEL <- plyr::mapvalues(proteo_query$Var1,meta_exo_ezh2_prot$Proteomics_sample,meta_exo_ezh2_prot$VG_label)
  proteo_query$proteomics_intensity <- as.numeric(proteo_query$proteomics_intensity )
  proteo_query <- proteo_query %>% mutate(proteomics_intensity = if_else(is.na(proteomics_intensity), 0, proteomics_intensity))
  proteo_query$`Tissue Type` <- plyr::mapvalues(proteo_query$Var1,meta_exo_ezh2_prot$Proteomics_sample,meta_exo_ezh2_prot$CellType)
  proteo_query$CPD <- plyr::mapvalues(proteo_query$Var1,meta_exo_ezh2_prot$Proteomics_sample,meta_exo_ezh2_prot$CPD)
  proteo_query$Batch <- plyr::mapvalues(proteo_query$Var1,meta_exo_ezh2_prot$Proteomics_sample,meta_exo_ezh2_prot$proteomics_batch)
  
  
  
  huvec_ymo <- proteo_query[proteo_query$VG_LABEL %in% c("HUVEC_YOUNG","HUVEC_MEDIUM","HUVEC_OLD"),] 
  p_ymo <- ggpubr::ggboxplot(huvec_ymo,x="VG_LABEL",y="proteomics_intensity",add = "jitter",
                             palette="npg",ylab="log2(Protein abundance)", xlab="")+
    theme_bw()+
    theme(text = element_text(size=12),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  huvec_rest_oe <- proteo_query[proteo_query$VG_LABEL %in% c("REST_CTRL","REST_OE"),] 
  huvec_rest_oe$VG_LABEL <- factor(huvec_rest_oe$VG_LABEL, levels = c("REST_CTRL","REST_OE"))
  p_rest <- ggpubr::ggboxplot(huvec_rest_oe,x="VG_LABEL",y="proteomics_intensity",add = "jitter",
                              palette="npg",ylab="log2(Protein abundance)", xlab="")+
    theme_bw()+
    theme(text = element_text(size=12),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  mvec <- proteo_query[proteo_query$VG_LABEL %in% c("MVEC_NEONATAL","MVEC_MIDDLE","MVEC_OLD"),] 
  # huvec_rest_oe$VG_LABEL <- factor(huvec_rest_oe$VG_LABEL, levels = c("REST_CTRL","REST_OE"))
  p_mvec <- ggpubr::ggboxplot(mvec,x="VG_LABEL",y="proteomics_intensity",add = "jitter",
                              palette="npg",ylab="log2(Protein abundance)", xlab="")+
    theme_bw()+
    theme(text = element_text(size=12),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  gsk <- proteo_query[proteo_query$VG_LABEL %in% c("CTRL_GSK","GSK"),] 
  gsk$CPD <- gsk$CPD %>% as.character() %>% as.numeric()
  p_gsk <- ggpubr::ggline(gsk,x="CPD",y="proteomics_intensity",add = "box",group = "VG_LABEL",color = "VG_LABEL",
                          palette="npg",ylab="log2(Protein abundance)", xlab="CPD",numeric.x.axis = T)+
    theme_bw()+
    theme(text = element_text(size=12),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  exo <- proteo_query[proteo_query$VG_LABEL %in% c("CTRL_EXO","EXO"),] 
  exo$CPD <- exo$CPD %>% as.character() %>% as.numeric()
  p_exo <- ggpubr::ggline(exo,x="CPD",y="proteomics_intensity",add = "box",group = "VG_LABEL",color = "VG_LABEL",
                          palette="npg",ylab="log2(Protein abundance)", xlab="CPD",numeric.x.axis = T)+
    theme_bw()+
    theme(text = element_text(size=12),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  kera_rapa <- proteo_query[proteo_query$VG_LABEL %in% c("FSK_ctr","FSK_Rapa","FSK_Forsk","FSK_R_F"),] 
  p_kera_rapa <- ggpubr::ggboxplot(kera_rapa,x="VG_LABEL",y="proteomics_intensity",add = "jitter",
                                   palette="npg",ylab="log2(Protein abundance)", xlab="")+
    theme_bw()+
    theme(text = element_text(size=12),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggarrange(p_ymo, p_rest, p_mvec, p_gsk,p_exo,p_kera_rapa,
            labels = c("HUVEC-YMO", "REST-OE", "MVEC","GSK","EXO","Keratinocytes"),
            ncol = 3, nrow = 2)
}
  
########################################################### 
#Plot ChIP seq dataset for HUVEC YMO REST and EGR1
########################################################### 

  
plot_static_browser_view <- function(query_gene){
  
  pp <- getDefaultPlotParams(plot.type=1)
  pp$leftmargin <- 0.15
  pp$topmargin <- 15
  pp$bottommargin <- 15
  pp$ideogramheight <- 5
  pp$data1inmargin <- 10
  pp$data1outmargin <- 0
  
  myGeneSymbols <- AnnotationDbi::select(org.Hs.eg.db,
                                         keys = query_gene,
                                         columns = c("SYMBOL","ENTREZID"),
                                         keytype = "SYMBOL")
  
  myGeneSymbolsTx <- AnnotationDbi::select(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                           keys = myGeneSymbols$ENTREZID,
                                           columns = c("GENEID", "TXID", "TXCHROM", "TXSTART", "TXEND"),
                                           keytype = "GENEID")
  
  region <- merge(myGeneSymbols, myGeneSymbolsTx, by.x = "ENTREZID", by.y = "GENEID")
  region$length <- region$TXEND-region$TXSTART
  region <- region[which.max(abs(region$length)),]
  region <- toGRanges(paste0(c(region$TXCHROM,":",region$TXSTART,"-",region$TXEND),collapse = ""))
  
  region <- region+5000
  
  
  kp <- plotKaryotype(zoom = region, cex=1, plot.params = pp,plot.type = 2)
  kpAddBaseNumbers(kp, tick.dist = 10000, minor.tick.dist = 2000,
                   add.units = TRUE, cex=1, tick.len = 3)
  kpAddMainTitle(kp, paste0("HUVEC-YMO ChIP seq profiles for ",query_gene), cex=2)
  
  genes.data <- makeGenesDataFromTxDb(txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                      karyoplot=kp,
                                      plot.transcripts = TRUE,
                                      plot.transcripts.structure = TRUE)
  
  genes.data <- addGeneNames(genes.data)
  genes.data <- mergeTranscripts(genes.data)
  
  kpPlotGenes(kp, data=genes.data, r0=0, r1=0.02, gene.name.cex = 1)
  # kpPlotRegions(kp, epic_gr,col=colByChr(epic_gr$methylated,colors = cols), r0=0.2, r1=0.25,avoid.overlapping = F,clipping = F)
  
  kpPlotRegions(kp, data = epic_38_gr,col = epic_38_gr$jhs_color, r0=0.03, r1=0.06)
  kpAddLabels(kp, labels = "HUMAN-EWAS", r0=0.04, r1=0.05,cex=1)
  
  kpPlotRegions(kp, data = epic_38_gr,col = epic_38_gr$huvec_color, r0=0.07, r1=0.10)
  kpAddLabels(kp, labels = "HUVEC-Oldv/sYoung", r0=0.08, r1=0.09,cex=1)
  
  kpPlotRegions(kp, promoters, col="red", r0=0.11, r1=0.13)
  kpAddLabels(kp, labels = "Promoters", r0=0.12, r1=0.13,cex=1)
  
  kpPlotRegions(kp, hmm.model, col=hmm.model$RGB, r0=0.14, r1=0.16)
  kpAddLabels(kp, labels = "Chromatin\nState (HMM)", r0=0.14, r1=0.15,cex=1)
  
  s3_folder <- "/share/vgupta/workspaces/huvec-analysis/HUVEC_YMO_deeptools_VG_run_all_chipseq/input_normalized_bigwigs/"
  # s3_folder <- "~/Google Drive/Shared drives/Raj Lab/VipulGupta/HUVEC_YMO/CHIPseq/Diffbind_all_samples/bw/"
  
  rest_files <- c(OLD_REST="OLD_REST.bw",MIDDLE_REST="MEDIUM_REST.bw",YOUNG_REST="YOUNG_REST.bw")
  egr1_files <- c(OLD_EGR1="OLD_EGR1.bw",MIDDLE_EGR1="MEDIUM_EGR1.bw",YOUNG_EGR1="YOUNG_EGR1.bw")
  ezh2_files <- c(OLD_EZH2="OLD_EZH2.bw",MIDDLE_EZH2="MEDIUM_EZH2.bw",YOUNG_EZH2="YOUNG_EZH2.bw")
  jarid2_files <- c(OLD_JARID2="OLD_JARID2.bw",MIDDLE_JARID2="MEDIUM_JARID2.bw",YOUNG_JARID2="YOUNG_JARID2.bw")
  k4me3_files <- c(OLD_K4me3="OLD_K4me3.bw",MIDDLE_K4me3="MEDIUM_K4me3.bw",YOUNG_K4me3="YOUNG_K4me3.bw")
  k27me3_files <- c(OLD_K27me3="OLD_K27me3.bw",MIDDLE_K27me3="MEDIUM_K27me3.bw",YOUNG_K27me3="YOUNG_K27me3.bw")
  
  
  
  
  computed.ymax.k27me3 <-0
  for(i in names(c(k27me3_files))){
    temp <- max(import(paste0(s3_folder,k27me3_files[i]),which=region)$score)
    computed.ymax.k27me3 <- ifelse(computed.ymax.k27me3 > temp,computed.ymax.k27me3,temp)
    computed.ymax.k27me3 <- ifelse(round(computed.ymax.k27me3,0) > 1.5,round(computed.ymax.k27me3,0),2)
  }
  computed.ymax.rest <-0
  for(i in names(c(rest_files))){
    temp <- max(import(paste0(paste0(s3_folder,rest_files[i])),which=region)$score)
    computed.ymax.rest <- ifelse(computed.ymax.rest > temp,computed.ymax.rest,temp)
    computed.ymax.rest <- ifelse(round(computed.ymax.rest,0) > 1.5,round(computed.ymax.rest,0),2)
  }
  computed.ymax.egr1 <-0
  for(i in names(c(egr1_files))){
    temp <- max(import(paste0(paste0(s3_folder, egr1_files[i])),which=region)$score)
    computed.ymax.egr1 <- ifelse(computed.ymax.egr1 > temp,computed.ymax.egr1,temp)
    computed.ymax.egr1 <- ifelse(round(computed.ymax.egr1,0) > 1.5,round(computed.ymax.egr1,0),2)
  }
  computed.ymax.k4me3 <-0
  for(i in names(k4me3_files)){
    temp <- max(import(paste0(s3_folder, k4me3_files[i]),which=region)$score)
    computed.ymax.k4me3 <- ifelse(computed.ymax.k4me3 > temp,computed.ymax.k4me3,temp)
    computed.ymax.k4me3 <- ifelse(round(computed.ymax.k4me3,0) > 1.5,round(computed.ymax.k4me3,0),2)
  }
  
  computed.ymax.ezh2 <-0
  for(i in names(c(ezh2_files))){
    temp <- max(import(paste0(paste0(s3_folder, ezh2_files[i])),which=region)$score)
    computed.ymax.ezh2 <- ifelse(computed.ymax.ezh2 > temp,computed.ymax.ezh2,temp)
    computed.ymax.ezh2 <- ifelse(round(computed.ymax.ezh2,0) > 1.5,round(computed.ymax.ezh2,0),2)
  }
  computed.ymax.jarid2 <-0
  for(i in names(c(jarid2_files))){
    temp <- max(import(paste0(paste0(s3_folder, jarid2_files[i])),which=region)$score)
    computed.ymax.jarid2 <- ifelse(computed.ymax.jarid2 > temp,computed.ymax.jarid2,temp)
    computed.ymax.jarid2 <- ifelse(round(computed.ymax.jarid2,0) > 1.5,round(computed.ymax.jarid2,0),2)
  }
  
  total.tracks <- length(k27me3_files)+length(k4me3_files)+length(ezh2_files)+length(jarid2_files)+length(rest_files)+length(egr1_files)
  
  #K27me3 Track
  out.at <- autotrack(1:3, total.tracks,  r0=0.17,r1=1,margin = 0.2)
  for(i in seq_len(length(k27me3_files))) {
    bigwig.file <- paste0(s3_folder, k27me3_files[i])
    at <- autotrack(i, 3, r0=out.at$r0, r1=out.at$r1,margin = 0.1)
    # kp <- kpPlotBigWig(kp, data=bigwig.file, ymax="visible.region",
    # kp <- kpPlotBigWig(kp, data=bigwig.file, ymax = "global",
    kp <- kpPlotBigWig(kp, data=bigwig.file, ymax=computed.ymax.k27me3,
                       r0=at$r0, r1=at$r1, col = "#AAFFAA")
    # computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
    # computed.ymax <- ifelse(computed.ymax > ceiling(kp$latest.plot$computed.values$ymax),computed.ymax,ceiling(kp$latest.plot$computed.values$ymax))
    kpAxis(kp, ymin=0, ymax=computed.ymax.k27me3, numticks = 2, r0=at$r0, r1=at$r1,cex=0.5)
    # kpAxis(kp, ymin=0, ymax=31, numticks = 2, r0=at$r0, r1=at$r1,cex=0.5)
    kpAddLabels(kp, labels = names(k27me3_files)[i], r0=at$r0, r1=at$r1,cex=1, label.margin = 0.035)
  }
  
  #k4me3 Track
  out.at <- autotrack(4:6, total.tracks,  r0=0.17,r1=1,margin = 0.2)
  for(i in seq_len(length(k4me3_files))) {
    bigwig.file <- paste0(s3_folder, k4me3_files[i])
    at <- autotrack(i, 3, r0=out.at$r0, r1=out.at$r1,margin = 0.1)
    # kp <- kpPlotBigWig(kp, data=bigwig.file, ymax="visible.region",
    # kp <- kpPlotBigWig(kp, data=bigwig.file, ymax = "global",
    kp <- kpPlotBigWig(kp, data=bigwig.file, ymax=computed.ymax.k4me3,
                       r0=at$r0, r1=at$r1, col = "#AAAAFF")
    # computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
    # computed.ymax <- ifelse(computed.ymax > ceiling(kp$latest.plot$computed.values$ymax),computed.ymax,ceiling(kp$latest.plot$computed.values$ymax))
    kpAxis(kp, ymin=0, ymax=computed.ymax.k4me3, numticks = 2, r0=at$r0, r1=at$r1,cex=0.5)
    # kpAxis(kp, ymin=0, ymax=31, numticks = 2, r0=at$r0, r1=at$r1,cex=0.5)
    kpAddLabels(kp, labels = names(k4me3_files)[i], r0=at$r0, r1=at$r1, cex=1, label.margin = 0.035)
  }
  #ezh2 Track
  out.at <- autotrack(7:9, total.tracks, r0=0.17,r1=1,margin = 0.2)
  for(i in seq_len(length(ezh2_files))) {
    bigwig.file <- paste0(s3_folder, ezh2_files[i])
    at <- autotrack(i, 3, r0=out.at$r0, r1=out.at$r1,margin = 0.1)
    # kp <- kpPlotBigWig(kp, data=bigwig.file, ymax="visible.region",
    # kp <- kpPlotBigWig(kp, data=bigwig.file, ymax = "global",
    kp <- kpPlotBigWig(kp, data=bigwig.file, ymax=computed.ymax.ezh2,
                       r0=at$r0, r1=at$r1, col = "#FFAAAA")
    # computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
    # computed.ymax <- ifelse(computed.ymax > ceiling(kp$latest.plot$computed.values$ymax),computed.ymax,ceiling(kp$latest.plot$computed.values$ymax))
    kpAxis(kp, ymin=0, ymax=computed.ymax.ezh2, numticks = 2, r0=at$r0, r1=at$r1,cex=0.5)
    # kpAxis(kp, ymin=0, ymax=31, numticks = 2, r0=at$r0, r1=at$r1,cex=0.5)
    kpAddLabels(kp, labels = names(ezh2_files)[i], r0=at$r0, r1=at$r1,cex=1, label.margin = 0.035)
  }
  
  #jarid2 Track
  out.at <- autotrack(10:12, total.tracks, r0=0.17,r1=1,margin = 0.2)
  for(i in seq_len(length(jarid2_files))) {
    bigwig.file <- paste0(s3_folder, jarid2_files[i])
    at <- autotrack(i, 3, r0=out.at$r0, r1=out.at$r1,margin = 0.1)
    # kp <- kpPlotBigWig(kp, data=bigwig.file, ymax="visible.region",
    # kp <- kpPlotBigWig(kp, data=bigwig.file, ymax = "global",
    kp <- kpPlotBigWig(kp, data=bigwig.file, ymax=computed.ymax.jarid2,
                       r0=at$r0, r1=at$r1, col = "#E4D00A")
    # computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
    # computed.ymax <- ifelse(computed.ymax > ceiling(kp$latest.plot$computed.values$ymax),computed.ymax,ceiling(kp$latest.plot$computed.values$ymax))
    kpAxis(kp, ymin=0, ymax=computed.ymax.jarid2, numticks = 2, r0=at$r0, r1=at$r1,cex=0.5)
    # kpAxis(kp, ymin=0, ymax=31, numticks = 2, r0=at$r0, r1=at$r1,cex=0.5)
    kpAddLabels(kp, labels = names(jarid2_files)[i], r0=at$r0, r1=at$r1,cex=1, label.margin = 0.035)
  }
  
  # REST track
  out.at <- autotrack(13:15, total.tracks,  r0=0.17,r1=1,margin = 0.2)
  for(i in seq_len(length(rest_files))) {
    bigwig.file <- paste0(s3_folder, rest_files[i])
    at <- autotrack(i, 3, r0=out.at$r0, r1=out.at$r1,margin = 0.1)
    # kp <- kpPlotBigWig(kp, data=bigwig.file, ymax="visible.region",
    # kp <- kpPlotBigWig(kp, data=bigwig.file, ymax = "global",
    kp <- kpPlotBigWig(kp, data=bigwig.file, ymax=computed.ymax.rest,
                       r0=at$r0, r1=at$r1, col = "cadetblue2")
    # computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
    # computed.ymax <- ifelse(computed.ymax > ceiling(kp$latest.plot$computed.values$ymax),computed.ymax,ceiling(kp$latest.plot$computed.values$ymax))
    kpAxis(kp, ymin=0, ymax=computed.ymax.rest, numticks = 2, r0=at$r0, r1=at$r1,cex=0.5)
    # kpAxis(kp, ymin=0, ymax=31, numticks = 2, r0=at$r0, r1=at$r1,cex=0.5)
    kpAddLabels(kp, labels = names(rest_files)[i], r0=at$r0, r1=at$r1,cex=1, label.margin = 0.035)
  }
  
  #EGR1 track
  out.at <- autotrack(16:18, total.tracks, r0=0.17,r1=1,margin = 0.2)
  
  for(i in seq_len(length(egr1_files))) {
    bigwig.file <- paste0(s3_folder, egr1_files[i])
    at <- autotrack(i, 3, r0=out.at$r0, r1=out.at$r1,margin = 0.1)
    # kp <- kpPlotBigWig(kp, data=bigwig.file, ymax="visible.region",
    # kp <- kpPlotBigWig(kp, data=bigwig.file, ymax = "global",
    kp <- kpPlotBigWig(kp, data=bigwig.file, ymax=computed.ymax.egr1,
                       r0=at$r0, r1=at$r1, col = "darkolivegreen1")
    # computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
    # computed.ymax <- ifelse(computed.ymax > ceiling(kp$latest.plot$computed.values$ymax),computed.ymax,ceiling(kp$latest.plot$computed.values$ymax))
    kpAxis(kp, ymin=0, ymax=computed.ymax.egr1, numticks = 2, r0=at$r0, r1=at$r1,cex=0.5)
    # kpAxis(kp, ymin=0, ymax=31, numticks = 2, r0=at$r0, r1=at$r1,cex=0.5)
    kpAddLabels(kp, labels = names(egr1_files)[i], r0=at$r0, r1=at$r1, cex=1, label.margin = 0.035)
  }
}






########################################################### 
### SHINY UI ###
########################################################### 
ui <- bootstrapPage(
  tags$head(includeHTML("gtag.html")),
  navbarPage(theme = shinytheme("flatly"), collapsible = TRUE,
             HTML('<a style="text-decoration:none;cursor:default;color:#FFFFFF;" class="active" href="#">HUVEC Young/Middle/Old Expression</a>'), id="nav",
             windowTitle = "HUVEC-y-m-o",
             
            
             tabPanel("PLOTS",
                      
                      sidebarLayout(
                        sidebarPanel(
                          span(tags$i(h5("Enter Gene names Check the expression")), style="color:#045a8d"),
                          selectizeInput('input_gene2', 'Gene name', genes,multiple = T,selected = c("TP53","TMEM240","DLX1"),options = list(delimiter = " ", create = T)),
                          submitButton("Submit"),
                          span(tags$i(h5("Click below to download the full datasets")), style="color:#045a8d"),
                          downloadButton('exp_mat',label="Expression"),
                          downloadButton('meth_mat',label="Methylation"),
                          width = 3
                          
                        ),
                        
                        mainPanel(
                          tabsetPanel(
                            # tabPanel("Bi-plot", plotlyOutput("biplot")),
                            tabPanel("Gene Level Expression",
                                     fluidRow(
                                       box(title="HUVEC Y/M/O (Generated@Altos)",solidHeader = TRUE,height = 900, width = 6, status = "primary",plotOutput("biplot")),
                                       box(title="HUVEC 10timepoint (legacy dataset)",solidHeader = TRUE,height = 900, width = 6, status = "primary",plotOutput("biplot2")),
                                       # box(title="Selected Genes normZ scores",solidHeader = TRUE,width=3,height = 900, status = "primary",dataTableOutput("brush_info"))
                                       )
                                     # column(width = 8,plotOutput("biplot",brush = brushOpts(id = "plot1_brush"))),
                                     # column(width = 4,h4("Selected Genes normZ scores"),dataTableOutput("brush_info"))
                                     ),
                            tabPanel("Transcript Level Expression", uiOutput('mytabs2')),
                            tabPanel("Methylation", uiOutput('mytabs')),
                            tabPanel("Proteomics",uiOutput('prot_plot')),
                            tabPanel("ChIPseq",uiOutput('chipseq_tab')),
                            tabPanel("Interventions",tabsetPanel(
                              tabPanel("EZH2-inibitor",plotOutput('ezh2_plot')),
                            tabPanel("Exosome",plotOutput('exo_plot')),
                            tabPanel("REST-OE",plotOutput('rest_plot')),
                            tabPanel("Rapa+Fors",plotOutput('kera_rapa_fors')))),
                            tabPanel("MVEC",plotOutput('mvec_plot'))
                            
                          )
                        )
                      )
             )
             )
)



server <- function(input, output, session) {
  
  output$mytabs = renderUI({
    nTabs = input$input_gene2
    
    myTabs = lapply(nTabs, function(x){
      tabPanel(x,renderPlot({
      meth_age_plot(x)
        },height = 800))})
    
    do.call(tabsetPanel, myTabs)
  })
  
  output$mytabs2 = renderUI({
    nTabs = input$input_gene2
    
    myTabs = lapply(nTabs, function(x){
      tabPanel(x,renderPlot({
        transcript_age_plot(x)
      },height = 600))})
    
    do.call(tabsetPanel, myTabs)
  })
  
  
  output$chipseq_tab = renderUI({
    nTabs = input$input_gene2
    
    myTabs = lapply(nTabs, function(x){
      tabPanel(x,renderPlot({
        plot_static_browser_view(x)
      },height = 600))})
    
    do.call(tabsetPanel, myTabs)
  })

 output$biplot <- renderPlot({
      exp_line_plot(input$input_gene2)
    },height = 600)
    # ,height = 1000, width = 1000)
    
 output$biplot2 <- renderPlot({
   exp_line_plot_huvec10(input$input_gene2)
 },height = 600)
 
 
 output$ezh2_plot <- renderPlot({
   ezh2_inihibition_plot(input$input_gene2)
 },height = 600)
 output$exo_plot <- renderPlot({
   exo_inihibition_plot(input$input_gene2)
 },height = 600)
 output$rest_plot <- renderPlot({
   rest_oe_plot(input$input_gene2)
 },height = 600)
 
 output$kera_rapa_fors <- renderPlot({
   kera_rapa_fors_plot(input$input_gene2)
 },height = 600)
 
 output$mvec_plot <- renderPlot({
   mvec_plot(input$input_gene2)
 },height = 600)

    
    output$exp_mat <- downloadHandler(filename = function() { 
      paste("HUVEC_YMO_Expression-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      # test <- reshape2::dcast(df.m[,c('Var1','Var2',"samples","age",'counts')],Var1 ~ Var2+samples+age)
      # colnames(test) <- unlist(lapply(names(test),function(x) strsplit(x, split = "\\hisat2/")[[1]][2]))
      write.csv(g.count.matrix, file)
    })
    
    output$meth_mat <- downloadHandler(filename = function() { 
      paste("HUVEC_YMO_Methylation-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(dat0, file)
    })
    
    output$chipseq <- renderUI({
      HTML("<b> Experiment in progress. Have patience ^-^ </b>")
    })
    
    
    output$prot_plot = renderUI({
      nTabs = input$input_gene2
      
      myTabs = lapply(nTabs, function(x){
        tabPanel(x,renderPlot({
          prot_plot_fun(x)
        },height = 600))})
      
      do.call(tabsetPanel, myTabs)
    })
    
    # output$prot_plot <- renderPlot({
    #   prot_plot_fun(input$input_gene2)
    # },height = 600)
    # # ,height = 1000, width = 1000)
    # 
    
}

shinyApp(ui, server)
