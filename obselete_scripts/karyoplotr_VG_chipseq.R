rm(list = ls())
library(AnnotationHub)
library(karyoploteR)
library(BiocFileCache)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(plyr)
library(ggplot2)
library(org.Hs.eg.db)
library(tibble)
library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)

# pvalue_cut <- 0.05/866836
pvalue_cut <- 1.0e-07
epic <- data.table::fread("../../../../EPIC/EPICAnnotation_age_DNAmTL.txt") %>%
  column_to_rownames("CGid") %>%
  mutate(age_associated=(P.age.JHS < pvalue_cut)) %>%
  mutate(methylated=ifelse((Z.age.JHS>0 & P.age.JHS < pvalue_cut),"HYPER",ifelse((Z.age.JHS<0 & P.age.JHS < pvalue_cut),"HYPO","NOCHANGE"))) 

cols <- c(NOCHANGE="lightgrey", HYPER="red","HYPO"="blue")
epic$jhs_color <- plyr::mapvalues(epic$methylated,names(cols),cols)
epic_gr <- makeGRangesFromDataFrame(epic,seqnames.field = "seqnames",start.field ="start",strand.field = "strand",keep.extra.columns = T)


pvalue_cut <- 1.0e-03
epic_huvec_DE <- data.table::fread("../../../METHYLATION/DE_meth_old_young.tsv")%>%
  column_to_rownames("V1") %>%
  mutate(age_associated=(adj.P.Val < pvalue_cut)) %>%
  mutate(methylated=ifelse((logFC>0.2 & adj.P.Val < pvalue_cut),"HYPER",ifelse((logFC<-0.2 & adj.P.Val < pvalue_cut),"HYPO","NOCHANGE"))) 
cols <- c(NOCHANGE="lightgrey", HYPER="red","HYPO"="blue")
epic_huvec_DE$huvec_color <- plyr::mapvalues(epic_huvec_DE$methylated,names(cols),cols)


epic_38 <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38) %>% as.data.frame()
epic_38 <- epic_38[,c("CHR_hg38","Start_hg38","End_hg38","Strand_hg38")]
epic_38 <- merge(epic_38,epic[,"jhs_color",drop=F],by=0) %>% column_to_rownames("Row.names")
epic_38 <- merge(epic_38,epic_huvec_DE[,"huvec_color",drop=F],by=0) %>% column_to_rownames("Row.names")
epic_38 <- epic_38[complete.cases(epic_38),]
epic_38_gr <- makeGRangesFromDataFrame(epic_38,seqnames.field = "CHR_hg38",start.field ="Start_hg38",end.field ="End_hg38",strand.field = "Strand_hg38",keep.extra.columns = T)



  
# base.url <- "~/Google Drive/Shared drives/Raj Lab/VipulGupta/HUVEC_YMO/CHIPseq/EGR1andREST_repeatfastqs/replicate_merged_bigwig/"
# setwd(base.url)



# download_hmm <- "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E005_15_coreMarks_hg38lift_dense.bed.gz"
hmm.file <- data.table::fread("../../../../E005_15_coreMarks_hg38lift_dense.bed.gz",skip = 1)
hmm.model <- makeGRangesFromDataFrame(hmm.file,seqnames.field = "V1",start.field ="V2",end.field ="V3",keep.extra.columns = T )
hmm.model$RGB <- sapply(strsplit(hmm.model$V9, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue=255))




promoters <- promoters(genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
# , upstream = 1000, downstream = 1000)


# base.url="app_db/"
# histone.marks <- c(
#   OLD_REST="OLD_REST.bw",
#   MIDDLE_REST="MEDIUM_REST.bw",
#   YOUNG_REST="YOUNG_REST.bw"
#   # ,YOUNG_input="3_0BGS_01FXPHE_FSK_Input_hg38_i49_uniqnorm_signal.bw"
# )
# 
# DNA.binding <- c(
#   OLD_EGR1="OLD_EGR1.bw",
#   MIDDLE_EGR1="MEDIUM_EGR1.bw",
#   YOUNG_EGR1="YOUNG_EGR1.bw"
# )


save(list=c("epic_38_gr","hmm.model","promoters"),file = "../app_db/browser_data.Rdata")



query_gene <- "LHFPL4"
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
  
  # s3_folder <- "/share/vgupta/workspaces/huvec-analysis/HUVEC_YMP_CHIPseq_bw/"
  s3_folder <- "~/Google Drive/Shared drives/Raj Lab/VipulGupta/HUVEC_YMO/CHIPseq/Diffbind_all_samples/bw/"
  
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
    computed.ymax.k27me3 <- ifelse(round(computed.ymax.k27me3,0) > 1,round(computed.ymax.k27me3,0),1)
  }
  computed.ymax.rest <-0
  for(i in names(c(rest_files))){
    temp <- max(import(paste0(paste0(s3_folder,rest_files[i])),which=region)$score)
    computed.ymax.rest <- ifelse(computed.ymax.rest > temp,computed.ymax.rest,temp)
    computed.ymax.rest <- ifelse(round(computed.ymax.rest,0) > 1,round(computed.ymax.rest,0),1)
  }
  computed.ymax.egr1 <-0
  for(i in names(c(egr1_files))){
    temp <- max(import(paste0(paste0(s3_folder, egr1_files[i])),which=region)$score)
    computed.ymax.egr1 <- ifelse(computed.ymax.egr1 > temp,computed.ymax.egr1,temp)
    computed.ymax.egr1 <- ifelse(round(computed.ymax.egr1,0) > 1,round(computed.ymax.egr1,0),1)
  }
  computed.ymax.k4me3 <-0
  for(i in names(k4me3_files)){
    temp <- max(import(paste0(s3_folder, k4me3_files[i]),which=region)$score)
    computed.ymax.k4me3 <- ifelse(computed.ymax.k4me3 > temp,computed.ymax.k4me3,temp)
    computed.ymax.k4me3 <- ifelse(round(computed.ymax.k4me3,0) > 1,round(computed.ymax.k4me3,0),1)
  }
  
  computed.ymax.ezh2 <-0
  for(i in names(c(ezh2_files))){
    temp <- max(import(paste0(paste0(s3_folder, ezh2_files[i])),which=region)$score)
    computed.ymax.ezh2 <- ifelse(computed.ymax.ezh2 > temp,computed.ymax.ezh2,temp)
    computed.ymax.ezh2 <- ifelse(round(computed.ymax.ezh2,0) > 1,round(computed.ymax.ezh2,0),1)
  }
  computed.ymax.jarid2 <-0
  for(i in names(c(jarid2_files))){
    temp <- max(import(paste0(paste0(s3_folder, jarid2_files[i])),which=region)$score)
    computed.ymax.jarid2 <- ifelse(computed.ymax.jarid2 > temp,computed.ymax.jarid2,temp)
    computed.ymax.jarid2 <- ifelse(round(computed.ymax.jarid2,0) > 1,round(computed.ymax.jarid2,0),1)
  }
  
  total.tracks <- length(k27me3_files)+length(k4me3_files)+length(ezh2_files)+length(jarid2_files)+length(rest_files)+length(egr1_files)
  
  #K27me3 Track
  out.at <- autotrack(1:3, total.tracks,  r0=0.17,r1=1,margin = 0.1)
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
  out.at <- autotrack(4:6, total.tracks,  r0=0.17,r1=1,margin = 0.1)
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
  out.at <- autotrack(7:9, total.tracks, r0=0.17,r1=1,margin = 0.1)
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
  out.at <- autotrack(10:12, total.tracks, r0=0.17,r1=1,margin = 0.1)
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
  out.at <- autotrack(13:15, total.tracks,  r0=0.17,r1=1,margin = 0.1)
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
  out.at <- autotrack(16:18, total.tracks, r0=0.17,r1=1,margin = 0.1)
  
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

