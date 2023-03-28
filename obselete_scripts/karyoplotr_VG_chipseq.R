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
  mutate(methylated=ifelse((logFC>0 & adj.P.Val < pvalue_cut),"HYPER",ifelse((logFC<0 & adj.P.Val < pvalue_cut),"HYPO","NOCHANGE"))) 
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


base.url="app_db/"
histone.marks <- c(
  OLD_REST="OLD_REST.bw",
  MIDDLE_REST="MEDIUM_REST.bw",
  YOUNG_REST="YOUNG_REST.bw"
  # ,YOUNG_input="3_0BGS_01FXPHE_FSK_Input_hg38_i49_uniqnorm_signal.bw"
)

DNA.binding <- c(
  OLD_EGR1="OLD_EGR1.bw",
  MIDDLE_EGR1="MEDIUM_EGR1.bw",
  YOUNG_EGR1="YOUNG_EGR1.bw"
)


save(list=c("epic_38_gr","hmm.model","promoters","histone.marks","DNA.binding","base.url"),file = "../app_db/browser_data.Rdata")



query_gene <- "LHFPL4"
plot_static_browser_view <- function(query_gene){
  
  pp <- getDefaultPlotParams(plot.type=1)
  pp$leftmargin <- 0.15
  pp$topmargin <- 15
  pp$bottommargin <- 15
  pp$ideogramheight <- 5
  pp$data1inmargin <- 10
  pp$data1outmargin <- 0
  
  myGeneSymbols <- select(org.Hs.eg.db,
                          keys = query_gene,
                          columns = c("SYMBOL","ENTREZID"),
                          keytype = "SYMBOL")
  
  myGeneSymbolsTx <- select(TxDb.Hsapiens.UCSC.hg38.knownGene,
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
  
  kpPlotGenes(kp, data=genes.data, r0=0, r1=0.05, gene.name.cex = 1)
  # kpPlotRegions(kp, epic_gr,col=colByChr(epic_gr$methylated,colors = cols), r0=0.2, r1=0.25,avoid.overlapping = F,clipping = F)
  
  kpPlotRegions(kp, data = epic_38_gr,col = epic_38_gr$jhs_color, r0=0.06, r1=0.10)
  kpAddLabels(kp, labels = "JHS-EPIC", r0=0.06, r1=0.08,cex=1)
  
  kpPlotRegions(kp, data = epic_38_gr,col = epic_38_gr$huvec_color, r0=0.11, r1=0.15)
  kpAddLabels(kp, labels = "HUVEC-Oldv/sYoung", r0=0.11, r1=0.13,cex=1)
  
  kpPlotRegions(kp, promoters, col="red", r0=0.16, r1=0.20)
  kpAddLabels(kp, labels = "Promoters", r0=0.16, r1=0.18,cex=1)
  
  kpPlotRegions(kp, hmm.model, col=hmm.model$RGB, r0=0.21, r1=0.25)
  kpAddLabels(kp, labels = "Chromatin\nState (HMM)", r0=0.21, r1=0.23,cex=1)
  
  
  
  computed.ymax.rest <-0
  for(i in names(c(histone.marks))){
    temp <- max(import(paste0(paste0(base.url, histone.marks[i])),which=region)$score)
    computed.ymax.rest <- ifelse(computed.ymax.rest > temp,computed.ymax.rest,temp)
    computed.ymax.rest <- round(computed.ymax.rest,0)
  }
  computed.ymax.egr1 <-0
  for(i in names(c(DNA.binding))){
    temp <- max(import(paste0(paste0(base.url, DNA.binding[i])),which=region)$score)
    computed.ymax.egr1 <- ifelse(computed.ymax.egr1 > temp,computed.ymax.egr1,temp)
    computed.ymax.egr1 <- round(computed.ymax.egr1,0)
  }
  
  
  #Histone marks
  total.tracks <- length(histone.marks)+length(DNA.binding)
  out.at <- autotrack(1:length(histone.marks), total.tracks, margin = 0.24, r0=0.35)
  
  for(i in seq_len(length(histone.marks))) {
    bigwig.file <- paste0(base.url, histone.marks[i])
    at <- autotrack(i, length(histone.marks), r0=out.at$r0, r1=out.at$r1, margin = 0.1)
    # kp <- kpPlotBigWig(kp, data=bigwig.file, ymax="visible.region",
    # kp <- kpPlotBigWig(kp, data=bigwig.file, ymax = "global",
    kp <- kpPlotBigWig(kp, data=bigwig.file, ymax=computed.ymax.rest,
                       r0=at$r0, r1=at$r1, col = "cadetblue2")
    # computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
    # computed.ymax <- ifelse(computed.ymax > ceiling(kp$latest.plot$computed.values$ymax),computed.ymax,ceiling(kp$latest.plot$computed.values$ymax))
    kpAxis(kp, ymin=0, ymax=computed.ymax.rest, numticks = 2, r0=at$r0, r1=at$r1,cex=0.5)
    # kpAxis(kp, ymin=0, ymax=31, numticks = 2, r0=at$r0, r1=at$r1,cex=0.5)
    kpAddLabels(kp, labels = names(histone.marks)[i], r0=at$r0, r1=at$r1, 
                cex=1, label.margin = 0.035)
  }
  
  #DNA binding proteins
  out.at <- autotrack((length(histone.marks)+1):(total.tracks), total.tracks, margin = 0.36, r0=0.45)
  
  for(i in seq_len(length(DNA.binding))) {
    bigwig.file <- paste0(base.url, DNA.binding[i])
    at <- autotrack(i, length(DNA.binding), r0=out.at$r0, r1=out.at$r1, margin = 0.1)
    # kp <- kpPlotBigWig(kp, data=bigwig.file, ymax="visible.region",
    # kp <- kpPlotBigWig(kp, data=bigwig.file, ymax = "global",
    kp <- kpPlotBigWig(kp, data=bigwig.file, ymax=computed.ymax.egr1,
                       r0=at$r0, r1=at$r1, col = "darkolivegreen1")
    # computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
    # computed.ymax <- ifelse(computed.ymax > ceiling(kp$latest.plot$computed.values$ymax),computed.ymax,ceiling(kp$latest.plot$computed.values$ymax))
    kpAxis(kp, ymin=0, ymax=computed.ymax.egr1, numticks = 2, r0=at$r0, r1=at$r1,cex=0.5)
    # kpAxis(kp, ymin=0, ymax=31, numticks = 2, r0=at$r0, r1=at$r1,cex=0.5)
    kpAddLabels(kp, labels = names(DNA.binding)[i], r0=at$r0, r1=at$r1, 
                cex=1, label.margin = 0.035)
  }
}

