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

# renv::snapshot()

########################################################### 
#Read the Meta data file
########################################################### 


load(file = "app_db/app_data.Rdata")

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
  
  proteo_query <- proteo[proteo$Genes %in% query_gene,] %>% as.data.frame()
  proteo_query <- reshape2::melt(t(proteo_query[,names(proteo_query)[-(61:65)]]),value.name = "proteomics_intensity")
  
  proteo_query$VG_LABEL <- plyr::mapvalues(proteo_query$Var1,p_meta$`MS Hub ID`,p_meta$VG_LABEL)
  proteo_query$proteomics_intensity <- as.numeric(proteo_query$proteomics_intensity )
  proteo_query <- proteo_query %>% mutate(proteomics_intensity = if_else(is.na(proteomics_intensity), 0, proteomics_intensity))
  proteo_query$`Tissue Type` <- plyr::mapvalues(proteo_query$Var1,p_meta$`MS Hub ID`,p_meta$`Tissue Type`)
  
  ggpubr::ggboxplot(proteo_query %>% as.data.frame(),x="VG_LABEL",y="proteomics_intensity",color = "Tissue Type",add = "jitter",title = query_gene,
                    palette="npg",ylab="log2(Protein abundance)", xlab="Samples")+
    theme_bw()+
    theme(text = element_text(size=12),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
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
                            tabPanel("ChIPseq",uiOutput('chipseq'))
                            
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

 output$biplot <- renderPlot({
      exp_line_plot(input$input_gene2)
    },height = 600)
    # ,height = 1000, width = 1000)
    
 output$biplot2 <- renderPlot({
   exp_line_plot_huvec10(input$input_gene2)
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
