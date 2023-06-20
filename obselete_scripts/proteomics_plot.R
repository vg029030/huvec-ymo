query_gene <- c("RBM3")
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
