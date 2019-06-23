load("results/preprocessing_2.RData")
load("results/annot.RData")
library(dplyr)
#library(biomaRt)
library(DESeq2)
library(tidyverse)
library(gplots)
library(RColorBrewer)

resCond_vV <- results(ddsObj, name="condition_Resistant_vs_Sensitive")
#load("results/resCond_vV.data")

resCond_vV  <- as.data.frame(resCond_vV) %>% 
    rownames_to_column("GeneID") %>% 
    left_join( annotL, "GeneID") 



#resCond_vV.0.05 <- resCond_vV.nona[(resCond_vV.nona$padj < 0.05),]
#resCond_vV.0.05 <- resCond_vV.0.05[!is.na(resCond_vV.0.05$GeneID),]
#dim(resCond_vV.0.05) 
#https://stephenturner.github.io/deseq-to-fgsea/
res2t <- resCond_vV

res2 <- res2t %>% 
  dplyr::select(Symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(Symbol) %>% 
  summarize(stat=mean(stat))
res2

library(fgsea)
ranks <- deframe(res2)
head(ranks, 20)

#pathways.hallmark <- gmtPathways("msigdb.v6.2.symbols.gmt")
pathways.hallmark <- gmtPathways("h.all.v6.2.symbols.gmt")

pathways.hallmark %>% 
  head() %>% 
  lapply(head)
  
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

###all pathways
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
  
 library("DT")
# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) #%>% 
#DT::datatable()


     
tiff("IBN000-gsea_hallmark.tif",width=1000,height=1000)   
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
theme_minimal()
dev.off()


     
pdf("IBN000-gsea_hallmark.pdf",width=10,height=10)   
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
theme_minimal()
dev.off()




tiff("IBN000-gsea_G2M.tif",width=1000,height=1000)  
plotEnrichment(pathways.hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks) +
labs(title="HALLMARK_G2M_CHECKPOINT")
dev.off()


tiff("IBN000-gsea_OXIDATIVE_PHOSPHORYLATION.tif",width=1000,height=1000)  
plotEnrichment(pathways.hallmark[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], ranks) +
labs(title="HALLMARK_OXIDATIVE_PHOSPHORYLATION")
dev.off()


tiff("IBN000-gsea_E2F.tif",width=1000,height=1000)  
plotEnrichment(pathways.hallmark[["HALLMARK_E2F_TARGETS"]], ranks) +
labs(title="HALLMARK_E2F_TARGETS")

dev.off()


tiff("IBN000-gsea_PI3K.tif",width=1000,height=1000)  
plotEnrichment(pathways.hallmark[["HALLMARK_PI3K_AKT_MTOR_SIGNALING"]], ranks) +
labs(title="HALLMARK_PI3K_AKT_MTOR_SIGNALING")

dev.off()


tiff("IBN000-gsea_G2M.tif",width=1000,height=1000)  
plotEnrichment(pathways.hallmark[["HALLMARK_G2M_CHECKPOINT"]], ranks) +
labs(title="HALLMARK_G2M_CHECKPOINT")
dev.off()



tiff("IBN000-gsea_HALLMARK_MYC_TARGETS_V2.tif",width=1000,height=1000)  
plotEnrichment(pathways.hallmark[["HALLMARK_MYC_TARGETS_V2"]], ranks) +
labs(title="HALLMARK_MYC_TARGETS_V2")
dev.off()




tiff("IBN000-gsea_HALLMARK_MTORC1_SIGNALING.tif",width=1000,height=1000)  
plotEnrichment(pathways.hallmark[["HALLMARK_MTORC1_SIGNALING"]], ranks) +
labs(title="HALLMARK_MTORC1_SIGNALING")
dev.off()



tiff("IBN000-gsea_HALLMARK_EMT.tif",width=1000,height=1000)  
plotEnrichment(pathways.hallmark[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]], ranks) +
labs(title="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
dev.off()


tiff("IBN000-gsea_HALLMARK_TNFA_SIGNALING_VIA_NFKB.tif",width=1000,height=1000)  
plotEnrichment(pathways.hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]], ranks) +
labs(title="HALLMARK_TNFA_SIGNALING_VIA_NFKB")
dev.off()


tiff("IBN000-gsea_HALLMARK_HEDGEHOG_SIGNALING.tif",width=1000,height=1000)  
plotEnrichment(pathways.hallmark[["HALLMARK_HEDGEHOG_SIGNALING"]], ranks) +
labs(title="HALLMARK_HEDGEHOG_SIGNALING")
dev.off()


[1] "HALLMARK_E2F_TARGETS"                      
 [2] "HALLMARK_G2M_CHECKPOINT"                   
 [3] "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
 [4] "HALLMARK_MYC_TARGETS_V2"                   
 [5] "HALLMARK_ANGIOGENESIS"                     
 [6] "HALLMARK_UV_RESPONSE_UP"                   
 [7] "HALLMARK_TNFA_SIGNALING_VIA_NFKB"          
 [8] "HALLMARK_IL6_JAK_STAT3_SIGNALING"          
 [9] "HALLMARK_ESTROGEN_RESPONSE_LATE"           
[10] "HALLMARK_SPERMATOGENESIS"                  
[11] "HALLMARK_COAGULATION"                      
[12] "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"        
[13] "HALLMARK_INFLAMMATORY_RESPONSE"            
[14] "HALLMARK_HEDGEHOG_SIGNALING"               
[15] "HALLMARK_MTORC1_SIGNALING"                 
[16] "HALLMARK_GLYCOLYSIS"                       
[17] "HALLMARK_ESTROGEN_RESPONSE_EARLY"          
[18] "HALLMARK_MITOTIC_SPINDLE"                  
[19] "HALLMARK_CHOLESTEROL_HOMEOSTASIS"          
[20] "HALLMARK_APOPTOSIS"                        
[21] "HALLMARK_MYC_TARGETS_V1"             




path <- fgseaRes[order(pval), ][2,]$pathway

tiff("IBN000-gsea_path.tif",width=1000,height=1000)  
plotEnrichment(pathways.hallmark[[path]],
               ranks) + labs(title=path)
dev.off()

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

tiff("IBN000-gsea_tab.tif",width=1000,height=1000) 
plotGseaTable(pathways.hallmark[topPathways], ranks, fgseaRes, 
              gseaParam = 0.5)
dev.off()









############################
pathways.kegg <- gmtPathways("c2.cp.kegg.v6.2.symbols.gmt")

pathways.kegg %>% 
  head() %>% 
  lapply(head)
  
fgseaRes <- fgsea(pathways=pathways.kegg, stats=ranks, nperm=1000)

###all pathways
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
  
 library("DT")
# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) #%>% 
  #DT::datatable()
            
     
tiff("IBN000-gsea_kegg.tif",width=1000,height=2000)   
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="kegg pathways NES from GSEA") + 
theme_minimal()
dev.off()


     
pdf("IBN000-gsea_kegg.pdf",width=10,height=20)   
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="kegg pathways NES from GSEA") + 
theme_minimal()
dev.off()

"KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY"


tiff("IBN000-gsea_KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY.tif",width=1000,height=1000)  
plotEnrichment(pathways.kegg[["KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY"]], ranks) +
labs(title="KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY")

dev.off()


