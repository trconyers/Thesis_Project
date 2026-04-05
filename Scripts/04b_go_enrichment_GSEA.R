########################### Selection_ORA_Transcriptomics ###########################
library(GOstats)
library(org.Dm.eg.db)
library(clusterProfiler)
library(Category)
library(enrichplot)

compareCluster_generic <- compareCluster(geneCluster = as.list.Bimap(org.Dm.egGO2EG)[1:3], fun = enrichGO, OrgDb = org.Dm.eg.db)
dotplot_generic <- dotplot(compareCluster_generic)
set_theme(new = dotplot_generic@theme); rm(dotplot_generic, compareCluster_generic)
new_theme <- get_theme()
new_theme$text@family <- "serif"
new_theme$text@size <- 13.5
new_theme$axis.text.x@size <- 13.5
new_theme$axis.text.y@size <- 13.5
new_theme$text@lineheight <- 0.8
new_theme$plot.title@hjust <- 0.5
set_theme(new = new_theme); rm(new_theme)

Immun_genes.T <- Immundat.T$Elements
Immun_genes.T <- SYMBOL2Gene(unique(unlist(strsplit(x = str_flatten_comma(Immun_genes.T), split = ", "))))
Immun_lds.T <- c(is.longevity(Immun_genes.T), is.immunity(Immun_genes.T), is.stress(Immun_genes.T))
Immun_Clusts.T <- list(Lifespan = is.longevity(Immun_genes.T), Defense = is.immunity(Immun_genes.T), Stress_Response = is.stress(Immun_genes.T), Other = setdiff(Immun_genes.T, Immun_lds.T))
Immun_enrich.T <- enrichGO(gene = Immun_genes.T, OrgDb = org.Dm.eg.db)
Immun_enrich.T_comp <- compareCluster(
  geneCluster = Immun_Clusts.T,
  fun = enrichGO,
  OrgDb = org.Dm.eg.db,
  ont = "ALL",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)
ImmunityBP_ora.T <- hyperGTest(new(Class = "GOHyperGParams",
                                   ontology = "BP",
                                   geneIds = Immun_enrich.T@gene,
                                   universeGeneIds = Immun_enrich.T@universe,
                                   annotation = "org.Dm.eg.db",
                                   pvalueCutoff = 0.05,
                                   testDirection = "over",
                                   conditional = TRUE))
ImmunityBP_ora.T_summ <- summary(ImmunityBP_ora.T)
ImmunityCC_ora.T <- hyperGTest(new(Class = "GOHyperGParams",
                                   ontology = "CC",
                                   geneIds = Immun_enrich.T@gene,
                                   universeGeneIds = Immun_enrich.T@universe,
                                   annotation = "org.Dm.eg.db",
                                   pvalueCutoff = 0.05,
                                   testDirection = "over",
                                   conditional = TRUE))
ImmunityCC_ora.T_summ <- summary(ImmunityCC_ora.T)
ImmunityMF_ora.T <- hyperGTest(new(Class = "GOHyperGParams",
                                   ontology = "MF",
                                   geneIds = Immun_enrich.T@gene,
                                   universeGeneIds = Immun_enrich.T@universe,
                                   annotation = "org.Dm.eg.db",
                                   pvalueCutoff = 0.05,
                                   testDirection = "over",
                                   conditional = TRUE))
ImmunityMF_ora.T_summ <- summary(ImmunityMF_ora.T)
colnames(ImmunityBP_ora.T_summ)[1] <- colnames(ImmunityCC_ora.T_summ)[1] <- colnames(ImmunityMF_ora.T_summ)[1] <- "GOID"
Immunity_ora.T_summ <- rbind.data.frame(ImmunityBP_ora.T_summ, ImmunityCC_ora.T_summ, ImmunityMF_ora.T_summ) %>%  mutate(P.adj = p.adjust(Pvalue, method = "BH")) %>%  filter(P.adj < 0.05)
rm(ImmunityBP_ora.T_summ, ImmunityCC_ora.T_summ, ImmunityMF_ora.T_summ)
rownames(Immunity_ora.T_summ) <- Immunity_ora.T_summ$GOID
Immun_enrich.T_comp@compareClusterResult <- Immun_enrich.T_comp@compareClusterResult[Immun_enrich.T_comp@compareClusterResult$ID %fin% Immunity_ora.T_summ$GOID,]
Immun_enrich.T_comp@compareClusterResult$p.adjust <- Immunity_ora.T_summ[Immun_enrich.T_comp@compareClusterResult$ID, "P.adj"]
Immun_enrich.T_comp@compareClusterResult$EnrichRatio <- (Immun_enrich.T_comp@compareClusterResult$Count)/(as.numeric(gsub("/.*$", "", Immun_enrich.T_comp@compareClusterResult$BgRatio)))
Immun_dot.T <- clusterProfiler::simplify(Immun_enrich.T_comp)
Immun_dotplot.T <- dotplot(Immun_dot.T, label_format = 32, by = "EnrichRatio", showCategory = 9) + get_theme() + labs(x = "Function", title = "Overrepresented GO Categories in Immunity-Selection Overlaps (Transcriptomics)")
Immun_dotplot.T
savePlotAsImage(
  file = "Graphs/Immun_GOs_Transcriptom.png",
  format = "png",
  width = 1200,
  height = 742
)

Stress_genes.T <- Stressdat.T$Elements[Stressdat.T$Degree>=Maxs(Stressdat.T$Degree) - 1]
Stress_genes.T <- unname(SYMBOL2Gene(unique(unlist(strsplit(x = str_flatten_comma(Stress_genes.T), split = ", ")))))
Stress_lds.T <- c(is.longevity(Stress_genes.T), is.immunity(Stress_genes.T), is.stress(Stress_genes.T))
Stress_Clusts.T <- list(Lifespan = is.longevity(Stress_genes.T), Defense = is.immunity(Stress_genes.T), Stress_Response = is.stress(Stress_genes.T), Other = setdiff(Stress_genes.T, Stress_lds.T))
Stress_enrich.T <- enrichGO(gene = Stress_genes.T, OrgDb = org.Dm.eg.db)
Stress_enrich.T_comp <- compareCluster(
  geneCluster = Stress_Clusts.T,
  fun = enrichGO,
  OrgDb = org.Dm.eg.db,
  ont = "ALL",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)
StressBP_ora.T <- hyperGTest(new(Class = "GOHyperGParams",
                                 ontology = "BP",
                                 geneIds = Stress_enrich.T@gene,
                                 universeGeneIds = Stress_enrich.T@universe,
                                 annotation = "org.Dm.eg.db",
                                 pvalueCutoff = 0.05,
                                 testDirection = "over",
                                 conditional = TRUE))
StressBP_ora.T_summ <- summary(StressBP_ora.T)
StressCC_ora.T <- hyperGTest(new(Class = "GOHyperGParams",
                                 ontology = "CC",
                                 geneIds = Stress_enrich.T@gene,
                                 universeGeneIds = Stress_enrich.T@universe,
                                 annotation = "org.Dm.eg.db",
                                 pvalueCutoff = 0.05,
                                 testDirection = "over",
                                 conditional = TRUE))
StressCC_ora.T_summ <- summary(StressCC_ora.T)
StressMF_ora.T <- hyperGTest(new(Class = "GOHyperGParams",
                                 ontology = "MF",
                                 geneIds = Stress_enrich.T@gene,
                                 universeGeneIds = Stress_enrich.T@universe,
                                 annotation = "org.Dm.eg.db",
                                 pvalueCutoff = 0.05,
                                 testDirection = "over",
                                 conditional = TRUE))
StressMF_ora.T_summ <- summary(StressMF_ora.T)
colnames(StressBP_ora.T_summ)[1] <- colnames(StressCC_ora.T_summ)[1] <- colnames(StressMF_ora.T_summ)[1] <- "GOID"
Stress_ora.T_summ <- rbind.data.frame(StressBP_ora.T_summ, StressCC_ora.T_summ, StressMF_ora.T_summ) %>%  mutate(P.adj = p.adjust(Pvalue, method = "BH")) %>%  filter(P.adj < 0.05)
rm(StressBP_ora.T_summ, StressCC_ora.T_summ, StressMF_ora.T_summ)
rownames(Stress_ora.T_summ) <- Stress_ora.T_summ$GOID
Stress_enrich.T_comp@compareClusterResult <- Stress_enrich.T_comp@compareClusterResult[Stress_enrich.T_comp@compareClusterResult$ID %fin% Stress_ora.T_summ$GOID,]
Stress_enrich.T_comp@compareClusterResult$p.adjust <- Stress_ora.T_summ[Stress_enrich.T_comp@compareClusterResult$ID, "P.adj"]
Stress_enrich.T_comp@compareClusterResult$EnrichRatio <- (Stress_enrich.T_comp@compareClusterResult$Count)/(as.numeric(gsub("/.*$", "", Stress_enrich.T_comp@compareClusterResult$BgRatio)))
Stress_dot.T <- clusterProfiler::simplify(Stress_enrich.T_comp)
Stress_dotplot.T <- dotplot(Stress_dot.T, label_format = 32, by = "EnrichRatio", showCategory = 4) + get_theme() + labs(x = "Function", title = "Overrepresented GO Categories in Stress-Selection Overlaps (Transcriptomics)")
Stress_dotplot.T
savePlotAsImage(
  file = "Graphs/Stress_GOs_Transcriptom.png",
  format = "png",
  width = 1200,
  height = 742
)

######################## Longev_Immun/Stress_ORA_Transcriptomics #########################
LI_genes.T <- unname(SYMBOL2Gene(unique(unlist(strsplit(x = Selectiondat.T$Elements[Selectiondat.T$Intersections=="Longevity & Immunity"], split = ", ")))))
LI_lds.T <- c(is.longevity(LI_genes.T), is.immunity(LI_genes.T), is.stress(LI_genes.T))
LI_Clusts.T <- list(Lifespan = is.longevity(LI_genes.T), Defense = is.immunity(LI_genes.T), Stress_Response = is.stress(LI_genes.T), Other = setdiff(LI_genes.T, LI_lds.T))
LI_enrich.T <- enrichGO(gene = LI_genes.T, OrgDb = org.Dm.eg.db)
LI_enrich.T_comp <- compareCluster(
  geneCluster = LI_Clusts.T,
  fun = enrichGO,
  OrgDb = org.Dm.eg.db,
  ont = "ALL",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)
LI_BP_ora.T <- hyperGTest(new(Class = "GOHyperGParams",
                              ontology = "BP",
                              geneIds = LI_enrich.T@gene,
                              universeGeneIds = LI_enrich.T@universe,
                              annotation = "org.Dm.eg.db",
                              pvalueCutoff = 0.05,
                              testDirection = "over",
                              conditional = TRUE))
LI_BP_ora.T_summ <- summary(LI_BP_ora.T)
LI_CC_ora.T <- hyperGTest(new(Class = "GOHyperGParams",
                              ontology = "CC",
                              geneIds = LI_enrich.T@gene,
                              universeGeneIds = LI_enrich.T@universe,
                              annotation = "org.Dm.eg.db",
                              pvalueCutoff = 0.05,
                              testDirection = "over",
                              conditional = TRUE))
LI_CC_ora.T_summ <- summary(LI_CC_ora.T)
LI_MF_ora.T <- hyperGTest(new(Class = "GOHyperGParams",
                              ontology = "MF",
                              geneIds = LI_enrich.T@gene,
                              universeGeneIds = LI_enrich.T@universe,
                              annotation = "org.Dm.eg.db",
                              pvalueCutoff = 0.05,
                              testDirection = "over",
                              conditional = TRUE))
LI_MF_ora.T_summ <- summary(LI_MF_ora.T)
colnames(LI_BP_ora.T_summ)[1] <- colnames(LI_CC_ora.T_summ)[1] <- colnames(LI_MF_ora.T_summ)[1] <- "GOID"
LI_ora.T_summ <- rbind.data.frame(LI_BP_ora.T_summ, LI_CC_ora.T_summ, LI_MF_ora.T_summ) %>%  mutate(P.adj = p.adjust(Pvalue, method = "BH")) %>%  filter(P.adj < 0.05)
rm(LI_BP_ora.T_summ, LI_CC_ora.T_summ, LI_MF_ora.T_summ)
rownames(LI_ora.T_summ) <- LI_ora.T_summ$GOID
LI_enrich.T_comp@compareClusterResult <- LI_enrich.T_comp@compareClusterResult[LI_enrich.T_comp@compareClusterResult$ID %fin% LI_ora.T_summ$GOID,]
LI_enrich.T_comp@compareClusterResult$p.adjust <- LI_ora.T_summ[LI_enrich.T_comp@compareClusterResult$ID, "P.adj"]
LI_enrich.T_comp@compareClusterResult$EnrichRatio <- (LI_enrich.T_comp@compareClusterResult$Count)/(as.numeric(gsub("/.*$", "", LI_enrich.T_comp@compareClusterResult$BgRatio)))
LI_dot.T <- clusterProfiler::simplify(LI_enrich.T_comp)
LI_dotplot.T <- dotplot(LI_dot.T, label_format = 25, showCategory = 4, by = "EnrichRatio") + get_theme() + labs(x = "Function", title = "Overrepresented GO Categories in Longevity/Immunity-Selection Overlap (Transcriptomics)")
LI_dotplot.T
savePlotAsImage(
  file = "Graphs/Longev_Immun_GOs_Transcriptom.png",
  format = "png",
  width = 1200,
  height = 742
)

LS_genes.T <- unname(SYMBOL2Gene(unique(unlist(strsplit(x = Selectiondat.T$Elements[Selectiondat.T$Intersections=="Longevity & Stress"], split = ", ")))))
LS_lds.T <- c(is.longevity(LS_genes.T), is.immunity(LS_genes.T), is.stress(LS_genes.T))
LS_Clusts.T <- list(Lifespan = is.longevity(LS_genes.T), Defense = is.immunity(LS_genes.T), Stress_Response = is.stress(LS_genes.T), Other = setdiff(LS_genes.T, LS_lds.T))
LS_enrich.T <- enrichGO(gene = LS_genes.T, OrgDb = org.Dm.eg.db)
LS_enrich.T_comp <- compareCluster(
  geneCluster = LS_Clusts.T,
  fun = enrichGO,
  OrgDb = org.Dm.eg.db,
  ont = "ALL",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)
LS_BP_ora.T <- hyperGTest(new(Class = "GOHyperGParams",
                              ontology = "BP",
                              geneIds = LS_enrich.T@gene,
                              universeGeneIds = LS_enrich.T@universe,
                              annotation = "org.Dm.eg.db",
                              pvalueCutoff = 0.05,
                              testDirection = "over",
                              conditional = TRUE))
LS_BP_ora.T_summ <- summary(LS_BP_ora.T)
LS_CC_ora.T <- hyperGTest(new(Class = "GOHyperGParams",
                              ontology = "CC",
                              geneIds = LS_enrich.T@gene,
                              universeGeneIds = LS_enrich.T@universe,
                              annotation = "org.Dm.eg.db",
                              pvalueCutoff = 0.05,
                              testDirection = "over",
                              conditional = TRUE))
LS_CC_ora.T_summ <- summary(LS_CC_ora.T)
LS_MF_ora.T <- hyperGTest(new(Class = "GOHyperGParams",
                              ontology = "MF",
                              geneIds = LS_enrich.T@gene,
                              universeGeneIds = LS_enrich.T@universe,
                              annotation = "org.Dm.eg.db",
                              pvalueCutoff = 0.05,
                              testDirection = "over",
                              conditional = TRUE))
LS_MF_ora.T_summ <- summary(LS_MF_ora.T)
colnames(LS_BP_ora.T_summ)[1] <- colnames(LS_CC_ora.T_summ)[1] <- colnames(LS_MF_ora.T_summ)[1] <- "GOID"
LS_ora.T_summ <- rbind.data.frame(LS_BP_ora.T_summ, LS_CC_ora.T_summ, LS_MF_ora.T_summ) %>%  mutate(P.adj = p.adjust(Pvalue, method = "BH")) %>%  filter(P.adj < 0.05)
rm(LS_BP_ora.T_summ, LS_CC_ora.T_summ, LS_MF_ora.T_summ)
rownames(LS_ora.T_summ) <- LS_ora.T_summ$GOID
LS_enrich.T_comp@compareClusterResult <- LS_enrich.T_comp@compareClusterResult[LS_enrich.T_comp@compareClusterResult$ID %fin% LS_ora.T_summ$GOID,]
LS_enrich.T_comp@compareClusterResult$p.adjust <- LS_ora.T_summ[LS_enrich.T_comp@compareClusterResult$ID, "P.adj"]
LS_enrich.T_comp@compareClusterResult$EnrichRatio <- (LS_enrich.T_comp@compareClusterResult$Count)/(as.numeric(gsub("/.*$", "", LS_enrich.T_comp@compareClusterResult$BgRatio)))
LS_dot.T <- clusterProfiler::simplify(LS_enrich.T_comp)
LS_dotplot.T <- dotplot(LS_dot.T, label_format = 25, by = "EnrichRatio", showCategory = 3) + get_theme() + labs(x = "Function", title = "Overrepresented GO Categories in Longevity/Stress-Selection Overlap (Transcriptomics)")
LS_dotplot.T
savePlotAsImage(
  file = "Graphs/Longev_Stress_GOs_Transcriptom.png",
  format = "png",
  width = 1200,
  height = 742
)

dotplot_list.T <- cowplot::plot_grid(LI_dotplot.T, LS_dotplot.T, labels = c("A", "B"), ncol = 2)
dotplot_list.T
savePlotAsImage(
  file = "Graphs/Longev_Overlaps_GOs_Transcriptom.png",
  format = "png",
  width = 1800,
  height = 742
)

########################### Gene_Table edit for GSEA ###########################
library(readxl)
library(writexl)

Gene_Table <- as.data.frame(read_excel("Data/Gene_Table.xlsx"))
Gene_Table.T_GSEA <- Gene_Table[Gene_Table$Sequencing=="Transcriptomic",]
FB_invs <- as.data.frame(read_tsv("Scripts/misc/withdrawns.txt"))
Gene_Table.T_GSEA <-Gene_Table.T_GSEA[!Gene_Table.T_GSEA$`Candidate FB IDs` %fin% FB_invs$FBIDs,]
rm(FB_invs)

Multi_Pops.T_GSEA <- Gene_Table.T_GSEA[, 1:3]
rownames(Multi_Pops.T_GSEA) <- rownames(Gene_Table.T_GSEA)
Multi_Pops.T_GSEA <- unique(Multi_Pops.T_GSEA)
Gene_Table.T_GSEA <- Gene_Table.T_GSEA[rownames(Multi_Pops.T_GSEA), ]
rm(Multi_Pops.T_GSEA)

Gene_Table_EG.T_GSEA <- Gene_Table.T_GSEA
Gene_Table_EG.T_GSEA$`Candidate FB IDs` <- FB2EG(Gene_Table.T_GSEA$`Candidate FB IDs`)
Gene_Table_EG.T_GSEA$`Direction (Transcriptomics)` <- factor(Gene_Table_EG.T_GSEA$`Direction (Transcriptomics)`)
levels(Gene_Table_EG.T_GSEA$`Direction (Transcriptomics)`) <- c(-1,1)
Gene_Table_EG.T_GSEA$`Direction (Transcriptomics)` <- as.list(as.integer(unfactor(Gene_Table_EG.T_GSEA$`Direction (Transcriptomics)`)))
names(Gene_Table_EG.T_GSEA$`Direction (Transcriptomics)`) <- Gene_Table_EG.T_GSEA$`Candidate FB IDs`

all_pops.T_GSEA <- Unique(Gene_Table.T_GSEA$`Populations Used`)
Gene_Table_List.T_GSEA <- list()
for (i in cli_progress_along(seq(all_pops.T_GSEA))) {
  Gene_Table_List.T_GSEA[[i]] <- unlist2(Gene_Table_EG.T_GSEA$`Direction (Transcriptomics)`[Gene_Table_EG.T_GSEA$`Populations Used` == all_pops.T_GSEA[i]])
}
names(Gene_Table_List.T_GSEA) <- all_pops.T_GSEA
Gene_Table_List.T_GSEA <- mixedSorts(Gene_Table_List.T_GSEA)

####################### GSEA_SelectionType_Lists #######################
Longevity <- names(Gene_Table_EG.T_GSEA$`Direction (Transcriptomics)`)[Gene_Table_EG.T_GSEA$`Selection Type` == "Longevity"]
Longevity.dups <- Gene_Table_EG.T_GSEA$`Direction (Transcriptomics)`[Gene_Table_EG.T_GSEA$`Selection Type` == "Longevity"][Longevity %fin% Longevity[duplicated(Longevity)]]
Longevity.dups <- list_flip(list_flip(Longevity.dups))
Longevity <- setdiff(Longevity, names(Longevity.dups))
Longevity.dups <- map(.x = Longevity.dups, .f = as.integer)
Longevity.dups <- map(.x = Longevity.dups, .f = sum)
Longev_Table.T_GSEA <- sort(unlist(c(Gene_Table_EG.T_GSEA$`Direction (Transcriptomics)`[Gene_Table_EG.T_GSEA$`Selection Type` == "Longevity"][Longevity], Longevity.dups)), decreasing = TRUE)
rm(Longevity, Longevity.dups)

Immunity <- names(Gene_Table_EG.T_GSEA$`Direction (Transcriptomics)`)[Gene_Table_EG.T_GSEA$`Selection Type` == "Immunity"]
Immunity.dups <- Gene_Table_EG.T_GSEA$`Direction (Transcriptomics)`[Gene_Table_EG.T_GSEA$`Selection Type` == "Immunity"][Immunity %fin% Immunity[duplicated(Immunity)]]
Immunity.dups <- list_flip(list_flip(Immunity.dups))
Immunity <- setdiff(Immunity, names(Immunity.dups))
Immunity.dups <- map(.x = Immunity.dups, .f = as.integer)
Immunity.dups <- map(.x = Immunity.dups, .f = sum)
Immun_Table.T_GSEA <- sort(unlist(c(Gene_Table_EG.T_GSEA$`Direction (Transcriptomics)`[Gene_Table_EG.T_GSEA$`Selection Type` == "Immunity"][Immunity], Immunity.dups)), decreasing = TRUE)
rm(Immunity, Immunity.dups)

Stress <- names(Gene_Table_EG.T_GSEA$`Direction (Transcriptomics)`)[Gene_Table_EG.T_GSEA$`Selection Type` == "Stress"]
Stress.dups <- Gene_Table_EG.T_GSEA$`Direction (Transcriptomics)`[Gene_Table_EG.T_GSEA$`Selection Type` == "Stress"][Stress %fin% Stress[duplicated(Stress)]]
Stress.dups <- list_flip(list_flip(Stress.dups))
Stress <- setdiff(Stress, names(Stress.dups))
Stress.dups <- map(.x = Stress.dups, .f = as.integer)
Stress.dups <- map(.x = Stress.dups, .f = sum)
Stress_Table.T_GSEA <- sort(unlist(c(Gene_Table_EG.T_GSEA$`Direction (Transcriptomics)`[Gene_Table_EG.T_GSEA$`Selection Type` == "Stress"][Stress], Stress.dups)), decreasing = TRUE)
rm(Stress, Stress.dups)

Gene_Lists.T_GSEA <- list(names(Longev_Table.T_GSEA), names(Immun_Table.T_GSEA), names(Stress_Table.T_GSEA))
names(Gene_Lists.T_GSEA) <- c("Longevity", "Immunity", "Stress")

########################### Selection_GSEA ###########################
library(GOstats)
library(Category)
library(clusterProfiler)
library(org.Dm.eg.db)
library(enrichplot)
library(GOSemSim)

heatplot_generic <- heatplot(new("gseaResult"))
set_theme(new = heatplot_generic@theme); rm(heatplot_generic)
new_theme <- get_theme()
new_theme$text@family <- "serif"
new_theme$text@size <- 13.5
new_theme$plot.title@hjust <- 0.5
new_theme$axis.text.x@angle <- 70
new_theme$axis.text.y@size <- 12.75
new_theme$axis.text.y@lineheight <- 0.85
set_theme(new = new_theme); rm(new_theme)

Longev_GSEA <- gseGO(gene = Longev_Table.T_GSEA, OrgDb = org.Dm.eg.db, ont = "ALL", minGSSize = 15, maxGSSize = Inf, eps = 0)
Longev_GSEA@result <- Longev_GSEA@result[!Longev_GSEA@result$ID %fin% BPCCMF,]
Longev_GSEA <- clusterProfiler::simplify(Longev_GSEA)
Longev_GSEA <- setReadable(Longev_GSEA, "org.Dm.eg.db", "ENTREZID")
Longev_heat <- heatplot(Longev_GSEA, foldChange = Longev_GSEA@geneList, showCategory = 15, label_format = 32)
Longev_heat <- Longev_heat@plot_env$p + Longev_heat@scales$scales[[2]] + labs(y = "Category", title = "GSEA of Candidate Genes From Longevity-Selection Studies (Transcriptomics)")
d <- godata("org.Dm.eg.db", keytype = "SYMBOL", ont = "BP", computeIC = FALSE)
geneSims <- mgeneSim(unique(Longev_heat@data$Gene), semData = d, measure = "Wang")
dissimilarity <- as.dist(1 - geneSims)
hcluster <- hclust(dissimilarity)
dend <- as.dendrogram(hcluster, hang = 0)
cut_hclust <- dendextend::cutree(dend, h = dendextend::heights_per_k.dendrogram(dend)["99"])
sim_clust <- split(x = names(cut_hclust), f = cut_hclust)
sim_clust <- map(.x = sim_clust, .f = nameVector)
sim_clust <- mixedSorts(map(.x = sim_clust, .f = function(x) str_remove(string = x, pattern = "CG")))
sim_clust <- map(.x = sim_clust, .f = rev)
sim_clust <- map(.x = sim_clust, .f = names)
sim_clust <- map(.x = sim_clust, .f = function(x) tcount(Longev_heat@data$Gene)[x])
sim_clust <- map(.x = sim_clust, .f = function(x) x[which.max(x)])
keep <- unlist(map(.x = sim_clust, .f = names))
Longev_heat@data <- Longev_heat@data[Longev_heat@data$Gene %fin% keep,]
Longev_heat
savePlotAsImage(
  file = "Graphs/Longev_GSEA_Transcriptom.png",
  format = "png",
  width = 1312,
  height = 742
)

Immun_GSEA <- gseGO(gene = Immun_Table.T_GSEA, OrgDb = org.Dm.eg.db, ont = "ALL", minGSSize = 15, maxGSSize = Inf, eps = 0)
Immun_GSEA@result <- Immun_GSEA@result[!Immun_GSEA@result$ID %fin% BPCCMF,]
Immun_GSEA <- clusterProfiler::simplify(Immun_GSEA)
Immun_GSEA <- setReadable(Immun_GSEA, "org.Dm.eg.db", "ENTREZID")
Immun_heat <- heatplot(Immun_GSEA, foldChange = Immun_GSEA@geneList, showCategory = 15, label_format = 32)
Immun_heat <- Immun_heat@plot_env$p + Immun_heat@scales$scales[[2]] + labs(y = "Category", title = "GSEA of Candidate Genes From Immunity-Selection Studies (Transcriptomics)")
geneSims <- mgeneSim(unique(Immun_heat@data$Gene), semData = d, measure = "Wang")
dissimilarity <- as.dist(1 - geneSims)
hcluster <- hclust(dissimilarity)
dend <- as.dendrogram(hcluster, hang = 0)
cut_hclust <- dendextend::cutree(dend, h = dendextend::heights_per_k.dendrogram(dend)["100"])
sim_clust <- split(x = names(cut_hclust), f = cut_hclust)
sim_clust <- map(.x = sim_clust, .f = nameVector)
sim_clust <- mixedSorts(map(.x = sim_clust, .f = function(x) str_remove(string = x, pattern = "CG")))
sim_clust <- map(.x = sim_clust, .f = rev)
sim_clust <- map(.x = sim_clust, .f = names)
sim_clust <- map(.x = sim_clust, .f = function(x) tcount(Immun_heat@data$Gene)[x])
sim_clust <- map(.x = sim_clust, .f = function(x) x[which.max(x)])
keep <- unlist(map(.x = sim_clust, .f = names))
Immun_heat@data <- Immun_heat@data[Immun_heat@data$Gene %fin% keep,]
Immun_heat
savePlotAsImage(
  file = "Graphs/Immun_GSEA_Transcriptom.png",
  format = "png",
  width = 1312,
  height = 742
)

Stress_GSEA <- gseGO(gene = Stress_Table.T_GSEA, OrgDb = org.Dm.eg.db, ont = "ALL", minGSSize = 15, maxGSSize = Inf, eps = 0)
Stress_GSEA@result <- Stress_GSEA@result[!Stress_GSEA@result$ID %fin% BPCCMF,]
Stress_GSEA <- clusterProfiler::simplify(Stress_GSEA)
Stress_GSEA <- setReadable(Stress_GSEA, "org.Dm.eg.db", "ENTREZID")
Stress_heat <- heatplot(Stress_GSEA, foldChange = Stress_GSEA@geneList, showCategory = 15, label_format = 32)
Stress_heat <- Stress_heat@plot_env$p + Stress_heat@scales$scales[[2]] + labs(y = "Category", title = "GSEA of Candidate Genes From Stress-Selection Studies (Transcriptomics)")
geneSims <- mgeneSim(unique(Stress_heat@data$Gene), semData = d, measure = "Wang")
dissimilarity <- as.dist(1 - geneSims)
hcluster <- hclust(dissimilarity)
dend <- as.dendrogram(hcluster, hang = 0)
cut_hclust <- dendextend::cutree(dend, h = dendextend::heights_per_k.dendrogram(dend)["59"])
sim_clust <- split(x = names(cut_hclust), f = cut_hclust)
sim_clust <- map(.x = sim_clust, .f = nameVector)
sim_clust <- mixedSorts(map(.x = sim_clust, .f = function(x) str_remove(string = x, pattern = "CG")))
sim_clust <- map(.x = sim_clust, .f = rev)
sim_clust <- map(.x = sim_clust, .f = names)
sim_clust <- map(.x = sim_clust, .f = function(x) tcount(Stress_heat@data$Gene)[x])
sim_clust <- map(.x = sim_clust, .f = function(x) x[which.max(x)])
keep <- unlist(map(.x = sim_clust, .f = names))
Stress_heat@data <- Stress_heat@data[Stress_heat@data$Gene %fin% keep,]
Stress_heat
savePlotAsImage(
  file = "Graphs/Stress_GSEA_Transcriptom.png",
  format = "png",
  width = 1312,
  height = 742
)
rm(d, geneSims, dissimilarity, hcluster, dend, cut_hclust, sim_clust, keep)

######################## Longev_Immun/Stress_GSEA #########################
Longev_Immun.T_GSEA <- Gene_Lists.T_GSEA[c("Longevity", "Immunity")]
Longev_Immun_Summ.T_GSEA <- SuperExactTest:::summary.msets(supertest(x = Longev_Immun.T_GSEA, n = 13986, degree = seq(Longev_Immun.T_GSEA)[-1]))
Longev_Immundat.T_GSEA <- Longev_Immun_Summ.T_GSEA$Table
View(Longev_Immundat.T_GSEA)
Longev_Immun_genes.T_GSEA <- unlist(stri_extract_all_words(Longev_Immundat.T_GSEA$Elements))
Longev_Immun.dups <- list(Longevity = Longev_Table.T_GSEA[Longev_Immun_genes.T_GSEA], Immunity = Immun_Table.T_GSEA[Longev_Immun_genes.T_GSEA])
Longev_Immun.dups <- purrr::list_transpose(Longev_Immun.dups)
Longev_Immun.dups <- map(.x = Longev_Immun.dups, .f = mean)
Longev_Immun_genes.T_GSEA <- sort(unlist(Longev_Immun.dups), decreasing = TRUE); rm(Longev_Immun.dups)
Longev_Immun_GSEA <- gseGO(gene = Longev_Immun_genes.T_GSEA, OrgDb = org.Dm.eg.db, ont = "ALL", minGSSize = 15, maxGSSize = Inf, eps = 0)
Longev_Immun_GSEA@result <- Longev_Immun_GSEA@result[!Longev_Immun_GSEA@result$ID %fin% BPCCMF,]
Longev_Immun_GSEA <- clusterProfiler::simplify(Longev_Immun_GSEA)
Longev_Immun_GSEA <- setReadable(Longev_Immun_GSEA, "org.Dm.eg.db", "ENTREZID")
Longev_Immun_heat <- heatplot(Longev_Immun_GSEA, foldChange = Longev_Immun_GSEA@geneList, showCategory = 15, label_format = 29)
Longev_Immun_heat <- Longev_Immun_heat@plot_env$p + Longev_Immun_heat@scales$scales[[2]] + labs(y = "Category", title = "GSEA of Longevity/Immunity-Selection Overlap (Transcriptomics)")
Longev_Immun_heat@data <- Longev_Immun_heat@data[Longev_Immun_heat@data$Gene %fin% names(tcount(Longev_Immun_heat@data$Gene)[!(tcount(Longev_Immun_heat@data$Gene)==1 & str_detect(string = names(tcount(Longev_Immun_heat@data$Gene)), pattern = "CG"))]),]
Longev_Immun_heat
savePlotAsImage(
  file = "Graphs/Longev_Immun_GSEA_Transcriptom.png",
  format = "png",
  width = 1200,
  height = 742
)

Longev_Stress.T_GSEA <- Gene_Lists.T_GSEA[c("Longevity", "Stress")]
Longev_Stress_Summ.T_GSEA <- SuperExactTest:::summary.msets(supertest(x = Longev_Stress.T_GSEA, n = 13986, degree = seq(Longev_Stress.T_GSEA)[-1]))
Longev_Stressdat.T_GSEA <- Longev_Stress_Summ.T_GSEA$Table
View(Longev_Stressdat.T_GSEA)
Longev_Stress_genes.T_GSEA <- unlist(stri_extract_all_words(Longev_Stressdat.T_GSEA$Elements))
Longev_Stress.dups <- list(Longevity = Longev_Table.T_GSEA[Longev_Stress_genes.T_GSEA], Stress = Stress_Table.T_GSEA[Longev_Stress_genes.T_GSEA])
Longev_Stress.dups <- purrr::list_transpose(Longev_Stress.dups)
Longev_Stress.dups <- map(.x = Longev_Stress.dups, .f = mean)
Longev_Stress_genes.T_GSEA <- sort(unlist(Longev_Stress.dups), decreasing = TRUE); rm(Longev_Stress.dups)
Longev_Stress_GSEA <- gseGO(gene = Longev_Stress_genes.T_GSEA, OrgDb = org.Dm.eg.db, ont = "ALL", minGSSize = 15, maxGSSize = Inf, eps = 0)
save.image()
