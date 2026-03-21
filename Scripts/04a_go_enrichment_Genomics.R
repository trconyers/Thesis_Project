########################### Selection_ORA ###########################
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

Longev_genes.G <- Longevdat.G$Elements[Longevdat.G$Degree>=Maxs(Longevdat.G$Degree) - 1]
Longev_genes.G <- unname(FB2EG(unique(unlist(strsplit(x = str_flatten_comma(Longev_genes.G), split = ", ")))))
Longev_lds.G <- c(is.longevity(Longev_genes.G), is.immunity(Longev_genes.G), is.stress(Longev_genes.G))
Longev_Clusts.G <- list(Lifespan = is.longevity(Longev_genes.G), Defense = is.immunity(Longev_genes.G), Stress_Response = is.stress(Longev_genes.G), Other = setdiff(Longev_genes.G, Longev_lds.G))
Longev_enrich.G <- enrichGO(gene = Longev_genes.G, OrgDb = org.Dm.eg.db)
Longev_enrich.G_comp <- compareCluster(
  geneCluster = Longev_Clusts.G,
  fun = enrichGO,
  OrgDb = org.Dm.eg.db,
  ont = "ALL",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)
LongevityBP_ora.G <- hyperGTest(new(Class = "GOHyperGParams",
                                    ontology = "BP",
                                    geneIds = Longev_enrich.G@gene,
                                    universeGeneIds = Longev_enrich.G@universe,
                                    annotation = "org.Dm.eg.db",
                                    pvalueCutoff = 0.05,
                                    testDirection = "over",
                                    conditional = TRUE))
LongevityBP_ora.G_summ <- summary(LongevityBP_ora.G)
LongevityCC_ora.G <- hyperGTest(new(Class = "GOHyperGParams",
                                    ontology = "CC",
                                    geneIds = Longev_enrich.G@gene,
                                    universeGeneIds = Longev_enrich.G@universe,
                                    annotation = "org.Dm.eg.db",
                                    pvalueCutoff = 0.05,
                                    testDirection = "over",
                                    conditional = TRUE))
LongevityCC_ora.G_summ <- summary(LongevityCC_ora.G)
LongevityMF_ora.G <- hyperGTest(new(Class = "GOHyperGParams",
                                    ontology = "MF",
                                    geneIds = Longev_enrich.G@gene,
                                    universeGeneIds = Longev_enrich.G@universe,
                                    annotation = "org.Dm.eg.db",
                                    pvalueCutoff = 0.05,
                                    testDirection = "over",
                                    conditional = TRUE))
LongevityMF_ora.G_summ <- summary(LongevityMF_ora.G)
colnames(LongevityBP_ora.G_summ)[1] <- colnames(LongevityCC_ora.G_summ)[1] <- colnames(LongevityMF_ora.G_summ)[1] <- "GOID"
Longevity_ora.G_summ <- rbind.data.frame(LongevityBP_ora.G_summ, LongevityCC_ora.G_summ, LongevityMF_ora.G_summ) %>%  mutate(P.adj = p.adjust(Pvalue, method = "BH")) %>%  filter(P.adj < 0.05)
rm(LongevityBP_ora.G_summ, LongevityCC_ora.G_summ, LongevityMF_ora.G_summ)
rownames(Longevity_ora.G_summ) <- Longevity_ora.G_summ$GOID
Longev_enrich.G_comp@compareClusterResult <- Longev_enrich.G_comp@compareClusterResult[Longev_enrich.G_comp@compareClusterResult$ID %fin% Longevity_ora.G_summ$GOID,]
Longev_enrich.G_comp@compareClusterResult$p.adjust <- Longevity_ora.G_summ[Longev_enrich.G_comp@compareClusterResult$ID, "P.adj"]
Longev_enrich.G_comp@compareClusterResult$EnrichRatio <- (Longev_enrich.G_comp@compareClusterResult$Count)/(as.numeric(gsub("/.*$", "", Longev_enrich.G_comp@compareClusterResult$BgRatio)))
Longev_dot.G <- clusterProfiler::simplify(Longev_enrich.G_comp)
Longev_dotplot.G <- dotplot(Longev_dot.G, label_format = 32, by = "EnrichRatio", showCategory = 4) + get_theme() + labs(x = "Function", title = "Overrepresented GO Categories in Longevity-Selection Overlaps (Genomics)")
Longev_dotplot.G
savePlotAsImage(
  file = "Graphs/Longev_GOs_Genom.png",
  format = "png",
  width = 1200,
  height = 742
)

Stress_genes.G <- Stressdat.G$Elements[Stressdat.G$Degree>=Maxs(Stressdat.G$Degree) - 1]
Stress_genes.G <- unname(FB2EG(unique(unlist(strsplit(x = str_flatten_comma(Stress_genes.G), split = ", ")))))
Stress_lds.G <- c(is.longevity(Stress_genes.G), is.immunity(Stress_genes.G), is.stress(Stress_genes.G))
Stress_Clusts.G <- list(Lifespan = is.longevity(Stress_genes.G), Defense = is.immunity(Stress_genes.G), Stress_Response = is.stress(Stress_genes.G), Other = setdiff(Stress_genes.G, Stress_lds.G))
Stress_enrich.G <- enrichGO(gene = Stress_genes.G, OrgDb = org.Dm.eg.db)
Stress_enrich.G_comp <- compareCluster(
  geneCluster = Stress_Clusts.G,
  fun = enrichGO,
  OrgDb = org.Dm.eg.db,
  ont = "ALL",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)
StressBP_ora.G <- hyperGTest(new(Class = "GOHyperGParams",
                                 ontology = "BP",
                                 geneIds = Stress_enrich.G@gene,
                                 universeGeneIds = Stress_enrich.G@universe,
                                 annotation = "org.Dm.eg.db",
                                 pvalueCutoff = 0.05,
                                 testDirection = "over",
                                 conditional = TRUE))
StressBP_ora.G_summ <- summary(StressBP_ora.G)
StressCC_ora.G <- hyperGTest(new(Class = "GOHyperGParams",
                                 ontology = "CC",
                                 geneIds = Stress_enrich.G@gene,
                                 universeGeneIds = Stress_enrich.G@universe,
                                 annotation = "org.Dm.eg.db",
                                 pvalueCutoff = 0.05,
                                 testDirection = "over",
                                 conditional = TRUE))
StressCC_ora.G_summ <- summary(StressCC_ora.G)
StressMF_ora.G <- hyperGTest(new(Class = "GOHyperGParams",
                                 ontology = "MF",
                                 geneIds = Stress_enrich.G@gene,
                                 universeGeneIds = Stress_enrich.G@universe,
                                 annotation = "org.Dm.eg.db",
                                 pvalueCutoff = 0.05,
                                 testDirection = "over",
                                 conditional = TRUE))
StressMF_ora.G_summ <- summary(StressMF_ora.G)
colnames(StressBP_ora.G_summ)[1] <- colnames(StressCC_ora.G_summ)[1] <- colnames(StressMF_ora.G_summ)[1] <- "GOID"
Stress_ora.G_summ <- rbind.data.frame(StressBP_ora.G_summ, StressCC_ora.G_summ, StressMF_ora.G_summ) %>%  mutate(P.adj = p.adjust(Pvalue, method = "BH")) %>%  filter(P.adj < 0.05)
rm(StressBP_ora.G_summ, StressCC_ora.G_summ, StressMF_ora.G_summ)
rownames(Stress_ora.G_summ) <- Stress_ora.G_summ$GOID
Stress_enrich.G_comp@compareClusterResult <- Stress_enrich.G_comp@compareClusterResult[Stress_enrich.G_comp@compareClusterResult$ID %fin% Stress_ora.G_summ$GOID,]
Stress_enrich.G_comp@compareClusterResult$p.adjust <- Stress_ora.G_summ[Stress_enrich.G_comp@compareClusterResult$ID, "P.adj"]
Stress_enrich.G_comp@compareClusterResult$EnrichRatio <- (Stress_enrich.G_comp@compareClusterResult$Count)/(as.numeric(gsub("/.*$", "", Stress_enrich.G_comp@compareClusterResult$BgRatio)))
Stress_dot.G <- clusterProfiler::simplify(Stress_enrich.G_comp)
Stress_dotplot.G <- dotplot(Stress_dot.G, label_format = 32, by = "EnrichRatio", showCategory = 4) + get_theme() + labs(x = "Function", title = "Overrepresented GO Categories in Stress-Selection Overlaps (Genomics)")
Stress_dotplot.G
savePlotAsImage(
  file = "Graphs/Stress_GOs_Genom.png",
  format = "png",
  width = 1200,
  height = 742
)

######################## Longev_Immun/Stress_ORA #########################
LI_genes.G <- unname(FB2EG(unique(unlist(strsplit(x = Selectiondat.G$Elements[Selectiondat.G$Intersections=="Longevity & Immunity"], split = ", ")))))
LI_lds.G <- c(is.longevity(LI_genes.G), is.immunity(LI_genes.G), is.stress(LI_genes.G))
LI_Clusts.G <- list(Lifespan = is.longevity(LI_genes.G), Defense = is.immunity(LI_genes.G), Stress_Response = is.stress(LI_genes.G), Other = setdiff(LI_genes.G, LI_lds.G))
LI_enrich.G <- enrichGO(gene = LI_genes.G, OrgDb = org.Dm.eg.db)
LI_enrich.G_comp <- compareCluster(
  geneCluster = LI_Clusts.G,
  fun = enrichGO,
  OrgDb = org.Dm.eg.db,
  ont = "ALL",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)
LI_BP_ora.G <- hyperGTest(new(Class = "GOHyperGParams",
                              ontology = "BP",
                              geneIds = LI_enrich.G@gene,
                              universeGeneIds = LI_enrich.G@universe,
                              annotation = "org.Dm.eg.db",
                              pvalueCutoff = 0.05,
                              testDirection = "over",
                              conditional = TRUE))
LI_BP_ora.G_summ <- summary(LI_BP_ora.G)
LI_CC_ora.G <- hyperGTest(new(Class = "GOHyperGParams",
                              ontology = "CC",
                              geneIds = LI_enrich.G@gene,
                              universeGeneIds = LI_enrich.G@universe,
                              annotation = "org.Dm.eg.db",
                              pvalueCutoff = 0.05,
                              testDirection = "over",
                              conditional = TRUE))
LI_CC_ora.G_summ <- summary(LI_CC_ora.G)
LI_MF_ora.G <- hyperGTest(new(Class = "GOHyperGParams",
                              ontology = "MF",
                              geneIds = LI_enrich.G@gene,
                              universeGeneIds = LI_enrich.G@universe,
                              annotation = "org.Dm.eg.db",
                              pvalueCutoff = 0.05,
                              testDirection = "over",
                              conditional = TRUE))
LI_MF_ora.G_summ <- summary(LI_MF_ora.G)
colnames(LI_BP_ora.G_summ)[1] <- colnames(LI_CC_ora.G_summ)[1] <- colnames(LI_MF_ora.G_summ)[1] <- "GOID"
LI_ora.G_summ <- rbind.data.frame(LI_BP_ora.G_summ, LI_CC_ora.G_summ, LI_MF_ora.G_summ) %>%  mutate(P.adj = p.adjust(Pvalue, method = "BH")) %>%  filter(P.adj < 0.05)
rm(LI_BP_ora.G_summ, LI_CC_ora.G_summ, LI_MF_ora.G_summ)
rownames(LI_ora.G_summ) <- LI_ora.G_summ$GOID
LI_enrich.G_comp@compareClusterResult <- LI_enrich.G_comp@compareClusterResult[LI_enrich.G_comp@compareClusterResult$ID %fin% LI_ora.G_summ$GOID,]
LI_enrich.G_comp@compareClusterResult$p.adjust <- LI_ora.G_summ[LI_enrich.G_comp@compareClusterResult$ID, "P.adj"]
LI_enrich.G_comp@compareClusterResult$EnrichRatio <- (LI_enrich.G_comp@compareClusterResult$Count)/(as.numeric(gsub("/.*$", "", LI_enrich.G_comp@compareClusterResult$BgRatio)))
LI_dot.G <- clusterProfiler::simplify(LI_enrich.G_comp)
LI_dotplot.G <- dotplot(LI_dot.G, label_format = 25, showCategory = 3, by = "EnrichRatio") + get_theme() + labs(x = "Function", title = "Overrepresented GO Categories in Longevity/Immunity-Selection Overlap (Genomics)")
LI_dotplot.G
savePlotAsImage(
  file = "Graphs/Longev_Immun_GOs_Genom.png",
  format = "png",
  width = 1200,
  height = 742
)

LS_genes.G <- unname(FB2EG(unique(unlist(strsplit(x = Selectiondat.G$Elements[Selectiondat.G$Intersections=="Longevity & Stress"], split = ", ")))))
LS_lds.G <- c(is.longevity(LS_genes.G), is.immunity(LS_genes.G), is.stress(LS_genes.G))
LS_Clusts.G <- list(Lifespan = is.longevity(LS_genes.G), Defense = is.immunity(LS_genes.G), Stress_Response = is.stress(LS_genes.G), Other = setdiff(LS_genes.G, LS_lds.G))
LS_enrich.G <- enrichGO(gene = LS_genes.G, OrgDb = org.Dm.eg.db)
LS_enrich.G_comp <- compareCluster(
  geneCluster = LS_Clusts.G,
  fun = enrichGO,
  OrgDb = org.Dm.eg.db,
  ont = "ALL",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)
LS_BP_ora.G <- hyperGTest(new(Class = "GOHyperGParams",
                              ontology = "BP",
                              geneIds = LS_enrich.G@gene,
                              universeGeneIds = LS_enrich.G@universe,
                              annotation = "org.Dm.eg.db",
                              pvalueCutoff = 0.05,
                              testDirection = "over",
                              conditional = TRUE))
LS_BP_ora.G_summ <- summary(LS_BP_ora.G)
LS_CC_ora.G <- hyperGTest(new(Class = "GOHyperGParams",
                              ontology = "CC",
                              geneIds = LS_enrich.G@gene,
                              universeGeneIds = LS_enrich.G@universe,
                              annotation = "org.Dm.eg.db",
                              pvalueCutoff = 0.05,
                              testDirection = "over",
                              conditional = TRUE))
LS_CC_ora.G_summ <- summary(LS_CC_ora.G)
LS_MF_ora.G <- hyperGTest(new(Class = "GOHyperGParams",
                              ontology = "MF",
                              geneIds = LS_enrich.G@gene,
                              universeGeneIds = LS_enrich.G@universe,
                              annotation = "org.Dm.eg.db",
                              pvalueCutoff = 0.05,
                              testDirection = "over",
                              conditional = TRUE))
LS_MF_ora.G_summ <- summary(LS_MF_ora.G)
colnames(LS_BP_ora.G_summ)[1] <- colnames(LS_CC_ora.G_summ)[1] <- colnames(LS_MF_ora.G_summ)[1] <- "GOID"
LS_ora.G_summ <- rbind.data.frame(LS_BP_ora.G_summ, LS_CC_ora.G_summ, LS_MF_ora.G_summ) %>%  mutate(P.adj = p.adjust(Pvalue, method = "BH")) %>%  filter(P.adj < 0.05)
rm(LS_BP_ora.G_summ, LS_CC_ora.G_summ, LS_MF_ora.G_summ)
rownames(LS_ora.G_summ) <- LS_ora.G_summ$GOID
LS_enrich.G_comp@compareClusterResult <- LS_enrich.G_comp@compareClusterResult[LS_enrich.G_comp@compareClusterResult$ID %fin% LS_ora.G_summ$GOID,]
LS_enrich.G_comp@compareClusterResult$p.adjust <- LS_ora.G_summ[LS_enrich.G_comp@compareClusterResult$ID, "P.adj"]
LS_enrich.G_comp@compareClusterResult$EnrichRatio <- (LS_enrich.G_comp@compareClusterResult$Count)/(as.numeric(gsub("/.*$", "", LS_enrich.G_comp@compareClusterResult$BgRatio)))
LS_dot.G <- clusterProfiler::simplify(LS_enrich.G_comp)
LS_dotplot.G <- dotplot(LS_dot.G, label_format = 25, showCategory = 3, by = "EnrichRatio") + get_theme() + labs(x = "Function", title = "Overrepresented GO Categories in Longevity/Stress-Selection Overlap (Genomics)")
LS_dotplot.G
savePlotAsImage(
  file = "Graphs/Longev_Stress_GOs_Genom.png",
  format = "png",
  width = 1200,
  height = 742
)

dotplot_list.G <- cowplot::plot_grid(LI_dotplot.G, LS_dotplot.G, labels = c("A", "B"), ncol = 2)
dotplot_list.G
savePlotAsImage(
  file = "Graphs/Longev_Overlaps_GOs_Genom.png",
  format = "png",
  width = 1800,
  height = 742
)
save.image()
