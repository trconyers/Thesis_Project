########################### Gene_Table edit ###########################
library(readxl)
library(writexl)

Gene_Table <- as.data.frame(read_excel("Data/Gene_Table.xlsx"))
Gene_Table.G <- Gene_Table[Gene_Table$Sequencing=="Genomic",]
FB_invs <- as.data.frame(read_tsv("misc/withdrawns.txt"))
Gene_Table.G <- Gene_Table.G[!Gene_Table.G$`Candidate FB IDs` %fin% FB_invs$FBIDs,]
Gene_Table.G$`Populations Used`[Gene_Table.G$`Populations Used` %fin% c("B vs. O (new)", "B vs. O (old)")] <- "B-type vs. O-type"
rm(FB_invs)

Multi_Pops.G <- Gene_Table.G[, 1:3]
rownames(Multi_Pops.G) <- rownames(Gene_Table.G)
Multi_Pops.G <- unique(Multi_Pops.G)
Gene_Table.G <- Gene_Table.G[rownames(Multi_Pops.G), ]
rm(Multi_Pops.G)

Gene_Table_EG.G <- Gene_Table.G
Gene_Table_EG.G$`Candidate FB IDs` <- FB2EG(Gene_Table.G$`Candidate FB IDs`)

all_pops.G <- Unique(Gene_Table.G$`Populations Used`)
Gene_Table_List.G <- list()
for (i in cli_progress_along(seq(all_pops.G))) {
  Gene_Table_List.G[[i]] <- Gene_Table_EG.G$`Candidate FB IDs`[Gene_Table_EG.G$`Populations Used` == all_pops.G[i]]
}
names(Gene_Table_List.G) <- all_pops.G

########################### Selection Type Overlaps ###########################
library(rstudioapi)
library(SuperExactTest)

eulerr::eulerr_options(eulerr.plot_opts)

Longev_Table.G <- Unique(Gene_Table_EG.G$`Candidate FB IDs`[Gene_Table_EG.G$`Selection Type` == "Longevity"])
Immun_Table.G <- Unique(Gene_Table_EG.G$`Candidate FB IDs`[Gene_Table_EG.G$`Selection Type` == "Immunity"])
Stress_Table.G <- Unique(Gene_Table_EG.G$`Candidate FB IDs`[Gene_Table_EG.G$`Selection Type` == "Stress"])

Gene_Lists.G <- list(Longev_Table.G, Immun_Table.G, Stress_Table.G)
names(Gene_Lists.G) <- c("Longevity", "Immunity", "Stress")

Longev_Immun.G <- Gene_Lists.G[c("Longevity", "Immunity")]
Longev_Immun_Venn.G <- venn.list(Longev_Immun.G)
Longev_Immun_Intersections.G <- supertest(x = Longev_Immun.G, n = 17871, degree = seq(Longev_Immun.G)[-1])
Longev_Immun_Summ.G <- SuperExactTest:::summary.msets(Longev_Immun_Intersections.G)
plot.venn(
  Longev_Immun_Venn.G,
  labels = list(fontsize = 19),
  quantities = list(fontsize = 19, labels = quantlabel.venn(Longev_Immun_Intersections.G)),
  main = list(
    label = "Candidate Genes from Longevity & Immunity Selection (Genomics)", fontsize = 21.5)
)
savePlotAsImage(
  file = "Graphs/Longev_Immun_Genom.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Longev_Immundat.G <- Longev_Immun_Summ.G$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH"))
View(Longev_Immundat.G)
Longev_Immundat.G$Elements <- str_flatten_comma(EG2FB(unlist(strsplit(x = Longev_Immundat.G$Elements, split = ", "))))
Longev_Immundat.G <- Longev_Immundat.G[,-2]
write_excel_csv(Longev_Immundat.G,
                file = "Data/Longev_Immun_Genom.csv",
                append = FALSE,
                escape = "none")

Gene_Venn.G <- venn.list(Gene_Lists.G[c("Longevity", "Immunity", "Stress")])
Selection_Intersections.G <- supertest(x = Gene_Lists.G, n = 17871, degree = seq(Gene_Lists.G)[-1])
Selection_Summ.G <- SuperExactTest:::summary.msets(Selection_Intersections.G)
plot.venn(
  Gene_Venn.G,
  labels = list(fontsize = 19),
  quantities = list(fontsize = 19, labels = quantlabel.venn(Selection_Intersections.G)),
  main = list(
    label = "Candidate Genes by Selection Type (Genomics)",
    fontsize = 22.5
  )
)
savePlotAsImage(
  file = "Graphs/sel_type_Genom.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Selectiondat.G <- Selection_Summ.G$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH")) %>%  filter(P.adj < 0.05)
View(Selectiondat.G)
Selectiondat.G <- sort_by.data.frame(x = Selectiondat.G, y = -Selectiondat.G$FE)
Selectiondat.G <- sort_by.data.frame(x = Selectiondat.G, y = -Selectiondat.G$Degree)
Selectiondat.G$Elements <- unlist(map(.x = strsplit(x = Selectiondat.G$Elements, split = ", "), .f = compose(str_flatten_comma, EG2FB)))
write_excel_csv(Selectiondat.G,
                file = "Data/sel_type_Genom.csv",
                append = FALSE,
                escape = "none")

gp <- get.gpar()
gp$fontfamily <- "serif"
LIS_obj.G <- supertest(Gene_Lists.G, n=17871, degree = seq(Gene_Lists.G)[-1])
LIS_summ.G <- SuperExactTest:::summary.msets(LIS_obj.G)
LIS_obj.G_adj <- LIS_obj.G; LIS_obj.G_adj$P.value <- p.adjust(LIS_obj.G_adj$P.value, method = "BH")
inject(plot(
    LIS_obj.G_adj,
    title = "Longevity, Immunity, & Stress Gene Overlaps (Genomics)",
    minMinusLog10PValue = 1,
    maxMinusLog10PValue = 4,
    !!!SuperExactTest.plot_opts
)); rm(LIS_obj.G_adj)
grab <- grid.grab()
grab$gp <- inject(gpar(!!!gp))
grab[[5]][[1]][["y"]] <- grab[[5]][[1]][["y"]] * 1.5
grid.newpage()
grid.draw(grab)
rm(gp,grab)
savePlotAsImage(
  file = "Graphs/LIS_Genom.svg",
  format = "svg",
  width = 1200,
  height = 742
)

########################### Population Overlaps ###########################
library(readxl)
library(rstudioapi)

Longev_Pops.G <- Unique(Gene_Table.G$`Populations Used`[Gene_Table.G$`Selection Type` == "Longevity"])
Immun_Pops.G <- Unique(Gene_Table.G$`Populations Used`[Gene_Table.G$`Selection Type` == "Immunity"])
Stress_Pops.G <- Unique(Gene_Table.G$`Populations Used`[Gene_Table.G$`Selection Type` == "Stress"])
Gene_Table_List.G_L <-  Gene_Table_List.G[Longev_Pops.G]
Gene_Table_List.G_I <-  Gene_Table_List.G[Immun_Pops.G]
Gene_Table_List.G_S <-  Gene_Table_List.G[Stress_Pops.G]
Longev_Pops_Venn.G <- euler.list(Gene_Table_List.G_L)
Immun_Pops_Venn.G <- euler.list(Gene_Table_List.G_I)
Stress_Pops_Venn.G <- euler.list(Gene_Table_List.G_S)

Longev_Intersections.G <- supertest(x = Gene_Table_List.G_L, n = 17871, degree = seq(Gene_Table_List.G_L)[-1])
Longev_Summ.G <- SuperExactTest:::summary.msets(Longev_Intersections.G)
plot.euler(
  Longev_Pops_Venn.G,
  legend = list(fontsize = 20),
  quantities = list(fontsize = 18, labels = quantlabel(Longev_Intersections.G)),
  main = list(
    label = "Candidate Genes from Genomic Longevity-Selection Studies",
    fontsize = 22.5,
    x = unit(0.675, "npc"),
    y = unit(0, "npc")
  )
)
savePlotAsImage(
  file = "Graphs/lifespan_Genom.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Immun_Intersections.G <- supertest(x = Gene_Table_List.G_I, n = 17871, degree = seq(Gene_Table_List.G_I)[-1])
Immun_Summ.G <- SuperExactTest:::summary.msets(Immun_Intersections.G)
plot.euler(
  Immun_Pops_Venn.G,
  legend = list(fontsize = 20),
  quantities = list(fontsize = 18, labels = quantlabel(Immun_Intersections.G)),
  main = list(
    label = "Candidate Genes from Genomic Immunity-Selection Studies",
    fontsize = 22.5,
    x = unit(0.645, "npc")
  )
)
savePlotAsImage(
  file = "Graphs/immune_Genom.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Stress_Intersections.G <- supertest(x = Gene_Table_List.G_S, n = 17871, degree = seq(Gene_Table_List.G_S)[-1])
Stress_Summ.G <- SuperExactTest:::summary.msets(Stress_Intersections.G)
plot.euler(
  Stress_Pops_Venn.G,
  legend = list(fontsize = 20),
  quantities = list(fontsize = 18, labels = quantlabel(Stress_Intersections.G)),
  main = list(
    label = "Candidate Genes from Genomic Stress-Selection Studies",
    fontsize = 22.5,
    x = unit(0.73, "npc")
  )
)
savePlotAsImage(
  file = "Graphs/stress_Genom.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Longevdat.G <- Longev_Summ.G$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH")) %>%  filter(P.adj < 0.05); View(Longevdat.G)
Longevdat.G <- sort_by.data.frame(x = Longevdat.G, y = -Longevdat.G$FE)
Longevdat.G <- sort_by.data.frame(x = Longevdat.G, y = -Longevdat.G$Degree)
Longevdat.G <- Longevdat.G[Longevdat.G$Observed.Overlap>0,]
Longevdat.G$Elements <- unlist(map(.x = strsplit(x = Longevdat.G$Elements, split = ", "), .f = compose(str_flatten_comma, EG2FB)))
write_excel_csv(Longevdat.G,
                file = "Data/longevity_pops_Genom.csv",
                append = FALSE,
                escape = "none")

Immundat.G <- Immun_Summ.G$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH")); View(Immundat.G)
Immundat.G <- sort_by.data.frame(x = Immundat.G, y = -Immundat.G$FE)
Immundat.G <- sort_by.data.frame(x = Immundat.G, y = -Immundat.G$Degree)
Immundat.G$Elements <- "NA"
write_excel_csv(Immundat.G,
                file = "Data/immunity_pops_Genom.csv",
                append = FALSE,
                escape = "none")

Stressdat.G <- Stress_Summ.G$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH")) %>%  filter(P.adj < 0.05); View(Stressdat.G)
Stressdat.G <- sort_by.data.frame(x = Stressdat.G, y = -Stressdat.G$FE)
Stressdat.G <- sort_by.data.frame(x = Stressdat.G, y = -Stressdat.G$Degree)
Stressdat.G <- Stressdat.G[Stressdat.G$Observed.Overlap>0,]
Stressdat.G$Elements <- unlist(map(.x = strsplit(x = Stressdat.G$Elements, split = ", "), .f = compose(str_flatten_comma, EG2FB)))
write_excel_csv(Stressdat.G,
                file = "Data/stress_pops_Genom.csv",
                append = FALSE,
                escape = "none")
save.image()
