############################################################### Full analysis ###############################################################
########################### Gene_Table edit ###########################
library(readxl)
library(writexl)

Gene_Table <- as.data.frame(read_excel("Data/Gene_Table.xlsx"))
Gene_Table.T <- Gene_Table[Gene_Table$Sequencing=="Transcriptomic",]
FB_invs <- as.data.frame(read_tsv("Scripts/misc/withdrawns.txt"))
Gene_Table.T <-Gene_Table.T[!Gene_Table.T$`Candidate FB IDs` %fin% FB_invs$FBIDs,]
rm(FB_invs)

Multi_Pops.T <- Gene_Table.T[, 1:3]
rownames(Multi_Pops.T) <- rownames(Gene_Table.T)
Multi_Pops.T <- unique(Multi_Pops.T)
Gene_Table.T <- Gene_Table.T[rownames(Multi_Pops.T), ]
rm(Multi_Pops.T)

Gene_Table_EG.T <- Gene_Table.T
Gene_Table_EG.T$`Candidate FB IDs` <- FB2EG(Gene_Table.T$`Candidate FB IDs`)

all_pops.T <- Unique(Gene_Table.T$`Populations Used`)
Gene_Table_List.T <- list()
for (i in cli_progress_along(seq(all_pops.T))) {
  Gene_Table_List.T[[i]] <- Gene_Table_EG.T$`Candidate FB IDs`[Gene_Table_EG.T$`Populations Used` == all_pops.T[i]]
}
names(Gene_Table_List.T) <- all_pops.T

########################### Selection Type Overlaps ###########################
library(rstudioapi)
library(SuperExactTest)

eulerr::eulerr_options(eulerr.plot_opts)

Longev_Table.T <- Unique(Gene_Table_EG.T$`Candidate FB IDs`[Gene_Table_EG.T$`Selection Type` == "Longevity"])
Immun_Table.T <- Unique(Gene_Table_EG.T$`Candidate FB IDs`[Gene_Table_EG.T$`Selection Type` == "Immunity"])
Stress_Table.T <- Unique(Gene_Table_EG.T$`Candidate FB IDs`[Gene_Table_EG.T$`Selection Type` == "Stress"])

Gene_Lists.T <- list(Longev_Table.T, Immun_Table.T, Stress_Table.T)
names(Gene_Lists.T) <- c("Longevity", "Immunity", "Stress")

Longev_Immun.T <- Gene_Lists.T[c("Longevity", "Immunity")]
Longev_Immun_Venn.T <- venn.list(Longev_Immun.T)
Longev_Immun_Intersections.T <- supertest(x = Longev_Immun.T, n = 13986, degree = seq(Longev_Immun.T)[-1])
Longev_Immun_Summ.T <- SuperExactTest:::summary.msets(Longev_Immun_Intersections.T)
plot.venn(
  Longev_Immun_Venn.T,
  labels = list(fontsize = 19),
  quantities = list(fontsize = 19, labels = quantlabel.venn(Longev_Immun_Intersections.T)),
  main = list(
    label = "Candidate Genes from Longevity & Immunity Selection (Transcriptomics)", fontsize = 20)
)
savePlotAsImage(
  file = "Graphs/Longev_Immun_Transcriptom.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Longev_Immundat.T <- Longev_Immun_Summ.T$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH"))
View(Longev_Immundat.T)
Longev_Immundat.T$Elements <- str_flatten_comma(EG2FB(unlist(strsplit(x = Longev_Immundat.T$Elements, split = ", "))))
write_excel_csv(Longev_Immundat.T,
                file = "Data/Longev_Immun_Transcriptom.csv",
                append = FALSE,
                escape = "none")

Gene_Venn.T <- venn.list(Gene_Lists.T[c("Longevity", "Immunity", "Stress")])
Selection_Intersections.T <- supertest(x = Gene_Lists.T, n = 13986, degree = seq(Gene_Lists.T)[-1])
Selection_Summ.T <- SuperExactTest:::summary.msets(Selection_Intersections.T)
plot.venn(
  Gene_Venn.T,
  labels = list(fontsize = 19),
  quantities = list(fontsize = 19, labels = quantlabel.venn(Selection_Intersections.T)),
  main = list(
    label = "Candidate Genes by Selection Type (Transcriptomics)",
    fontsize = 20
  )
)
savePlotAsImage(
  file = "Graphs/sel_type_Transcriptom.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Selectiondat.T <- Selection_Summ.T$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH")) %>%  filter(P.adj < 0.05)
View(Selectiondat.T)
Selectiondat.T <- sort_by.data.frame(x = Selectiondat.T, y = -Selectiondat.T$FE)
Selectiondat.T <- sort_by.data.frame(x = Selectiondat.T, y = -Selectiondat.T$Degree)
Selectiondat.T$Elements <- unlist(map(.x = strsplit(x = Selectiondat.T$Elements, split = ", "), .f = compose(str_flatten_comma, EG2FB)))
write_excel_csv(Selectiondat.T,
                file = "Data/sel_type_Transcriptom.csv",
                append = FALSE,
                escape = "none")

gp <- get.gpar()
gp$fontfamily <- "serif"
LIS_obj.T <- supertest(Gene_Lists.T, n=13986, degree = seq(Gene_Lists.T)[-1])
LIS_summ.T <- SuperExactTest:::summary.msets(LIS_obj.T)
LIS_obj.T_adj <- LIS_obj.T; LIS_obj.T_adj$P.value <- p.adjust(LIS_obj.T_adj$P.value, method = "BH")
inject(plot(
    LIS_obj.T_adj,
    title = "Longevity, Immunity, & Stress Gene Overlaps (Transcriptomics)",
    minMinusLog10PValue = 1,
    maxMinusLog10PValue = 4,
    !!!SuperExactTest.plot_opts
)); rm(LIS_obj.T_adj)
grab <- grid.grab()
grab$gp <- inject(gpar(!!!gp))
grab[[5]][[1]][["y"]] <- grab[[5]][[1]][["y"]] * 1.5
grid.newpage()
grid.draw(grab)
rm(gp,grab)
savePlotAsImage(
  file = "Graphs/LIS_Transcriptom.svg",
  format = "svg",
  width = 1200,
  height = 742
)

########################### Population Overlaps ###########################
library(readxl)
library(rstudioapi)

Longev_Pops.T <- Unique(Gene_Table.T$`Populations Used`[Gene_Table.T$`Selection Type` == "Longevity"])
Immun_Pops.T <- Unique(Gene_Table.T$`Populations Used`[Gene_Table.T$`Selection Type` == "Immunity"])
Stress_Pops.T <- Unique(Gene_Table.T$`Populations Used`[Gene_Table.T$`Selection Type` == "Stress"])
Gene_Table_List.T_L <-  Gene_Table_List.T[Longev_Pops.T]
Gene_Table_List.T_I <-  Gene_Table_List.T[Immun_Pops.T]
Gene_Table_List.T_S <-  Gene_Table_List.T[Stress_Pops.T]
Longev_Pops_Venn.T <- euler.list(Gene_Table_List.T_L)
Immun_Pops_Venn.T <- euler.list(Gene_Table_List.T_I)
Stress_Pops_Venn.T <- euler.list(Gene_Table_List.T_S)

Longev_Intersections.T <- supertest(x = Gene_Table_List.T_L, n = 13986, degree = seq(Gene_Table_List.T_L)[-1])
Longev_Summ.T <- SuperExactTest:::summary.msets(Longev_Intersections.T)
plot.euler(
  Longev_Pops_Venn.T,
  legend = list(fontsize = 20),
  quantities = list(fontsize = 18, labels = quantlabel(Longev_Intersections.T)),
  main = list(
    label = "Candidate Genes from Transcriptomic Longevity-Selection Studies",
    fontsize = 20,
    x = unit(0.64, "npc")
  )
)
savePlotAsImage(
  file = "Graphs/lifespan_Transcriptom.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Immun_Intersections.T <- supertest(x = Gene_Table_List.T_I, n = 13986, degree = seq(Gene_Table_List.T_I)[-1])
Immun_Summ.T <- SuperExactTest:::summary.msets(Immun_Intersections.T)
plot.euler(
  Immun_Pops_Venn.T,
  legend = list(fontsize = 20),
  quantities = list(fontsize = 18, labels = quantlabel(Immun_Intersections.T)),
  main = list(
    label = "Candidate Genes from Transcriptomic Immunity-Selection Studies",
    fontsize = 20,
    x = unit(0.64, "npc"),
    y = unit(1, "npc")
  )
)
savePlotAsImage(
  file = "Graphs/immune_Transcriptom.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Stress_Intersections.T <- supertest(x = Gene_Table_List.T_S, n = 13986, degree = seq(Gene_Table_List.T_S)[-1])
Stress_Summ.T <- SuperExactTest:::summary.msets(Stress_Intersections.T)
plot.euler(
  Stress_Pops_Venn.T,
  legend = list(fontsize = 20),
  quantities = list(fontsize = 18, labels = quantlabel(Stress_Intersections.T)),
  main = list(
    label = "Candidate Genes from Transcriptomic Stress-Selection Studies",
    fontsize = 20,
    x = unit(0.61, "npc"),
    y = unit(0, "npc")
  )
)
savePlotAsImage(
  file = "Graphs/stress_Transcriptom.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Longevdat.T <- Longev_Summ.T$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH"))
View(Longevdat.T)
Longevdat.T <- sort_by.data.frame(x = Longevdat.T, y = -Longevdat.T$FE)
Longevdat.T <- sort_by.data.frame(x = Longevdat.T, y = -Longevdat.T$Degree)
Longevdat.T <- Longevdat.T[Longevdat.T$Observed.Overlap>0,]
Longevdat.T$Elements <- "non-significant"
write_excel_csv(Longevdat.T,
                file = "Data/longevity_pops_Transcriptom.csv",
                append = FALSE,
                escape = "none")

Immundat.T <- Immun_Summ.T$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH")) %>%  filter(P.adj < 0.05)
View(Immundat.T)
Immundat.T <- sort_by.data.frame(x = Immundat.T, y = -Immundat.T$FE)
Immundat.T <- sort_by.data.frame(x = Immundat.T, y = -Immundat.T$Degree)
Immundat.T <- Immundat.T[Immundat.T$Observed.Overlap>0,]
Immundat.T$Elements <- unlist(map(.x = strsplit(x = Immundat.T$Elements, split = ", "), .f = compose(str_flatten_comma, EG2FB)))
write_excel_csv(Immundat.T,
                file = "Data/immunity_pops_Transcriptom.csv",
                append = FALSE,
                escape = "none")

Stressdat.T <- Stress_Summ.T$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH")) %>%  filter(P.adj < 0.05)
View(Stressdat.T)
Stressdat.T <- sort_by.data.frame(x = Stressdat.T, y = -Stressdat.T$FE)
Stressdat.T <- sort_by.data.frame(x = Stressdat.T, y = -Stressdat.T$Degree)
Stressdat.T <- Stressdat.T[Stressdat.T$Observed.Overlap>0,]
Stressdat.T$Elements <- unlist(map(.x = strsplit(x = Stressdat.T$Elements, split = ", "), .f = compose(str_flatten_comma, EG2FB)))
write_excel_csv(Stressdat.T,
                file = "Data/stress_pops_Transcriptom.csv",
                append = FALSE,
                escape = "none")

############################################################### Up-regulated ###############################################################
########################### Gene_Table edit ###########################
library(readxl)
library(writexl)

Gene_Table <- as.data.frame(read_excel("Data/Gene_Table.xlsx"))
Gene_Table.T_Up <- Gene_Table[Gene_Table$Sequencing=="Transcriptomic" & Gene_Table$`Direction (Transcriptomics)`=="Up",]
FB_invs <- as.data.frame(read_tsv("Scripts/misc/withdrawns.txt"))
Gene_Table.T_Up <-Gene_Table.T_Up[!Gene_Table.T_Up$`Candidate FB IDs` %fin% FB_invs$FBIDs,]
rm(FB_invs)

Multi_Pops.T_Up <- Gene_Table.T_Up[, 1:3]
rownames(Multi_Pops.T_Up) <- rownames(Gene_Table.T_Up)
Multi_Pops.T_Up <- unique(Multi_Pops.T_Up)
Gene_Table.T_Up <- Gene_Table.T_Up[rownames(Multi_Pops.T_Up), ]
rm(Multi_Pops.T_Up)

Gene_Table_EG.T_Up <- Gene_Table.T_Up
Gene_Table_EG.T_Up$`Candidate FB IDs` <- FB2EG(Gene_Table.T_Up$`Candidate FB IDs`)

all_pops.T_Up <- Unique(Gene_Table.T_Up$`Populations Used`)
Gene_Table_List.T_Up <- list()
for (i in cli_progress_along(seq(all_pops.T_Up))) {
  Gene_Table_List.T_Up[[i]] <- Gene_Table_EG.T_Up$`Candidate FB IDs`[Gene_Table_EG.T_Up$`Populations Used` == all_pops.T_Up[i]]
}
names(Gene_Table_List.T_Up) <- all_pops.T_Up

########################### Selection Type Overlaps ###########################
library(rstudioapi)
library(SuperExactTest)

eulerr::eulerr_options(eulerr.plot_opts)

Longev_Table.T_Up <- Unique(Gene_Table_EG.T_Up$`Candidate FB IDs`[Gene_Table_EG.T_Up$`Selection Type` == "Longevity"])
Immun_Table.T_Up <- Unique(Gene_Table_EG.T_Up$`Candidate FB IDs`[Gene_Table_EG.T_Up$`Selection Type` == "Immunity"])
Stress_Table.T_Up <- Unique(Gene_Table_EG.T_Up$`Candidate FB IDs`[Gene_Table_EG.T_Up$`Selection Type` == "Stress"])

Gene_Lists.T_Up <- list(Longev_Table.T_Up, Immun_Table.T_Up, Stress_Table.T_Up)
names(Gene_Lists.T_Up) <- c("Longevity", "Immunity", "Stress")

Longev_Immun.T_Up <- Gene_Lists.T_Up[c("Longevity", "Immunity")]
Longev_Immun_Venn.T_Up <- venn.list(Longev_Immun.T_Up)
Longev_Immun_Intersections.T_Up <- supertest(x = Longev_Immun.T_Up, n = 13986, degree = seq(Longev_Immun.T_Up)[-1])
Longev_Immun_Summ.T_Up <- SuperExactTest:::summary.msets(Longev_Immun_Intersections.T_Up)
plot.venn(
  Longev_Immun_Venn.T_Up,
  labels = list(fontsize = 19),
  quantities = list(fontsize = 19, labels = quantlabel.venn(Longev_Immun_Intersections.T_Up)),
  main = list(
    label = "Up-regulated Candidate Genes from Longevity & Immunity Selection", fontsize = 20)
)
savePlotAsImage(
  file = "Graphs/Longev_Immun_Transcriptom_Up.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Longev_Immundat.T_Up <- Longev_Immun_Summ.T_Up$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH"))
View(Longev_Immundat.T_Up)
Longev_Immundat.T_Up$Elements <- "non-significant"
write_excel_csv(Longev_Immundat.T_Up,
                file = "Data/Longev_Immun_Transcriptom_Up.csv",
                append = FALSE,
                escape = "none")

Gene_Venn.T_Up <- venn.list(Gene_Lists.T_Up[c("Longevity", "Immunity", "Stress")])
Selection_Intersections.T_Up <- supertest(x = Gene_Lists.T_Up, n = 13986, degree = seq(Gene_Lists.T_Up)[-1])
Selection_Summ.T_Up <- SuperExactTest:::summary.msets(Selection_Intersections.T_Up)
plot.venn(
  Gene_Venn.T_Up,
  labels = list(fontsize = 19),
  quantities = list(fontsize = 19, labels = quantlabel.venn(Selection_Intersections.T_Up)),
  main = list(
    label = "Up-regulated Candidate Genes by Selection Type",
    fontsize = 20
  )
)
savePlotAsImage(
  file = "Graphs/sel_type_Transcriptom_Up.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Selectiondat.T_Up <- Selection_Summ.T_Up$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH"))
View(Selectiondat.T_Up)
Selectiondat.T_Up <- sort_by.data.frame(x = Selectiondat.T_Up, y = -Selectiondat.T_Up$FE)
Selectiondat.T_Up <- sort_by.data.frame(x = Selectiondat.T_Up, y = -Selectiondat.T_Up$Degree)
Selectiondat.T_Up$Elements <- "non-significant"
write_excel_csv(Selectiondat.T_Up,
                file = "Data/sel_type_Transcriptom_Up.csv",
                append = FALSE,
                escape = "none")

gp <- get.gpar()
gp$fontfamily <- "serif"
LIS_obj.T_Up <- supertest(Gene_Lists.T_Up, n=13986, degree = seq(Gene_Lists.T_Up)[-1])
LIS_summ.T_Up <- SuperExactTest:::summary.msets(LIS_obj.T_Up)
LIS_obj.T_Up_adj <- LIS_obj.T_Up; LIS_obj.T_Up_adj$P.value <- p.adjust(LIS_obj.T_Up_adj$P.value, method = "BH")
inject(plot(
    LIS_obj.T_Up_adj,
    title = "Longevity, Immunity, & Stress Up-regulated Gene Overlaps",
    minMinusLog10PValue = 1,
    maxMinusLog10PValue = 4,
    !!!SuperExactTest.plot_opts
)); rm(LIS_obj.T_Up_adj)
grab <- grid.grab()
grab$gp <- inject(gpar(!!!gp))
grab[[5]][[1]][["y"]] <- grab[[5]][[1]][["y"]] * 1.5
grid.newpage()
grid.draw(grab)
rm(gp,grab)
savePlotAsImage(
  file = "Graphs/LIS_Transcriptom_Up.svg",
  format = "svg",
  width = 1200,
  height = 742
)

########################### Population Overlaps ###########################
library(readxl)
library(rstudioapi)

Longev_Pops.T_Up <- Unique(Gene_Table.T_Up$`Populations Used`[Gene_Table.T_Up$`Selection Type` == "Longevity"])
Immun_Pops.T_Up <- Unique(Gene_Table.T_Up$`Populations Used`[Gene_Table.T_Up$`Selection Type` == "Immunity"])
Stress_Pops.T_Up <- Unique(Gene_Table.T_Up$`Populations Used`[Gene_Table.T_Up$`Selection Type` == "Stress"])
Gene_Table_List.T_L.Up <-  Gene_Table_List.T_Up[Longev_Pops.T_Up]
Gene_Table_List.T_I.Up <-  Gene_Table_List.T_Up[Immun_Pops.T_Up]
Gene_Table_List.T_S.Up <-  Gene_Table_List.T_Up[Stress_Pops.T_Up]
Longev_Pops_Venn.T_Up <- euler.list(Gene_Table_List.T_L.Up)
Immun_Pops_Venn.T_Up <- euler.list(Gene_Table_List.T_I.Up)
Stress_Pops_Venn.T_Up <- euler.list(Gene_Table_List.T_S.Up)

Longev_Intersections.T_Up <- supertest(x = Gene_Table_List.T_L.Up, n = 13986, degree = seq(Gene_Table_List.T_L.Up)[-1])
Longev_Summ.T_Up <- SuperExactTest:::summary.msets(Longev_Intersections.T_Up)
plot.euler(
  Longev_Pops_Venn.T_Up,
  legend = list(fontsize = 20),
  quantities = list(fontsize = 18, labels = quantlabel(Longev_Intersections.T_Up)),
  main = list(
    label = "Up-regulated Candidate Genes from Longevity-Selection Studies",
    fontsize = 20,
    x = unit(0.64, "npc")
  )
)
savePlotAsImage(
  file = "Graphs/lifespan_Transcriptom_Up.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Immun_Intersections.T_Up <- supertest(x = Gene_Table_List.T_I.Up, n = 13986, degree = seq(Gene_Table_List.T_I.Up)[-1])
Immun_Summ.T_Up <- SuperExactTest:::summary.msets(Immun_Intersections.T_Up)
plot.euler(
  Immun_Pops_Venn.T_Up,
  legend = list(fontsize = 20),
  quantities = list(fontsize = 18, labels = quantlabel(Immun_Intersections.T_Up)),
  main = list(
    label = "Up-regulated Candidate Genes from Immunity-Selection Studies",
    fontsize = 20,
    x = unit(0.64, "npc"),
    y = unit(1, "npc")
  )
)
savePlotAsImage(
  file = "Graphs/immune_Transcriptom_Up.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Stress_Intersections.T_Up <- supertest(x = Gene_Table_List.T_S.Up, n = 13986, degree = seq(Gene_Table_List.T_S.Up)[-1])
Stress_Summ.T_Up <- SuperExactTest:::summary.msets(Stress_Intersections.T_Up)
plot.euler(
  Stress_Pops_Venn.T_Up,
  legend = list(fontsize = 20),
  quantities = list(fontsize = 18, labels = quantlabel(Stress_Intersections.T_Up)),
  main = list(
    label = "Up-regulated Candidate Genes from Stress-Selection Studies",
    fontsize = 20,
    x = unit(0.63, "npc"),
    y = unit(0.15, "npc")
  )
)
savePlotAsImage(
  file = "Graphs/stress_Transcriptom_Up.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Longevdat.T_Up <- Longev_Summ.T_Up$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH"))
View(Longevdat.T_Up)
Longevdat.T_Up <- sort_by.data.frame(x = Longevdat.T_Up, y = -Longevdat.T_Up$FE)
Longevdat.T_Up <- sort_by.data.frame(x = Longevdat.T_Up, y = -Longevdat.T_Up$Degree)
Longevdat.T_Up <- Longevdat.T_Up[Longevdat.T_Up$Observed.Overlap>0,]
Longevdat.T_Up$Elements <- "non-significant"
write_excel_csv(Longevdat.T_Up,
                file = "Data/longevity_pops_Transcriptom_Up.csv",
                append = FALSE,
                escape = "none")

Immundat.T_Up <- Immun_Summ.T_Up$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH")) %>%  filter(P.adj < 0.05)
View(Immundat.T_Up)
Immundat.T_Up <- sort_by.data.frame(x = Immundat.T_Up, y = -Immundat.T_Up$FE)
Immundat.T_Up <- sort_by.data.frame(x = Immundat.T_Up, y = -Immundat.T_Up$Degree)
Immundat.T_Up <- Immundat.T_Up[Immundat.T_Up$Observed.Overlap>0,]
Immundat.T_Up$Elements <- unlist(map(.x = strsplit(x = Immundat.T_Up$Elements, split = ", "), .f = compose(str_flatten_comma, EG2FB)))
write_excel_csv(Immundat.T_Up,
                file = "Data/immunity_pops_Transcriptom_Up.csv",
                append = FALSE,
                escape = "none")

Stressdat.T_Up <- Stress_Summ.T_Up$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH")) %>%  filter(P.adj < 0.05)
View(Stressdat.T_Up)
Stressdat.T_Up <- sort_by.data.frame(x = Stressdat.T_Up, y = -Stressdat.T_Up$FE)
Stressdat.T_Up <- sort_by.data.frame(x = Stressdat.T_Up, y = -Stressdat.T_Up$Degree)
Stressdat.T_Up <- Stressdat.T_Up[Stressdat.T_Up$Observed.Overlap>0,]
Stressdat.T_Up$Elements <- unlist(map(.x = strsplit(x = Stressdat.T_Up$Elements, split = ", "), .f = compose(str_flatten_comma, EG2FB)))
write_excel_csv(Stressdat.T_Up,
                file = "Data/stress_pops_Transcriptom_Up.csv",
                append = FALSE,
                escape = "none")

############################################################### Down-regulated ###############################################################
########################### Gene_Table edit ###########################
library(readxl)
library(writexl)

Gene_Table <- as.data.frame(read_excel("Data/Gene_Table.xlsx"))
Gene_Table.T_Down <- Gene_Table[Gene_Table$Sequencing=="Transcriptomic" & Gene_Table$`Direction (Transcriptomics)`=="Down",]
FB_invs <- as.data.frame(read_tsv("Scripts/misc/withdrawns.txt"))
Gene_Table.T_Down <-Gene_Table.T_Down[!Gene_Table.T_Down$`Candidate FB IDs` %fin% FB_invs$FBIDs,]
rm(FB_invs)

Multi_Pops.T_Down <- Gene_Table.T_Down[, 1:3]
rownames(Multi_Pops.T_Down) <- rownames(Gene_Table.T_Down)
Multi_Pops.T_Down <- unique(Multi_Pops.T_Down)
Gene_Table.T_Down <- Gene_Table.T_Down[rownames(Multi_Pops.T_Down), ]
rm(Multi_Pops.T_Down)

Gene_Table_EG.T_Down <- Gene_Table.T_Down
Gene_Table_EG.T_Down$`Candidate FB IDs` <- FB2EG(Gene_Table.T_Down$`Candidate FB IDs`)

all_pops.T_Down <- Unique(Gene_Table.T_Down$`Populations Used`)
Gene_Table_List.T_Down <- list()
for (i in cli_progress_along(seq(all_pops.T_Down))) {
  Gene_Table_List.T_Down[[i]] <- Gene_Table_EG.T_Down$`Candidate FB IDs`[Gene_Table_EG.T_Down$`Populations Used` == all_pops.T_Down[i]]
}
names(Gene_Table_List.T_Down) <- all_pops.T_Down

########################### Selection Type Overlaps ###########################
library(rstudioapi)
library(SuperExactTest)

eulerr::eulerr_options(eulerr.plot_opts)

Longev_Table.T_Down <- Unique(Gene_Table_EG.T_Down$`Candidate FB IDs`[Gene_Table_EG.T_Down$`Selection Type` == "Longevity"])
Immun_Table.T_Down <- Unique(Gene_Table_EG.T_Down$`Candidate FB IDs`[Gene_Table_EG.T_Down$`Selection Type` == "Immunity"])
Stress_Table.T_Down <- Unique(Gene_Table_EG.T_Down$`Candidate FB IDs`[Gene_Table_EG.T_Down$`Selection Type` == "Stress"])

Gene_Lists.T_Down <- list(Longev_Table.T_Down, Immun_Table.T_Down, Stress_Table.T_Down)
names(Gene_Lists.T_Down) <- c("Longevity", "Immunity", "Stress")

Longev_Immun.T_Down <- Gene_Lists.T_Down[c("Longevity", "Immunity")]
Longev_Immun_Venn.T_Down <- venn.list(Longev_Immun.T_Down)
Longev_Immun_Intersections.T_Down <- supertest(x = Longev_Immun.T_Down, n = 13986, degree = seq(Longev_Immun.T_Down)[-1])
Longev_Immun_Summ.T_Down <- SuperExactTest:::summary.msets(Longev_Immun_Intersections.T_Down)
plot.venn(
  Longev_Immun_Venn.T_Down,
  labels = list(fontsize = 19),
  quantities = list(fontsize = 19, labels = quantlabel.venn(Longev_Immun_Intersections.T_Down)),
  main = list(
    label = "Down-regulated Candidate Genes from Longevity & Immunity Selection", fontsize = 20)
)
savePlotAsImage(
  file = "Graphs/Longev_Immun_Transcriptom_Down.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Longev_Immundat.T_Down <- Longev_Immun_Summ.T_Down$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH"))
View(Longev_Immundat.T_Down)
Longev_Immundat.T_Down$Elements <- unlist(map(.x = strsplit(x = Longev_Immundat.T_Down$Elements, split = ", "), .f = compose(str_flatten_comma, EG2FB)))
write_excel_csv(Longev_Immundat.T_Down,
                file = "Data/Longev_Immun_Transcriptom_Down.csv",
                append = FALSE,
                escape = "none")

Gene_Venn.T_Down <- venn.list(Gene_Lists.T_Down[c("Longevity", "Immunity", "Stress")])
Selection_Intersections.T_Down <- supertest(x = Gene_Lists.T_Down, n = 13986, degree = seq(Gene_Lists.T_Down)[-1])
Selection_Summ.T_Down <- SuperExactTest:::summary.msets(Selection_Intersections.T_Down)
plot.venn(
  Gene_Venn.T_Down,
  labels = list(fontsize = 19),
  quantities = list(fontsize = 19, labels = quantlabel.venn(Selection_Intersections.T_Down)),
  main = list(
    label = "Down-regulated Candidate Genes by Selection Type",
    fontsize = 20
  )
)
savePlotAsImage(
  file = "Graphs/sel_type_Transcriptom_Down.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Selectiondat.T_Down <- Selection_Summ.T_Down$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH")) %>%  filter(P.adj < 0.05)
View(Selectiondat.T_Down)
Selectiondat.T_Down <- sort_by.data.frame(x = Selectiondat.T_Down, y = -Selectiondat.T_Down$FE)
Selectiondat.T_Down <- sort_by.data.frame(x = Selectiondat.T_Down, y = -Selectiondat.T_Down$Degree)
Selectiondat.T_Down$Elements <- unlist(map(.x = strsplit(x = Selectiondat.T_Down$Elements, split = ", "), .f = compose(str_flatten_comma, EG2FB)))
write_excel_csv(Selectiondat.T_Down,
                file = "Data/sel_type_Transcriptom_Down.csv",
                append = FALSE,
                escape = "none")

gp <- get.gpar()
gp$fontfamily <- "serif"
LIS_obj.T_Down <- supertest(Gene_Lists.T_Down, n=13986, degree = seq(Gene_Lists.T_Down)[-1])
LIS_summ.T_Down <- SuperExactTest:::summary.msets(LIS_obj.T_Down)
LIS_obj.T_Down_adj <- LIS_obj.T_Down; LIS_obj.T_Down_adj$P.value <- p.adjust(LIS_obj.T_Down_adj$P.value, method = "BH")
inject(plot(
    LIS_obj.T_Down_adj,
    title = "Longevity, Immunity, & Stress Down-regulated Gene Overlaps",
    minMinusLog10PValue = 1,
    maxMinusLog10PValue = 4,
    !!!SuperExactTest.plot_opts
)); rm(LIS_obj.T_Down_adj)
grab <- grid.grab()
grab$gp <- inject(gpar(!!!gp))
grab[[5]][[1]][["y"]] <- grab[[5]][[1]][["y"]] * 1.5
grid.newpage()
grid.draw(grab)
rm(gp,grab)
savePlotAsImage(
  file = "Graphs/LIS_Transcriptom_Down.svg",
  format = "svg",
  width = 1200,
  height = 742
)

########################### Population Overlaps ###########################
library(readxl)
library(rstudioapi)

Longev_Pops.T_Down <- Unique(Gene_Table.T_Down$`Populations Used`[Gene_Table.T_Down$`Selection Type` == "Longevity"])
Immun_Pops.T_Down <- Unique(Gene_Table.T_Down$`Populations Used`[Gene_Table.T_Down$`Selection Type` == "Immunity"])
Stress_Pops.T_Down <- Unique(Gene_Table.T_Down$`Populations Used`[Gene_Table.T_Down$`Selection Type` == "Stress"])
Gene_Table_List.T_L.Down <-  Gene_Table_List.T_Down[Longev_Pops.T_Down]
Gene_Table_List.T_I.Down <-  Gene_Table_List.T_Down[Immun_Pops.T_Down]
Gene_Table_List.T_S.Down <-  Gene_Table_List.T_Down[Stress_Pops.T_Down]
Longev_Pops_Venn.T_Down <- euler.list(Gene_Table_List.T_L.Down)
Immun_Pops_Venn.T_Down <- euler.list(Gene_Table_List.T_I.Down)
Stress_Pops_Venn.T_Down <- euler.list(Gene_Table_List.T_S.Down)

Longev_Intersections.T_Down <- supertest(x = Gene_Table_List.T_L.Down, n = 13986, degree = seq(Gene_Table_List.T_L.Down)[-1])
Longev_Summ.T_Down <- SuperExactTest:::summary.msets(Longev_Intersections.T_Down)
plot.euler(
  Longev_Pops_Venn.T_Down,
  legend = list(fontsize = 20),
  quantities = list(fontsize = 18, labels = quantlabel(Longev_Intersections.T_Down)),
  main = list(
    label = "Down-regulated Candidate Genes from Longevity-Selection Studies",
    fontsize = 20,
    x = unit(0.64, "npc")
  )
)
savePlotAsImage(
  file = "Graphs/lifespan_Transcriptom_Down.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Immun_Intersections.T_Down <- supertest(x = Gene_Table_List.T_I.Down, n = 13986, degree = seq(Gene_Table_List.T_I.Down)[-1])
Immun_Summ.T_Down <- SuperExactTest:::summary.msets(Immun_Intersections.T_Down)
plot.euler(
  Immun_Pops_Venn.T_Down,
  legend = list(fontsize = 20),
  quantities = list(fontsize = 18, labels = quantlabel(Immun_Intersections.T_Down)),
  main = list(
    label = "Down-regulated Candidate Genes from Immunity-Selection Studies",
    fontsize = 20,
    x = unit(0.64, "npc"),
    y = unit(1, "npc")
  )
)
savePlotAsImage(
  file = "Graphs/immune_Transcriptom_Down.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Stress_Intersections.T_Down <- supertest(x = Gene_Table_List.T_S.Down, n = 13986, degree = seq(Gene_Table_List.T_S.Down)[-1])
Stress_Summ.T_Down <- SuperExactTest:::summary.msets(Stress_Intersections.T_Down)
plot.euler(
  Stress_Pops_Venn.T_Down,
  legend = list(fontsize = 20),
  quantities = list(fontsize = 18, labels = quantlabel(Stress_Intersections.T_Down)),
  main = list(
    label = "Down-regulated Candidate Genes from Stress-Selection Studies",
    fontsize = 20,
    x = unit(0.6, "npc"),
    y = unit(0.1, "npc")
  )
)
savePlotAsImage(
  file = "Graphs/stress_Transcriptom_Down.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Longevdat.T_Down <- Longev_Summ.T_Down$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH")) %>%  filter(P.adj < 0.05)
View(Longevdat.T_Down)
Longevdat.T_Down <- sort_by.data.frame(x = Longevdat.T_Down, y = -Longevdat.T_Down$FE)
Longevdat.T_Down <- sort_by.data.frame(x = Longevdat.T_Down, y = -Longevdat.T_Down$Degree)
Longevdat.T_Down <- Longevdat.T_Down[Longevdat.T_Down$Observed.Overlap>0,]
Longevdat.T_Down$Elements <- unlist(map(.x = strsplit(x = Longevdat.T_Down$Elements, split = ", "), .f = compose(str_flatten_comma, EG2FB)))
write_excel_csv(Longevdat.T_Down,
                file = "Data/longevity_pops_Transcriptom_Down.csv",
                append = FALSE,
                escape = "none")

Immundat.T_Down <- Immun_Summ.T_Down$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH"))
View(Immundat.T_Down)
Immundat.T_Down <- sort_by.data.frame(x = Immundat.T_Down, y = -Immundat.T_Down$FE)
Immundat.T_Down <- sort_by.data.frame(x = Immundat.T_Down, y = -Immundat.T_Down$Degree)
Immundat.T_Down <- Immundat.T_Down[Immundat.T_Down$Observed.Overlap>0,]
Immundat.T_Down$Elements <- "non-significant"
write_excel_csv(Immundat.T_Down,
                file = "Data/immunity_pops_Transcriptom_Down.csv",
                append = FALSE,
                escape = "none")

Stressdat.T_Down <- Stress_Summ.T_Down$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH")) %>%  filter(P.adj < 0.05)
View(Stressdat.T_Down)
Stressdat.T_Down <- sort_by.data.frame(x = Stressdat.T_Down, y = -Stressdat.T_Down$FE)
Stressdat.T_Down <- sort_by.data.frame(x = Stressdat.T_Down, y = -Stressdat.T_Down$Degree)
Stressdat.T_Down <- Stressdat.T_Down[Stressdat.T_Down$Observed.Overlap>0,]
Stressdat.T_Down$Elements <- unlist(map(.x = strsplit(x = Stressdat.T_Down$Elements, split = ", "), .f = compose(str_flatten_comma, EG2FB)))
write_excel_csv(Stressdat.T_Down,
                file = "Data/stress_pops_Transcriptom_Down.csv",
                append = FALSE,
                escape = "none")
save.image()
