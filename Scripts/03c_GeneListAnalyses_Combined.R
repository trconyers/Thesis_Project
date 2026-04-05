########################### Gene_Table edit ###########################
library(readxl)
library(writexl)

Gene_Table <- as.data.frame(read_excel("Data/Gene_Table.xlsx"))
Gene_Table <- Gene_Table[Gene_Table$`Candidate FB IDs` %fin% protein_coding,]
FB_invs <- as.data.frame(read_tsv("Scripts/misc/withdrawns.txt"))
Gene_Table <-Gene_Table[!Gene_Table$`Candidate FB IDs` %fin% FB_invs$FBIDs,]
Gene_Table$`Populations Used`[Gene_Table$`Populations Used` %fin% c("B vs. O (new)", "B vs. O (old)")] <- "B-type vs. O-type"
rm(FB_invs)

Multi_Pops <- Gene_Table[, 1:3]
rownames(Multi_Pops) <- rownames(Gene_Table)
Multi_Pops <- unique(Multi_Pops)
Gene_Table <- Gene_Table[rownames(Multi_Pops), ]
rm(Multi_Pops)

Gene_Table_EG <- Gene_Table
Gene_Table_EG$`Candidate FB IDs` <- FB2EG(Gene_Table$`Candidate FB IDs`)

all_pops <- Unique(Gene_Table$`Populations Used`)
Gene_Table_List <- list()
for (i in cli_progress_along(seq(all_pops))) {
  Gene_Table_List[[i]] <- Gene_Table_EG$`Candidate FB IDs`[Gene_Table_EG$`Populations Used` == all_pops[i]]
}
names(Gene_Table_List) <- all_pops

########################### Selection Type Overlaps ###########################
library(rstudioapi)
library(SuperExactTest)

eulerr::eulerr_options(eulerr.plot_opts)

Longev_Table <- Unique(Gene_Table_EG$`Candidate FB IDs`[Gene_Table_EG$`Selection Type` == "Longevity"])
Immun_Table <- Unique(Gene_Table_EG$`Candidate FB IDs`[Gene_Table_EG$`Selection Type` == "Immunity"])
Stress_Table <- Unique(Gene_Table_EG$`Candidate FB IDs`[Gene_Table_EG$`Selection Type` == "Stress"])

Gene_Lists <- list(Longev_Table, Immun_Table, Stress_Table)
names(Gene_Lists) <- c("Longevity", "Immunity", "Stress")

Longev_Immun <- Gene_Lists[c("Longevity", "Immunity")]
Longev_Immun_Venn <- venn.list(Longev_Immun)
Longev_Immun_Intersections <- supertest(x = Longev_Immun, n = 13986, degree = seq(Longev_Immun)[-1])
Longev_Immun_Summ <- SuperExactTest:::summary.msets(Longev_Immun_Intersections)
plot.venn(
  Longev_Immun_Venn,
  labels = list(fontsize = 19),
  quantities = list(fontsize = 19, labels = quantlabel.venn(Longev_Immun_Intersections)),
  main = list(
    label = "Candidate Coding Genes from Longevity & Immunity Selection", fontsize = 22.5)
)
savePlotAsImage(
  file = "Graphs/Longev_Immun.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Longev_Immundat <- Longev_Immun_Summ$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH"))
View(Longev_Immundat)
Longev_Immundat$Elements <- str_flatten_comma(Gene2SYMBOL(unlist(strsplit(x = Longev_Immundat$Elements, split = ", "))))
write_excel_csv(Longev_Immundat,
                file = "Data/Longev_Immun.csv",
                append = FALSE,
                escape = "none")

Gene_Venn <- venn.list(Gene_Lists[c("Longevity", "Immunity", "Stress")])
Selection_Intersections <- supertest(x = Gene_Lists, n = 13986, degree = seq(Gene_Lists)[-1])
Selection_Summ <- SuperExactTest:::summary.msets(Selection_Intersections)
plot.venn(
  Gene_Venn,
  labels = list(fontsize = 19),
  quantities = list(fontsize = 19, labels = quantlabel.venn(Selection_Intersections)),
  main = list(
    label = "Candidate Coding Genes by Selection Type",
    fontsize = 22.5
  )
)
savePlotAsImage(
  file = "Graphs/sel_type.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Selectiondat <- Selection_Summ$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH")) %>%  filter(P.adj < 0.05)
View(Selectiondat)
Selectiondat <- sort_by.data.frame(x = Selectiondat, y = -Selectiondat$FE)
Selectiondat <- sort_by.data.frame(x = Selectiondat, y = -Selectiondat$Degree)
Selectiondat$Elements <- unlist(map(.x = strsplit(x = Selectiondat$Elements, split = ", "), .f = compose(str_flatten_comma, Gene2SYMBOL)))
write_excel_csv(Selectiondat,
                file = "Data/sel_type.csv",
                append = FALSE,
                escape = "none")

gp <- get.gpar()
gp$fontfamily <- "serif"
LIS_obj <- supertest(Gene_Lists, n=13986, degree = seq(Gene_Lists)[-1])
LIS_summ <- SuperExactTest:::summary.msets(LIS_obj)
LIS_obj_adj <- LIS_obj; LIS_obj_adj$P.value <- p.adjust(LIS_obj_adj$P.value, method = "BH")
inject(plot(
    LIS_obj_adj,
    title = "Longevity, Immunity, & Stress Gene Overlaps",
    minMinusLog10PValue = 1,
    maxMinusLog10PValue = 4,
    !!!SuperExactTest.plot_opts
)); rm(LIS_obj_adj)
grab <- grid.grab()
grab$gp <- inject(gpar(!!!gp))
grab[[5]][[1]][["y"]] <- grab[[5]][[1]][["y"]] * 1.5
grid.newpage()
grid.draw(grab)
rm(gp,grab)
savePlotAsImage(
  file = "Graphs/LIS.svg",
  format = "svg",
  width = 1200,
  height = 742
)

########################### Population Overlaps ###########################
library(readxl)
library(rstudioapi)

Longev_Pops <- Unique(Gene_Table$`Populations Used`[Gene_Table$`Selection Type` == "Longevity"])
Immun_Pops <- Unique(Gene_Table$`Populations Used`[Gene_Table$`Selection Type` == "Immunity"])
Stress_Pops <- Unique(Gene_Table$`Populations Used`[Gene_Table$`Selection Type` == "Stress"])
Gene_Table_List_L <-  Gene_Table_List[Longev_Pops]
Gene_Table_List_I <-  Gene_Table_List[Immun_Pops]
Gene_Table_List_S <-  Gene_Table_List[Stress_Pops]
Longev_Pops_Venn <- euler.list(Gene_Table_List_L)
Immun_Pops_Venn <- euler.list(Gene_Table_List_I)
Stress_Pops_Venn <- euler.list(Gene_Table_List_S)

Longev_Intersections <- supertest(x = Gene_Table_List_L, n = 13986, degree = seq(Gene_Table_List_L)[-1])
Longev_Summ <- SuperExactTest:::summary.msets(Longev_Intersections)
plot.euler(
  Longev_Pops_Venn,
  legend = list(fontsize = 20),
  quantities = list(fontsize = 18, labels = quantlabel(Longev_Intersections)),
  main = list(
    label = "Candidate Coding Genes from Longevity-Selection Studies",
    fontsize = 22.5,
    x = unit(0.66, "npc")
  )
)
savePlotAsImage(
  file = "Graphs/lifespan.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Immun_Intersections <- supertest(x = Gene_Table_List_I, n = 13986, degree = seq(Gene_Table_List_I)[-1])
Immun_Summ <- SuperExactTest:::summary.msets(Immun_Intersections)
plot.euler(
  Immun_Pops_Venn,
  legend = list(fontsize = 20),
  quantities = list(fontsize = 18, labels = quantlabel(Immun_Intersections)),
  main = list(
    label = "Candidate Coding Genes from Immunity-Selection Studies",
    fontsize = 22.5,
    x = unit(0.64, "npc")
  )
)
savePlotAsImage(
  file = "Graphs/immune.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Stress_Intersections <- supertest(x = Gene_Table_List_S, n = 13986, degree = seq(Gene_Table_List_S)[-1])
Stress_Summ <- SuperExactTest:::summary.msets(Stress_Intersections)
plot.euler(
  Stress_Pops_Venn,
  legend = list(fontsize = 20),
  quantities = list(fontsize = 18, labels = quantlabel(Stress_Intersections)),
  main = list(
    label = "Candidate Coding Genes from Stress-Selection Studies",
    fontsize = 22.5,
    x = unit(0.73, "npc"),
    y = unit(1, "npc")
  )
)
savePlotAsImage(
  file = "Graphs/stress.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Longevdat <- Longev_Summ$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH")) %>%  filter(P.adj < 0.05)
View(Longevdat)
Longevdat <- sort_by.data.frame(x = Longevdat, y = -Longevdat$FE)
Longevdat <- sort_by.data.frame(x = Longevdat, y = -Longevdat$Degree)
Longevdat <- Longevdat[Longevdat$Observed.Overlap>0,]
Longevdat$Elements <- unlist(map(.x = strsplit(x = Longevdat$Elements, split = ", "), .f = compose(str_flatten_comma, Gene2SYMBOL)))
write_excel_csv(Longevdat,
                file = "Data/longevity_pops.csv",
                append = FALSE,
                escape = "none")

Immundat <- Immun_Summ$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH")) %>%  filter(P.adj < 0.05)
View(Immundat)
Immundat <- sort_by.data.frame(x = Immundat, y = -Immundat$FE)
Immundat <- sort_by.data.frame(x = Immundat, y = -Immundat$Degree)
Immundat <- Immundat[Immundat$Observed.Overlap>0,]
Immundat$Elements <- unlist(map(.x = strsplit(x = Immundat$Elements, split = ", "), .f = compose(str_flatten_comma, Gene2SYMBOL)))
write_excel_csv(Immundat,
                file = "Data/immunity_pops.csv",
                append = FALSE,
                escape = "none")

Stressdat <- Stress_Summ$Table  %>%  mutate(P.adj = p.adjust(P.value, method = "BH")) %>%  filter(P.adj < 0.05)
View(Stressdat)
Stressdat <- sort_by.data.frame(x = Stressdat, y = -Stressdat$FE)
Stressdat <- sort_by.data.frame(x = Stressdat, y = -Stressdat$Degree)
Stressdat <- Stressdat[Stressdat$Observed.Overlap>0,]
Stressdat$Elements <- unlist(map(.x = strsplit(x = Stressdat$Elements, split = ", "), .f = compose(str_flatten_comma, Gene2SYMBOL)))
write_excel_csv(Stressdat,
                file = "Data/stress_pops.csv",
                append = FALSE,
                escape = "none")
save.image()
