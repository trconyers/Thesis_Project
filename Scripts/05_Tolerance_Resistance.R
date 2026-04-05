library(readxl)
library(writexl)
library(rstudioapi)
library(SuperExactTest)

Gene.Table <- as.data.frame(read_excel("Data/Gene Table.xlsx"))
Gene.Table$`Populations Used`[Gene.Table$`Populations Used` %fin% c("B vs. O (new)", "B vs. O (old)")] <- "B-type vs. O-type"
Multi_Pops <- Gene.Table[, 1:3]
rownames(Multi_Pops) <- rownames(Gene.Table)
Multi_Pops <- unique(Multi_Pops)
Gene.Table <- Gene.Table[rownames(Multi_Pops), ]
rm(Multi_Pops)
Gene.Table_EG <- Gene.Table
Gene.Table_EG$`Candidate FB IDs` <- FB2EG(Gene.Table$`Candidate FB IDs`)

Longev.Table <- Unique(Gene.Table_EG$`Candidate FB IDs`[Gene.Table_EG$`Selection Type` == "Longevity"])
Immun.Table <- Unique(Gene.Table_EG$`Candidate FB IDs`[Gene.Table_EG$`Selection Type` == "Immunity"])
Stress.Table <- Unique(Gene.Table_EG$`Candidate FB IDs`[Gene.Table_EG$`Selection Type` == "Stress"])

Howick_data <- as.data.frame(read_tsv("Immunity/Tolerance.v.Resistance/Howick_2017/Howick_SNPs.tsv"))
Howick_Tolerance <- Howick_data[Howick_data$immun_type=="tolerance",]
tolerance <- gene_mapper(chromosome = Howick_Tolerance$chromosome, start = Howick_Tolerance$start, end = Howick_Tolerance$end)
Howick_Resistance <- Howick_data[Howick_data$immun_type=="resistance",]
resistance <- gene_mapper(chromosome = Howick_Resistance$chromosome, start = Howick_Resistance$start, end = Howick_Resistance$end)

Tol_Res_list <- list(
  Longevity_Immunity = SuperExactTest::intersect(Longev.Table, Immun.Table),
  Longevity_Stress = SuperExactTest::intersect(Longev.Table, Stress.Table),
  Longevity = Longev.Table,
  Immunity = Immun.Table,
  Stress = Stress.Table,
  Tolerance = tolerance,
  Resistance = resistance
)
Tol.Res_obj <- supertest(x = Tol_Res_list, n = 17871, degree = seq(Tol_Res_list)[-1])

TolRes_plotter <- function(y) {
  keep <- str_detect(string = deBarcode(names(y$overlap.sizes), y$set.names), pattern = "ance")
  y$overlap.sizes <- y$overlap.sizes[keep]
  y$overlap.expected <- y$overlap.expected[keep]
  y$P.value <- y$P.value[keep]
  y$P.value <- p.adjust(y$P.value, method = "BH")
  .SuperExactTest.plot_opts <- SuperExactTest.plot_opts[-3]
  .SuperExactTest.plot_opts$margin[[3]] <- 3.5
  .SuperExactTest.plot_opts$color.scale.pos[[2]] <- 0.7
  .SuperExactTest.plot_opts$color.scale.pos[[3]] <- 0.95
  inject(
    SuperExactTest:::plot.msets.landscape(
      x = y,
      degree = 2,
      sort.by = "p-value",
      minMinusLog10PValue = 1,
      maxMinusLog10PValue = 4,
      title = "Tolerance and Resistance Overlaps",!!!.SuperExactTest.plot_opts
    )
  )
}

gp <- get.gpar()
gp$fontfamily <- "serif"
TolRes_plotter(Tol.Res_obj)
grab <- grid.grab()
grab$gp <- inject(gpar(!!!gp))
grab[[5]][[1]][["y"]] <- grab[[5]][[1]][["y"]] * 1.5
grid.newpage()
grid.draw(grab)
rm(gp,grab)
savePlotAsImage(
  file = "Graphs/def_type.svg",
  format = "svg",
  width = 1200,
  height = 742
)

Defense_Intersections <- SuperExactTest:::summary.msets(Tol.Res_obj)
keep <- names(str_match(string = deBarcode(Defense_Intersections$Barcode, Defense_Intersections$set.names), pattern = "ance"))
Defensedat <- Defense_Intersections$Table[keep,] %>%  mutate(P.adj = p.adjust(P.value, method = "BH"))
rm(keep)
View(Defensedat)
Defensedat <- sort_by.data.frame(x = Defensedat, y = -Defensedat$FE)
Defensedat <- sort_by.data.frame(x = Defensedat, y = Defensedat$P.value)
Defensedat <- Defensedat[Defensedat$Degree==2,]
Defensedat <- Defensedat[Defensedat$Observed.Overlap>0,-2]
Defensedat$Elements <- unlist(map(.x = strsplit(x = Defensedat$Elements, split = ", "), .f = compose(str_flatten_comma, Gene2SYMBOL)))
Defensedat$Elements[Defensedat$P.adj>=0.05] <- "non-significant"
write_excel_csv(Defensedat,
                file = "Data/def_type.csv",
                append = FALSE,
                escape = "none")
save.image()
