quantlabel <- function(Intersections) {
  otab <- Intersections$P.value
  names(otab) <- deBarcode(
    barcode = names(otab),
    setnames = Intersections$set.names,
    collapse = "&"
  )
  otab <- p.adjust(otab, method = "BH")
  otab2 <- otab[which(otab < 0.05)]
  sigsets <- names(otab2)
  foo <- plot.euler(euler.list(Intersections$x), quantities = TRUE)
  quantlabels <- rep("", length(foo$data$fitted.values))
  quantlabels[foo$data$centers$id] <- as.character(foo$data$centers$quantities)
  allsets <- rownames(foo$data$centers)
  sigs <- foo$data$centers$id[allsets %fin% sigsets]
  if (n_distinct(symnum(otab[allsets[allsets %fin% sigsets]], corr = FALSE, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", ""))) > 1) {
    quantlabels[sigs] <- paste2(quantlabels[sigs], symnum(otab[allsets[allsets %fin% sigsets]], corr = FALSE, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "")))
  } else {
    quantlabels[sigs] <- paste2(quantlabels[sigs], "*")
  }
  return(quantlabels)
}
quantlabel.venn <- function(Intersections) {
  otab <- Intersections$P.value
  names(otab) <- deBarcode(
    barcode = names(otab),
    setnames = Intersections$set.names,
    collapse = "&"
  )
  otab <- p.adjust(otab, method = "BH")
  otab2 <- otab[which(otab < 0.05)]
  sigsets <- names(otab2)
  foo <- plot.venn(venn.list(Intersections$x), quantities = TRUE)
  quantlabels <- rep("", length(foo$data$fitted.values))
  quantlabels[foo$data$centers$id] <- as.character(foo$data$centers$quantities)
  allsets <- rownames(foo$data$centers)
  sigs <- foo$data$centers$id[allsets %fin% sigsets]
  if (n_distinct(symnum(otab[allsets[allsets %fin% sigsets]], corr = FALSE, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", ""))) > 1) {
    quantlabels[sigs] <- paste2(quantlabels[sigs], symnum(otab[allsets[allsets %fin% sigsets]], corr = FALSE, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "")))
  } else {
    quantlabels[sigs] <- paste2(quantlabels[sigs], "*")
  }
  return(quantlabels)
}

is.longevity <- function(x) {
  library(org.Dm.eg.db)
  if (Any(startsWith(x, "FBgn"))) {
    x <- FB2EG(x)
  }
  .EG2GO_list <- goProfiles::GOTermsList(x, orgPkg = "org.Dm.eg.db")
  .EG2GO_list <- map(.x = .EG2GO_list, .f = unname)
  .EG2GO_list <- uniques(.EG2GO_list)
  .EG2GO_list <- map(.x = .EG2GO_list, .f = function(x) x[!x %fin% BPCCMF])
  .EG2GO_list <- rmNULL(.EG2GO_list)
  cats_List <- as.List(.EG2GO_list)
  longevity <- c("GO:0008340", "GO:0090398", "GO:0012501", "GO:0040007")
  longev_GOs <- unlist2(map(.x = getGOOFFSPRING(longevity), .f = "Offspring"))
  longev_GOs <- unique(c(longevity, longev_GOs))
  longev_GOs <- setdiff(longev_GOs, c("GO:0045087", unlist2(map(.x = getGOOFFSPRING("GO:0045087"), .f = "Offspring")), "GO:0006950", unlist2(map(.x = getGOOFFSPRING("GO:0006950"), .f = "Offspring"))))
  gene_res <- cats_List@partitioning@NAMES[unique(togroup(x = cats_List@partitioning, j = which(cats_List@unlistData %fin% longev_GOs)))]
  KEGG_genes <- (KEGGREST::keggGet("dme04213")[[1]]$GENE)[(1:56)*2]
  KEGG_genes <- SYMBOL2Gene(str_split_i(string = KEGG_genes, pattern = ";", i = 1))
  gene_res <- union(gene_res, intersect(x, KEGG_genes))
  genage_models <- as.data.frame(read_csv("Data/genage_models_export.csv"))
  genage_models <- genage_models[genage_models$`Lifespan Effect` %fin% c("increase", "decrease"),]
  gene_res <- union(gene_res, intersect(x, genage_models$`Entrez Gene ID`))
  answer <- unique(mixedsort(gene_res))
  return(answer)
}
is.immunity <- function(x) {
  library(org.Dm.eg.db)
  if (Any(startsWith(x, "FBgn"))) {
    x <- FB2EG(x)
  }
  .EG2GO_list <- goProfiles::GOTermsList(x, orgPkg = "org.Dm.eg.db")
  .EG2GO_list <- map(.x = .EG2GO_list, .f = unname)
  .EG2GO_list <- uniques(.EG2GO_list)
  .EG2GO_list <- map(.x = .EG2GO_list, .f = function(x) x[!x %fin% BPCCMF])
  .EG2GO_list <- rmNULL(.EG2GO_list)
  cats_List <- as.List(.EG2GO_list)
  immun_GOs <- unlist2(map(.x = getGOOFFSPRING(c("GO:0006955", "GO:0006952")), .f = "Offspring"))
  immun_GOs <- unique(c("GO:0045087", immun_GOs))
  immun_GOs <- intersect(c("GO:0002376", "GO:0006952", "GO:0051707", unlist2(map(.x = getGOOFFSPRING(c("GO:0002376", "GO:0006952", "GO:0051707")), .f = "Offspring"))), c(immun_GOs, unlist2(map(.x = getGOANCESTORS("GO:0045087"), .f = "Ancestor"))))
  gene_res <- cats_List@partitioning@NAMES[unique(togroup(x = cats_List@partitioning, j = which(cats_List@unlistData %fin% immun_GOs)))]
  answer <- unique(mixedsort(gene_res))
  return(answer)
}
is.stress <- function(x) {
  library(org.Dm.eg.db)
  if (Any(startsWith(x, "FBgn"))) {
    x <- FB2EG(x)
  }
  .EG2GO_list <- goProfiles::GOTermsList(x, orgPkg = "org.Dm.eg.db")
  .EG2GO_list <- map(.x = .EG2GO_list, .f = unname)
  .EG2GO_list <- uniques(.EG2GO_list)
  .EG2GO_list <- map(.x = .EG2GO_list, .f = function(x) x[!x %fin% BPCCMF])
  .EG2GO_list <- rmNULL(.EG2GO_list)
  cats_List <- as.List(.EG2GO_list)
  stress_GOs <- unlist2(map(.x = getGOOFFSPRING("GO:0006950"), .f = "Offspring"))
  stress_GOs <- unique(c("GO:0006950", stress_GOs))
  stress_GOs <- setdiff(stress_GOs, c(intersect(c("GO:0006952", unlist2(map(.x = getGOOFFSPRING("GO:0006952"), .f = "Offspring"))), c("GO:0045087", unlist2(map(.x = getGOOFFSPRING(c("GO:0006955", "GO:0006952")), .f = "Offspring")), unlist2(map(.x = getGOANCESTORS("GO:0045087"), .f = "Ancestor")))), "GO:0008340", "GO:0090398", "GO:0012501", "GO:0040007", unlist2(map(.x = getGOOFFSPRING(c("GO:0008340", "GO:0090398", "GO:0012501", "GO:0040007")), .f = "Offspring"))))
  gene_res <- cats_List@partitioning@NAMES[unique(togroup(x = cats_List@partitioning, j = which(cats_List@unlistData %fin% stress_GOs)))]
  answer <- unique(mixedsort(gene_res))
  return(answer)
}

save.image()
