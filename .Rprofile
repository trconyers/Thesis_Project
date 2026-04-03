if (interactive()) {
  suppressMessages(require(usethis))
}

# Configure BiocManager to use Posit Package Manager:
options(BioC_mirror = "https://packagemanager.posit.co/bioconductor/latest")
options(BIOCONDUCTOR_CONFIG_FILE = "https://packagemanager.posit.co/bioconductor/latest/config.yaml")

# Configure a CRAN snapshot compatible with Bioconductor 3.22:
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/latest"))
options(width = 79, digits = 7)
options(show.signif.stars = FALSE)
setHook(packageEvent("grDevices", "onLoad"), function(...)
  grDevices::ps.options(horizontal = FALSE))
set.seed(1234)

#######ENVIRONMENT########
.bioc <- c("annotate", "clusterProfiler", "GEOquery", "GO.db", "GOstats", "topGO", "Biobase", "limma")
.init_pkgs <- c("base", "boot", "class", "cluster", "codetools", "compiler", "datasets", "foreign", "ggplot2", "graphics", "grDevices", "grid", "KernSmooth", "lattice", "MASS", "Matrix", "methods", "mgcv", "nlme", "nnet", "parallel", "rpart", "spatial", "splines", "stats", "stats4", "survival", "tcltk", "tools", "usethis", "utils", "XML")
.my_pkgs <- c("AcidBase", "AcidDevTools", "AcidGenerics", "AcidGenomes", "AcidTest", "admisc", "annotate", "assertthat", "basejump", "Biobase", "BiocGenerics", "BiocManager", "BiocVersion", "bit", "blob", "cachem", "callr", "checkmate", "cli", "clipr", "collapse", "colorspace", "cpp11", "crayon", "data.table", "DBI", "devtools", "dplyr", "easy.utils", "fansi", "fastmap", "fastmatch", "fgsea", "generics", "GenomicFeatures", "GetoptLong", "goalie", "gplots", "graph", "gtools", "hms", "installr", "jamba", "lifecycle", "memoise", "mulea", "pbapply", "pillar", "pipette", "pkgbuild", "pkgconfig", "pkgdown", "pkgload", "plogr", "plyr", "prettyunits", "progress", "purrr", "R.methodsS3", "R.oo", "R.utils", "R6", "RCurl", "readr", "readxl", "remotes", "Rfast", "Rfast2", "rlang", "roxygen2", "RSQLite", "rstudioapi", "rtracklayer", "S4Vectors", "seqinr", "statmod", "stringi", "stringr", "SuperExactTest", "tibble", "tidyselect", "tzdb", "utf8", "vctrs", "vroom", "xml2", "XVector")
Sys.setenv(RCMDCHECK_ERROR_ON = "error")
options(defaultPackages = c(.bioc, .init_pkgs, .my_pkgs))
# Might need to do this:
# xfun::pkg_attach(options("defaultPackages")$defaultPackages)

#####CUSTOM FUNCTIONS#####
as.List <- partial(.f = as, Class = "CompressedList")
as.na <- function(x) {
  dplyr::na_if(x = x, y = x)
}
as.string <- function(x)  toString(x)
dataframe <- function(data = NA_integer_, nrow = 1, ncol = 1) {
  compose(as.data.frame.matrix, matrix)(data, nrow, ncol)
}
EG2FB <- function(x) {
  library(org.Dm.eg.db)
  Entrez <- as.list.Bimap(org.Dm.egFLYBASE)
  unlist2(Entrez[as.character(x)])
}
FB2EG <- function(x) {
  library(org.Dm.eg.db)
  FlyBase <- as.list.Bimap(org.Dm.egFLYBASE2EG)
  unlist2(FlyBase[x])
}
Gene2GOs <- function(x) {
  library(org.Dm.eg.db)
  if (Any(startsWith(x, "FBgn"))) {
    x <- FB2EG(x)
  }
  Gene2GOList <- goProfiles::GOTermsList(x, orgPkg = "org.Dm.eg.db")
  Gene2GOList <- map(
    .x = Gene2GOList,
    .f = function(x) x[!x %fin% BPCCMF]
  )
  Gene2GOList <- rmNULL(Gene2GOList)
  res <- unique(unlist2(Gene2GOList))
  return(res)
}
GOs2Gene <- function(y) {
  library(org.Dm.eg.db)
  GOList <- AnnotationDbi::mget(y, org.Dm.egGO2ALLEGS, ifnotfound = NA)
  res <- unique(unlist2(GOList))
  return(res)
}
getClassMethods <- function(x) {
  mths <- attr(methods(class = class(x), all.names = TRUE), "info")
  rownames(mths) <- str_remove(string = rownames(mths), pattern = "-method")
  notsym <- str_starts(string = mths$generic, pattern = "[:letter:]|\\.")
  mths <- mths[notsym, ]
  mths <- sort_by.data.frame(x = mths, y = rownames(mths))
  mths <- sort_by.data.frame(x = mths, y = mths$from)
  View(mths)
  invisible(mths)
}
getFDef <- function(fn, env = character()) {
  if (!missing(env)) {
    getAnywhere(as.character(substitute(fn)))$objs[[paste2("package:", env)]]
  }  else {
    result <- getAnywhere(as.character(substitute(fn)))
    uniq <- which(!duplicated(str_split_i(
      string = result$where,
      pattern = ":",
      i = 2
    )))
    if (length(uniq) > 1) {
      return(result$objs[uniq])
    } else {
      return(result$objs[[uniq]])
    }
  }
}
getGOOFFSPRING <- function(x) {
  if (!is.character(x)) {
    stop("need a character argument")
  }
  if (length(x) == 0) {
    return(list())
  }
  loadNamespace("GO.db")
  MF_offspring <- mget(x, envir = GO.db::GOMFOFFSPRING, ifnotfound = NA)
  BP_offspring <- mget(x, envir = GO.db::GOBPOFFSPRING, ifnotfound = NA)
  CC_offspring <- mget(x, envir = GO.db::GOCCOFFSPRING, ifnotfound = NA)
  lapply(setNames(seq_along(x), x), function(i) {
    xi_offspring <- MF_offspring[[i]]
    if (!identical(xi_offspring, NA)) {
      return(list(Ontology = "MF", Offspring = xi_offspring))
    }
    xi_offspring <- BP_offspring[[i]]
    if (!identical(xi_offspring, NA)) {
      return(list(Ontology = "BP", Offspring = xi_offspring))
    }
    xi_offspring <- CC_offspring[[i]]
    if (!identical(xi_offspring, NA)) {
      return(list(Ontology = "CC", Offspring = xi_offspring))
    }
    list()
  })
}
getGOSYNONYMS <- function(x) {
  library(GO.db)
  GOlist <- as.list.Bimap(GOSYNONYM)
  SYNONYM <- map(.x = GOlist, .f = GOID)
  SYNONYM_t <- list_flip(SYNONYM)
  SYNONYMS <- c(SYNONYM, SYNONYM_t)
  res <- SYNONYMS[x]
  res <- rmNULL(res)
  return(res)
}
GOSynGenes <- function(x, background) {
  GeneGOs <- Gene2GOs(x)
  GeneGOs <- GeneGOs[!GeneGOs %fin% BPCCMF]
  GeneGOs <- unique(c(unlist(getGOSYNONYMS(GeneGOs)), GeneGOs))
  GOGenes <- na.omit(unique(unlist(GOs2Gene(GeneGOs))))
  if (Any(startsWith(background, "FBgn"))) {
    background <- FB2EG(background)
  }
  res <- GOGenes[GOGenes %fin% background & (!GOGenes %fin% x)]
  return(res)
}
importEnv <- function(envname, parent.frame = .GlobalEnv) {
  Env <- as.list.environment(asNamespace(envname))
  Biobase::multiassign(x = names(Env),
                       value = Env,
                       envir = parent.frame)
}
list.flatten <- function(x)  map(.x = x, .f = unlist2)
list_flip <- function(list) {
  revlist <- reverseSplit(list)
  revlist <- map(.x = revlist, .f = mixedsort)
  sorted <- unsplit(revlist, names(revlist))
  revlist <- mixedsort_by(x = revlist, y = sorted)
  return(revlist)
}
map <- function(.x, .f) {
  .fn <- function(.x, .f, ..., .progress = TRUE) {
    purrr:::map_("list", .x, .f, ..., .progress = .progress)
  }
  .fn(.x, .f)
}
mixedrank <- function(x, na.last = TRUE) {
  oo <- mixedorder(x, na.last = na.last)
  ans <- integer(length(oo))
  ans[oo] <- seq_along(oo)
  ans <- ans[selfmatch(x)]
  ans
}
mixedsort_by <- function(x, y, ...) {
  if (length(dim(x)) != 2) {
    x[mixedorder(y, ...)]
  } else {
    if (inherits(y, "formula")) y <- .formula2varlist(y, x)
    if (!is.list(y)) y <- list(y)
    o <- do.call(mixedorder, c(unname(y), list(...)))
    x[o, , drop = FALSE]
  }
}
na.omit <- collapse::na_rm
open_file <- function(path) {
  usethis:::create_directory(fs::path_dir(path))
  fs::file_create(path)
  usethis:::ui_bullets(c(`_` = "Modify {.path {usethis:::pth(path)}}."))
  rstudioapi::navigateToFile(path)
  invisible(path)
}
sample <- function(x, size = 1, replace, prob) {
  .fn <- function(x, size, replace = TRUE, prob = NULL)  {
    if (length(x) == 1L && is.numeric(x) && is.finite(x) && x >= 1) {
      if (missing(size)) {
        size <- x
      }
      sample.int(x, size, replace, prob)
    } else {
      if (missing(size)) {
        size <- length(x)
      }
      x[sample.int(length(x), size, replace, prob)]
    }
  }
  .fn(x, size)
}
str_match <- function (string, pattern, negate = FALSE) {
  string[tidyplus::str_detect2(string, pattern, negate)]
}
which.is.na <- collapse::whichNA
restartR <- function()  .Options$restart()
.First <- function()  cat("\n Welcome to RStudio!\n\n")
.Last <- function()  cat("\n Goodbye!\n\n")
