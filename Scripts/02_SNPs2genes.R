library(GenomicRanges)
library(GenomicFeatures)
library(Rfast)
library(IRanges)
library(readr)

gene_mapper <- function(chromosome, start, end) {
  if (All(end - start == 0)) {
    query <- GPos(
      seqnames = Rle(paste0("chr", chromosome)),
      pos = IRanges(start, end),
      strand = Rle("*")
    )
  } else {
    coords <- cbind(start, end)
    fixed <- rowSort(coords)
    start <- as.integer(fixed[, 1])
    end <- as.integer(fixed[, 2])
    query <- GRanges(
      seqnames = Rle(paste0("chr", chromosome)),
      ranges = IRanges(start, end),
      strand = Rle("*")
    )
  }
  subject <- genes(dmel_txdb)
  subject <- resize(subject, width = (width(subject) + 2000), fix = "center")
  subject <- trim(subject)
  start(subject[as.character(subject@seqnames) == "chrM" & start(subject) < 0]) <- 1
  if (All(width(query) <= Mins(width(subject)))) {
    hits <- mergeByOverlaps(query, subject, type = "within")
  } else {
    hits <- mergeByOverlaps(query, subject, type = "any")
  }
  hits[, 1] <- as.character(hits[, 1])
  hits[, 2] <- as.character(hits[, 2])
  res <- Unique(hits$gene_id)
  return(res)
}

##################### Map coords to genes. #####################
Burke_data <- as.data.frame(read_tsv("Longevity/Genomics/Burke_2010/Burke_SNPs.tsv"))
Burke_data <- Burke_data[Burke_data$L10FET.ACOvCO>4, ]
Burke_genes <- gene_mapper(chromosome = Burke_data$chromosome, start = Burke_data$start, end = Burke_data$end)
Burke <- data.frame(check.names = FALSE,
                    `Candidate FB IDs` = Burke_genes,
                    `Selection Type` = "Longevity",
                    `Populations Used` = "A-type vs. C-type",
                    `Source Articles` = "Burke 2010",
                    `Year` = 2010,
                    `Sequencing` = "Genomic",
                    `Direction (Transcriptomics)` = NA_character_
)
###
Graves_data <- as.data.frame(read_tsv("Longevity/Genomics/Graves_2017/Graves_SNPs.tsv"))
Graves_data <- Graves_data[which(Graves_data$ACO_CO_pval < 5.67e-145 | Graves_data$AO_nCO_pval < 1.64e-166), ]
Graves_genes <- gene_mapper(chromosome = Graves_data$chromosome, start = Graves_data$start, end = Graves_data$end)
Graves <- data.frame(check.names = FALSE,
                     `Candidate FB IDs` = Graves_genes,
                     `Selection Type` = "Longevity",
                     `Populations Used` = "A-type vs. C-type",
                     `Source Articles` = "Graves 2017",
                     `Year` = 2017,
                     `Sequencing` = "Genomic",
                     `Direction (Transcriptomics)` = NA_character_
)
###
Carnes_SNP_data <- as.data.frame(read_tsv("Longevity/Genomics/Carnes_2015/Carnes_SNPs.tsv"))
Carnes_SNP_genes <- gene_mapper(chromosome = Carnes_SNP_data$chromosome, start = Carnes_SNP_data$start, end = Carnes_SNP_data$end)
Carnes_SNPs <- data.frame(check.names = FALSE,
                     `Candidate FB IDs` = Carnes_SNP_genes,
                     `Selection Type` = "Longevity",
                     `Populations Used` = "B vs. O (old)",
                     `Source Articles` = "Carnes 2015",
                     `Year` = 2015,
                     `Sequencing` = "Genomic",
                     `Direction (Transcriptomics)` = NA_character_
)
###
Remolina_data <- as.data.frame(read_tsv("Longevity/Genomics/Remolina_2012/Remolina_SNPs.tsv"))
Remolina_genes <- gene_mapper(chromosome = Remolina_data$chromosome, start = Remolina_data$start, end = Remolina_data$end)
Remolina <- data.frame(check.names = FALSE,
                       `Candidate FB IDs` = Remolina_genes,
                       `Selection Type` = "Longevity",
                       `Populations Used` = "C vs. S (Remolina)",
                       `Source Articles` = "Remolina 2012",
                       `Year` = 2012,
                       `Sequencing` = "Genomic",
                       `Direction (Transcriptomics)` = NA_character_
)
###
Hoedjes_data <- as.data.frame(read_tsv("Longevity/Genomics/Hoedjes_2019/Hoedjes_SNPs.tsv"))
Hoedjes_genes <- gene_mapper(chromosome = Hoedjes_data$chromosome, start = Hoedjes_data$start, end = Hoedjes_data$end)
Hoedjes <- data.frame(check.names = FALSE,
                      `Candidate FB IDs` = Hoedjes_genes,
                      `Selection Type` = "Longevity",
                      `Populations Used` = "E vs. P",
                      `Source Articles` = "Hoedjes 2019",
                      `Year` = 2019,
                      `Sequencing` = "Genomic",
                      `Direction (Transcriptomics)` = NA_character_
)
###
Fabian_data <- as.data.frame(read_tsv("Longevity/Genomics/Fabian_2018/Fabian_SNPs.tsv"))
Fabian_genes <- gene_mapper(chromosome = Fabian_data$chromosome, start = Fabian_data$start, end = Fabian_data$end)
Fabian <- data.frame(check.names = FALSE,
                     `Candidate FB IDs` = Fabian_genes,
                     `Selection Type` = "Longevity",
                     `Populations Used` = "R vs. L",
                     `Source Articles` = "Fabian 2018",
                     `Year` = 2018,
                     `Sequencing` = "Genomic",
                     `Direction (Transcriptomics)` = NA_character_
)
###
Jalvingh_data <- as.data.frame(read_tsv("Immunity/Genomics/Jalvingh_2014/Jalvingh_SNPs.tsv"))
Jalvingh_genes <- gene_mapper(chromosome = Jalvingh_data$chromosome, start = Jalvingh_data$start, end = Jalvingh_data$end)
Jalvingh <- data.frame(check.names = FALSE,
                       `Candidate FB IDs` = Jalvingh_genes,
                       `Selection Type` = "Immunity",
                       `Populations Used` = "C vs. S (Kraaijeveld)",
                       `Source Articles` = "Jalvingh 2014",
                       `Year` = 2014,
                       `Sequencing` = "Genomic",
                       `Direction (Transcriptomics)` = NA_character_
)
###
Shahrestani_data <- as.data.frame(read_tsv("Immunity/Genomics/Shahrestani_2021/Shahrestani_SNPs.tsv"))
Shahrestani_genes <- gene_mapper(chromosome = Shahrestani_data$chromosome, start = Shahrestani_data$start, end = Shahrestani_data$end)
Shahrestani <- data.frame(check.names = FALSE,
                          `Candidate FB IDs` = Shahrestani_genes,
                          `Selection Type` = "Immunity",
                          `Populations Used` = "C vs. S (Shahrestani)",
                          `Source Articles` = "Shahrestani 2021",
                          `Year` = 2021,
                          `Sequencing` = "Genomic",
                          `Direction (Transcriptomics)` = NA_character_
)
###
Martins_data <- as.data.frame(read_tsv("Immunity/Genomics/Martins_2014/Martins_SNPs.tsv"))
Martins_genes <- gene_mapper(chromosome = Martins_data$chromosome, start = Martins_data$start, end = Martins_data$end)
Martins <- data.frame(check.names = FALSE,
                      `Candidate FB IDs` = Martins_genes,
                      `Selection Type` = "Immunity",
                      `Populations Used` = "Sys_Oral",
                      `Source Articles` = "Martins 2014",
                      `Year` = 2014,
                      `Sequencing` = "Genomic",
                      `Direction (Transcriptomics)` = NA_character_
)
###
Kang_data <- as.data.frame(read_tsv("Stress/Genomics/Kang_2016/Kang_SNPs.tsv"))
Kang_genes <- gene_mapper(chromosome = Kang_data$chromosome, start = Kang_data$start, end = Kang_data$end)
Kang <- data.frame(check.names = FALSE,
                   `Candidate FB IDs` = Kang_genes,
                   `Selection Type` = "Stress",
                   `Populations Used` = "C vs. S (Kang)",
                   `Source Articles` = "Kang 2016",
                   `Year` = 2016,
                   `Sequencing` = "Genomic",
                   `Direction (Transcriptomics)` = NA_character_
)
###
Kawecki_data <- as.data.frame(read_tsv("Stress/Genomics/Kawecki_2021/Kawecki_SNPs.tsv"))
Kawecki_genes <- gene_mapper(chromosome = Kawecki_data$chromosome, start = Kawecki_data$start, end = Kawecki_data$end)
Kawecki <- data.frame(check.names = FALSE,
                      `Candidate FB IDs` = Kawecki_genes,
                      `Selection Type` = "Stress",
                      `Populations Used` = "C vs. S (Kolss)",
                      `Source Articles` = "Kawecki 2021",
                      `Year` = 2021,
                      `Sequencing` = "Genomic",
                      `Direction (Transcriptomics)` = NA_character_
)
###
TS_data <- as.data.frame(read_tsv("Stress/Genomics/Telonis-Scott_2012/TS_SNPs.tsv"))
TS_genes <- gene_mapper(chromosome = TS_data$chromosome, start = TS_data$start, end = TS_data$end)
TS <- data.frame(check.names = FALSE,
                 `Candidate FB IDs` = TS_genes,
                 `Selection Type` = "Stress",
                 `Populations Used` = "C vs. S (TS)",
                 `Source Articles` = "Telonis-Scott 2012",
                 `Year` = 2012,
                 `Sequencing` = "Genomic",
                 `Direction (Transcriptomics)` = NA_character_
)
###
Griffin_data <- as.data.frame(read_tsv("Stress/Genomics/Griffin_2017/Griffin_SNPs.tsv"))
Griffin_data <- split.data.frame(x = Griffin_data, f = Griffin_data$`#`)
Griffin_genes <- unlist(map(.x = Griffin_data, .f = function(x) {
  gene_mapper(chromosome = x$chromosome, start = x$start, end = x$end)
  }
))
Griffin_data <- rbindList(Griffin_data, keepListNames = FALSE)
Griffin_genes <- unique(Griffin_genes[tcount(Griffin_genes) > 1])
Griffin_rejs <- read_tsv(file = "Stress/Genomics/Griffin_2017/Putative_lab_adaptation_genes.txt", col_names = FALSE)[[1]]
Griffin_genes <- setdiff(Griffin_genes, Griffin_rejs); rm(Griffin_rejs)
Griffin <- data.frame(check.names = FALSE,
                      `Candidate FB IDs` = Griffin_genes,
                      `Selection Type` = "Stress",
                      `Populations Used` = "Control vs. Desiccation (Griffin)",
                      `Source Articles` = "Griffin 2017",
                      `Year` = 2017,
                      `Sequencing` = "Genomic",
                      `Direction (Transcriptomics)` = NA_character_
)
###
Zhou_data <- as.data.frame(read_tsv("Stress/Genomics/Zhou_2011/Zhou_SNPs.tsv"))
Zhou_genes <- gene_mapper(chromosome = Zhou_data$chromosome, start = Zhou_data$start, end = Zhou_data$end)
Zhou <- data.frame(check.names = FALSE,
                   `Candidate FB IDs` = Zhou_genes,
                   `Selection Type` = "Stress",
                   `Populations Used` = "Control vs. Hypoxia",
                   `Source Articles` = "Zhou 2011",
                   `Year` = 2011,
                   `Sequencing` = "Genomic",
                   `Direction (Transcriptomics)` = NA_character_
)
###
Michalak_data <- as.data.frame(read_tsv("Stress/Genomics/Michalak_2019/Michalak_SNPs.tsv"))
Michalak_CS <- Michalak_data[Michalak_data$SELECTION_TYPE=="cs",]
Michalak_genes.CS <- gene_mapper(chromosome = Michalak_CS$chromosome, start = Michalak_CS$start, end = Michalak_CS$end)
Michalak.CS <- data.frame(check.names = FALSE,
                          `Candidate FB IDs` = Michalak_genes.CS,
                          `Selection Type` = "Stress",
                          `Populations Used` = "UC vs. CS",
                          `Source Articles` = "Michalak 2019",
                          `Year` = 2019,
                          `Sequencing` = "Genomic",
                          `Direction (Transcriptomics)` = NA_character_
)
Michalak_DS <- Michalak_data[Michalak_data$SELECTION_TYPE=="ds",]
Michalak_genes.DS <- gene_mapper(chromosome = Michalak_DS$chromosome, start = Michalak_DS$start, end = Michalak_DS$end)
Michalak.DS <- data.frame(check.names = FALSE,
                          `Candidate FB IDs` = Michalak_genes.DS,
                          `Selection Type` = "Stress",
                          `Populations Used` = "UC vs. D",
                          `Source Articles` = "Michalak 2019",
                          `Year` = 2019,
                          `Sequencing` = "Genomic",
                          `Direction (Transcriptomics)` = NA_character_
)
Michalak_HEAT <- Michalak_data[Michalak_data$SELECTION_TYPE %in% c("hs","ks"),]
Michalak_genes.HEAT <- gene_mapper(chromosome = Michalak_HEAT$chromosome, start = Michalak_HEAT$start, end = Michalak_HEAT$end)
Michalak.HEAT <- data.frame(check.names = FALSE,
                            `Candidate FB IDs` = Michalak_genes.HEAT,
                            `Selection Type` = "Stress",
                            `Populations Used` = "UC vs. Heat",
                            `Source Articles` = "Michalak 2019",
                            `Year` = 2019,
                            `Sequencing` = "Genomic",
                            `Direction (Transcriptomics)` = NA_character_
)
Michalak_SS <- Michalak_data[Michalak_data$SELECTION_TYPE=="ss",]
Michalak_genes.SS <- gene_mapper(chromosome = Michalak_SS$chromosome, start = Michalak_SS$start, end = Michalak_SS$end)
Michalak.SS <- data.frame(check.names = FALSE,
                          `Candidate FB IDs` = Michalak_genes.SS,
                          `Selection Type` = "Stress",
                          `Populations Used` = "UC vs. ST",
                          `Source Articles` = "Michalak 2019",
                          `Year` = 2019,
                          `Sequencing` = "Genomic",
                          `Direction (Transcriptomics)` = NA_character_
)
Michalak <- rbind.data.frame(Michalak.CS, Michalak.DS, Michalak.HEAT, Michalak.SS)
rm(Michalak.CS, Michalak.DS, Michalak.HEAT, Michalak.SS)

############### Compile gene lists and add to Gene_Table ###############
GenesFromSNPs <- rbind.data.frame(Burke, Graves, Carnes_SNPs, Remolina, Hoedjes, Fabian, Jalvingh, Shahrestani, Martins, Kang, Kawecki, TS, Griffin, Zhou, Michalak)
Gene_Table <- as.data.frame(read_excel("Data/Gene_Table.xlsx"))
Gene_Table <- rbind.data.frame(Gene_Table, GenesFromSNPs)
write_xlsx(Gene_Table, path = "Data/Gene_Table.xlsx", col_names = TRUE)
save.image()
