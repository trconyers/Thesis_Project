##################### Make BSgenome package for r6.61. #####################
library(BSgenomeForge)
library(easy.utils)
library(S4Vectors)
library(readxl)
library(stringr)
library(Biostrings)

wd <- getwd()
UCSC_cnvrt <- function(chr) {
  Convert <- as.data.frame(read_excel("Data/UCSC.xlsx"))
  x <- factor(chr)
  map <- as.list(Convert$FlyBase)
  names(map) <- Convert$UCSC
  res <- as.character(replaceEntries(x, map))
  unname(Vectorize(str_remove)(res, "chr"))
}

download.file(url = "https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.61_FB2024_06/fasta/dmel-all-chromosome-r6.61.fasta.gz", destfile = "misc/dmel-all-chromosome-r6.61.fa.gz")
dna <- readDNAStringSet("misc/dmel-all-chromosome-r6.61.fa.gz")
names(dna) <- str_split_i(string = names(dna), pattern = " ", i = 1)
names(dna) <- UCSC_cnvrt(names(dna))
names(dna) <- paste0("chr", names(dna))
writeXStringSet(x = dna, filepath = "misc/dmel-all-chromosome-r6.61.fa.gz")
rm(dna)

fastaTo2bit("misc/dmel-all-chromosome-r6.61.fa.gz",
            "misc/seqs/dmel-r6.61.2bit")
destdir <- tempdir()
forgeBSgenomeDataPkg(
  "misc/BSgenome.Dmelanogaster.FlyBase.dm6-seed",
  seqs_srcdir = "misc/seqs",
  destdir = destdir,
  replace = TRUE,
  verbose = TRUE
)

setwd(destdir)
pkg <- pkgbuild::build(binary = TRUE, manual = TRUE)
setwd(wd)

install.packages(pkg, repos = NULL, type = "win.binary")

##################### Make TxDb for r6.61. #####################
library(purrr)
library(rtracklayer)
library(Seqinfo)
library(IRanges)
library(BSgenome.Dmelanogaster.FlyBase.dm6)
dmel_genome <- BSgenome.Dmelanogaster.FlyBase.dm6

gtfFile <- "misc/dmel-all-r6.61.gtf.gz"
download.file(url = "https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.61_FB2024_06/gtf/dmel-all-r6.61.gtf.gz", destfile = gtfFile)

mcolnames <- function(x) compose(colnames, mcols)(x)
`mcolnames<-` <- function(x, value) {
  mcols(x) <- `colnames<-`(x = mcols(x), value = value)
  return(x)
}

GTFgr <- import.gff2(gtfFile)
seqlevels(GTFgr) <- UCSC_cnvrt(seqlevels(GTFgr))
seqlevels(GTFgr) <- paste0("chr", seqlevels(GTFgr))
mcolnames(GTFgr) <- str_replace(string = mcolnames(GTFgr), pattern = "symbol", replacement = "name")

transcripts <- GTFgr[GTFgr$type!="gene"]
names(transcripts) <- transcripts$transcript_id

genes <- GTFgr[GTFgr$type=="gene"]
names(genes) <- genes$gene_id
genes$transcript_id <- Ensembl_Canon[names(genes),2]
genes$transcript_name <- transcripts[genes$transcript_id]$transcript_name
names(genes) <- genes$transcript_id

dmel_gr <- c(transcripts,genes)
dmel_gr <- sort_by(x = dmel_gr, y = dmel_gr$transcript_id)
dmel_gr@seqinfo <- seqinfo(dmel_genome)[levels(dmel_gr@seqnames)]
problems <- which(dmel_gr$type=="exon" & names(dmel_gr) %in% c("FBtr0100857", "FBtr0100863", "FBtr0433500", "FBtr0433501"))
end(dmel_gr[problems[1:2]]) <- end(dmel_gr[problems[1:2]]) + 2
start(dmel_gr[problems[3:4]]) <- start(dmel_gr[problems[3:4]]) - 2
rm(GTFgr,transcripts,genes,problems)

dmel_txdb <- txdbmaker::makeTxDbFromGRanges(dmel_gr, taxonomyId = 7227)
save.image()
