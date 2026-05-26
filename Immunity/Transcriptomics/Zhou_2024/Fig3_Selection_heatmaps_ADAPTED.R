### doi:10.17863/CAM.101651 ###
## https://www.repository.cam.ac.uk/bitstreams/6601be17-a344-4101-99be-7827375bb901/download ##

library(RColorBrewer)
library(gplots)
library(edgeR)

#for plotting heatmaps from RNAseq analysis using edgeR

#load read counts
#log fold expression
#create matrix of expression and scale expression using z-score normalization
#plot scaled expression

####extract counts per million gene expression data####

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#import reads data
with_symbol <- read.table("Selected_FB_RNAseq_counts.txt", header = TRUE)
with_symbol$validated_id <- replaceEntries(x = with_symbol$validated_id, map = MAP)
with_symbol <- with_symbol[2:39]
str(with_symbol)
#import design matrix
design_matrix <- read.table("Selected_FB_RNAseq_design_matrix.csv", 
                            sep = ",", dec = ".", header = TRUE)
str(design_matrix)

#factor and relevel design matrix
design_matrix$selection <- factor(design_matrix$selection,
                                  levels = c("Control", "Nsref"))
design_matrix$selection <- relevel(design_matrix$selection, ref = "Control")

design_matrix$treatment <- factor(design_matrix$treatment,
                                  levels = c("INF", "CTRL"))
design_matrix$treatment <- relevel(design_matrix$treatment, ref = "CTRL")
design_matrix$time <- factor(design_matrix$time, levels = c(2,12,24))

design_matrix$group <- interaction(design_matrix$selection, 
                                   design_matrix$treatment, 
                                   design_matrix$time)
design_matrix$group <- factor(design_matrix$group, levels = unique(design_matrix$group))
design_matrix
str(design_matrix)



# make DGEList object 
DGE_all <- DGEList(counts = with_symbol[,3:38], 
                   group = design_matrix$group, 
                   genes = with_symbol[,1:2])
FB_invs <- as.data.frame(read_tsv("E:/thesis/Scripts/misc/withdrawns.txt"))
DGE_all <- DGE_all[!(DGE_all$genes$validated_id %fin% FB_invs$FBIDs),]; rm(FB_invs)

# filtering & normalization 
keep <- filterByExpr(DGE_all, group = design_matrix$group)
y_all <- DGE_all[keep, , keep.lib.sizes = FALSE]

y_all <- calcNormFactors(y_all)
y_all$samples


group <- design_matrix$group
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
design

y_all <- estimateDisp(y_all, design, robust = TRUE)
y_all$common.dispersion

fit_QLF <- glmQLFit(y_all, design, robust = TRUE)


# get counts per million per gene 
cpm_log <- cpm(y_all, log = TRUE, prior.count = 0.25)
cpm_symbol <- cbind(y_all$genes, cpm_log)
head(cpm_symbol)




#### get top 30 DEGs under uninfected condition ####
SelectionEffects_CTRL <- makeContrasts( ((Nsref.CTRL.2 + Nsref.CTRL.12 + Nsref.CTRL.24)/3 ) - ((Control.CTRL.2 + Control.CTRL.12 + Control.CTRL.24)/3 ), levels = group)


qlf_CTRL <- glmQLFTest(fit_QLF, con = SelectionEffects_CTRL)

#filter log fold change for greater than 1
ctrl_avg_results <- as.data.frame(topTags(qlf_CTRL, Inf))
ctrl_avg_deg <- ctrl_avg_results[abs(ctrl_avg_results$logFC) > 1 & 
                                   ctrl_avg_results$FDR < 0.05, ]
rownames(ctrl_avg_deg) <- ctrl_avg_deg$validated_id
Zhou.24_Uninf_genes <- ctrl_avg_deg$validated_id

#### get top 30 DEGs under infected conditions ####
SelectionEffects_Inf <- makeContrasts( ((Nsref.INF.2 + Nsref.INF.12 + Nsref.INF.24)/3 ) - ((Control.INF.2 + Control.INF.12 + Control.INF.24)/3 ), levels = group)


qlf_Inf <- glmQLFTest(fit_QLF, con = SelectionEffects_Inf)


#filter log fold change for greater than 1
inf_avg_results <- as.data.frame(topTags(qlf_Inf, Inf))
inf_avg_deg <- inf_avg_results[abs(inf_avg_results$logFC) > 1 & 
                                 inf_avg_results$FDR <0.05, ]
rownames(inf_avg_deg) <- inf_avg_deg$validated_id
Zhou.24_Inf_genes <- inf_avg_deg$validated_id

#genes present as significant in both groups
both_sig <- intersect(Zhou.24_Inf_genes, Zhou.24_Uninf_genes)

avg_deg <- rbind.data.frame(ctrl_avg_deg[(!ctrl_avg_deg$validated_id %in% both_sig),], inf_avg_deg[(!inf_avg_deg$validated_id %in% both_sig),])
both_sig_deg <- cbind.data.frame(ctrl_avg_deg[both_sig,1:2], ((ctrl_avg_deg[both_sig,3:7] + inf_avg_deg[both_sig,3:7])/2))
avg_deg <- rbind.data.frame(avg_deg, both_sig_deg)

Direction <- as.factor(avg_deg$logFC>0)
levels(Direction) <- c("Down", "Up")
avg_deg <- cbind.data.frame(avg_deg, Direction = unfactor(Direction))

writexl::write_xlsx(x = avg_deg, path = "avg_DEGs.xlsx")
