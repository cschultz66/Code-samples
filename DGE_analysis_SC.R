options(stringsAsFactors = FALSE)
library(limma)
library(edgeR)
library(rtracklayer)
library(Glimma)
library(magrittr)
library(gplots)
library(WGCNA)
enableWGCNAThreads()
library(org.Sc.eg.db)
library(GOstats)
library(AnnotationHub)
library(KEGGREST)

#### QC AND NORMALIZATION ####

# Read in "Targets.txt" file containing sample info and add columns
# based on the Genotype X Day groups:

targets <- readTargets("results/Targets_Final.txt")
targets$group <- paste(targets$Geno, targets$Day, sep = ".")
targets$GpF <- factor(targets$group, levels = c("WT.0", "WT.3", "Mut.0", "Mut.3"))
targets$col <- as.numeric(targets$GpF)

# Record read fates in each sample for QC

read.fate <- targets[,c("QCfiltered","unmapped","multimapped","not.in.gene",
                        "ambiguous","in.a.gene")] / targets$Total * 100

jpeg("results/Demo_readfate.jpeg", width=12, height=6, units="in", res=300, 
     quality=100)
barplot( t(read.fate), beside = TRUE, col = 1:6, las = 2, ylim = c( 0,110), legend.text = FALSE,
         ylab = "Percent of total reads", main = "Read fates per sample", cex.names = 0.8,
         names.arg = targets$Label)
legend ("topright", legend = colnames( read.fate), col = 1:6, fill = 1:6, ncol = 3 )
dev.off()

# Make a design matrix with one column/coefficient per group:

design <- model.matrix(~0 + d.filt$samples$group)
colnames(design) <- levels(d.filt$samples$group)
rownames(design) <- d.filt$samples$Label  

# Make contrast matrix to get the 4 pairwise comparisons and the interaction:

cont.matrix <- makeContrasts(Mut_3vs0 = Mut.3 - Mut.0, 
                             WT_3vs0 = WT.3 - WT.0,
                             t0_MvsWT = Mut.0 - WT.0, 
                             t3_MvsWT = Mut.3 - WT.3,
                             Interact = (Mut.3 - Mut.0) - (WT.3 - WT.0), 
                             levels=design)


#### Read in Count Data #### 

# Read in the raw counts and plot their overall distributions:

d <- readDGE(targets, path = "results/featureCounts/", columns = c(1,7), 
             labels = targets$Label, comment.char = "#", header = TRUE)

# Get minimal gene annotation to get propagated through 

plain.egs <-  gsub("GeneID:","", rownames(d$counts))
info.orgDB <- select(org.Sc.eg.db, keys = plain.egs, keytype="ENTREZID", 
                     columns = c("SYMBOL", "GENENAME"))
rownames(info.orgDB) <- rownames(d$counts)
d$genes <- info.orgDB

# Look at distributions of raw read counts (transformed)

x11()
plotDensities( log2( d$counts + 0.1 ), group = d$samples$GpF, col = 1:4, 
               main = "log2(raw counts + 0.1)", legend = "topright")


# Plot the total library sizes in millions for each sample:

x11()
barplot( colSums(d$counts) / 1e6, ylab = "total number of reads (million)", las = 2, 
         col = d$samples$col, main = "Library sizes", cex.axis = 0.8 )

#### Filtering #### 

# Filter out genes without at least 100 counts in 3 samples

d.filt <- d[filterByExpr(d, design = design, min.count = 100, large.n = 3), , keep.lib.sizes=FALSE]

# Calculate TMM normalization factors 

d.filt <- calcNormFactors(d.filt)

# Draw MDS plot to visualize sample grouping

glMDSPlot(d.filt, top = 5000, labels = targets$Label,
          groups = targets[,c("Geno","Day","Rep","group")],
          folder = "./results/Rstats/glimma-plots", html = "MDS-Plot_trendValues", 
          launch = TRUE)

# Estimate dispersions
d.filt <- estimateDisp(d.filt, design, robust = TRUE)
# setting robust = TRUE will prevent outlier genes with very large or small 
# dispersions from affecting the estimates

fit.edgeR <- glmQLFit(d.filt, design, robust = TRUE)

# Perform quasi-likelihood F-test for each contrast of interest

eR.Mut_3vs0 <- glmQLFTest(fit.edgeR, contrast = cont.matrix[ , 1])
eR.WT_3vs0 <- glmQLFTest(fit.edgeR, contrast = cont.matrix[ , 2])
eR.t0_MvsWT <- glmQLFTest(fit.edgeR, contrast = cont.matrix[ , 3])
eR.t3_MvsWT <- glmQLFTest(fit.edgeR, contrast = cont.matrix[ , 4])
eR.Interact <- glmQLFTest(fit.edgeR, contrast = cont.matrix[ , 5])

edgeR.coded <- new("TestResults", cbind(Mut_3vs0 = decideTestsDGE(eR.Mut_3vs0)[ , 1],
                                        WT_3vs0 = decideTestsDGE(eR.WT_3vs0)[ , 1],
                                        t0_MvsWT = decideTestsDGE(eR.t0_MvsWT)[ , 1], 
                                        t3_MvsWT = decideTestsDGE(eR.t3_MvsWT)[ , 1],
                                        Interact = decideTestsDGE(eR.Interact)[ , 1]))

rownames(edgeR.coded) <- rownames(eR.Mut_3vs0$table)
colnames(edgeR.coded) <- colnames(cont.matrix)
attr(edgeR.coded, "labels") <- c("Down", "NotSig", "Up")

# Get all the genes that are significant in either pairwise time comparison or interaction test

G.all <- rownames(edgeR.coded)[edgeR.coded[,1] != 0 |
                                 edgeR.coded[,2] != 0 |
                                 edgeR.coded[,5] != 0]

#### HEATMAPS ####

# Make a scaled heatmap of the G.all genes. Select the normalized 
# expression values and create a color scale:

logCPM <- cpm(d.filt, prior.count = 3, log = TRUE)
heatdata <- logCPM[G.all, ]
heatdata.scaled <- heatdata %>% t() %>% scale() %>% t()
color.scale <- colorpanel(100, "blue", "white", "red")

# Calculate gene dendrogram and split into clusters

cluster.all <- hclust(dist(heatdata.scaled))
cutree.clust.treeF <- cutreeDynamic(dendro = cluster.all, method = "tree", 
                                    minClusterSize = 25, deepSplit = F)
treeF.rowcols <- labels2colors(cutree.clust.treeF)

# Write out a high-res heatmap to the current working directory:

jpeg("results/Rstats/Demo_Heatmap.jpeg", width=5, height=10, units="in", res=300, 
     quality=100)
out.heatall <- heatmap.2(heatdata.scaled, col = color.scale, 
                         density.info = "none", trace = "none", labRow = "",
                         margins = c(7,2), key.xlab = "SD from mean",
                         Colv = FALSE, dendrogram = "row", Rowv = as.dendrogram(cluster.all),
                         keysize = 1, main = "scaled heatmap",
                         RowSideColors=treeF.rowcols)
dev.off()


#### ANNOTATION ####

# Get additional gene annotation like chromosome, GO and KEGG

# Get chromosome location from gtf file:

gtf0 <- import("data/genome/GCF_000001635.26_GRCm38.noAlts.noPatches.gtf")
# Subset down to just the exon rows:
gtf0 <- gtf0[gtf0$type == "exon"]

# Remove rows with duplicate gene_ids:

gtf0 <- gtf0[!duplicated(gtf0$gene_id)]
gtf.gene.df <- as.data.frame(gtf0)

#keep seqnames, start, strand, gbkey and gene_id

gtf.gene.df <- gtf.gene.df[,c(1,2,5,14,22)]

# Add annotation data to filtered list of genes

d.filt$genes <- cbind(d.filt$genes, gtf.gene.df[rownames(d.filt$genes),1:4])
rm(gtf0)

# Form lists of differentially expressed GO terms and KEGG pathways

G.ent <- gsub("GeneID:","", G.all)
GOTERM <- goana(G.ent, species = "Sc", trend = TRUE)
deGO <- GOTERM[GOTERM$DE > 0,]

KEGG <- kegga(G.ent, species = "Sc", trend = TRUE)
deKEGG <- KEGG[KEGG$DE > 0,]

#### SAVING RESULTS ####

write.table(cbind(KEGGID = rownames(deKEGG), deKEGG), 
            file = "results/Rstats/KEGG_overrepTests.txt", row.names = FALSE, sep = "\t")

write.table(cbind(GOID = rownames(deGO), deGO), 
            file = "results/Rstats/GO_overrepTests.txt", row.names = FALSE, sep = "\t")

write.table(eR.Mut_3vs0, file = "results/Rstats/eR.Mut_3vs0.txt", row.names = FALSE, sep = "\t")
write.table(er.WT_3vs0, file = "results/Rstats/er.WT_3vs0.txt", row.names = FALSE, sep = "\t")
write.table(eR.t0_MvsWT, file = "results/Rstats/eR.t0_MvsWT.txt", row.names = FALSE, sep = "\t")
write.table(eR.t3_MvsWT, file = "results/Rstats/eR.t3_MvsWT.txt", row.names = FALSE, sep = "\t")
write.table(eR.Interact, file = "results/Rstats/eR.Interact.txt", row.names = FALSE, sep = "\t")

save.image("results/Rstats/DGE_analysis_SC.RData")