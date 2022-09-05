## 1.normalize htseq count matrix using DESeq2
htseq_count_iyBomTerr1.2 <- read.delim("~/htseq_count_iyBomTerr1.2.txt", header=FALSE)
#iyBomTerr1.2 13251 probes, Bter_1.0 11803 probes

library(readr)
library(DESeq2)

#construct count matrix
#column name = "gene_id" + sample names
colnames(htseq_count_iyBomTerr1.2)<- c("gene_id","CLO_C02_Q", "CLO_C02_W", "CLO_C27_Q", "CLO_C27_W", "CLO_C38_Q", "CLO_C38_W", "CLO_C67_Q", "CLO_C67_W",
                                       "CON_C06_Q", "CON_C06_W", "CON_C34_Q", "CON_C34_W", "CON_C48_Q", "CON_C48_W", "CON_C61_Q", "CON_C61_W", 
                                       "IMI_C07_Q", "IMI_C07_W", "IMI_C32_Q", "IMI_C32_W", "IMI_C45_Q", "IMI_C45_W", "IMI_C55_Q", "IMI_C55_W")
#remove the summary statistics (last 4 rows) from htseq count matrix
htseq_count_iyBT_no_stats <- htseq_count_iyBomTerr1.2[-((nrow(htseq_count_iyBomTerr1.2)-4):nrow(htseq_count_iyBomTerr1.2)),]
#remove row name
rownames(htseq_count_iyBT_no_stats) <- htseq_count_iyBT_no_stats[,1]
htseq_count_iyBT_matrix <- htseq_count_iyBT_no_stats[,-1]

#construct coldata (sample info) for the matrix
htseq_count_iyBt_coldata <- data.frame(neonicotinoid = c("CLO", "CLO", "CLO", "CLO", "CLO", "CLO", "CLO", "CLO",
                                                         "CON", "CON", "CON", "CON", "CON", "CON", "CON", "CON", 
                                                         "IMI", "IMI", "IMI", "IMI", "IMI", "IMI", "IMI", "IMI"),
                                       caste = c("Q","W","Q","W","Q","W","Q","W",
                                                 "Q","W","Q","W","Q","W","Q","W",
                                                 "Q","W","Q","W","Q","W","Q","W"),
                                       row.names = colnames(htseq_count_iyBT_matrix))
htseq_count_iyBt_coldata$neonicotinoid <- as.factor(htseq_count_iyBt_coldata$neonicotinoid)
htseq_count_iyBt_coldata$caste <- as.factor(htseq_count_iyBt_coldata$caste)
all(rownames(htseq_count_iyBt_coldata) == colnames(htseq_count_iyBT_matrix))

#data input for DESeq2
count_matrix_DESeq_iyBT <- DESeqDataSetFromMatrix(countData = htseq_count_iyBT_matrix,
                                                  colData = htseq_count_iyBt_coldata,
                                                  design = ~neonicotinoid,
                                                  tidy = FALSE,
                                                  ignoreRank = FALSE)
#DESeq2 calculation
count_matrix_DESeq_iyBT_calculation <- DESeq(count_matrix_DESeq_iyBT)
count_matrix_DESeq_iyBT_result <- results(count_matrix_DESeq_iyBT_calculation)
count_matrix_DESeq_iyBT_result

#variance calculation
count_matrix_DESeq_iyBT_var_stab <- varianceStabilizingTransformation(count_matrix_DESeq_iyBT_calculation)
count_matrix_DESeq_iyBT_var_stab_result <- getVarianceStabilizedData(count_matrix_DESeq_iyBT_calculation)
#calculate variance for each gene (each row represents a gene)
library(genefilter)
count_matrix_DESeq_iyBT_var_stab_result_rowvar <- rowVars(count_matrix_DESeq_iyBT_var_stab_result)
summary(count_matrix_DESeq_iyBT_var_stab_result_rowvar)

#data normalization
#remove genes with expression sample-wise variation < 95 quantile to reduce noise
q95_rowvar_iyBT <- quantile(rowVars(count_matrix_DESeq_iyBT_var_stab_result), 0.95)
var_normalized_iyBT <- count_matrix_DESeq_iyBT_var_stab_result[ count_matrix_DESeq_iyBT_var_stab_result_rowvar > q95_rowvar_iyBT, ]

##2.WGCNA using WGCNA
#WGCNA input
#transpose the expression matrix, genes as columns, samples as rows
WGCNA.input_iyBT <- t(var_normalized_iyBT)
#663 gene probes

library(WGCNA)
#allow multi-core calculation to speed up
allowWGCNAThreads(nThreads = 8)

#check if there is any deviation value
sampleTree_iyBT <- hclust(dist(WGCNA.input_iyBT), method ="average");
plot(sampleTree_iyBT)

#test soft thresholds
#maximum power is 30
powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2))
soft.threshold_iyBT <- pickSoftThreshold(WGCNA.input_iyBT,
                                         powerVector = powers,
                                         verbose = 5)
#Scale-free topology fit index as a function of the soft-threshold power
#Pick a soft threshold power near the curve of the plot
plot(soft.threshold_iyBT$fitIndices[, 1],
     -sign(soft.threshold_iyBT$fitIndices[, 3]) * soft.threshold_iyBT$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence"))
text(soft.threshold_iyBT$fitIndices[, 1],
     -sign(soft.threshold_iyBT$fitIndices[, 3]) * soft.threshold_iyBT$fitIndices[, 2],
     labels = powers, cex=0.9, col = "red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="yellow")
#no power is higher than 0.9
abline(h=0.85,col="red")
#power=26 is over 0.85
plot(soft.threshold_iyBT$fitIndices[, 1],
     soft.threshold_iyBT$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity"))
text(soft.threshold_iyBT$fitIndices[, 1],
     soft.threshold_iyBT$fitIndices[, 5],
     labels = powers,
     cex = 0.9, col = "red")
#mean connectivity should not be too low
abline(h=100,col="red")

#one-step network building
WGCNA.network_iyBT_signed <- blockwiseModules(WGCNA.input_iyBT, #input
                                              # == Adjacency Function ==
                                              power = 26, # power
                                              networkType = "signed", #signed or unsigned
                                              # == Tree and Block Options ==
                                              deepSplit = 2,
                                              pamRespectsDendro = F,
                                              # detectCutHeight = 0.75,
                                              minModuleSize = 30,
                                              maxBlockSize = 4000,
                                              # == Module Adjustments ==
                                              reassignThreshold = 0,
                                              mergeCutHeight = 0.25,
                                              # == TOM == Archive the run results in TOM file (saves time)
                                              saveTOMs = T,
                                              saveTOMFileBase = "iyBT",
                                              # == Output Options
                                              numericLabels = T,
                                              verbose = 3)

#Convert labels to colors for plotting
mergedColors_iyBT_signed = labels2colors(WGCNA.network_iyBT_signed$colors)
#Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  WGCNA.network_iyBT_signed$dendrograms[[1]],
  mergedColors_iyBT_signed[WGCNA.network_iyBT$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )
table(labels2colors(WGCNA.network_iyBT_signed$colors))
#change default colors for better display
customColorOrder = c("blue", "red", "yellow")
mergedColors_iyBT_signed_updated = labels2colors(WGCNA.network_iyBT_signed$colors, colorSeq = customColorOrder)
plotDendroAndColors(
  WGCNA.network_iyBT_signed$dendrograms[[1]],
  mergedColors_iyBT_signed_updated,
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

## 3.calculate module-wise relations
#show genes in the module
module_df_iyBT_signed <- data.frame(gene_id = names(WGCNA.network_iyBT_signed$colors),
                             colors = labels2colors(WGCNA.network_iyBT_signed$colors, colorSeq = customColorOrder))
library(readr)
write_delim(module_df_iyBT_signed, file = "gene_modules_iyBT.txt", delim = "\t")
#export genes in different modules
module_blue <- as.data.frame(module_df_iyBT[which(module_df_iyBT_signed$colors=="blue"),])
module_red <- as.data.frame(module_df_iyBT[which(module_df_iyBT_signed$colors=="red"),])
module_yellow <- as.data.frame(module_df_iyBT[which(module_df_iyBT_signed$colors=="yellow"),])

# sined network: Get Module Eigengenes with updated color
MEs0_iyBT_signed_updated <- moduleEigengenes(WGCNA.input_iyBT, mergedColors_iyBT_signed_updated)$eigengenes
# Reorder modules so similar modules are next to each other
MEs0_iyBT_signed_updated <- orderMEs(MEs0_iyBT_signed_updated)
module_order_iyBT_updated <- names(MEs0_iyBT_signed_updated) %>% gsub("ME","", .)

#plot adjacency between modules
plotEigengeneNetworks(MEs0_iyBT_signed, 
                      "Eigengene adjacency heatmap", #title 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2),
                      plotAdjacency = T,
                      printAdjacency = T,
                      plotDendrograms = T, 
                      xLabelsAngle = 90)

#visualize co-expression network
dissTOM_iyBT_signed <- 1-TOMsimilarityFromExpr(WGCNA.input_iyBT, power = 26, networkType = "signed")
#transform the power of dissTOM, making medium strong links more obvious in the heatmap
plotTOM_iyBT_signed <- dissTOM_iyBT_signed^7
#set diagonal as NA
diag(plotTOM_iyBT_signed) <- NA
gene_tree_iyBT_signed <- WGCNA.network_iyBT_signed$dendrograms[[1]]
#make sure the colors of cells range from lemonchiffon to red
library(gplots)
heatmap_color <- colorpanel(250, "red", "orange", "lemonchiffon")
TOMplot(plotTOM_iyBT_signed, gene_tree_iyBT_signed, mergedColors_iyBT_signed, col = heatmap_color, main = "Network heatmap plot of all genes")
#patches on diagonal are modules, reds are high overlap, whites are low overlap

## 4.modules-treatments correlation
#build sample treatments info matrix, treatments represented by numbers
num_sample <- c("CLO_C02_Q", "CLO_C02_W", "CLO_C27_Q", "CLO_C27_W", "CLO_C38_Q", "CLO_C38_W", "CLO_C67_Q", "CLO_C67_W",
                "CON_C06_Q", "CON_C06_W", "CON_C34_Q", "CON_C34_W", "CON_C48_Q", "CON_C48_W", "CON_C61_Q", "CON_C61_W", 
                "IMI_C07_Q", "IMI_C07_W", "IMI_C32_Q", "IMI_C32_W", "IMI_C45_Q", "IMI_C45_W", "IMI_C55_Q", "IMI_C55_W")
#0=CONTROL, 1=CON, 2=IMI
num_neonicotinoid <- c(1,1,1,1,1,1,1,1,
                       0,0,0,0,0,0,0,0,
                       2,2,2,2,2,2,2,2)
#0=QUEEN, 1=WORKER
num_caste <- c(0,1,0,1,0,1,0,1,
               0,1,0,1,0,1,0,1,
               0,1,0,1,0,1,0,1)
datTraits_num <- data.frame(num_sample, num_neonicotinoid, num_caste)
rownames(datTraits_num) <- datTraits_num$num_sample
datTraits_num <- datTraits_num[, 2:3]
colnames(datTraits_num) <- c("neonicotinoid", "caste")
#pearson's correlation
module_trait_cor_num_signed <- cor(MEs0_iyBT_signed_updated, datTraits_num, use = "p") #color updated
nsample <- nrow(WGCNA.input_iyBT)
module_trait_pvalue_num_signed <- corPvalueStudent(module_trait_cor_num_signed, nsample)
#extract text to show on the correlation heatmap
text_Matrix_num_signed <- paste(signif(module_trait_cor_num_signed, 2), "\n(",
                                signif(module_trait_pvalue_num_signed, 1), ")", sep = "")
dim(text_Matrix_num_signed) <- dim(module_trait_cor_num_signed)
#adjust margin for a suitable plot display
sizeGrWindow(20,20)
par(mar = c(6, 6.5, 4.1, 2.1))
labeledHeatmap(Matrix = module_trait_cor_num_signed,
               xLabels = colnames(datTraits_num),
               yLabels = names(MEs0_iyBT_signed_updated),
               ySymbols = c("blue", "red", "yellow", "grey"),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = text_Matrix_num_signed,
               cex.text = 0.9,
               setStdMargins = FALSE,
               main = paste("Module-treatment relationships"),
               font.lab.x = 2,
               font.lab.y = 2)

## 5.export the network to cytoscape 
modules <- c("blue","brown","turquoise")
probes <- row.names(var_normalized_iyBT)
inmodule <- is.finite(match(mergedColors_iyBT_signed,modules))
probes_in_module <- probes[inmodule]
#select corresponding topological overlap
modTOM <- TOM_iyBT_signed[inmodule,inmodule]
dimnames(modTOM) <- list(probes_in_module,probes_in_module)
#export the network into edge and node list files
edge_list_0.4 <- exportNetworkToCytoscape(modTOM,
                                          edgeFile = paste("CytoscapeInput-edges-0.4", paste(modules, collapse="-"), ".txt", sep=""),
                                          nodeFile = paste("CytoscapeInput-nodes-0.4", paste(modules, collapse="-"), ".txt", sep=""),
                                          weighted = TRUE,
                                          threshold = 0.4,
                                          nodeNames = probes_in_module,
                                          nodeAttr = mergedColors_iyBT_signed[inmodule])

## 6.identify hub genes using dgha
library(dhga)
#construc gene con-expression network at threshold level 0.4
dhga.input_iyBT <- as.data.frame(var_normalized_iyBT)
dhga_adj_matrix_iyBT <- Adjacency(dhga.input_iyBT, beta = 26, threshold = 0.4)

#differential hub status of genes in the network, threshold 0.0001
dhga_diffhub_iyBT <- DiffHub(dhga.input_iyBT[,-(9:16)], dhga.input_iyBT[,9:16], 
                             16,8,80,26, alpha = 0.0001, plot = TRUE)

#random sampling is used in DiffHub, to increase reliability, repeat the calculation for 10 times
H <- data.frame(first = 0, second = 0, third =0, fourth =0, fifth=0, sixth=0, seventh=0, eighth=0, ninth=0, tenth=0)
H <- data.frame(gene=rep(0,10), row.names = c(1,2,3,4,5,6,7,8,9,10))
#HK = housekeeping hub genes, HC = hub genes unique to control, HT = hub genes unique to neonicotinoid treatment
HK <- data.frame(blank = rep(0,300))
HC <- data.frame(blank = rep(0,150))
HT <- data.frame(blank = rep(0,100))
temp_H <- data.frame(gene = 0)
temp_C <- data.frame(gene = 0)
temp_T <- data.frame(gene = 0)

for (i in 1:10) {
  dhga_diffhub_iyBT <- DiffHub(dhga.input_iyBT[,-(9:16)], dhga.input_iyBT[,9:16], 
                               16,8,80,26, alpha = 0.0001, plot = FALSE)
  temp_H <- as.data.frame(dhga_diffhub_iyBT[which(dhga_diffhub_iyBT$HubStatus=="Housekeeping Hub"),1])
  temp_C <- as.data.frame(dhga_diffhub_iyBT[which(dhga_diffhub_iyBT$HubStatus=="Unique Hub to Normal"),1])
  temp_T <- as.data.frame(dhga_diffhub_iyBT[which(dhga_diffhub_iyBT$HubStatus=="Unique Hub to Stress"),1])
  colnames(temp_H) <- "gene"
  colnames(temp_C) <- "gene"
  colnames(temp_T) <- "gene"
  temp_H_sup <- data.frame(gene = rep(0,(300-nrow(temp_H))))
  temp_C_sup <- data.frame(gene = rep(0,(300-nrow(temp_C))))
  temp_T_sup <- data.frame(gene = rep(0,(300-nrow(temp_T))))
  temp_H <- rbind(temp_H,temp_H_sup)
  temp_C <- rbind(temp_C,temp_C_sup)
  temp_T <- rbind(temp_T,temp_T_sup)
  HK <- cbind(HK, temp_H)
  HC <- cbind(HC, temp_C)
  HT <- cbind(HT, temp_T)
}

#identify common hub genes in each attempt
HK <- HK[,2:11]
HC <- HC[,2:11]
HT <- HT[,2:11]
colnames(HK) <- c(1,2,3,4,5,6,7,8,9,10)
colnames(Hc) <- c(1,2,3,4,5,6,7,8,9,10)
colnames(HT) <- c(1,2,3,4,5,6,7,8,9,10)
HK_common <- Reduce(intersect, list(HK[,1],HK[,2],HK[,3],HK[,4],HK[,5],HK[,6],HK[,7],HK[,8],HK[,9],HK[,10]))
HC_common <- Reduce(intersect, list(HC[,1],HC[,2],HC[,3],HC[,4],HC[,5],HC[,6],HC[,7],HC[,8],HC[,9],HC[,10]))
HT_common <- Reduce(intersect, list(HT[,1],HT[,2],HT[,3],HT[,4],HT[,5],HT[,6],HT[,7],HT[,8],HT[,9],HT[,10]))
#remove the last row
HK_common <- as.data.frame(HK_common[1:length(HK_common)-1])
HC_common <- as.data.frame(HC_common[1:length(HC_common)-1])
HT_common <- as.data.frame(HT_common[1:length(HT_common)-1])
colnames(HK_common) <- colnames(HC_common) <- colnames(HT_common) <- "gene"

#identify module membership of hub genes
HK_red <- intersect(HK_common$gene,module_red$gene_id)
HC_red <- intersect(HC_common$gene,module_red$gene_id)
HT_red <- intersect(HT_common$gene,module_red$gene_id)
HK_blue <- intersect(HK_common$gene,module_blue$gene_id)
HC_blue <- intersect(HC_common$gene,module_blue$gene_id)
HT_blue <- intersect(HT_common$gene,module_blue$gene_id)
HK_yellow <- intersect(HK_common$gene,module_yellow$gene_id)
HC_yellow <- intersect(HC_common$gene,module_yellow$gene_id)
HT_yellow <- intersect(HT_common$gene,module_yellow$gene_id)

#reduce number of hub genes if needed
#identify hub genes based on gene connection significance values, choose the top 10 as hub
dhga_hub_con_sig_iyBT <- hub.pval.cutoff(dhga.input_iyBT, beta = 26, m = 24, s = 80, n = 10)
dhga_con_sig_iyBT <- pvalue.hub(dhga.input_iyBT, beta = 26, m = 24, s = 80, plot = TRUE)
dhga_con_sig_value_iyBT <- pvalue.hub(dhga.input_iyBT, beta = 26, m = 24, s = 80, plot = FALSE)

#identify hub genes based on weighted gene score, top 20 as hub genes
dhga_hub_wgs_iyBT <- hub.wgs(dhga.input_iyBT, beta = 26, n = 20)
dhga_wgs_iyBT <- WeightedGeneScore(dhga.input_iyBT, beta = 26, plot=TRUE)
head(dhga_wgs_iyBT)

#house keeping genes
dhga_hub_common_iyBT <- dhga_diffhub_iyBT[which(dhga_diffhub_iyBT$HubStatus=="Housekeeping Hub"),]
hub_common_wgs <- data.frame(Genes=0,WGS=0)
temp <- data.frame(Genes=0,WGS=0)
for (i in 1:length(dhga_hub_common_iyBT$Genes)) {
  temp$Genes <- dhga_hub_common_iyBT$Genes[i]
  temp$WGS <- dhga_wgs_iyBT[which(rownames(dhga_wgs_iyBT)==dhga_hub_common_iyBT$Genes[i]),]
  hub_common_wgs <- rbind(hub_common_wgs,temp)
}
hub_common_wgs<-hub_common_wgs[-1,]

#stress hub genes
dhga_hub_stress_iyBT <- dhga_diffhub_iyBT[which(dhga_diffhub_iyBT$HubStatus=="Unique Hub to Stress"),]
hub_stress_wgs <- data.frame(Genes=0,WGS=0)
temp <- data.frame(Genes=0,WGS=0)
for (i in 1:length(dhga_hub_stress_iyBT$Genes)) {
  temp$Genes <- dhga_hub_stress_iyBT$Genes[i]
  temp$WGS <- dhga_wgs_iyBT[which(rownames(dhga_wgs_iyBT)==dhga_hub_stress_iyBT$Genes[i]),]
  hub_stress_wgs <- rbind(hub_stress_wgs,temp)
}
hub_stress_wgs<-hub_stress_wgs[-1,]
sort(hub_stress_wgs$WGS, decreasing = TRUE)

#control hub genes
dhga_hub_control_iyBT <- dhga_diffhub_iyBT[which(dhga_diffhub_iyBT$HubStatus=="Unique Hub to Normal"),]
hub_control_wgs <- data.frame(Genes=0,WGS=0)
temp <- data.frame(Genes=0,WGS=0)
for (i in 1:length(dhga_hub_control_iyBT$Genes)) {
  temp$Genes <- dhga_hub_control_iyBT$Genes[i]
  temp$WGS <- dhga_wgs_iyBT[which(rownames(dhga_wgs_iyBT)==dhga_hub_control_iyBT$Genes[i]),]
  hub_control_wgs <- rbind(hub_control_wgs,temp)
}
hub_control_wgs<-hub_control_wgs[-1,]


## 7.gene ontology enrichment analysis using TopGO
library(topGO)
#read in the Bter1.0 GO terms
Bt_geneID2GO <- readMappings(file = "~/Bt_GO_term.txt")
str(head(Bt_geneID2GO))

#GO enrichment analysis for modules
#blue module
gene_universe <- rbind(module_blue, module_red, module_yellow)
blue_list <- factor(as.integer(gene_universe$gene_id %in% module_blue$gene_id))
names(blue_list) <- gene_universe$gene_id
#blue module cellular component
blue_go_CC <- new("topGOdata",
                  description = "GO_test_BF",
                  ontology    = "CC",
                  allGenes    = blue_list,
                  annot       = annFUN.gene2GO,
                  gene2GO     = Bt_geneID2GO,
                  nodeSize    = 10) # Modify to reduce/increase stringency.
#fisher test to test significance
blue_resultFisher_CC <- runTest(blue_go_CC, algorithm = "classic", statistic = "fisher")
#Kolmogorov–Smirnov test to calculate significance
blue_resultKS_CC <- runTest(blue_go_CC, algorithm = "classic", statistic = "ks")
#Kolmogorov–Smirnov test with elim method to calculate significance
blue_resultKS.elim_cc <- runTest(blue_go_CC, algorithm = "elim", statistic = "ks")
#show top 10 result, ranking in classic Kolmogorov–Smirnov test
blue_allRes_CC <- GenTable(blue_go_CC, classicFisher = blue_resultFisher_cc,
                           classicKS = blue_resultKS_cc, elimKS = blue_resultKS.elim_cc,
                           orderBy = "classicKS", ranksOf = "classicKS", topNodes = 10)
#blue module molecular function
blue_go_MF <- new("topGOdata",
                  description = "GO_test_BF",
                  ontology    = "MF",
                  allGenes    = blue_list,
                  annot       = annFUN.gene2GO,
                  gene2GO     = Bt_geneID2GO,
                  nodeSize    = 10) # Modify to reduce/increase stringency.
blue_resultFisher_MF <- runTest(blue_go_MF, algorithm = "classic", statistic = "fisher")
blue_resultKS_MF <- runTest(blue_go_MF, algorithm = "classic", statistic = "ks")
blue_resultKS.elim_MF <- runTest(blue_go_MF, algorithm = "elim", statistic = "ks")
blue_allRes_MF <- GenTable(blue_go_MF, classicFisher = blue_resultFisher_MF,
                           classicKS = blue_resultKS_MF, elimKS = blue_resultKS.elim_MF,
                           orderBy = "classicKS", ranksOf = "classicKS", topNodes = 10)
#blue module biological process
blue_go_BP <- new("topGOdata",
                  description = "GO_test_BF",
                  ontology    = "BP",
                  allGenes    = blue_list,
                  annot       = annFUN.gene2GO,
                  gene2GO     = Bt_geneID2GO,
                  nodeSize    = 10) # Modify to reduce/increase stringency.
blue_resultFisher_BP <- runTest(blue_go_BP, algorithm = "classic", statistic = "fisher") #3 term found
blue_resultKS_BP <- runTest(blue_go_BP, algorithm = "classic", statistic = "ks") #28 terms
blue_resultKS.elim_BP <- runTest(blue_go_BP, algorithm = "elim", statistic = "ks") #8 terms
blue_allRes_BP <- GenTable(blue_go_BP, classicFisher = blue_resultFisher_BP,
                           classicKS = blue_resultKS_BP, elimKS = blue_resultKS.elim_BP,
                           orderBy = "classicKS", ranksOf = "classicKS", topNodes = 10)

#blue module: plot GO enrichment analysis result
blue_GO_top10_all <- rbind(blue_allRes_CC[,c(2,3,7)],blue_allRes_BP[,c(2,3,7)],blue_allRes_MF[,c(2,3,7)])
blue_GO_top10_all$category <- c("CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC",
                                "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", 
                                "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF")
blue_GO_top10_all$Annotated <- as.numeric(blue_GO_top10_all$Annotated)
blue_GO_top10_all$classicKS <- as.numeric(blue_GO_top10_all$classicKS)
blue_GO_top10_all$category <- as.factor(blue_GO_top10_all$category)
library(ggplot2)
blue_GO_plot <- ggplot(blue_GO_top10_all, aes(x=Annotated, y=reorder(Term, as.numeric(category)), size = Annotated)) +
  geom_point(alpha=0.7, aes(col = classicKS)) +
  ggtitle("GO enrichment analysis \n blue module") + 
  xlab("number of genes") + 
  ylab("GO terms") + 
  labs(col = "p-value") +
  labs(size = "gene numbers") +
  theme(axis.text = element_text(face = "bold")) +
  scale_color_gradient(low="blue", high="red") +
  geom_tile(aes(width = Inf, fill = category), alpha = 0.4) + 
  scale_fill_manual(values = c("azure2", "beige", "darkseagreen2")) #set background colors for the three categories

#repeated for module red and yellow

#GO enrichment analysis for hub genes
#neonicotinoid treatment unique hubs
HT_list <- factor(as.integer(gene_universe$gene_id %in% HT_common$gene))
names(HT_list) <- gene_universe$gene_id
#Biological process
HT_go_BP <- new("topGOdata",
                description = "GO_test_BF",
                ontology    = "BP",
                allGenes    = HT_list,
                annot       = annFUN.gene2GO,
                gene2GO     = Bt_geneID2GO,
                nodeSize    = 10) # Modify to reduce/increase stringency.
HT_resultFisher_BP <- runTest(HT_go_BP, algorithm = "classic", statistic = "fisher") #83 term found
HT_resultKS_BP <- runTest(HT_go_BP, algorithm = "classic", statistic = "ks") #3 terms
HT_resultKS.elim_BP <- runTest(HT_go_BP, algorithm = "elim", statistic = "ks") #1 terms
HT_allRes_BP <- GenTable(HT_go_BP, classicFisher = HT_resultFisher_BP,
                         classicKS = HT_resultKS_BP, elimKS = HT_resultKS.elim_BP,
                         orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)
#cellular component
HT_go_CC <- new("topGOdata",
                description = "GO_test_BF",
                ontology    = "CC",
                allGenes    = HT_list,
                annot       = annFUN.gene2GO,
                gene2GO     = Bt_geneID2GO,
                nodeSize    = 10) # Modify to reduce/increase stringency.
HT_resultFisher_CC <- runTest(HT_go_CC, algorithm = "classic", statistic = "fisher") #83 term found
HT_resultKS_CC <- runTest(HT_go_CC, algorithm = "classic", statistic = "ks") #3 terms
HT_resultKS.elim_CC <- runTest(HT_go_CC, algorithm = "elim", statistic = "ks") #1 terms
HT_allRes_CC <- GenTable(HT_go_CC, classicFisher = HT_resultFisher_CC,
                         classicKS = HT_resultKS_CC, elimKS = HT_resultKS.elim_CC,
                         orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)
#molecular function
HT_go_MF <- new("topGOdata",
                description = "GO_test_BF",
                ontology    = "MF",
                allGenes    = HT_list,
                annot       = annFUN.gene2GO,
                gene2GO     = Bt_geneID2GO,
                nodeSize    = 10) # Modify to reduce/increase stringency.
HT_resultFisher_MF <- runTest(HT_go_MF, algorithm = "classic", statistic = "fisher") #83 term found
HT_resultKS_MF <- runTest(HT_go_MF, algorithm = "classic", statistic = "ks") #3 terms
HT_resultKS.elim_MF <- runTest(HT_go_MF, algorithm = "elim", statistic = "ks") #1 terms
HT_allRes_MF <- GenTable(HT_go_MF, classicFisher = HT_resultFisher_MF,
                         classicKS = HT_resultKS_MF, elimKS = HT_resultKS.elim_MF,
                         orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)

#repeat for HK = housekeeping hubs, and HC = control unique hubs

#retrive gene ID for a GO id
#get the GO id list for housekeeping hub genes cellular component GO enrichment analysis
HK_CC_TERM <- genesInTerm(HK_go_CC)
#retrive genes related to term "mitochondrion"
HK_CC_TERM["GO:0005739"]
#can retrive for any GO terms interested






