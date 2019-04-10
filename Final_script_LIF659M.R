# Final script for LIF659M
#--------------------------------------------------------#
# Sigurdur Karvelsson, MARCH 2019

# Install packages (if needed):
install.packages('dendextend')
install.packages('colorspace')
install.packages('factoextra')
install.packages('NbClust')
install.packages('gplots')
install.packages('tidyr')
install.packages('survival')
install.packages('survminer')
install.packages('ggforce')
install.packages('igraph')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AnnotationHub", version = "3.8")
BiocManager::install("clusterProfiler", version = "3.8")
BiocManager::install("org.Hs.eg.db", version = "3.8")
BiocManager::install("topGO", version = "3.8")
BiocManager::install("DOSE", version = "3.8")
BiocManager::install("ReactomePA", version = "3.8")

# Load packages:
library(dendextend)
library(colorspace)
library(factoextra)
library(NbClust)
library(gplots)
library(tidyr)
library(survminer)
library(survival)
library(BiocManager)
library(clusterProfiler)
library(org.Hs.eg.db)
library(topGO)
library(DOSE)
library(ReactomePA)
library(ggforce)
library(igraph)
#--------------------------------------------------------#
# PARAMETERS:
# How many clusters? Range is 2-8.
k_used = 4
patients = 'all' # as opposed to 'TNBC'
clustering_method = 'hierarchical' # as opposed to 'k_means'
genes = 'top_variance' # as opposed to 'top_variance'
save_images = TRUE
# wd = folder with data, wd2 = folder for exporting figures
wd <- "C:/Users/sigur/Dropbox/R_Dropbox"
wd2 <- "C:/Users/sigur/Dropbox/R_Dropbox/OUTPUT_FOLDER"
setwd(wd)
source("./STKfunctions/loadAll.r")
loadAll()
source("./FJfunctions/loadAll.r")
loadAll()
set.seed(222)
# Color palette used throughout:
col4 <- ColBrew('JCO',k_used)
col5 <- ColBrew('JCO')
#--------------------------------------------------------#
## 1. Load all the required data
# Clinical:
clinical_data <- read.table('data_bcr_clinical_data_patient.txt',header=T,sep='\t',fill=T)
# Extract TNBC patients:
id_ER <- which(clinical_data$ER_STATUS_BY_IHC == 'Negative')
id_PR <- which(clinical_data$PR_STATUS_BY_IHC == 'Negative')
id_HER2 <- union(which(clinical_data$HER2_IHC_SCORE == 1), which(clinical_data$HER2_IHC_SCORE == 0))
id_TNBC <- union(union(id_ER,id_PR),id_HER2)
TNBC_patients <- clinical_data[,c("PATIENT_ID","OS_MONTHS","OS_STATUS","RACE","DFS_MONTHS")]
TNBC_patients$PATIENT_ID <- gsub("-",".",TNBC_patients$PATIENT_ID)
TNBC_patients$PATIENT_ID <- paste0(TNBC_patients$PATIENT_ID,'.01')

if (patients == 'all'){
  clin_data <- TNBC_patients
} else if (patients == 'TNBC'){
  clin_data <- TNBC_patients[id_TNBC,]
}

# RNAseq data:
raw_rna <- read.table('data_RNA_Seq_v2_expression_median.txt',header=T)
id_used1 <- which(colnames(raw_rna) %in% clin_data$PATIENT_ID)
raw_rna <- cbind(raw_rna[,c(1,2)],raw_rna[,id_used1]) #only use the TNBC patients

# Load metabolic genes:
MES_genes <- read.table("REACTOME_metabolic_genes.txt",header=T)
MES_genes <- data.frame(MES_genes$hgnc_symbol)

# Load RECON2 genes:
recon_genes <- read.table('Recon_genes.txt',header=T)

# Extract rnaseq data from TNBC patients:
idx <- which(raw_rna$Hugo_Symbol %in% MES_genes$MES_genes.hgnc_symbol)
idx_recon <- which(raw_rna$Hugo_Symbol %in% recon_genes$recon_genes)
MMS_rna_data <- raw_rna[idx,]
#rownames(MMS_rna_data) <- MMS_rna_data$Hugo_Symbol
row_sub = apply(MMS_rna_data, 1, function(row) all(row !=0 ))
#MMS_rna_data <- MMS_rna_data[ - as.numeric(which(apply(MMS_rna_data, 2, var) == 0))]
MMS_rna_data <- MMS_rna_data[row_sub,]
rownames(MMS_rna_data) <- make.names(MMS_rna_data$Hugo_Symbol, unique=TRUE)

# Check which genes have highest variability - use top ones for clustering...
checker3 <- as.matrix(MMS_rna_data[,c(-1,-2)])
rownames(checker3) <- MMS_rna_data[,1]
coef.var <- rep(0,nrow(checker3))
for (i in 1:nrow(checker3)){
  coef.var[i] = sd(checker3[i,])/mean(checker3[i,])
}

# Check what genes we wanted to use:
if (genes == 'top_variance'){
  MMS_rna_data <- MMS_rna_data[which(coef.var>0.1),]
} else if (genes == 'all'){
  MMS_rna_data <- MMS_rna_data
}
#--------------------------------------------------------#
## 2. Show that the data is not dependent on race nor total survival:
MES_data_1 = prepareRNAdata(MMS_rna_data,clin_data,"survival",'Hugo')
MES_data_2 = prepareRNAdata(MMS_rna_data,clin_data,"race",'Hugo')
PCA_survival <- metPCA(MES_data_1,groups="three_groups",scale="auto",colors=col5,pcs=c(1,2),narm=T,minvalue=T,log=T) + ggtitle('Survival and expression of metabolic genes')
PCA_survival
if (save_images){
  setwd(wd2)
  ggsave('PCA_survival_MS_genes_MARCH2019_with_errors.png',width = 8, height = 6, units = "in")
  PCA_race <- metPCA(MES_data_2,groups="races",scale="auto",colors=col5,pcs=c(1,2),narm=T,minvalue=T,log=T) + ggtitle('Race and expression of metabolic genes')
  PCA_race
  ggsave('PCA_race_MS_genes_MARCH2019_with_errors.png',width = 8, height = 6, units = "in")
  setwd(wd)
}
PCA_race <- metPCA(MES_data_2,groups="races",scale="auto",colors=col5,pcs=c(1,2),narm=T,minvalue=T,log=T) + ggtitle('Survival and expression of metabolic genes')
PCA_race
#--------------------------------------------------------#
## 3. Cluster the data into 6 clusters:
MMS_z_data <- standardizeTCGA(MMS_rna_data)
Clustering_MMS <- Hierarchical_clustering(MMS_z_data,'Hugo','Spearman')
# Extract k clusters:
mycl <- cutree(Clustering_MMS,k=k_used)
df <- t(MMS_z_data[,3:ncol(MMS_z_data)])
df2 <- df[, apply(df, 2, function(x) !any(is.na(x)))] 
df2 <- data.frame(df2)
dd <- cbind(mycl,df2)
colnames(dd)[1] = 'clusters'
dd$clusters <- as.factor(dd$clusters)
input = rep(0,nrow(dd))
dd_fin <- cbind(dd$clusters,input,dd[,2:ncol(dd)])
colnames(dd_fin)[1] = 'clusters'
PCA <- metPCA(dd_fin,groups="clusters",colors=col4,pcs=c(1,2))
PCA <- PCA + ggtitle('PCA of hierarchical clustering')
PCA
if (save_images){
  setwd(wd2)
  ggsave('PCA_Hierarchical_clustering.png',width = 8, height = 6, units = "in")
  setwd(wd)
}

#--------------------------------------------------------#
## 4. Check optimal cluster amount:
# Silhouette method
fviz_nbclust(df2, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")
#Elbow method
fviz_nbclust(df2, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")
if (save_images){
  setwd(wd2)
  ggsave('Optimum_k_elbow_method.png',width = 8, height = 6, units = "in")
  setwd(wd)
}
km.res <- kmeans(df2, k_used, nstart = 50)
aggregate(df2, by=list(cluster=km.res$cluster), mean)
dd <- cbind(df2, cluster = km.res$cluster)
km.res$size
dd_outp <- cbind(dd$cluster,dd[,1:ncol(dd)-1])
colnames(dd_outp)[1] = 'clusters'
dd_outp$clusters <- as.factor(dd_outp$clusters)
input = rep(0,nrow(dd_outp))
dd_fin_kmeans <- cbind(dd_outp$clusters,input,dd_outp[,2:ncol(dd_outp)])
colnames(dd_fin_kmeans)[1] = 'clusters'

PCA2 <- metPCA(dd_fin_kmeans,groups="clusters",scale="auto",colors=col4,pcs=c(1,2))
PCA2 <- PCA2 + ggtitle('PCA of k-means clustering')
PCA2
if (save_images){
  setwd(wd2)
  ggsave('PCA_kmeans_clustering.png',width = 8, height = 6, units = "in")
  setwd(wd)
}
#--------------------------------------------------------#
## 5. Plot hierarchical clustering (k = 6)
iris_species <- levels(as.factor(mycl))
dend <- as.dendrogram(Clustering_MMS)
dend2 <- color_branches(dend, k=k_used,col =col4[unique(as.numeric(mycl)[order.dendrogram(dend)])])
dend3 <- color_branches(dend, k=k_used,col =col4)
labels_colors(dend2) <- col4[as.numeric(mycl)[order.dendrogram(dend)]] 
dend2 <- hang.dendrogram(dend2,hang_height=0.1)
dend2 <- set(dend2, "labels_cex", 0.5)
if (save_images){
  setwd(wd2)
  png("Hierarchical_clustering_k_used_APRIL2019.png",width = 600, height = 600, units = "px")
  plot(dend2, 
       main = "Hierarchical clustering", 
       horiz =  TRUE,  nodePar = list(cex = .007))
  legend("topleft", legend = iris_species, fill = col4,
         title = 'Clusters')
  dev.off()
  setwd(wd)
}

#--------------------------------------------------------#
## 6. Plot a heatmap:
# choose a color palette for the heat map
myheatcol <- bluered(100)
MMS_z_data2 = MMS_z_data[order(MMS_z_data[,3]),]
rownames(MMS_z_data2) <- make.names(MMS_z_data2$Hugo_Symbol, unique=TRUE)
checker <- MMS_z_data2
mat_used =as.matrix(t(checker[,c(-1,-2)]))
mat_used <- mat_used[,colSums(is.na(mat_used))<nrow(mat_used)]
checker2 <- as.matrix(checker[,c(-1,-2)])
cor_mat <- cor(mat_used,method = 'spearman') # for usage if one wants to check the correlation of geenes
if (save_images){
  setwd(wd2)
  png("Pairwise_spearman_correlation_of_metabolic_genes_APRIL2019.png",width = 800, height = 600, units = "px") 
  heatmap.2(cor_mat, main="Pairwise Spearman-correlation of metabolic genes",symm = T,
            col = myheatcol,trace="none",scale="none")
  dev.off()
} else {
  heatmap.2(cor_mat, main="Pairwise Spearman-correlation of metabolic genes",symm = T,
            col = myheatcol,trace="none",scale="none")
}
# Now plot genes vs patients:
dist.spear <- function(x) as.dist(1-cor(t(x),method='spearman'))
hclust.ward <- function(x) hclust(x, method="ward.D2")
if (save_images){
  png("TNBC_Patients_vs_MS_genes_APRIL2019.png",width = 800, height = 600, units = "px") 
  heatmap.2(checker2, main="Patients clustered by metabolic gene expression",
            hclustfun = hclust.ward,distfun = dist.spear,
            col = myheatcol,trace="none",scale="none", ColSideColors=col4[mycl])
  dev.off()
  setwd(wd)
} else {
  heatmap.2(checker2, main="Patients clustered by metabolic gene expression",
            hclustfun = hclust.ward,distfun = dist.spear,
            col = myheatcol,trace="none",scale="none", ColSideColors=col4[mycl])
}
#--------------------------------------------------------#
# 8. Extract the significant genes from the raw RNA seq data, based on clustering on Metabolic genes:
if (clustering_method == 'k_means'){
  dd_fin <- dd_fin_kmeans
} else if (clustering_method == 'hierarchical'){
  dd_fin <- dd_fin
}
rownames(raw_rna) <- rownames(raw_rna) <- make.names(as.character(raw_rna$Hugo_Symbol),TRUE)
total_rna_z_scored = data.frame(t(raw_rna))
total_rna_z_scored <- total_rna_z_scored[sapply(total_rna_z_scored, function(x) !any(is.na(x)))] 
total_rna_z_scored <- total_rna_z_scored[!(rowSums(is.na(total_rna_z_scored))),]

# remove from dataframe columns containing only 0:
for (i in 1:ncol(total_rna_z_scored)){
  total_rna_z_scored[,i] <- as.numeric(as.character(total_rna_z_scored[,i]))
}
total_rna_z_scored <- total_rna_z_scored[!(rowSums(is.na(total_rna_z_scored))),]
total_rna_z_scored <- total_rna_z_scored[, colSums(total_rna_z_scored,na.rm=T) != 0] 

dd_cluster1 <- dd[which(dd_fin$clusters==1),]
dd_cluster2 <- dd[which(dd_fin$clusters==2),]
dd_cluster3 <- dd[which(dd_fin$clusters==3),]
dd_cluster4 <- dd[which(dd_fin$clusters==4),]
dd_cluster5 <- dd[which(dd_fin$clusters==5),]
dd_cluster6 <- dd[which(dd_fin$clusters==6),]
dd_cluster7 <- dd[which(dd_fin$clusters==7),]
dd_cluster8 <- dd[which(dd_fin$clusters==8),]

jojo <- rownames(total_rna_z_scored)
jojo1 <- rownames(dd_cluster1)
jojo2 <- rownames(dd_cluster2)
jojo3 <- rownames(dd_cluster3)
jojo4 <- rownames(dd_cluster4)
jojo5 <- rownames(dd_cluster5)
jojo6 <- rownames(dd_cluster6)
jojo7 <- rownames(dd_cluster7)
jojo8 <- rownames(dd_cluster8)

dd_cluster1_rna <- total_rna_z_scored[which(jojo %in% jojo1),]
dd_cluster2_rna <- total_rna_z_scored[which(jojo %in% jojo2),]
dd_cluster3_rna <- total_rna_z_scored[which(jojo %in% jojo3),]
dd_cluster4_rna <- total_rna_z_scored[which(jojo %in% jojo4),]
dd_cluster5_rna <- total_rna_z_scored[which(jojo %in% jojo5),]
dd_cluster6_rna <- total_rna_z_scored[which(jojo %in% jojo6),]
dd_cluster7_rna <- total_rna_z_scored[which(jojo %in% jojo7),]
dd_cluster8_rna <- total_rna_z_scored[which(jojo %in% jojo8),]

if (k_used == 2){
  dd_cluster_rna_all <- rbind(dd_cluster1_rna,dd_cluster2_rna)
} else if (k_used == 3){
  dd_cluster_rna_all <- rbind(dd_cluster1_rna,dd_cluster2_rna,dd_cluster3_rna)
} else if (k_used == 4){
  dd_cluster_rna_all <- rbind(dd_cluster1_rna,dd_cluster2_rna,dd_cluster3_rna,dd_cluster4_rna)
} else if (k_used == 5){
  dd_cluster_rna_all <- rbind(dd_cluster1_rna,dd_cluster2_rna,dd_cluster3_rna,dd_cluster4_rna,dd_cluster5_rna)
} else if (k_used == 6){
  dd_cluster_rna_all <- rbind(dd_cluster1_rna,dd_cluster2_rna,dd_cluster3_rna,dd_cluster4_rna,dd_cluster5_rna,dd_cluster6_rna)
} else if (k_used == 7){
  dd_cluster_rna_all <- rbind(dd_cluster1_rna,dd_cluster2_rna,dd_cluster3_rna,dd_cluster4_rna,dd_cluster5_rna,dd_cluster6_rna,dd_cluster7_rna)
} else if (k_used == 8){
  dd_cluster_rna_all <- rbind(dd_cluster1_rna,dd_cluster2_rna,dd_cluster3_rna,dd_cluster4_rna,dd_cluster5_rna,dd_cluster6_rna,dd_cluster7_rna,dd_cluster8_rna)
}

cluster1 <- extractDEgenesTCGA(dd_cluster1_rna,
                               not_in(dd_cluster_rna_all,dd_cluster1_rna),
                               1.5,0.05,1)
cluster2 <- extractDEgenesTCGA(dd_cluster2_rna,
                               not_in(dd_cluster_rna_all,dd_cluster2_rna),
                               1.5,0.05,1)
cluster3 <- extractDEgenesTCGA(dd_cluster3_rna,
                               not_in(dd_cluster_rna_all,dd_cluster3_rna),
                               1.5,0.05,1)
cluster4 <- extractDEgenesTCGA(dd_cluster4_rna,
                               not_in(dd_cluster_rna_all,dd_cluster4_rna),
                               1.5,0.05,1)
cluster5 <- extractDEgenesTCGA(dd_cluster5_rna,
                               not_in(dd_cluster_rna_all,dd_cluster5_rna),
                               1.5,0.05,1)
cluster6 <- extractDEgenesTCGA(dd_cluster6_rna,
                               not_in(dd_cluster_rna_all,dd_cluster6_rna),
                               1.5,0.05,1)
cluster7 <- extractDEgenesTCGA(dd_cluster7_rna,
                               not_in(dd_cluster_rna_all,dd_cluster7_rna),
                               1.5,0.05,1)
cluster8 <- extractDEgenesTCGA(dd_cluster8_rna,
                               not_in(dd_cluster_rna_all,dd_cluster8_rna),
                               1.5,0.05,1)

# Export cluster results (optional):
if (save_images){
  setwd(wd2)
  write.table(cluster1, file = "cluster1_signifant_genes_260319.txt", sep = "\t",
              row.names = FALSE,quote=FALSE)
  write.table(cluster2, file = "cluster2_signifant_genes_260319.txt", sep = "\t",
              row.names = FALSE,quote=FALSE)
  write.table(cluster3, file = "cluster3_signifant_genes_260319.txt", sep = "\t",
              row.names = FALSE,quote=FALSE)
  write.table(cluster4, file = "cluster4_signifant_genes_260319.txt", sep = "\t",
              row.names = FALSE,quote=FALSE)
  write.table(cluster5, file = "cluster5_signifant_genes_260319.txt", sep = "\t",
              row.names = FALSE,quote=FALSE)
  write.table(cluster6, file = "cluster6_signifant_genes_260319.txt", sep = "\t",
              row.names = FALSE,quote=FALSE)
  write.table(cluster7, file = "cluster7_signifant_genes_260319.txt", sep = "\t",
              row.names = FALSE,quote=FALSE)
  write.table(cluster8, file = "cluster8_signifant_genes_260319.txt", sep = "\t",
              row.names = FALSE,quote=FALSE)
  setwd(wd)
}

#--------------------------------------------------------#
# 9. GO-analysis:
if (k_used == 2){
  cluster1_entrez <- mapIds(org.Hs.eg.db, as.character(cluster1$Group1), 'ENTREZID', 'SYMBOL')
  cluster2_entrez <- mapIds(org.Hs.eg.db, as.character(cluster2$Group1), 'ENTREZID', 'SYMBOL')
  gene_list <- list(as.character(unlist(cluster1_entrez,use.names = F)),
                    as.character(unlist(cluster2_entrez,use.names = F)))
  names(gene_list) <- c('CL1','CL2')
} else if (k_used == 3){
  cluster1_entrez <- mapIds(org.Hs.eg.db, as.character(cluster1$Group1), 'ENTREZID', 'SYMBOL')
  cluster2_entrez <- mapIds(org.Hs.eg.db, as.character(cluster2$Group1), 'ENTREZID', 'SYMBOL')
  cluster3_entrez <- mapIds(org.Hs.eg.db, as.character(cluster3$Group1), 'ENTREZID', 'SYMBOL')
  gene_list <- list(as.character(unlist(cluster1_entrez,use.names = F)),
                    as.character(unlist(cluster2_entrez,use.names = F)),
                    as.character(unlist(cluster3_entrez,use.names = F)))
  names(gene_list) <- c('CL1','CL2','CL3')
} else if (k_used == 4){
  cluster1_entrez <- mapIds(org.Hs.eg.db, as.character(cluster1$Group1), 'ENTREZID', 'SYMBOL')
  cluster2_entrez <- mapIds(org.Hs.eg.db, as.character(cluster2$Group1), 'ENTREZID', 'SYMBOL')
  cluster3_entrez <- mapIds(org.Hs.eg.db, as.character(cluster3$Group1), 'ENTREZID', 'SYMBOL')
  cluster4_entrez <- mapIds(org.Hs.eg.db, as.character(cluster4$Group1), 'ENTREZID', 'SYMBOL')
  gene_list <- list(as.character(unlist(cluster1_entrez,use.names = F)),
                    as.character(unlist(cluster2_entrez,use.names = F)),
                    as.character(unlist(cluster3_entrez,use.names = F)),
                    as.character(unlist(cluster4_entrez,use.names = F)))
  names(gene_list) <- c('CL1','CL2','CL3','CL4')
} else if (k_used == 5){
  cluster1_entrez <- mapIds(org.Hs.eg.db, as.character(cluster1$Group1), 'ENTREZID', 'SYMBOL')
  cluster2_entrez <- mapIds(org.Hs.eg.db, as.character(cluster2$Group1), 'ENTREZID', 'SYMBOL')
  cluster3_entrez <- mapIds(org.Hs.eg.db, as.character(cluster3$Group1), 'ENTREZID', 'SYMBOL')
  cluster4_entrez <- mapIds(org.Hs.eg.db, as.character(cluster4$Group1), 'ENTREZID', 'SYMBOL')
  cluster5_entrez <- mapIds(org.Hs.eg.db, as.character(cluster5$Group1), 'ENTREZID', 'SYMBOL')
  gene_list <- list(as.character(unlist(cluster1_entrez,use.names = F)),
                    as.character(unlist(cluster2_entrez,use.names = F)),
                    as.character(unlist(cluster3_entrez,use.names = F)),
                    as.character(unlist(cluster4_entrez,use.names = F)),
                    as.character(unlist(cluster5_entrez,use.names = F)))
  names(gene_list) <- c('CL1','CL2','CL3','CL4','CL5')
} else if (k_used == 6){
  cluster1_entrez <- mapIds(org.Hs.eg.db, as.character(cluster1$Group1), 'ENTREZID', 'SYMBOL')
  cluster2_entrez <- mapIds(org.Hs.eg.db, as.character(cluster2$Group1), 'ENTREZID', 'SYMBOL')
  cluster3_entrez <- mapIds(org.Hs.eg.db, as.character(cluster3$Group1), 'ENTREZID', 'SYMBOL')
  cluster4_entrez <- mapIds(org.Hs.eg.db, as.character(cluster4$Group1), 'ENTREZID', 'SYMBOL')
  cluster5_entrez <- mapIds(org.Hs.eg.db, as.character(cluster5$Group1), 'ENTREZID', 'SYMBOL')
  cluster6_entrez <- mapIds(org.Hs.eg.db, as.character(cluster6$Group1), 'ENTREZID', 'SYMBOL')
  gene_list <- list(as.character(unlist(cluster1_entrez,use.names = F)),
                    as.character(unlist(cluster2_entrez,use.names = F)),
                    as.character(unlist(cluster3_entrez,use.names = F)),
                    as.character(unlist(cluster4_entrez,use.names = F)),
                    as.character(unlist(cluster5_entrez,use.names = F)),
                    as.character(unlist(cluster6_entrez,use.names = F)))
  names(gene_list) <- c('CL1','CL2','CL3','CL4','CL5','CL6')
} else if (k_used == 7){
  cluster1_entrez <- mapIds(org.Hs.eg.db, as.character(cluster1$Group1), 'ENTREZID', 'SYMBOL')
  cluster2_entrez <- mapIds(org.Hs.eg.db, as.character(cluster2$Group1), 'ENTREZID', 'SYMBOL')
  cluster3_entrez <- mapIds(org.Hs.eg.db, as.character(cluster3$Group1), 'ENTREZID', 'SYMBOL')
  cluster4_entrez <- mapIds(org.Hs.eg.db, as.character(cluster4$Group1), 'ENTREZID', 'SYMBOL')
  cluster5_entrez <- mapIds(org.Hs.eg.db, as.character(cluster5$Group1), 'ENTREZID', 'SYMBOL')
  cluster6_entrez <- mapIds(org.Hs.eg.db, as.character(cluster6$Group1), 'ENTREZID', 'SYMBOL')
  cluster7_entrez <- mapIds(org.Hs.eg.db, as.character(cluster7$Group1), 'ENTREZID', 'SYMBOL')
  gene_list <- list(as.character(unlist(cluster1_entrez,use.names = F)),
                    as.character(unlist(cluster2_entrez,use.names = F)),
                    as.character(unlist(cluster3_entrez,use.names = F)),
                    as.character(unlist(cluster4_entrez,use.names = F)),
                    as.character(unlist(cluster5_entrez,use.names = F)),
                    as.character(unlist(cluster6_entrez,use.names = F)),
                    as.character(unlist(cluster7_entrez,use.names = F)))
  names(gene_list) <- c('CL1','CL2','CL3','CL4','CL5','CL6','CL7')
} else if (k_used == 8){
  cluster1_entrez <- mapIds(org.Hs.eg.db, as.character(cluster1$Group1), 'ENTREZID', 'SYMBOL')
  cluster2_entrez <- mapIds(org.Hs.eg.db, as.character(cluster2$Group1), 'ENTREZID', 'SYMBOL')
  cluster3_entrez <- mapIds(org.Hs.eg.db, as.character(cluster3$Group1), 'ENTREZID', 'SYMBOL')
  cluster4_entrez <- mapIds(org.Hs.eg.db, as.character(cluster4$Group1), 'ENTREZID', 'SYMBOL')
  cluster5_entrez <- mapIds(org.Hs.eg.db, as.character(cluster5$Group1), 'ENTREZID', 'SYMBOL')
  cluster6_entrez <- mapIds(org.Hs.eg.db, as.character(cluster6$Group1), 'ENTREZID', 'SYMBOL')
  cluster7_entrez <- mapIds(org.Hs.eg.db, as.character(cluster7$Group1), 'ENTREZID', 'SYMBOL')
  cluster8_entrez <- mapIds(org.Hs.eg.db, as.character(cluster8$Group1), 'ENTREZID', 'SYMBOL')
  gene_list <- list(as.character(unlist(cluster1_entrez,use.names = F)),
                    as.character(unlist(cluster2_entrez,use.names = F)),
                    as.character(unlist(cluster3_entrez,use.names = F)),
                    as.character(unlist(cluster4_entrez,use.names = F)),
                    as.character(unlist(cluster5_entrez,use.names = F)),
                    as.character(unlist(cluster6_entrez,use.names = F)),
                    as.character(unlist(cluster7_entrez,use.names = F)),
                    as.character(unlist(cluster8_entrez,use.names = F)))
  names(gene_list) <- c('CL1','CL2','CL3','CL4','CL5','CL6','CL7','CL8')
}

x=compareCluster(gene_list, fun='enrichGO',OrgDb='org.Hs.eg.db',ont='BP')
x2=compareCluster(gene_list, fun='enrichPathway')
x3=compareCluster(gene_list,fun="enrichKEGG")
dotplot(x, showCategory=5, includeAll=FALSE,title = 'GO (Biological process) enriched in Clusters')
if (save_images){
  setwd(wd2)
  ggsave('GO_BP_Cluster1_to_Clusterk_APRIL2019.png',width = 10, height = 8, units = "in")
  setwd(wd)
}
dotplot(x2, showCategory=5, includeAll=FALSE,title = 'Reactome pathways enriched in Clusters')
if (save_images){
  setwd(wd2)
  ggsave('Reactome_Cluster1_to_Clusterk_APRIL2019.png',width = 14, height = 8, units = "in")
  setwd(wd)
}
dotplot(x3, showCategory=5, includeAll=FALSE,title = 'KEGG pathways enriched in Clusters')
if (save_images){
  setwd(wd2)
  ggsave('KEGG_Cluster1_to_Clusterk_APRIL2019.png',width = 14, height = 8, units = "in")
  setwd(wd)
}

#--------------------------------------------------------#
# 10. Survival analysis:
# Identify what cluster every patient belongs to:
#idx_clind <- which(clin_data$PATIENT_ID %in% rownames(dd_fin))
index_vector <- c()
for (i in 1:dim(MMS_rna_data)[2]){
  idx <- which(as.character(clin_data$PATIENT_ID) == colnames(MMS_rna_data)[i])
  if (length(idx) != 0){
    index_vector <- c(index_vector,idx)
  }
}

idx_clind <- index_vector
idx_dd_fin <- index_vector2
id_living = which(clin_data$OS_STATUS == "LIVING")
id_dead = which(clin_data$OS_STATUS == "DECEASED")
status <- rep(0,dim(clin_data)[1])
status[id_living] = 1
status[id_dead] = 2
months <- suppressWarnings(as.numeric(as.character(clin_data$OS_MONTHS[idx_clind])))
months3 <- suppressWarnings(as.numeric(as.character(clin_data$DFS_MONTHS[idx_clind])))
months <- ceiling(months)
months3 <- ceiling(months3)
cor_data2 <- data.frame(cbind(clin_data$PATIENT_ID[idx_clind],status[idx_clind],months))

index_vector2 <- c()
for (i in 1:length(cor_data2$V1)){
  idx2 <- which(rownames(dd_fin) == cor_data2$V1[i])
  if (length(idx2) != 0){
    index_vector2 <- c(index_vector2,idx2)
  }
}
cor_data2 <- cbind(cor_data2,dd_fin$clusters[index_vector2])
colnames(cor_data2) = c('patient','status','survival','cluster')
cor_data2$survival <- suppressWarnings(as.numeric(as.character(cor_data2$survival)))
cor_data2$status <- as.numeric(as.character(cor_data2$status))

# Try to include columns with one cluster against others:
cor_data2 <- cbind(cor_data2,rep(0,dim(cor_data2)[1]))
colnames(cor_data2) = c('patient','status','survival','cluster','CLUSTER2')
id_cl1 <- which(cor_data2$cluster == 1)
id_cl_not1 <- which(cor_data2$cluster != 1)
cor_data2$CLUSTER2[id_cl1] = 1
cor_data2$CLUSTER2[id_cl_not1] = 2

res.cox <- coxph(Surv(survival, status) ~ cluster, data = cor_data2)
fit <- survfit(Surv(survival, status) ~ cluster, data = cor_data2)
p <- ggsurvplot(fit, data = cor_data2, risk.table=T,pval=T, conf.int=T,
                xlim = c(0,120),xlab = 'Time in months',
                break.time.by = 12, ggtheme = theme_light())
p <- ggpar(p,font.x = c(14, "bold", "black"),
           font.y = c(14, "bold", "black"))
p

fit_1 <- survfit(Surv(survival, status) ~ CLUSTER2, data = cor_data2)
p_1 <- ggsurvplot(fit_1, data = cor_data2, risk.table=T,pval=T, conf.int=T,
                xlim = c(0,120),xlab = 'Time in months',
                break.time.by = 12, 
                legend.labs = c('cluster=1','Rest'),
                ggtheme = theme_light())
p_1 <- ggpar(p_1,font.x = c(14, "bold", "black"),
           font.y = c(14, "bold", "black"))
p_1




if (save_images){
  setwd(wd2)
  ggsave("COXPH_clusters_MARCH2019.png", plot=print(p), width = 8, height = 6, dpi = 1000)
  ggsave("COXPH_cluster1_vs_rest_MARCH2019.png", plot=print(p_1), width = 8, height = 6, dpi = 1000)
  setwd(wd)
}

#--------------------------------------------------------#
# 11. Check what MMS genes are significantly different between groups:
MMS_cluster1 <- MMS_rna_data[,which(colnames(MMS_rna_data) %in% jojo1)]
MMS_cluster2 <- MMS_rna_data[,which(colnames(MMS_rna_data) %in% jojo2)]
MMS_cluster3 <- MMS_rna_data[,which(colnames(MMS_rna_data) %in% jojo3)]
MMS_cluster4 <- MMS_rna_data[,which(colnames(MMS_rna_data) %in% jojo4)]
MMS_cluster5 <- MMS_rna_data[,which(colnames(MMS_rna_data) %in% jojo5)]
MMS_cluster6 <- MMS_rna_data[,which(colnames(MMS_rna_data) %in% jojo6)]
MMS_cluster7 <- MMS_rna_data[,which(colnames(MMS_rna_data) %in% jojo7)]
MMS_cluster8 <- MMS_rna_data[,which(colnames(MMS_rna_data) %in% jojo8)]

# Check differential expression between the two clusters, and then group them further...
MMS_1 <- data.frame(t(MMS_cluster1))
MMS_2 <- data.frame(t(MMS_cluster2))
MMS_3 <- data.frame(t(MMS_cluster3))
MMS_4 <- data.frame(t(MMS_cluster4))
MMS_5 <- data.frame(t(MMS_cluster5))
MMS_6 <- data.frame(t(MMS_cluster6))
MMS_7 <- data.frame(t(MMS_cluster7))
MMS_8 <- data.frame(t(MMS_cluster8))
if (k_used == 2){
  MMS_all <- rbind(MMS_1,MMS_2)
} else if (k_used == 3){
  MMS_all <- rbind(MMS_1,MMS_2,MMS_3)
} else if (k_used == 4){
  MMS_all <- rbind(MMS_1,MMS_2,MMS_3,MMS_4)
} else if (k_used == 5){
  MMS_all <- rbind(MMS_1,MMS_2,MMS_3,MMS_4,MMS_5)
} else if (k_used == 6){
  MMS_all <- rbind(MMS_1,MMS_2,MMS_3,MMS_4,MMS_5,MMS_6)
} else if (k_used == 7){
  MMS_all <- rbind(MMS_1,MMS_2,MMS_3,MMS_4,MMS_5,MMS_6,MMS_7)
} else if (k_used == 8){
  MMS_all <- rbind(MMS_1,MMS_2,MMS_3,MMS_4,MMS_5,MMS_6,MMS_7,MMS_8)
}

cluster1_2 <- extractDEgenesTCGA(MMS_1,not_in(MMS_all,MMS_1),0.3,0.05,1)
cluster2_2 <- extractDEgenesTCGA(MMS_2,not_in(MMS_all,MMS_2),0.3,0.05,1)
cluster3_2 <- extractDEgenesTCGA(MMS_3,not_in(MMS_all,MMS_3),0.3,0.05,1)
cluster4_2 <- extractDEgenesTCGA(MMS_4,not_in(MMS_all,MMS_4),0.3,0.05,1)
cluster5_2 <- extractDEgenesTCGA(MMS_5,not_in(MMS_all,MMS_5),0.3,0.05,1)
cluster6_2 <- extractDEgenesTCGA(MMS_6,not_in(MMS_all,MMS_6),0.3,0.05,1)
cluster7_2 <- extractDEgenesTCGA(MMS_7,not_in(MMS_all,MMS_7),0.3,0.05,1)
cluster8_2 <- extractDEgenesTCGA(MMS_8,not_in(MMS_all,MMS_8),0.3,0.05,1)


# Export cluster results (optional):
if (save_images){
  setwd(wd2)
  write.table(cluster1_2, file = "cluster1_2_signifant_genes.txt", sep = "\t",
              row.names = FALSE,quote=FALSE)
  write.table(cluster2_2, file = "cluster2_2_signifant_genes.txt", sep = "\t",
              row.names = FALSE,quote=FALSE)
  write.table(cluster3_2, file = "cluster3_2_signifant_genes.txt", sep = "\t",
              row.names = FALSE,quote=FALSE)
  write.table(cluster4_2, file = "cluster4_2_signifant_genes.txt", sep = "\t",
              row.names = FALSE,quote=FALSE)
  write.table(cluster5_2, file = "cluster5_2_signifant_genes.txt", sep = "\t",
              row.names = FALSE,quote=FALSE)
  write.table(cluster6_2, file = "cluster6_2_signifant_genes.txt", sep = "\t",
              row.names = FALSE,quote=FALSE)
  write.table(cluster7_2, file = "cluster7_2_signifant_genes.txt", sep = "\t",
              row.names = FALSE,quote=FALSE)
  write.table(cluster8_2, file = "cluster8_2_signifant_genes.txt", sep = "\t",
              row.names = FALSE,quote=FALSE)
}

# Now perform a GO-analysis on only the MMS genes:
if (k_used == 2){
  cluster1_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster1_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster2_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster2_2$Group1), 'ENTREZID', 'SYMBOL')
  gene_list2 <- list(as.character(unlist(cluster1_entrez_2,use.names = F)),
                    as.character(unlist(cluster2_entrez_2,use.names = F)))
  names(gene_list2) <- c('CL1','CL2')
} else if (k_used == 3){
  cluster1_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster1_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster2_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster2_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster3_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster3_2$Group1), 'ENTREZID', 'SYMBOL')
  gene_list2 <- list(as.character(unlist(cluster1_entrez_2,use.names = F)),
                    as.character(unlist(cluster2_entrez_2,use.names = F)),
                    as.character(unlist(cluster3_entrez_2,use.names = F)))
  names(gene_list2) <- c('CL1','CL2','CL3')
} else if (k_used == 4){
  cluster1_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster1_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster2_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster2_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster3_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster3_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster4_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster4_2$Group1), 'ENTREZID', 'SYMBOL')
  gene_list2 <- list(as.character(unlist(cluster1_entrez_2,use.names = F)),
                    as.character(unlist(cluster2_entrez_2,use.names = F)),
                    as.character(unlist(cluster3_entrez_2,use.names = F)),
                    as.character(unlist(cluster4_entrez_2,use.names = F)))
  names(gene_list2) <- c('CL1','CL2','CL3','CL4')
} else if (k_used == 5){
  cluster1_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster1_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster2_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster2_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster3_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster3_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster4_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster4_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster5_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster5_2$Group1), 'ENTREZID', 'SYMBOL')
  gene_list2 <- list(as.character(unlist(cluster1_entrez_2,use.names = F)),
                    as.character(unlist(cluster2_entrez_2,use.names = F)),
                    as.character(unlist(cluster3_entrez_2,use.names = F)),
                    as.character(unlist(cluster4_entrez_2,use.names = F)),
                    as.character(unlist(cluster5_entrez_2,use.names = F)))
  names(gene_list2) <- c('CL1','CL2','CL3','CL4','CL5')
} else if (k_used == 6){
  cluster1_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster1_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster2_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster2_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster3_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster3_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster4_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster4_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster5_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster5_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster6_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster6_2$Group1), 'ENTREZID', 'SYMBOL')
  gene_list2 <- list(as.character(unlist(cluster1_entrez_2,use.names = F)),
                    as.character(unlist(cluster2_entrez_2,use.names = F)),
                    as.character(unlist(cluster3_entrez_2,use.names = F)),
                    as.character(unlist(cluster4_entrez_2,use.names = F)),
                    as.character(unlist(cluster5_entrez_2,use.names = F)),
                    as.character(unlist(cluster6_entrez_2,use.names = F)))
  names(gene_list2) <- c('CL1','CL2','CL3','CL4','CL5','CL6')
} else if (k_used == 7){
  cluster1_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster1_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster2_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster2_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster3_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster3_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster4_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster4_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster5_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster5_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster6_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster6_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster7_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster7_2$Group1), 'ENTREZID', 'SYMBOL')
  gene_list2 <- list(as.character(unlist(cluster1_entrez_2,use.names = F)),
                    as.character(unlist(cluster2_entrez_2,use.names = F)),
                    as.character(unlist(cluster3_entrez_2,use.names = F)),
                    as.character(unlist(cluster4_entrez_2,use.names = F)),
                    as.character(unlist(cluster5_entrez_2,use.names = F)),
                    as.character(unlist(cluster6_entrez_2,use.names = F)),
                    as.character(unlist(cluster7_entrez_2,use.names = F)))
  names(gene_list2) <- c('CL1','CL2','CL3','CL4','CL5','CL6','CL7')
} else if (k_used == 8){
  cluster1_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster1_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster2_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster2_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster3_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster3_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster4_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster4_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster5_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster5_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster6_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster6_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster7_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster7_2$Group1), 'ENTREZID', 'SYMBOL')
  cluster8_entrez_2 <- mapIds(org.Hs.eg.db, as.character(cluster8_2$Group1), 'ENTREZID', 'SYMBOL')
  gene_list2 <- list(as.character(unlist(cluster1_entrez_2,use.names = F)),
                    as.character(unlist(cluster2_entrez_2,use.names = F)),
                    as.character(unlist(cluster3_entrez_2,use.names = F)),
                    as.character(unlist(cluster4_entrez_2,use.names = F)),
                    as.character(unlist(cluster5_entrez_2,use.names = F)),
                    as.character(unlist(cluster6_entrez_2,use.names = F)),
                    as.character(unlist(cluster7_entrez_2,use.names = F)),
                    as.character(unlist(cluster8_entrez_2,use.names = F)))
  names(gene_list2) <- c('CL1','CL2','CL3','CL4','CL5','CL6','CL7','CL8')
}
x3=compareCluster(gene_list2, fun='enrichGO',OrgDb='org.Hs.eg.db',ont='BP')
x4=compareCluster(gene_list2, fun='enrichPathway')
x5=compareCluster(gene_list2,fun='enrichKEGG')
dotplot(x3, showCategory=5, includeAll=FALSE,title = 'GO (Biological process) enriched in Clusters (metabolic genes)')
if (save_images){
  setwd(wd2)
  ggsave('GO_BP_Cluster1_to_Clusterk_METABOLIC_MARCH2019.png',width = 10, height = 8, units = "in")
  setwd(wd)
}
dotplot(x4, showCategory=5, includeAll=FALSE,title = 'Reactome pathways enriched in Clusters (metabolic genes)')
if (save_images){
  setwd(wd2)
  ggsave('Reactome_Cluster1_to_Clusterk_METABOLIC_MARCH2019.png',width = 14, height = 8, units = "in")
  setwd(wd)
}
dotplot(x5, showCategory=5, includeAll=FALSE,title = 'KEGG pathways enriched in Clusters (metabolic genes)')
if (save_images){
  setwd(wd2)
  ggsave('KEGG_Cluster1_to_Clusterk_METABOLIC_APRIL2019.png',width = 14, height = 8, units = "in")
  setwd(wd)
}

#--------------------------------------------------------#
# 12. It seems that cluster 1 has a slightly lower survival change than the rest. Find the (near)-optimal signature from 
# those genes:
q_vector <- seq(0.2,0.8,by=0.1)

# Find the best single gene predictors in cluster of interest:
cand_sign_genes <- c()
for (i in 1:dim(cluster1_2)[1]){
  test2 <- data.frame(cluster1_2$Group1[i])
  colnames(test2) <- 'Group1'
  p_val_cand <- findBestSurvival(MMS_rna_data,test2,q_vector,cor_data2,print_plot=FALSE)
  print(p_val_cand)
  if (p_val_cand < 0.05){
    print('EUREKA')
    cand_sign_genes <- c(cand_sign_genes,as.character(cluster1_2$Group1[i]))
  } else {
    cand_sign_genes <- cand_sign_genes
  }
}
cand_sign_genes <- data.frame(cand_sign_genes)

# Now test which genes are positively and negatively correlated with worse surival:
harmful_genes <- c()
beneficial_genes <- c()
for (i in 1:dim(cand_sign_genes)[1]){
  test3 <- data.frame(cand_sign_genes[i,])
  effect <- findBestSurvival(MMS_rna_data,test3,q_vector,cor_data2,effect_type = TRUE)
  if (effect == 'Harmful!'){
    harmful_genes <- c(harmful_genes,as.character(cand_sign_genes$cand_sign_genes[i]))
  } else if (effect == 'Beneficial!'){
    beneficial_genes <- c(beneficial_genes,as.character(cand_sign_genes$cand_sign_genes[i]))
  }
}
harmful_genes <- data.frame(harmful_genes)
beneficial_genes <- data.frame(beneficial_genes)
colnames(harmful_genes) <- 'Group1'
colnames(beneficial_genes) <- 'Group1'

# Save the genes of interest....
setwd(wd2)
write.table(harmful_genes, file = "genes_correlated_with_worse_survival.txt", sep = "\t",
            row.names = FALSE,quote=FALSE)
write.table(beneficial_genes, file = "genes_correlated_with_better_survival.txt", sep = "\t",
            row.names = FALSE,quote=FALSE)
setwd(wd)
# Now check what pairs/triplicates are the best (or just the whole list of harmful genes:)
q_vector <- seq(0.2,0.8,by=0.1)
p_ultime <- findBestSurvival(MMS_rna_data,harmful_genes,q_vector,cor_data2,print_plot=TRUE,effect_type = FALSE)
p_ultime
p_ultime2 <- findBestSurvival(MMS_rna_data,beneficial_genes,q_vector,cor_data2,print_plot=TRUE,effect_type = FALSE)
p_ultime2

if (save_images){
  setwd(wd2)
  ggsave("COXPH_harmful_signature_APRIL2019.png", plot=print(p_ultime), width = 8, height = 6, dpi = 1000)
  setwd(wd)
}

