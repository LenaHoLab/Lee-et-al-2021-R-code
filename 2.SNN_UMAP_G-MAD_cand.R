rm(list=ls(all=TRUE));
options(stringsAsFactors=FALSE);

library(dplyr)
library(Seurat)
library(pheatmap)
library(ggplot2)

# Load data
folder="C:/Users/e0545037/Desktop/Baptiste/Bioinfo/VI_sORFs/Identifying_Mito_VISEPs/G-MAD/Results/Matrix_Camera/SKM/"
genes_to_use=c("Mitocarta","Random_V2","Candidates")
MitoVI=read.csv("C:/Users/e0545037/Desktop/Baptiste/Bioinfo/VI_sORFs/Identifying_Mito_VISEPs/Data/mitoDE_2+.csv",header = T)
MitoSEPs=read.csv("C:/Users/e0545037/Desktop/Baptiste/Bioinfo/Fasta_csv_files/MitoCarta.sorfs.csv",header = T)

# Prepare data
D=list()
for (G in genes_to_use) {
  D[[G]]=read.csv(paste(folder,G,"_Matrix.csv",sep=""),header = T,row.names = 1)
  D[[G]]=D[[G]][!duplicated(D[[G]][,1]),]
  rownames(D[[G]])=D[[G]][,1]
  D[[G]]=D[[G]][,-1]
}

Dataset=c()
for(d in 1:length(genes_to_use)){
  Dataset=c(Dataset,rep(genes_to_use[d],ncol(D[[d]])))
}

df=cbind(D[[1]],D[[2]],D[[3]])

df[which(is.na(df),arr.ind = T)]=0


# Start analysis
M=CreateSeuratObject(counts=df,project = "G-MAD_Matrix")
M$Category=Dataset

M=FindVariableFeatures(M,selection.method = "vst",nfeatures = 5000)
VariableFeaturePlot(M)
M=ScaleData(M,features = rownames(M))
M=Seurat::RunPCA(M,features=VariableFeatures(object = M))
DimPlot(M, reduction = "pca",label = T,pt.size = 1)
M=FindNeighbors(M, dims = 1:50)
M=FindClusters(M,resolution = 0.5)

colorRampPalette(c("blue","darkgoldenrod1"))(8)

M=Seurat::RunUMAP(M,dims=1:50)
DimPlot(M,reduction = "umap",pt.size = 1,label = F)+scale_color_manual(values = c(
  "#B68453", 
  "#FFB90F",    
  "#241ADC", 
  "#4834BA",  
  "#916975",
  "#DA9E31",
  "#6D4F98",
  "#0000FF"))


Idents(M)=M$Category
DimPlot(M,reduction = "umap",pt.size = 1,label = T)
Idents(M)=M$seurat_clusters


# calculation of percentage
clusters=as.numeric(M$seurat_clusters)
names(clusters)=colnames(df)

Percentages=data.frame(Cluster=0,nb=0,nb_mito=0,nb_VI=0,nb_cand=0,pct_mito_Rand=0)
for (i in unique(clusters)) {
  clust.i=clusters[clusters==i]
  nb_mito=table(names(clust.i)%in%colnames(D[[1]]))[2]
  nb_rand=table(names(clust.i)%in%colnames(D[[2]]))[2]
  nb_cand=table(names(clust.i)%in%colnames(D[[3]]))[2]
  nb_VI=table(MitoVI$IDENTIFIER%in%names(clust.i))[2]
  pct_mito=100*nb_mito/(nb_mito+nb_rand)
  nb_VI=table(MitoVI$IDENTIFIER%in%names(clust.i))[2]
  Percentages=rbind(Percentages,c(i-1,nb_mito+nb_rand,nb_mito,nb_VI,nb_cand,pct_mito))
}
Percentages=Percentages[-1,]

Percentages$pct_mito_Rand[Percentages$Cluster==(clusters["C15orf48"]-1)]

Candidates=read.csv("C:/Users/gmsbsjk/Desktop/e0545037/Bioinfo/VI_sORFs/Results/Selected_CCDS_lfc_0.5.V5.csv",header = T)
Candidates=unique(Candidates$IDENTIFIER)
df=data.frame(IDENTIFIER=Candidates,Score=0)

for (i in 1:length(Candidates)) {
  if (Candidates[i]%in%colnames(D[[3]])) {
    j=Candidates[i]
    clust.i=clusters[j]-1
    df[i,2]=Percentages$pct_mito_Rand[Percentages$Cluster==clust.i]
  } else {df[i,2]=NA}
}

write.csv(df,"C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/VI_sORFs/Identifying_Mito_VISEPs/G-MAD/Results/Cand_score.csv")













