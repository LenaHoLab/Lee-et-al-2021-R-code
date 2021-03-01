rm(list=ls(all=TRUE));
options(stringsAsFactors=FALSE);


# Load libraries
library(limma)
library(data.table)
set.seed(666)

# Load Pathways and Data
load.pathways <- function(dir='./data/utils data/pathway/', species='human'){
  pathway.files <- list.files(dir, pattern='data', full.names=T)
  pathway.file <- grep(pattern=species, x=pathway.files, value=T)
  pathways <- readRDS(pathway.file)
  
  # filter out gene sets with gene numbers less than 15 and larger than 1000
  pathways <- pathways[which(lengths(pathways) >= 15)]
  pathways <- pathways[which(lengths(pathways) <= 1000)]
  
  return(pathways)
}

# Load RNAseq of interest
Folder="C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/VI_sORFs/Identifying_Mito_VISEPs/G-MAD/"
data=paste0(Folder,"/Data/5.SKM_120862.csv")
pathway_folder=paste0(Folder,"pathway/")
pathways <- load.pathways(dir=pathway_folder, species="human")

# Remove duplicated gene IDs
arrayData=read.csv(data,header = T)
arrayData=arrayData[!duplicated(arrayData$EntrezID),]
arrayData=arrayData[!is.na(arrayData$EntrezID),]
rownames(arrayData)=arrayData$EntrezID
arrayData=arrayData[,-1]
gene_names=arrayData$IDENTIFIER
arrayData=arrayData[,-1]

# Load gene list to run (Need identifier)

geness=gene_names

# Load genes of interest (Here Mitocarta)
ToRun=read.csv("C:/Users/gmsbsjk/Desktop/Baptiste/Bioinfo/Fasta_csv_files/Human.MitoCarta2.0.csv")
ToRun=unique(ToRun$IDENTIFIER)
ToRun=ToRun[ToRun%in%gene_names]

geness=geness[!geness%in%ToRun]
  
# Prepare Final Matrix
pathway_name=readRDS(paste0(Folder,"pathway/human_pathway_name.RDS"))
df=data.frame(row.names = pathway_name$path_id,Path_name=pathway_name$path_name)
df=subset(df,rownames(df)%in%names(pathways))


# Run Camera loop
for (i in 1:length(ToRun)) {
  gene_of_interest=ToRun[i]
  index <- ids2indices(pathways, rownames(arrayData), remove.empty=F)
  all.gene.num <- nrow(arrayData)
  
  target.gene.i=which(gene_names==gene_of_interest)[1]
  target.gene <- rownames(arrayData)[target.gene.i]
  message(i, '/',length(ToRun))
  gene.i <- as.numeric(arrayData[target.gene.i, ])
  
  design <- model.matrix( ~ gene.i)
  results <- camera(arrayData, index, design, inter.gene.cor=0.01, allow.neg.cor=T, use.ranks=F, sort=F)
  results$PValue <- results$PValue / 2
  
  results$path_id <- rownames(results)
  results$logp <- round(-log10(results$PValue), digits=3)
  results[which(results$Direction == 'Down'),'logp'] <- -results[which(results$Direction == 'Down'),'logp']
  df=cbind(df,results[rownames(df),"logp"])
  colnames(df)[i+1]=gene_of_interest
}



write.csv(df,paste0(Folder,"Results/Matrix_Camera/Random_V2_Matrix.csv"))







