install_packages <- function() {
  BiocManager::install("clusterProfiler")
  BiocManager::install("pathview")
  BiocManager::install("enrichplot")
  install.packages('samr')
  BiocManager::install("impute")
}

library(readxl)
library(knitr)
library(dplyr)
library(ggforce)
library(ggplot2)
library(scales)
library(reshape2)  # for melt
library(cowplot)#Load in file
library(NanoStringNCTools)
library("RColorBrewer")
library(gplots)
library(edgeR)
library(plyr)
library(DESeq2)

group_func <- function(type = "segment",count_file_name = "Filter_3%_target_count_matrix.xlsx"){
  data <- read_excel("Annotation_File.xlsx")
  target_demoData <- read_excel(count_file_name)
  if (type == "mouse"){
    group <- data$mouse_number
    group <- mapvalues(group,from = c(1,2,3,4,5,6,7,8,9,10,11), to = c("IM100","IM101","IM102","IM103","IM104","IM66","IM67","IM69","IM71","IM98","IM99"))
    tocheck <- data$mouse_number
  }
  if (type == "region") {
    group <- data$region_number
    group <- mapvalues(group,from = c(1,2,3,4,5,6), to = c("ACAd","ACAv","AI","CPu","MOp","MOs"))
    tocheck <- data$region_number
  }
  if (type == "segment") {
    group <- data$segment_number
    group <- mapvalues(group,from = c(1,2), to = c("NeuN","pSyn"))
    tocheck <- data$segment_number
  }
  group <-factor(group)
  check1 <- data$name
  check2 <- colnames(target_demoData[,-1])
  group1 <- c()
  for (j in 1:length(check2)) {
    for (i in 1:length(check1)) {
      c1 <- check1[i]
      c2 <- check2[j]
      if (c1 == c2) {
        group1 <- c(group1,tocheck[i])
      }
    }
  }
  if (type == "mouse"){
    group1 <- mapvalues(group1,from = c(1,2,3,4,5,6,7,8,9,10,11), to = c("IM100","IM101","IM102","IM103","IM104","IM66","IM67","IM69","IM71","IM98","IM99"))
  }
  if (type == "region") {
    group1 <- mapvalues(group1,from = c(1,2,3,4,5,6), to = c("ACAd","ACAv","AI","CPu","MOp","MOs"))
  }
  if (type == "segment") {
    group1 <- mapvalues(group1,from = c(1,2), to = c("NeuN","pSyn"))
  }
  group1 <- factor(group1)
  return(group1)
}

specify <- function (type,log_transformed) {
  names <- colnames(log_transformed)
  list <- c()
  new_group <- c()
  for (i in 1: length(to_sep)) {
    thing <- to_sep[i]
    if (thing == "NeuN") {
      list <- c(list, names[i])
    }
    else {
      new_group <- c(new_group,group[i])
    }
  }
  
  log_transformed <- log_transformed[ , !names(log_transformed) %in% list]
  if (type == "mouse") {
    new_group <- mapvalues(new_group,from = c(1,2,3,4,5,6,7,8,9,10,11), to = c("IM100","IM101","IM102","IM103","IM104","IM66","IM67","IM69","IM71","IM98","IM99"))
  }
  if (type == "region") {
    new_group <- mapvalues(new_group,from = c(1,2,3,4,5,6), to = c("ACAd","ACAv","AI","CPu","MOp","MOs"))
  }
  new_group <- factor(new_group)
  colnames(log_transformed) <- new_group
  
  return(log_transformed,group)
  
}

VST_N_P <- function(group, log_transformed){
  #Make it the correct file type for VST
  d0 <- DGEList(as.matrix(log_transformed))
  #Plot it in a way that retains the multiple dimensions. Helps you get a look at the data for patterns before you run this.
  plotMDS(d0,col=as.numeric(group))
  print("Plot 1 done.")
  #Get everything into right type to put into a DESeqDataSet object.
  metadata <- DataFrame(group)
  rownames(metadata) <- colnames(log_transformed)
  vst <- data.matrix(log_transformed)
  metadata <-DataFrame(group)
  #Put into the object
  vst <- DESeqDataSetFromMatrix(round(vst), colData = metadata, design = ~group)
  print("DESeq object created.")
  #Preform VST
  y <- varianceStabilizingTransformation(vst) 
  print("VST performed.")
  #Get results in a form that can be put into limma
  y <- assay(y)
  #Limma
  mm <-model.matrix(~0 + group)
  colnames(mm) <- make.names(colnames(mm))
  fit <- lmFit(y, mm)
  #Compare the two groups - Change "pSyn" and "NeuN" if group names differ.
  #This will result in a negative fold change representing down regulation and positive, up regulation.
  contrast.matrix <- makeContrasts(grouppSyn - groupNeuN, levels = colnames(coef(fit)))
  #ebays to preform DE as needed by limma
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  #Visualize data in q volcano plot
  volcanoplot(fit2,style = "P")
  print("Plot 2 done.")
  #Find positive fold change
  fold <- data.frame(fit2$coefficients)
  write.csv(fold, "All_Genes.csv")
  print("All Genes and Fold Changes")
  names <- rownames(fold)
  fold$genes <- names
  fold <- data.frame(fold)
  fold <- fold[order(-fold$grouppSyn...groupNeuN),]
  top <- fold[1:25,]
  write.csv(top,"Top_Positive_VST.csv")
  print("Up regulated genes saved.")
  #Find Negative fold change
  fold <- data.frame(fit2$coefficients)
  names <- rownames(fold)
  fold$genes <- names
  fold <- data.frame(fold)
  fold <- fold[order(fold$grouppSyn...groupNeuN),]
  top <- fold[1:25,]
  write.csv(top,"Top_Negative_VST.csv")
  print("Down regulated genes saved.")
}
#Change line 46 if group names differ from NeuN an pSyn
#All warnings/notes are okay (There should be three messages by the end).
#Perform_VST("Filter_3%_target_count_matrix.xlsx", group)

GSEA <- function () {
  library(clusterProfiler)
  library(enrichplot)
  organism = "org.Mm.eg.db"
  #BiocManager::install(organism, character.only = TRUE)
  library(organism, character.only = TRUE)
  genes = read.csv("the_current_top_genes.csv", header=TRUE)
  gene_list <- genes$grouppSyn...groupNeuN
  gene_list <- genes$Estimate
  names(gene_list) <- genes$Gene
  gene_list = sort(gene_list, decreasing = TRUE)
  gse <- gseGO(geneList=gene_list, keyType = "SYMBOL",OrgDb = organism)
  require(DOSE)
  dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
}

group <- group_func()
Perform_VST(group)
GSEA()

perform_samSeq <- function(count_file_name){
  library(samr)
  target_demoData <- read_excel(count_file_name)
  print("File read in.")
  row <- target_demoData$TargetName
  rownames(target_demoData) <- row
  data <- target_demoData[,-1]
  rownames(data) <- row
  data <- as.matrix(data)
  mode(data) <- "integer"
  samfit <- SAMseq(data, group, resp.type = "Two class unpaired") 
  
}
#perform_samSeq("Filter_3%_target_count_matrix.xlsx")


Perform_on_specifics <- function(group,count_file_name = "Filter_3%_target_count_matrix.xlsx"){
  #Read in file
  target_demoData <- read_excel(count_file_name)
  print("File read in.")
  #log10 transform the data to get a normal distribution
  log_transformed <- log(target_demoData[,-1])
  row <- target_demoData$TargetName
  rownames(log_transformed) <- row
  print("Log transformed.")
  names <- colnames(log_transformed)
  list <- c()
  new_group <- c()
  for (i in 1: length(to_sep)) {
    thing <- to_sep[i]
    if (thing == "NeuN") {
      list <- c(list, names[i])
    }
    else {
      new_group <- c(new_group,group[i])
    }
  }
  log_transformed <- log_transformed[ , !names(log_transformed) %in% list]
  #MICE #new_group <- mapvalues(new_group,from = c(1,2,3,4,5,6,7,8,9,10,11), to = c("IM100","IM101","IM102","IM103","IM104","IM66","IM67","IM69","IM71","IM98","IM99"))
  #Region #new_group <- mapvalues(new_group,from = c(1,2,3,4,5,6), to = c("ACAd","ACAv","AI","CPu","MOp","MOs"))
  new_group <- factor(new_group)
  colnames(log_transformed) <- new_group
  #Make it the correct file type for VST
  d0 <- DGEList(as.matrix(log_transformed))
  #Plot it in a way that retains the multiple dimensions. Helps you get a look at the data for patterns before you run this.
  plotMDS(d0,col=as.numeric(new_group))
}
#group <- group_func("region")
#to_sep <- group_func()
#Perform_on_specifics(group)










