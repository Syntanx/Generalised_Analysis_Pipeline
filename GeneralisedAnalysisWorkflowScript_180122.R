# Code Backend #
# Part One - Import, clean-up, exploratory-stats graphs, test-stats, dysregulation lists and volcano plots

#### Boilerplate Section ####

# Load Libraries #
require(ggplot2)
require(reshape2)
require(scales)
require(gplots)
# require(ggbiplot)
# require(ggfortify)
require(stringr)
require(sva)
require(hash)

#### functions ####
sleep_time <- 1

checkMake <- function(subDir="",mainDir="."){
  if ( (length(str_detect(string = subDir,pattern = "[[:upper:]]{1}:[\\\\/]+"))>0) & (mainDir==".") ) {mainDir <- subDir; subDir <- ""}
  
  path <- file.path(mainDir,subDir)
    
  if (dir.exists(path)){
    # setwd(file.path(mainDir, subDir))
    cat("Folder Exists\n")
  } else {
    dir.create(path)
    # setwd(file.path(mainDir, subDir))
    cat(paste(path," was created.\n",sep = ""))
  }
  return(paste0(path,"/"))
}

testExclusivity <- function(grouping_logical){
  if (any(apply(X = as.data.frame(grouping_logical),sum)==2)){
    cat(paste0("Grouping ",group," appears to not be exclusive.\nSome samples appear to match multiple groups"))
    stop()
  } else {return(TRUE)}
}
condenseGroup <- function(group, data=data){
  grouping_logical <- sapply(X = get(group), function(x){grepl(pattern = x, x = get(latestData)[["Sample"]])})
  testExclusivity(grouping_logical)
  for (x in get(group)){testData[[group]][grouping_logical[[x]]] <- x}
}
groupsFun <- function(groupsVector,data=get(latestData)){
  tempData <- data
  for(group in groupsVector){condenseGroup(group = group, data = data)}
}

medNorm = function(data_in,col_index){
  data = data_in[,col_index]
  colmed = apply(data,2,median)
  colmult = median(colmed)/colmed
  data_norm = sweep(x = data,MARGIN = 2,STATS = colmult,FUN = "*")
  return(data_norm)
  cat(apply(data_norm,2,median))
}

extractFromCols <- function(input, cols = colnames(latestData), rename = NULL){
  
  colIndexes <- sapply(X = input,FUN = function(x){grep(pattern = x,x = cols)})
  colNames <- cols[colIndexes]
  return(list(Indexes <- colIndexes, Names <- colNames))
}

#### Boilerplate Continued ####
# Files #
ProjectDirectory <- checkMake("Y:/2017/TB/Clemens/Paper_with_Alex/Alex/RevampedScript/")  <- choose.dir(default = getwd()) # Where to find the data and where to write the output
OutputDirectory <-   # You can specify, or else will make and/or write to the directory "Output" inside the ProjectDirectory
DataFile <- "Reruns_proteinGroups.txt"
DataFileParameters <- hash(delimiter = ",", skiplines = 0,header = T) # Delimiter: Comma (","), Tab ("\t"), or other. (";"). NULL => default (comma).
dataFileExtraParam <- hash(contamColSearch <- "contaminant", revColSearch <- "reverse", bySiteColSearch <- "by\\.Site" )

# Major Options #
QuantType <- "iBAQ"   # Can be anything to signify which column names to use for Quantitation. c("LFQ","intensity","iBAQ")
MedianNormalise <- TRUE  # Apply MedNorm funciton to median normalise quantitation values. c(TRUE,FALSE)
IncludePresenceAbsence <- TRUE  # Should analytes be additionally analysed by presence/absence or only by t.test (or anova... # How general can I make this?)
StatisticalTest <- "ttest"  # c("ttest","anova")
StatisticalParameters <- hash(statistical_test <- StatisticalTest, pVal_Cut_off <- 0.05, multiple_testing_correction <- MultipleTestingCorrection, min_fold_change <- MinFoldChange)  # Any parameters to be passed to the statistical test.
# MinFoldChange <- c(0,"fc")  # Require a minimum fold-change for consideration in statistical test. c(x,"sd") means x*the standard deviation of fold change for that sample. c(x,"fc") is the default and means x is the fold-change
# MultipleTestingCorrection <- NULL  # What method of multiple testing correction should be employed? c(NULL,"BH","Bonferroni") BH - Benjamini-Hochberg
FastaExtractionPattern <- hash(accession_ID <- NULL, gene_name <- NULL, description <- NULL)  # If doing proteomics and using a Uniprot Fasta, give the regular-expression to be used to identify you protein information (and gene names).

# Experimental Section #
Group1 <- extractFromCols(c("C","R"))
Group2 <- extractFromCols(c("CW","Cyto","CD"))

CompareBetween <- Group1
CompareAmong <- Group2

# Figure Dimension defaults
Scale = 1
unit_string <- "cm"  # c("cm","mm","in","px")
width <- 30*Scale
height <- 20*Scale
res <- 300 #dpi


#### Script ####

# Set-Up #
originalDirectory <- getwd()  # Can set the working directory back to this at the end if desired. Probably not though.
setwd(checkMake(ProjectDirectory))
outputDir <- ifelse(test = is.null(OutputDirectory),no = OutputDirectory,yes = checkMake("Output"))
input <- read.table(file = DataFile,sep = DataFileParameters$delimiter,skip = DataFileParameters$skiplines,header = DataFileParameters$header)

int_type <- QuantType

# Assign contam, rev and bySite cols
tempFun2 <- function(vars){
  sapply(X = vars, FUN = tempFun1)
}
tempFun1 <- function(var){
  cols <- colnames(input)
  temp <- grep(pattern = dataFileExtraParam[var],x = cols)
  assign(paste0(var,"Col"),hash(Index <- temp, Name <- cols[temp]))
}
tempFun2(c("contam","rev","bySite"))



# path_proteinGroups = DataFile
# output = checkMake(mainDir = outputDir,subDir = NULL)

# Load Tables #
# proteinGroups  = read.csv(path_proteinGroups,header = T,sep = "\t")
# summaryFile = read.csv(path_summaryFile,header = T,sep = "\t")

# pat = "[CR][\._]T[1-3][\._][1-3]"
# str_match

# Clean out Contaminants and Reverse Hits
cleanData <- function(data, contam=T,rev=T,bySite=T, MQ=T){
  if(MQ){return(cleanData_MQ(data, contam=contam,rev=rev,bySite=bySite))}
}
cleanData_MQ <- function(data,contam=T,rev=T,bySite=T){
  contamIndex = grep(pattern = "\\+",as.factor(data[contamCol]))
  revIndex = grep(pattern = "\\+",as.character(data[revCol]))
  bySite = grep(pattern = "\\+",data[bySiteCol])
  cleanIndex = unique(c(contamIndex,revIndex,bySite))
  
  return(proteinGroups[-cleanIndex,])
}

data_clean <- cleanData(input); latestData <- "data_clean"

# Intensity and iBAQ Columns
intensity_index_clean = grep(pattern = "Intensity.",x = colnames(get(latestData)))
iBAQ_index_clean = grep(pattern = "iBAQ.",x = colnames(get(latestData)))
LFQ_index_clean = grep(pattern = "LFQ",x = colnames(get(latestData)))
allIntIndex <- union(union(intensity_index_clean,iBAQ_index_clean),LFQ_index_clean)
firstIntCol <- min(allIntIndex)
excludeIndex <- setdiff(c(intensity_index_clean, iBAQ_index_clean, LFQ_index_clean),get(paste0(QuantType,"_index_clean")))

# Restrict data to chosen intensity columns
assign(x = paste(latestData,QuantType,sep = "_"),value = get(latestData)[-excludeIndex]); latestData <- paste(latestData,QuantType,sep = "_")

# Transform Tables #
data_temp <- get(latestData)
assign(paste(latestData,"long",sep = "_"), 
       value = melt(data = data_temp, id.vars = seq(1,firstIntCol,by=1), measure.vars = grep(pattern = QuantType,x = colnames(data_temp)), value.name = paste0(QuantType,"_Intensity", variable.name = "Sample")))
latestData <- paste(latestData,"long",sep = "_")

# Assign Groups to observations
Groups <- c("Group1","Group2")
assign(paste0(latestData,"_grouped"),
       value = groupsFun(groupsVector = Groups,data = get(latestData)))



# Basic Data Quality
# To look at some basics of the data. To begin, how many proteins were identified/quantified in the initial analysis across all files?
# ```{r echo=FALSE,message=FALSE}
# Identified
ggplot_numProteins_clean = ggplot(proteinGroups_clean_long[grepl(pattern = "Intensity",proteinGroups_clean_long$Sample) & (proteinGroups_clean_long$Intensity>0),])
ggplot_numProteins_clean + geom_bar() + aes(x = Sample) + ggtitle("Proteins Identified") + theme(axis.text.x = element_text(angle = 90,vjust = 0.5))

# iBAQ Quantified
ggplot_numProteins_quant_clean = ggplot(proteinGroups_clean_long[grepl(pattern = "iBAQ",proteinGroups_clean_long$Sample) & (proteinGroups_clean_long$Intensity>0),])
ggplot_numProteins_quant_clean + geom_bar() + aes(x = Sample) + ggtitle("Proteins Quantified (iBAQ)") + theme(axis.text.x = element_text(angle = 90,vjust = 0.5))

# Identified Peptides from SummaryFile
# peptides = summaryFile[19:36,]
# peptides_plot = ggplot(data = peptides,aes(x=Raw.file,y=Peptide.Sequences.Identified)) + geom_bar(stat = "identity") + labs(x = "Experiment", y = "Peptide Sequences Identified", title = "Peptide Sequences Identified by Experiment") + theme(title = element_text(family = "Times",face = "bold"), axis.text.x = element_text(vjust = 1,hjust = 1,angle = 45))
# peptides_plot
# ggsave(filename = "Sequences.jpg",plot = peptides_plot,scale = 2)
# ```
# 
# Hmm, one sample looks a bit low. Looking at protein abundance/normality (in a loading and not Gaussian sense):
#   
#   ```{r echo=FALSE,message=FALSE}
# ggplot_unnorm_proteinGroups_clean_LFQ = ggplot(data = proteinGroups_clean_long[proteinGroups_clean_long$LFQ==T,])
# ggplot_unnorm_proteinGroups_clean_LFQ + geom_boxplot() + aes(x = Sample, y = Intensity) + scale_y_continuous(trans=log2_trans()) #+ coord_trans(y="log10")
 
# Intensities before normalising.
# ggplot_unnorm_proteinGroups_clean_Int = ggplot(data = proteinGroups_clean_long[proteinGroups_clean_long$LFQ==F,])
# ggplot_unnorm_proteinGroups_clean_Int + geom_boxplot() + aes(x = Sample, y = Intensity) + scale_y_continuous(limits = c(0,4e9)) #+ scale_y_continuous(trans=log2_trans()) #+ coord_trans(y="log10")
boxplot(x = log(x = proteinGroups_clean[intensity_index_clean],base = 2),main = "Logged Intensity before Norm")
boxplot(x = log(x = proteinGroups_clean[iBAQ_index_clean],base = 2),main = "Logged iBAQ values before Norm",notch=T,boxlwd=1,pars = list(boxwex = 0.7))
# boxplot(x = proteinGroups_clean[iBAQ_index_clean])
# boxplot(x = proteinGroups_clean[LFQ_index_clean], ylim = c(0,2e8),main = "LFQ before Norm")
 
# Median Normalize Intensity
# ind = grep("Intensity\\.",colnames(proteinGroups_clean))
# data_Int = proteinGroups_clean[,intensity_index_clean]
# data_LFQ = proteinGroups_clean[,LFQ_index_clean]

  proteinGroups_clean_norm = proteinGroups_clean
  proteinGroups_clean_norm[,intensity_index_clean] = medNorm(proteinGroups_clean[,intensity_index_clean])
  proteinGroups_clean_norm[,iBAQ_index_clean] = medNorm(proteinGroups_clean[,iBAQ_index_clean])
  
# Intensities after normalising.
boxplot(x = log(x = proteinGroups_clean_norm[intensity_index_clean],base = 2),main = "Logged Intensity after medNorm")
boxplot(x = log(x = proteinGroups_clean_norm[iBAQ_index_clean],base = 2),main = "Logged iBAQ values after medNorm")
# boxplot(x = log(x = proteinGroups_clean_norm[LFQ_index_clean],base = 2))
# boxplot(x = proteinGroups_clean_norm[LFQ_index_clean], ylim = c(0,2e8),main = "LFQ after medNorm")
# Can't median normalise LFQs I think - zero's do not imply low abundance.

# 
# # Remelt data...
# proteinGroups_clean_long_norm = melt(data = proteinGroups_clean_norm, id.vars = c("id","Protein.IDs"), measure.vars = c(intensity_index_clean,LFQ_index_clean),value.name = "Intensity",variable.name = "Sample")
# proteinGroups_clean_long_norm = transform(proteinGroups_clean_long_norm, LFQ = grepl(pattern = "LFQ", x = proteinGroups_clean_long_norm$Sample))
# 
# # Intensities after normalising.
# ggplot_mednorm_proteinGroups_clean_Int = ggplot(data = proteinGroups_clean_long_norm[proteinGroups_clean_long_norm$LFQ==F,])
# ggplot_mednorm_proteinGroups_clean_Int + geom_boxplot() + aes(x = Sample, y = Intensity) +scale_y_continuous(limits=c(0,1e8)) #+ scale_y_continuous(trans=log2_trans()) #+ coord_trans(y="log10")
# ```
 
# Again, the same sample looks out. Plotting a heatmap (with hierarchical clustering):
   
#   ```{r echo=FALSE,message=FALSE}
####################################################
# Working with Normalised iBAQ values
proteinGroups_clean_norm_iBAQ = proteinGroups_clean_norm[,c(grep(pattern = "Protein.IDs",colnames(proteinGroups_clean)),iBAQ_index_clean)]
cor_proteinGroups_clean_norm_iBAQ_bySample = cor(proteinGroups_clean_norm_iBAQ[,-1])
hc_proteinGroups_clean_norm_iBAQ_bySample = hclust(d = dist(x = cor_proteinGroups_clean_norm_iBAQ_bySample,method = "euclidean"),method = "ward.D2")
hc_proteinGroups_clean_norm_iBAQ_by = hclust(d = dist(x = cor_proteinGroups_clean_norm_iBAQ_bySample,method = "euclidean"),method = "ward.D2")
hc_proteinGroups_clean_norm_iBAQ_byProtein = hclust(d = dist(x = proteinGroups_clean_norm_iBAQ[,-1],method = "euclidean"),method = "ward.D2")

NA_from_0_proteinGroups_clean_norm_iBAQ = proteinGroups_clean_norm_iBAQ
NA_from_0_proteinGroups_clean_norm_iBAQ[proteinGroups_clean_norm_iBAQ==0] <- NA

my.colours = colorpanel(n = 1000,low = "blue",mid = "yellow",high = "red")
# hm_proteinGroups_clean_norm_iBAQ = heatmap.2(x = as.matrix(NA_from_0_proteinGroups_clean_norm_iBAQ[,-1]),Rowv = hc_proteinGroups_clean_norm_iBAQ_byProtein,Colv = hc_proteinGroups_clean_norm_iBAQ_bySample$order,dendrogram = "column",scale = "row",na.rm = T,col=my.colours,trace = "none")

# # Working with LFQ values
# proteinGroups_clean_norm_LFQ = proteinGroups_clean_norm[,c(grep(pattern = "Protein.IDs",colnames(proteinGroups_clean)),LFQ_index_clean)]
# cor_proteinGroups_clean_norm_LFQ_bySample = cor(proteinGroups_clean_norm_LFQ[,-1])
# hc_proteinGroups_clean_norm_LFQ_bySample = hclust(d = dist(x = cor_proteinGroups_clean_norm_LFQ_bySample,method = "euclidean"),method = "ward.D2")
# hc_proteinGroups_clean_norm_LFQ_by = hclust(d = dist(x = cor_proteinGroups_clean_norm_LFQ_bySample,method = "euclidean"),method = "ward.D2")
# hc_proteinGroups_clean_norm_LFQ_byProtein = hclust(d = dist(x = proteinGroups_clean_norm_LFQ[,-1],method = "euclidean"),method = "ward.D2")
# 
# NA_from_0_proteinGroups_clean_norm_LFQ = proteinGroups_clean_norm_LFQ
# NA_from_0_proteinGroups_clean_norm_LFQ[proteinGroups_clean_norm_LFQ==0] <- NA
# # NA_from_0_proteinGroups_clean_Int = proteinGroups_clean_norm[,c(grep(pattern = "Protein.IDs",colnames(proteinGroups_clean)),intensity_index_clean)]
# # NA_from_0_proteinGroups_clean_Int[NA_from_0_proteinGroups_clean_Int==0] <- NA
# 
# my.colours = colorpanel(n = 1000,low = "blue",mid = "yellow",high = "red")
# 
# # Working with Normalised Intensity values
# proteinGroups_clean_norm_intensity = proteinGroups_clean_norm[,c(grep(pattern = "Protein.IDs",colnames(proteinGroups_clean)),intensity_index_clean)]
# cor_proteinGroups_clean_norm_intensity_bySample = cor(proteinGroups_clean_norm_intensity[,-1])
# hc_proteinGroups_clean_norm_intensity_bySample = hclust(d = dist(x = cor_proteinGroups_clean_norm_intensity_bySample,method = "euclidean"),method = "ward.D2")
# hc_proteinGroups_clean_norm_intensity_by = hclust(d = dist(x = cor_proteinGroups_clean_norm_intensity_bySample,method = "euclidean"),method = "ward.D2")
# hc_proteinGroups_clean_norm_intensity_byProtein = hclust(d = dist(x = proteinGroups_clean_norm_intensity[,-1],method = "euclidean"),method = "ward.D2")
# 
# NA_from_0_proteinGroups_clean_norm_intensity = proteinGroups_clean_norm_intensity
# NA_from_0_proteinGroups_clean_norm_intensity[proteinGroups_clean_norm_intensity==0] <- NA
# # NA_from_0_proteinGroups_clean_Int = proteinGroups_clean_norm[,c(grep(pattern = "Protein.IDs",colnames(proteinGroups_clean)),intensity_index_clean)]
# # NA_from_0_proteinGroups_clean_Int[NA_from_0_proteinGroups_clean_Int==0] <- NA
# 
# my.colours = colorpanel(n = 1000,low = "blue",mid = "yellow",high = "red")
# proteinGroups_clean_norm_intensity = heatmap.2(x = as.matrix(NA_from_0_proteinGroups_clean_norm_intensity[,-1]),Rowv = hc_proteinGroups_clean_norm_intensity_byProtein,Colv = hc_proteinGroups_clean_norm_intensity_bySample$order,dendrogram = "column",scale = "row",na.rm = T,col=my.colours,trace = "none")
# hm_proteinGroups_clean_LFQ
# ```
# 
# We see sample handling plays a major role in the clustering ("\_1", Alexander; "\_2/3", Elise) - this suggests that experimental variation is comparable with the biological variation we desire to measure: Either experimental error is large or biological variation is small (or a combination thereof). 
# 
# We again see sample R\_T1\_1 appears to be an outlier. This is made even more clear by examining the data through Principle Components Analysis (PCA). {Note to self: see (here)[http://www.vince.vu/software/] for ideas to make PCA plots prettier.}
# 
# ```{r echo=FALSE,message=FALSE}
# # prcomp
# pca_proteinGroups_clean_LFQ= prcomp(x = t(na.omit(NA_from_0_proteinGroups_clean_LFQ[,-1])),scale. = T)
# # autoplot(pca_proteinGroups_clean_LFQ)
# ggbiplot(pca_proteinGroups_clean_LFQ,var.scale = 0,obs.scale = 1,varname.abbrev = T,labels = colnames(proteinGroups_clean_LFQ[,-1]))
# 
# # variation explained
# pca = pca_proteinGroups_clean_LFQ
# variances = data.frame(variances=pca$sdev**2, pcomp=1:length(pca$sdev)) #gives variances of PCs; **2 means ^2
# varplot = ggplot(variances, aes(pcomp, variances)) +
#   geom_bar(stat="identity", fill="gray") #+ 
# #   geom_line() 
# varplot
# 
# # prcomp
# pca_proteinGroups_clean_LFQ= prcomp(x = t(na.omit(NA_from_0_proteinGroups_clean_LFQ[,-1])),scale. = T)
# # autoplot(pca_proteinGroups_clean_LFQ)
# ggbiplot(pca_proteinGroups_clean_LFQ,var.scale = 0,obs.scale = 1,varname.abbrev = T,labels = colnames(proteinGroups_clean_LFQ[,-1]))
# 
# # variation explained
# pca = pca_proteinGroups_clean_LFQ
# variances = data.frame(variances=pca$sdev**2, pcomp=1:length(pca$sdev)) #gives variances of PCs; **2 means ^2
# varplot = ggplot(variances, aes(pcomp, variances)) +
#   geom_bar(stat="identity", fill="gray") #+ 
# #   geom_line() 
# varplot
# 
# 

output_pca <- checkMake(mainDir = output,subDir = "PCA_plots")
# # PCA - Intensity
# # values = NA_from_0_proteinGroups_clean_LFQ[,-1]
# values = NA_from_0_proteinGroups_clean_norm_intensity[,-1]
# groups = str_match(string = colnames(values),pattern = "([CR])[123]")[,2]
# replicate = str_match(string = colnames(values),pattern = "[CR]([123])")[,2]
# batch = str_match(string = colnames(values),pattern = "C[CR][1-3]")
# batch[grep(x = colnames(values),pattern = "C[CR][1-3]")] <- 1
# batch[grep(x = colnames(values),pattern = "\\.[CR][1-3]")] <- 2
# partition = str_match(string = colnames(values),pattern = "_([CWcytoCD]{2,4})")[,2]
# df     = transform(as.data.frame(t(values)), ExperimentalGroup = groups, Replicate = replicate, CellPartition = partition, Batch = batch)
# # df[,"ExperimentalGroup"] <- as.character(groups)
# # row.names(df) <- str_replace(string = colnames(values),pattern = "LFQ.*\\.([CR])",replacement = "\\1")
# pca = prcomp(x = t(na.omit(values)),scale. = T)
# autoplot(pca, data = df, colour = "ExperimentalGroup", shape = "Batch",label = F, label.size = 3, label.hjust = -0.2,main = "Before ComBat Adjustment - Intensity",frame=F)
# 
# # ggbiplot(pca,var.scale = 0,obs.scale = 1,varname.abbrev = T,labels = colnames(proteinGroups_clean_LFQ[,-1]))
# # biplot(pca_proteinGroups_clean_LFQ)
# # plot(pca_proteinGroups_clean_LFQ)
# 
# # generate & show a bar graph & elbow plot of PC contributions to variance
# # pca = pca_proteinGroups_clean_LFQ
# variances = data.frame(variances=pca$sdev**2, pcomp=1:length(pca$sdev)) #gives variances of PCs; **2 means ^2
# 
# varplot = ggplot(variances, aes(pcomp, variances)) +
#   geom_bar(stat="identity", fill="gray") #+ 
# #   geom_line() 
# varplot + ggtitle("PC Contributions to Variance Coverage")


# PCA - iBAQ
factorRename <- function(string,search = c("C","R"),replace = c("Control","Rifampicin Treated")){
  len <- length(search)
  if (len != length(replace)) {cat("Replacement options not equal in length to existing")}
  temp <- string
  for (i in 1:length(search)) {
    temp <- str_replace_all(string = temp,pattern = search[i],replacement = replace[i])
  }
  out <- temp
  return(out)
  # subTable <- data.frame(Search=search, Replace=replace)
  # out <- apply(X = subTable, MARGIN = 1, FUN = function(x){
  #   str_replace_all(string = string, pattern = x[[1]],replacement = x[[2]])[2]})
  # out
}
# values = NA_from_0_proteinGroups_clean_LFQ[,-1]
values = NA_from_0_proteinGroups_clean_norm_iBAQ[,-1]
groups = str_match(string = colnames(values),pattern = "([CR])[123]")[,2]
groups <- factorRename(string = groups)
replicate = str_match(string = colnames(values),pattern = "[CR]([123])")[,2]
batch = str_match(string = colnames(values),pattern = "C[CR][1-3]")
batch[grep(x = colnames(values),pattern = "C[CR][1-3]")] <- 1
batch[grep(x = colnames(values),pattern = "\\.[CR][1-3]")] <- 2
partition = factorRename(string = str_match(string = colnames(values),pattern = "_([CWcytoCD]{2,4})")[,2], search = c("CW","CD","Cyto"),replace = c("Cell Wall","Cell Debris","Cytosol"))
df     = transform(as.data.frame(t(values)), Experimental_Group = groups, Replicate = replicate, Cell_Partition = partition, Batch = batch)
# df[,"ExperimentalGroup"] <- as.character(groups)
# row.names(df) <- str_replace(string = colnames(values),pattern = "LFQ.*\\.([CR])",replacement = "\\1")
pca = prcomp(x = t(na.omit(values)),scale. = T)

mult = 1
# tiff(filename = paste0(output_pca,"PCA_2-4",".tiff"),width = width*mult,height = height*mult,units = unit_string,res = res)
autoplot(pca, x = 2, y = 4, data = df, colour = "Experimental_Group", shape = "Cell_Partition",label = F, label.size = 3, label.hjust = -0.2,main = "Principal Components Analyis",frame=F)
# Sys.sleep(sleep_time)
# dev.off()
# Sys.sleep(sleep_time)

# ggbiplot(pca,var.scale = 0,obs.scale = 1,varname.abbrev = T,labels = colnames(proteinGroups_clean_LFQ[,-1]))
# biplot(pca_proteinGroups_clean_LFQ)
# plot(pca_proteinGroups_clean_LFQ)

# generate & show a bar graph & elbow plot of PC contributions to variance
# pca = pca_proteinGroups_clean_LFQ
variances = data.frame(variances=pca$sdev**2, pcomp=1:length(pca$sdev)) #gives variances of PCs; **2 means ^2

varplot = ggplot(variances, aes(pcomp, variances)) +
  geom_bar(stat="identity", fill="gray") #+ 
#   geom_line() 
varplot + ggtitle("PC Contributions to Variance Coverage")

# ```
# All looks fine - proceeding with iBAQ values from hereon:


## ComBat ##       # Not Currently being used
# But then I remembered about ComBat and used that. The algorithm requires proteins have non-zero variance in each batch wich resulted in the dropping of [620] proteins - but these probably weren't making their way into the t-tests anyway. But ultimately we need only decide if this way is better in toto.

# ```{r echo=FALSE,message=FALSE}
# Need to correct for batch effects in at least CW samples - CD and Cyto samples look pretty good.

# # ComBat Batch effects
# is.0 = function(arr){
#   return(as.numeric(arr)==0)
# }
# removeZeroVarianceProteins <- function(data,batch){
#   # This finds those proteins with only 0's (hence 0 variance) in across either all my or all Elise's samples
#   ind = apply(X = data[,-1][,batch=="2"],MARGIN = 1,FUN = function(x){all(is.0(x))}) |
#     apply(X = data[,-1][,batch=="1"],MARGIN = 1,FUN = function(x){all(is.0(x))})
#   dat = data[!ind,]
#   return(dat)
# }

# CW samples
# CW_clean_intensity_index <- grep(x = colnames(proteinGroups_clean_norm_intensity),pattern = "CW")
# CW_sub_intensity <- proteinGroups_clean_norm_intensity[,c(grep(pattern = "Protein.IDs",x = colnames(proteinGroups_clean_norm_intensity)), CW_clean_intensity_index)]

# data = proteinGroups_clean_norm
# # CW_batch = str_replace(string = str_replace_na(string = str_match(string = colnames(data[,-1]),pattern = "\\.(C)[CR][1-3]")[,2],replacement = 1),pattern = "C",replacement = 2)
# batch = batch
# # timePoints = str_match(string = colnames(data),pattern = ".*T([1-3]).*")[,2]
# partition = str_match(string = colnames(data[,-1]),pattern = "_([CWDyto]{2,4})")[,2]
# treatmentsGroup = sub(pattern = "R",replacement = "2",
#   x = sub(pattern = "C",replacement = "1",
#     x = str_match(string = colnames(data[,-1]),pattern = "([RC])[1-3].*")[,2])
#   )
# # # The following is being performed twice - removeZeroVarianceProteins() already does this inside of doCombat.
# ind = apply(X = data[,-1][,batch=="2"],MARGIN = 1,FUN = function(x){all(is.0(x))}) |
#   apply(X = data[,-1][,batch=="1"],MARGIN = 1,FUN = function(x){all(is.0(x))})
# proteinGroups_clean_norm_iBAQ_cut = proteinGroups_clean_norm_iBAQ[!ind,]
# 
# doCombat <- function(data, batch){
#   dat = removeZeroVarianceProteins(data = data,batch = batch)
#   
#   was0 = as.matrix(dat[,-1])==0
#   
#   model = model.matrix(~partition+treatmentsGroup)
#   data_ComBat = ComBat(dat = dat[,-1],batch = as.factor(batch),mod = model)
#   data_ComBat[was0] <- 0
#   
#   repCombatValue = NA # setting a value to zero effectively removes it from the analysis as I exclude proteins with any zero's within a time point.
#   relativeDifs = abs( data_ComBat - dat[,-1] )/dat[,-1] # How strongly was each value shifted?
#   data_ComBat[relativeDifs>0.5] <- repCombatValue # Set to 0 (exclude) values that shifted by more than half the original value.
#   return(data_ComBat)
# }

# # CW_sub_iBAQ_ComBat = doCombat(data = CW_sub_iBAQ,batch = CW_batch)
# # proteinGroups_clean_norm_iBAQ_cut[,CW_clean_iBAQ_index] <- CW_sub_iBAQ_ComBat
# proteinGroups_clean_norm_iBAQ_cut_ComValues = doCombat(data = proteinGroups_clean_norm_iBAQ_cut,batch = batch)
# proteinGroups_clean_norm_iBAQ_cut_ComBat <- proteinGroups_clean_norm_iBAQ_cut
# proteinGroups_clean_norm_iBAQ_cut_ComBat[,-1] <- proteinGroups_clean_norm_iBAQ_cut_ComValues
# 
# # Repeat PCA
# NA_from_0_proteinGroups_clean_norm_iBAQ_cut_ComBat = proteinGroups_clean_norm_iBAQ_cut_ComBat
# NA_from_0_proteinGroups_clean_norm_iBAQ_cut_ComBat[NA_from_0_proteinGroups_clean_norm_iBAQ_cut_ComBat==0] <- NA
# 
# values = NA_from_0_proteinGroups_clean_norm_iBAQ_cut_ComBat[,-1]
# groups = str_match(string = colnames(values),pattern = "([CR])[123]")[,2]
# replicate = str_match(string = colnames(values),pattern = "[CR]([123])")[,2]
# batch = str_match(string = colnames(values),pattern = "C[CR][1-3]")
# batch[grep(x = colnames(values),pattern = "C[CR][1-3]")] <- 1
# batch[grep(x = colnames(values),pattern = "\\.[CR][1-3]")] <- 2
# partition = str_match(string = colnames(values),pattern = "_([CWcytoCD]{2,4})")[,2]
# df     = transform(as.data.frame(t(values)), ExperimentalGroup = groups, Replicate = replicate, CellPartition = partition, Batch = batch)
# # df[,"ExperimentalGroup"] <- as.character(groups)
# # row.names(df) <- str_replace(string = colnames(values),pattern = "LFQ.*\\.([CR])",replacement = "\\1")
# pca = prcomp(x = t(na.omit(values)),scale. = T)
# autoplot(pca, data = df, colour = "ExperimentalGroup", shape = "Batch",label = F, label.size = 3, label.hjust = -0.2,main = "iBAQ After ComBat Adjustment",frame=F)
# 




# Write to file for Perseus
# dat[,-1] <- proteinGroups_clean_LFQ_ComBat[,-1]

# ComBat corrected:
# export2Perseus = proteinGroups[(proteinGroups$Protein.IDs %in% proteinGroups_clean_norm_iBAQ_cut_ComBat$Protein.IDs),]
# export2Perseus[,iBAQ_index_clean] <- proteinGroups_clean_norm_iBAQ_cut_ComBat[,-1] 

# Unadjusted values:
# export2Perseus = proteinGroups[(proteinGroups$Protein.IDs %in% proteinGroups_clean_norm_iBAQ$Protein.IDs),]
export2Perseus = proteinGroups[(proteinGroups$Protein.IDs %in% proteinGroups_clean_norm_iBAQ$Protein.IDs),]
export2Perseus[,iBAQ_index_clean] <- proteinGroups_clean_norm_iBAQ[,-1]
# write.table(x = export2Perseus,file = "C:/Users/Alexander/Desktop/RifAnalysisTemp/All2/ComBat_Adjusted_Full.txt",quote = F,sep = "\t",row.names = F)



# export2Perseus = proteinGroups_clean_norm # Use iBAQ instead of LFQ



# Now show the heatmaps
# Original one (all2uding outlier)
# hm_proteinGroups_clean_LFQ = heatmap.2(x = as.matrix(NA_from_0_proteinGroups_clean_LFQ[,-1]),Rowv = hc_proteinGroups_clean_LFQ_byProtein,Colv = hc_proteinGroups_clean_LFQ_bySample$order,dendrogram = "column",scale = "row",na.rm = T,col=my.colours,trace = "none",main = "Before ComBat")

# # ComBat Adjusted
# cor_proteinGroups_clean_LFQ_ComBat_bySample = cor(proteinGroups_clean_LFQ_ComBat[,-1])
# hc_proteinGroups_clean_LFQ_ComBat_bySample = hclust(d = dist(x = cor_proteinGroups_clean_LFQ_ComBat_bySample,method = "euclidean"),method = "ward.D2")
# # hc_proteinGroups_clean_LFQ_byProtein = hclust(d = dist(x = cor_proteinGroups_clean_LFQ_bySample,method = "euclidean"),method = "ward.D2")
# hc_proteinGroups_clean_LFQ_ComBat_byProtein = hclust(d = dist(x = proteinGroups_clean_LFQ_ComBat[,-1],method = "euclidean"),method = "ward.D2")
# 
# NA_from_0_proteinGroups_clean_LFQ_ComBat = proteinGroups_clean_LFQ
# NA_from_0_proteinGroups_clean_LFQ_ComBat[proteinGroups_clean_LFQ==0]<-NA
# 
# # hm_proteinGroups_clean_LFQ_ComBat = heatmap.2(x = as.matrix(NA_from_0_proteinGroups_clean_LFQ_ComBat[,-1]),Rowv = hc_proteinGroups_clean_LFQ_ComBat_byProtein,Colv = hc_proteinGroups_clean_LFQ_ComBat_bySample$order,dendrogram = "column",scale = "row",na.rm = T,col=my.colours,trace = "none",main = "After ComBat")
# 
# # PCA
# values = NA_from_0_proteinGroups_clean_LFQ[,-1]
# groups = str_match(string = colnames(values),pattern = "([CR])")[,2]
# handlers=str_match(string = colnames(values),pattern = "\\.([1-3])")[,2]
# time   = as.numeric(str_match(string = colnames(values),pattern = "T([1-3])")[,2])
# df     = transform(t(values),Group = groups, Handler = handlers, TimePoint = time)
# row.names(df) <- str_replace(string = colnames(values),pattern = "LFQ.*\\.([CR])",replacement = "\\1")
# pca = prcomp(x = t(na.omit(values)),scale. = T)
# autoplot(pca, data = df, colour = "Group", shape = "Handler",label = T, label.size = 3, label.hjust = -0.2,main = "Before ComBat Adjustment",xlim = c(-40,42),frame=T)
# # ggbiplot(pca,var.scale = 0,obs.scale = 1,varname.abbrev = T,labels = colnames(proteinGroups_clean_LFQ[,-1]))
# # biplot(pca_proteinGroups_clean_LFQ)
# # plot(pca_proteinGroups_clean_LFQ)
# 
# #generate & show a bar graph & elbow plot of PC contributions to variance
# pca = pca_proteinGroups_clean_LFQ
# variances = data.frame(variances=pca$sdev**2, pcomp=1:length(pca$sdev)) #gives variances of PCs; **2 means ^2
# varplot = ggplot(variances, aes(pcomp, variances)) +
#   geom_bar(stat="identity", fill="gray") #+ 
# #   geom_line() 
# varplot + ggtitle("PC Contributions to Variance Coverage")
# 
# 
# # pca_proteinGroups_clean_LFQ = princomp(x = t(na.omit(NA_from_0_proteinGroups_clean_LFQ[,-1])),cor = T,na.action = "na.omit")
# values = proteinGroups_clean_LFQ_ComBat
# groups = str_match(string = colnames(values),pattern = "([CR])")[,2]
# handlers=str_match(string = colnames(values),pattern = "\\.([1-3])")[,2]
# time   = as.numeric(str_match(string = colnames(values),pattern = "T([1-3])")[,2])
# df     = transform(t(values),Group = as.factor(groups), Handler = as.factor(handlers), TimePoint = as.factor(time))
# row.names(df) <- str_replace(string = colnames(values),pattern = "LFQ.*\\.([CR])",replacement = "\\1")
# pca = prcomp(x = t(na.omit(values)),scale. = T)
# df2 = as.data.frame(cbind(PC1=as.numeric(pca$x[,"PC1"]),PC2=as.numeric(pca$x[,"PC2"]),Group=groups,Handler=handlers,TimePoint=time,name=rownames(pca$x)))
# df2$name <- str_replace(string = colnames(values),pattern = "LFQ.*\\.([CR])",replacement = "\\1")
# # View(df2)
# acPCA = ggplot(pca$x, aes(PC1,PC2, color=df2$Group)) + geom_point() + scale_x_continuous() + scale_y_continuous() + ggtitle("Principle Component Analysis after Batch Correction") + scale_colour_discrete(name="Treatment Group") + theme(title = element_text(family = "times",face = "bold")) #+ geom_text(mapping = aes(label=df2$name), show.legend=F) 
# acPCA
# # ggsave(filename = "After ComBat PCA.jpg",plot = acPCA,scale = 2)
# 
# plot(pca$x,xlim = c(-41,41),ylim=c(-45,62))
# # identified <- identify(pca$x[,"PC1"],pca$x[,"PC2"],labels=df2$name,pos=T)

# df2$pos <- NA
# df2[identified$ind,]$pos <- identified$pos
# ggplot(pca$x, aes(as.numeric(PC1),as.numeric(PC2))) + geom_point() +
#   geom_point(data=subset(df2,!is.na(pos)),aes(color=Group)) +
#   geom_text(data=subset(df2,pos == 1),aes(label=name),vjust=1.1) + 
#   geom_text(data=subset(df2,pos == 2),aes(label=name),hjust=1.1) +
#   geom_text(data=subset(df2,pos == 3),aes(label=name),vjust=-0.5) + 
#   geom_text(data=subset(df2,pos == 4),aes(label=name),hjust=-0.1)

# plot(pca)
# autoplot(pca, data = df, colour = "Group", shape = "TimePoint",label = T, label.size = 3, label.hjust = -0.2,main = "Principle Component Analysis after Batch Correction",xlim = c(-40,42),frame=F,frame.type = "norm" )
# ggbiplot(pca,var.scale = 0,obs.scale = 1,varname.abbrev = T,labels = colnames(proteinGroups_clean_LFQ[,-1]))
# biplot(pca_proteinGroups_clean_LFQ)
# plot(pca_proteinGroups_clean_LFQ)

# #generate & show a bar graph & elbow plot of PC contributions to variance
# pca = pca_proteinGroups_clean_LFQ
# variances = data.frame(variances=pca$sdev**2, pcomp=1:length(pca$sdev)) #gives variances of PCs; **2 means ^2
# varplot = ggplot(variances, aes(pcomp, variances)) +
#   geom_bar(stat="identity", fill="gray") #+ 
# #   geom_line() 
# varplot + ggtitle("PC Contributions to Variance Coverage")

# ```

# Put all that t-test coding to work and make a custom volcano plot.

# ```{r,echo=FALSE}
library(dplyr)
# Split into time points.
input = export2Perseus  #Choose input data frame
pVal.CutOff = pVal_Cut_off
groups <- c("Control","Rifampicin Treated")

Text_Cols.names = c(
  "Majority.protein.IDs",
  "Protein.IDs",
  "Fasta.headers"
)

# Common_Cols = grep(pattern = ,colnames(input)) #Columns common to each (split) dataset.
quant_type = "iBAQ"  #could be "iBAQ" or "intensity"
CW_Cols = grep(pattern = paste0(quant_type,".*CW"),colnames(input),value=T)
Cyto_Cols = grep(pattern = paste0(quant_type,".*Cyto"),colnames(input),value=T)
CD_Cols = grep(pattern = paste0(quant_type,".*CD"),colnames(input),value=T)

data_CW = input[,c(CW_Cols,Text_Cols.names)]
data_Cyto = input[,c(Cyto_Cols,Text_Cols.names)]
data_CD = input[,c(CD_Cols,Text_Cols.names)]

# Function to see if all values to be compared are above zero
checkrows = function(data,req=3,groups=groups){
  
  cols_1 = grep(pattern = paste0("[C][1-3]"),colnames(data))
  cols_2 = grep(pattern = paste0("[R][1-3]"),colnames(data))
  
  keep_1 = apply(X = data[,cols_1],MARGIN = 1,FUN = function(row){(sum(row>0,na.rm = T) >= req) & !any(is.na(row))})
  keep_2 = apply(X = data[,cols_2],MARGIN = 1,FUN = function(row){(sum(row>0,na.rm = T) >= req) & !any(is.na(row))})
  keep  <- (keep_1 & keep_2)
  # keep = apply(X = data[,1:8],MARGIN = 1,FUN = function(row){(sum(row>0,na.rm = T) >= req) & !any(is.na(row))})
  return(keep)
}
# Function to remove rows not satisfying the above criterion
cutrows = function(data,req=3,groups=groups){
  outdata = data[checkrows(data[,1:8],req,groups),]
  return(outdata)
}

add_ttest_FC = function(data,groups=c("C","R")){
  cols_1 = grep(pattern = paste0("[",as.character(groups[1]),"][1-4]"),colnames(data))
  cols_2 = grep(pattern = paste0("[",as.character(groups[2]),"][1-4]"),colnames(data))
  
  pVals = apply(X = data,MARGIN = 1,FUN = function(row){
    row_n <- as.numeric(row)
    row_n[row_n==0] <- NA
    t.test(row_n[cols_1],row_n[cols_2],var.equal=T)$p.value})
  logFC = apply(X = data,MARGIN = 1,FUN = function(row){
    row_n <- as.numeric(row)
    calc = t.test(row_n[cols_1],row_n[cols_2])
    log2(calc$estimate[2]/calc$estimate[1])})
  return(cbind(data,pVals,logFC))
}

getFCthres = function(data,col="logFC"){
  FCthres = 2*abs(sd(data[,paste(col)]))   # Need to look at this. I feel doubling the deviation at the log level means with a slightly larger deviation I restrict "strong" to things with much larger fold changes cause you've gone from logFC of ~1 to 1.4 implying fold change of 2 to 2.6. Or 1.6 = 3 fold.
}

volcanoPlot = function(data,time,FCthres){
  #   FCthres = 2*abs(sd(data$logFC))
  data$threshold = as.factor(abs(data$logFC) > FCthres & data$pVals < pVal.CutOff)
  
  ##Construct the plot object
  g = ggplot(data=data, aes(x=logFC, y=-log10(pVals))) + # May need to define that my_palette. colour=my_palette
    geom_point(alpha=0.5, size=1, aes(color=threshold)) +
    theme(legend.position = "none",title = element_text(face = "bold",size = 12)) +
    # xlim(c(-max(abs(data$logFC)), max(abs(data$logFC)))) + ylim(c(0, 4)) +
    xlim(c(-2, 2)) + ylim(c(0, 4)) +
    xlab(expression(log[2]*"(Fold Change)")) + ylab(expression(-log[10]*"(p-Value)")) + ggtitle(paste("Volcano Plot for Time Point",time)) +
    geom_vline(xintercept = c(FCthres,-FCthres),linetype = "dashed") + 
    geom_hline(yintercept = -log10(pVal.CutOff),linetype = "dashed") +
    #       scale_color_brewer(type = "qual",palette = 1)
    scale_colour_manual(values = c("black","blue"))
  g
}

output_vplots <- output    # Images are all the wrong size/dimensions still
# Uncomment to override defaults from here on
# Scale = 1
# unit_string <- "cm"
# width <- 30*Scale #px
# height <- 20*Scale #px
# res <- 300 #ppi

datain = data_CW
data_CW_clean = cutrows(datain) # Remove any rows with zeroes
data_CW_added = add_ttest_FC(data_CW_clean) #Calculate and add columns with fold change (FC) and pVal from t.test()
CW_FCthres = getFCthres(data_CW_added)
# tiff(filename = paste0(output_vplots,"vPlot_","CW",".tiff"),width = width,height = height,units = unit_string,res = res)
# volcanoPlot(data_CW_added,time="CW",CW_FCthres)
# Sys.sleep(sleep_time)
# dev.off()
# Sys.sleep(sleep_time)
# ggsave(filename = "CW",device = "tiff",path = output_vplots,width = width,height = height,dpi = res,units = )
paste("CW pVal:",length(which(data_CW_added$pVals<pVal.CutOff)))
paste("CW pVal & FCthres:",length(which(data_CW_added$pVals<pVal.CutOff & abs(data_CW_added$logFC)>CW_FCthres)))
data_CW_pValC = filter(data_CW_added,pVals<pVal.CutOff)
data_CW_FC = filter(data_CW_pValC,abs(logFC)>CW_FCthres)

datain = data_Cyto
data_Cyto_clean = cutrows(datain) # Remove any rows with zeroes
data_Cyto_added = add_ttest_FC(data_Cyto_clean) #Calculate and add columns with fold change (FC) and pVal from t.test()
Cyto_FCthres = getFCthres(data_Cyto_added)
# tiff(filename = paste0(output_vplots,"vPlot_","Cyto",".tiff"),width = width,height = height,units = unit_string,res = 300)
# volcanoPlot(data_Cyto_added,time="Cyto",Cyto_FCthres)
# Sys.sleep(sleep_time)
# dev.off()
# Sys.sleep(sleep_time)
paste("Cyto pVal:",length(which(data_Cyto_added$pVals<pVal.CutOff)))
paste("Cyto pVal & FCthres:",length(which(data_Cyto_added$pVals<pVal.CutOff & abs(data_Cyto_added$logFC)>Cyto_FCthres)))
data_Cyto_pValC = filter(data_Cyto_added,pVals<pVal.CutOff)
data_Cyto_FC = filter(data_Cyto_pValC,abs(logFC)>Cyto_FCthres)

datain = data_CD
mult = 1
data_CD_clean = cutrows(datain) # Remove any rows with zeroes
data_CD_added = add_ttest_FC(data_CD_clean) #Calculate and add columns with fold change (FC) and pVal from t.test()
CD_FCthres = getFCthres(data_CD_added)
# tiff(filename = paste0(output_vplots,"vPlot_","CD",".tiff"),width = width*mult,height = height*mult,units = unit_string,res = 300)
# volcanoPlot(data_CD_added,time="CD",CD_FCthres)
# Sys.sleep(sleep_time)
# dev.off()
# Sys.sleep(sleep_time)
paste("CD pVal:",length(which(data_CD_added$pVals<pVal.CutOff)))
paste("CD pVal & FCthres:",length(which(data_CD_added$pVals<pVal.CutOff & abs(data_CD_added$logFC)>CD_FCthres)))
data_CD_pValC = filter(data_CD_added,pVals<pVal.CutOff)
data_CD_FC = filter(data_CD_pValC,abs(logFC)>CD_FCthres)

maxLength = max(nrow(data_CW_pValC),nrow(data_Cyto_pValC),nrow(data_CD_pValC))

### Present and Absent ###
doPresAbs <- F

if (doPresAbs){

getPresAbs = function(data,cols1,cols2) {
  data2 = data
  data2$Up <- !apply(data[,cols2]==0,1,any) & apply(data[,cols1]==0,1,all)
  data2$Down <- !apply(data[,cols1]==0,1,any) & apply(data[,cols2]==0,1,all)
  presAbs = data2[(data2$Up|data2$Down),]
  return(presAbs)
}

data = data_CW  
CW_pa_up = filter(
  .data = getPresAbs(
    data = data, 
    cols1 = grep(pattern = paste0("[C][1-3]"),colnames(data)),
    cols2 = grep(pattern = paste0("[R][1-3]"),colnames(data))),
  Down==F)$Protein.IDs
CW_pa_down = filter(
  .data = getPresAbs(
    data = data, 
    cols1 = grep(pattern = paste0("[C][1-3]"),colnames(data)),
    cols2 = grep(pattern = paste0("[R][1-3]"),colnames(data))),
  Down==T)$Protein.IDs
data = data_Cyto
Cyto_pa_up = filter(
  .data = getPresAbs(
    data = data, 
    cols1 = grep(pattern = paste0("[C][1-3]"),colnames(data)),
    cols2 = grep(pattern = paste0("[R][1-3]"),colnames(data))),
  Down==F)$Protein.IDs
Cyto_pa_down = filter(
  .data = getPresAbs(
    data = data, 
    cols1 = grep(pattern = paste0("[C][1-3]"),colnames(data)),
    cols2 = grep(pattern = paste0("[R][1-3]"),colnames(data))),
  Down==T)$Protein.IDs
data = data_CD
CD_pa_up = filter(
  .data = getPresAbs(
    data = data, 
    cols1 = grep(pattern = paste0("[C][1-3]"),colnames(data)),
    cols2 = grep(pattern = paste0("[R][1-3]"),colnames(data))),
  Down==F)$Protein.IDs
CD_pa_down = filter(
  .data = getPresAbs(
    data = data, 
    cols1 = grep(pattern = paste0("[C][1-3]"),colnames(data)),
    cols2 = grep(pattern = paste0("[R][1-3]"),colnames(data))),
  Down==T)$Protein.IDs


# # Filtering done in Excel :( Sorry future Alex!!
# CW_pa_down = read.csv("./PresAbs/CW_down.txt",header = F,as.is = T)
# Cyto_pa_down = read.csv("./PresAbs/Cyto_down.txt",header = F,as.is = T)
# CD_pa_down = read.csv("./PresAbs/CD_down.txt",header = F,as.is = T)
# CW_pa_up = read.csv("./PresAbs/CW_up.txt",header = F,as.is = T)
# Cyto_pa_up = read.csv("./PresAbs/Cyto_up.txt",header = F,as.is = T)
# CD_pa_up = read.csv("./PresAbs/CD_up.txt",header = F,as.is = T)
# CW_pa_combined = unique(c(CW_pa_down$V1,CW_pa_up$V1))
# Cyto_pa_combined = unique(c(Cyto_pa_down$V1,Cyto_pa_up$V1))
# CD_pa_combined = unique(c(CD_pa_down$V1,CD_pa_up$V1))

data_CW_pa = rbind(
  data.frame(Majority.protein.IDs = CW_pa_up,logFC = rep(x = 2,length.out = length(CW_pa_up)),pVals = rep(x = 0,length.out = length(CW_pa_up)) ),
  data.frame(Majority.protein.IDs = CW_pa_down,logFC = rep(x = -2,length.out = length(CW_pa_down)),pVals = rep(x = 0,length.out = length(CW_pa_down)) )
)
data_Cyto_pa = rbind(
  data.frame(Majority.protein.IDs = Cyto_pa_up,logFC = rep(x = 2,length.out = length(Cyto_pa_up)),pVals = rep(x = 0,length.out = length(Cyto_pa_up)) ),
  data.frame(Majority.protein.IDs = Cyto_pa_down,logFC = rep(x = -2,length.out = length(Cyto_pa_down)),pVals = rep(x = 0,length.out = length(Cyto_pa_down)) )
)
data_CD_pa = rbind(
  data.frame(Majority.protein.IDs = CD_pa_up,logFC = rep(x = 2,length.out = length(CD_pa_up)),pVals = rep(x = 0,length.out = length(CD_pa_up)) ),
  data.frame(Majority.protein.IDs = CD_pa_down,logFC = rep(x = -2,length.out = length(CD_pa_down)),pVals = rep(x = 0,length.out = length(CD_pa_down)) )
)
}

dysregulatedProteins.Fastas = list(
  CW <- as.character(data_CW_pValC$Fasta.headers),
  Cyto <- as.character(data_Cyto_pValC$Fasta.headers),
  CD <- as.character(data_CD_pValC$Fasta.headers)
)
dysregulatedProteins.Majority.protein.IDs = list( 
  CW = as.character(data_CW_pValC$Majority.protein.IDs),
  Cyto = as.character(data_Cyto_pValC$Majority.protein.IDs),
  CD = as.character(data_CD_pValC$Majority.protein.IDs)
)

searchTerm <- "[T]ransporter"
# View(
identified <- proteinGroups_clean %>% filter(grepl(searchTerm,x = Fasta.headers)) %>% select(c(contains(paste0(int_type,".")),contains("Fasta"))) %>% select(c(contains("CW"),contains("Fasta"))) %>% mutate(Sum = rowSums(.[1:8]), Count = Sum>0) %>% summarise(sum(Count))
# 73 in Cell Wall | 75 in Rif
d1 <- mutate(.data = data_CD_added,Compartment = "CD") #%>% filter(grepl(searchTerm,x = Fasta.headers) & pVals<pVal_Cut_off)
d2 <- mutate(.data = data_CW_added,Compartment = "CW")
d3 <- mutate(.data = data_Cyto_added,Compartment = "Cyto")

d4 <- bind_rows(d1,d2,d3) %>% filter(grepl(searchTerm,x = Fasta.headers) & pVals<pVal_Cut_off)
summarise(.data = d4, n_distinct(Majority.protein.IDs)) # 26 in CellWall | 2 in Rif (T2) | 17 in Rif (T2) with iBAQ and ComBat
# 4/7 Down in Cyto | 1/2 Down in Rif T2
# 8/10 Down in CW
# 12/17 Down in CD
# )
# ```
