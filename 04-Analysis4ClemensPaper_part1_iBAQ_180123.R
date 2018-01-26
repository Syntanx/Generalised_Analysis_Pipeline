# Code Backend #
# Part One - Import, clean-up, exploratory-stats graphs, test-stats, dysregulation lists and volcano plots

## Boilerplate Section ##
# Load Libraries #
require(ggplot2)
require(reshape2)
require(scales)
require(gplots)
# require(ggbiplot)
require(ggfortify)
require(stringr)
# require(sva)
require(dplyr)
require(tidyr)

checkMake <- function(mainDir,subDir){
  path <- file.path(mainDir,subDir)
  if (file.exists(path)){
    # setwd(file.path(mainDir, subDir))
    cat("Folder Exists\n")
  } else {
    dir.create(path)
    # setwd(file.path(mainDir, subDir))
    cat(paste(path," was created.\n",sep = ""))
  }
  return(paste0(path,"/"))
}

int_type <- "iBAQ"   # one of "iBAQ", "LFQ" and "Intensity"
plot_width <- 30
plot_height <- 20
plot_units <- "cm"
plot_resolution <- 300 #dpi
####################
setwd("Y:/2017/TB/Clemens/Paper_with_Alex/")
# Set Paths #
path = getwd()

path_proteinGroups  = paste(path, "Reruns_proteinGroups.txt",sep = "/")
# path_summaryFile = paste(path, "summary.txt",sep="/")

# output = "G:/2015/TB/Alexander/Proteomics/RifTimePoints1/RifPaper/Output/" # Where ggplots etc. ought to be saved.
output = checkMake(mainDir = path,subDir = "Alex/Output_iBAQ_6_2018_01_23/")
# Scale = 5 # Set generally for ggplots

# Load Tables #
proteinGroups  = read.csv(path_proteinGroups,header = T,sep = "\t")
# summaryFile = read.csv(path_summaryFile,header = T,sep = "\t")

# pat = "[CR][\._]T[1-3][\._][1-3]"
# str_match

# Clean out Contaminants and Reverse Hits
contamIndex = grep(pattern = "\\+",as.factor(proteinGroups$Potential.contaminant))
revIndex = grep(pattern = "\\+",as.character(proteinGroups$Reverse))
bySite = grep(pattern = "\\+",proteinGroups$Only.identified.by.site)
cleanIndex = unique(c(contamIndex,revIndex,bySite))

proteinGroups_clean = proteinGroups[-cleanIndex,]

# Intensity and iBAQ Columns
intensity_index_clean = grep(pattern = "Intensity.",x = colnames(proteinGroups_clean))
iBAQ_index_clean = grep(pattern = "iBAQ.",x = colnames(proteinGroups_clean))
LFQ_index_clean = grep(pattern = "LFQ",x = colnames(proteinGroups_clean))

# Transform Tables #
proteinGroups_clean_long = melt(data = proteinGroups_clean, id.vars = c("id","Protein.IDs"), measure.vars = c(intensity_index_clean,LFQ_index_clean,iBAQ_index_clean),value.name = "Intensity",variable.name = "Sample")
proteinGroups_clean_long = transform(proteinGroups_clean_long, LFQ = grepl(pattern = "LFQ", x = proteinGroups_clean_long$Sample), iBAQ = grepl(pattern = "iBAQ", x = proteinGroups_clean_long$Sample))

#```


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
medNorm = function(data_in,col_index){
  data = data_in[,col_index]
  colmed = apply(data,2,median)
  colmult = median(colmed)/colmed
  data_norm = sweep(x = data,MARGIN = 2,STATS = colmult,FUN = "*")
  return(data_norm)
  cat(apply(data_norm,2,median))
}
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
# PCA - Intensity
# values = NA_from_0_proteinGroups_clean_LFQ[,-1]
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
# values = NA_from_0_proteinGroups_clean_LFQ[,-1]
values = NA_from_0_proteinGroups_clean_norm_iBAQ[,-1]
groups = str_match(string = colnames(values),pattern = "([CR])[123]")[,2]
replicate = str_match(string = colnames(values),pattern = "[CR]([123])")[,2]
batch = str_match(string = colnames(values),pattern = "C[CR][1-3]")
batch[grep(x = colnames(values),pattern = "C[CR][1-3]")] <- 1
batch[grep(x = colnames(values),pattern = "\\.[CR][1-3]")] <- 2
partition = str_match(string = colnames(values),pattern = "_([CWcytoCD]{2,4})")[,2]
df     = transform(as.data.frame(t(values)), ExperimentalGroup = groups, Replicate = replicate, CellPartition = partition, Batch = batch)
# df[,"ExperimentalGroup"] <- as.character(groups)
# row.names(df) <- str_replace(string = colnames(values),pattern = "LFQ.*\\.([CR])",replacement = "\\1")
pca = prcomp(x = t(na.omit(values)),scale. = T)

mult = 1
tiff(filename = paste0(output_pca,"PCA_2-4",".tiff"),width = plot_width*mult,height = plot_height*mult,units = plot_units,res = plot_resolution)
autoplot(pca, x = 2, y = 4, data = df, colour = "ExperimentalGroup", shape = "CellPartition",label = F, label.size = 3, label.hjust = -0.2,main = "PCA using Median Normalised iBAQ Values",frame=F)
dev.off()

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
# 
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
# 
# # CW samples
# # CW_clean_intensity_index <- grep(x = colnames(proteinGroups_clean_norm_intensity),pattern = "CW")
# # CW_sub_intensity <- proteinGroups_clean_norm_intensity[,c(grep(pattern = "Protein.IDs",x = colnames(proteinGroups_clean_norm_intensity)), CW_clean_intensity_index)]
# 
# data = proteinGroups_clean_norm
# # CW_batch = str_replace(string = str_replace_na(string = str_match(string = colnames(data[,-1]),pattern = "\\.(C)[CR][1-3]")[,2],replacement = 1),pattern = "C",replacement = 2)
# batch = batch
# # timePoints = str_match(string = colnames(data),pattern = ".*T([1-3]).*")[,2]
# partition = str_match(string = colnames(data[,-1]),pattern = "_([CWDyto]{2,4})")[,2]
# treatmentsGroup = sub(pattern = "R",replacement = "2",
#   x = sub(pattern = "C",replacement = "1",
#     x = str_match(string = colnames(data[,-1]),pattern = "([RC])[1-3].*")[,2])
#   )
# # The following is being performed twice - removeZeroVarianceProteins() already does this inside of doCombat.
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
pVal.CutOff = 0.05
groups <- c("Control","Treated")

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
  data$threshold = as.factor(
    (abs(data$logFC) > FCthres & data$pVals < pVal.CutOff)*1 +
    (data$pVals < pVal.CutOff)*1
  )
  
  ##Construct the plot object
  g = ggplot(data=data, aes(x=logFC, y=-log10(pVals))) + # May need to define that my_palette. colour=my_palette
    geom_point(alpha=0.5, size=1, aes(color=threshold)) +
    theme_bw() +
    theme(panel.grid = element_line(linetype = "dashed")) +
    theme(legend.position = "none",title = element_text(face = "bold"),text = element_text(size = 14)) + #,size = 12
    # xlim(c(-max(abs(data$logFC)), max(abs(data$logFC)))) + ylim(c(0, 4)) +
    xlim(c(-6.5,6.5)) + ylim(c(0,7.1)) +
    # xlim(c(min(data$logFC), max(data$logFC))) + ylim(c(0, max(-(log10(data$pVals))))) +
    xlab(expression(log[2]*"(Fold Change)")) + ylab(expression(-log[10]*"(p-Value)")) + ggtitle(paste("Volcano Plot for",time, "Fraction")) +
    geom_vline(xintercept = c(FCthres,-FCthres),linetype = "dashed") + 
    geom_hline(yintercept = -log10(pVal.CutOff),linetype = "dashed") +
    #       scale_color_brewer(type = "qual",palette = 1)
    scale_colour_manual(values = c("black","red","blue"))
  g
}

output_vplots <- checkMake(mainDir = output,subDir = "Volcano_Plots")    # Images are all the wrong size/dimensions still
plot_scale = 1
plot_optional_text <- "_20_small_14_multicolour_bw_axes_tiff"   # set as "" for default naming
# unit_string <- "cm"
# width <- 30*Scale #px
# height <- 20*Scale #px
# res <- 300 #ppi

datain = data_CW
data_CW_clean = cutrows(datain) # Remove any rows with zeroes
data_CW_added = add_ttest_FC(data_CW_clean) #Calculate and add columns with fold change (FC) and pVal from t.test()
CW_FCthres = getFCthres(data_CW_added)
volPlot <- volcanoPlot(data_CW_added,time="Cell Wall",CW_FCthres)
volPlot
ggsave(filename = paste0(output_vplots,"vPlot_","CW",plot_optional_text,".tiff"),plot = volPlot, width = plot_width*plot_scale,height = plot_height*plot_scale,units = plot_units,dpi = plot_resolution)
paste("CW pVal:",length(which(data_CW_added$pVals<pVal.CutOff)))
paste("CW pVal & FCthres:",length(which(data_CW_added$pVals<pVal.CutOff & abs(data_CW_added$logFC)>CW_FCthres)))
data_CW_pValC = filter(data_CW_added,pVals<pVal.CutOff)
data_CW_FC = filter(data_CW_pValC,abs(logFC)>CW_FCthres)
data_CW_pValC_comp <- transform(data_CW_pValC, Compartment = rep("CW",times=nrow(data_CW_pValC)))

datain = data_Cyto
data_Cyto_clean = cutrows(datain) # Remove any rows with zeroes
data_Cyto_added = add_ttest_FC(data_Cyto_clean) #Calculate and add columns with fold change (FC) and pVal from t.test()
Cyto_FCthres = getFCthres(data_Cyto_added)
volPlot <- volcanoPlot(data_Cyto_added,time="Cytosolic",Cyto_FCthres)
volPlot
ggsave(filename = paste0(output_vplots,"vPlot_","Cyto",plot_optional_text,".tiff"),plot = volPlot, width = plot_width*plot_scale,height = plot_height*plot_scale,units = plot_units,dpi = plot_resolution)
paste("Cyto pVal:",length(which(data_Cyto_added$pVals<pVal.CutOff)))
paste("Cyto pVal & FCthres:",length(which(data_Cyto_added$pVals<pVal.CutOff & abs(data_Cyto_added$logFC)>Cyto_FCthres)))
data_Cyto_pValC = filter(data_Cyto_added,pVals<pVal.CutOff)
data_Cyto_FC = filter(data_Cyto_pValC,abs(logFC)>Cyto_FCthres)
data_Cyto_pValC_comp <- transform(data_Cyto_pValC, Compartment = rep("Cyto",times=nrow(data_Cyto_pValC)))

datain = data_CD
data_CD_clean = cutrows(datain) # Remove any rows with zeroes
data_CD_added = add_ttest_FC(data_CD_clean) #Calculate and add columns with fold change (FC) and pVal from t.test()
CD_FCthres = getFCthres(data_CD_added)
volPlot <- volcanoPlot(data_CD_added,time="Cell Debris",CD_FCthres)
volPlot
ggsave(filename = paste0(output_vplots,"vPlot_","CD",plot_optional_text,".tiff"),plot = volPlot, width = plot_width*plot_scale,height = plot_height*plot_scale,units = plot_units,dpi = plot_resolution)
paste("CD pVal:",length(which(data_CD_added$pVals<pVal.CutOff)))
paste("CD pVal & FCthres:",length(which(data_CD_added$pVals<pVal.CutOff & abs(data_CD_added$logFC)>CD_FCthres)))
data_CD_pValC = filter(data_CD_added,pVals<pVal.CutOff)
data_CD_FC = filter(data_CD_pValC,abs(logFC)>CD_FCthres)
data_CD_pValC_comp <- transform(data_CD_pValC, Compartment = rep("CD",times=nrow(data_CD_pValC)))

maxLength = max(nrow(data_CW_pValC),nrow(data_Cyto_pValC),nrow(data_CD_pValC))

data_all_pValC <- full_join(full_join(data_CD_pValC_comp, data_Cyto_pValC_comp), data_CW_pValC_comp)

assignStars <- function(numVec,three=0.001,two=0.01,one=0.05,ns=1,stars=c("NS","*","**","***")){
  stars[as.numeric(numVec<three) + as.numeric(numVec<two) + as.numeric(numVec<one) + as.numeric(numVec<ns)]
}

data_all_pValC <- transform(data_all_pValC, Stars = assignStars(pVals))
#### Present and Absent ####
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

#### Dysreg Out ####
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
# ```

