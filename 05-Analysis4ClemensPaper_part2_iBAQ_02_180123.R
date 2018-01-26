# ---
#   title: "RifAnalysis part 2"
# author: "Alexander Giddey"
# date: "December 14, 2015"
# output: html_document
# ---
#   
#   Generate profile plots for all proteins in a list (eg. all proteins we talk about in the paper)
# ```{r echo=FALSE}
require(tidyr)
library(gridExtra)
library(grid)
library(stringr)
library(ggplot2)
library(dplyr)
library(hash)

int_type = "iBAQ" # Put "Intensity" or "LFQ" or "iBAQ" (each with capitals)

proteome = read.table(file = "Y:/2017/Protein_databases_2017/Uniprot_db/uniprot-Msmeg_Ref_Proteome_6600_2017_01_19.tsv",header = T,sep = "\t",fill = F,quote = "\"")

checkMake <- function(mainDir,subDir){
  path <- file.path(mainDir,subDir)
  if (file.exists(path)){
    # setwd(file.path(mainDir, subDir))
    cat("Folder Exists")
  } else {
    dir.create(path)
    # setwd(file.path(mainDir, subDir))
    cat(paste(path," was created.",sep = ""))
  }
  return(paste0(path,"/"))
}

# proteinsFile = "G:/2015/TB/Alexander/Proteomics/RifTimePoints1/RifPaper/DysregulatedProteins.txt"
# proteinList = read.table(file = proteinsFile,header = F,fill = T,colClasses = "character",na.strings = "") # Reading missing values wrong. Doesn't matter for now though.
proteinList = dysregulatedProteins.Majority.protein.IDs

proteinsUniq = as.character(unique(c(proteinList[[1]],proteinList[[2]],proteinList[[3]])))
proteinsUniq = proteinsUniq[!is.na(proteinsUniq)]
searchData = export2Perseus
# out = character(0)
# for (i in 1:length(proteinsUniq)) {
#   out[i] = grep(proteinsUniq[i],searchData$Majority.protein.IDs,value=F)
# }
# View(out)
coll = lapply(X = proteinsUniq, FUN = function(x){grep(pattern = x, x = searchData$Majority.protein.IDs,value = F)}) # find rows in data with proteins in our list.
ordProt = searchData[as.numeric(coll),] #take rows in the order they were grep'ped above - i.e. same order as list that called them.
colIndex = c(
  grep(paste0(int_type,"\\."),colnames(ordProt)),
  grep("Maj.*protein",colnames(ordProt)),
  grep("Fasta",colnames(ordProt))
)
smallProt = ordProt[,colIndex]
# colnames(smallProt)[1:(grep("\\.IDs",colnames(smallProt))-1)] <- str_match(colnames(smallProt)[1:18],pattern = "LFQ.*([CR].*)\\.[1-3]")[,2]
# colnames(smallProt)
int_cols = grep(int_type,colnames(smallProt))
int_groups <- str_match(colnames(smallProt[int_cols]),"([CR])[1-3]_([CWDyto]{2,4})")
# int_groups = paste0( rep(c("C[1-3]_","R[1-3]_"),each=3), rep(c("CW","Cyto","CD"),2) )
treatment_group <- int_groups[,2]
compartment_group <- int_groups[,3]
treatment_group <- str_replace_all(string = treatment_group, pattern = "C", replacement = "Control")
treatment_group <- str_replace_all(string = treatment_group, pattern = "R", replacement = "Rifampicin")
compartment_group <- str_replace_all(string = compartment_group, pattern = "CW", replacement = "Cell Wall")
compartment_group <- str_replace_all(string = compartment_group, pattern = "CD", replacement = "Cellular Debris")
compartment_group <- str_replace_all(string = compartment_group, pattern = "Cyto", replacement = "Cytosolic Fraction")
grouping_group <- paste(int_groups[,2],int_groups[,3],sep = "_")          # Essential these are in 
finding_group <- paste(int_groups[,2],int_groups[,3],sep = "[1-3]_")      # the same order
# smallProt <- rbind(smallProt,grouping_group)


int_groups_cols = as.data.frame(sapply(unique(finding_group),function(group){grep(group,colnames(smallProt[,int_cols]))}))
colnames(int_groups_cols) <- unique(grouping_group)

prePlot = function(data,i){
  temp = data[i,int_cols]
  # num_blocks = length(int_groups)
  blocks = apply(int_groups_cols,2,function(col){as.numeric(temp[col])})# list(as.numeric(temp[1:3]),as.numeric(temp[4:6]),as.numeric(temp[7:9]),as.numeric(temp[10:12]),as.numeric(temp[13:15]),as.numeric(temp[16:18]))
  means = colMeans(blocks)
  SDs = as.numeric(apply(blocks,2,sd))
  emin= means-SDs
  emax= means+SDs
  bnames = unique(grouping_group)
  compartment = str_match(string = bnames,pattern = "_([CWDyto]{2,4})")[,2]
  compartment <- str_replace_all(string = compartment, pattern = "CW", replacement = "Cell Wall")
  compartment <- str_replace_all(string = compartment, pattern = "CD", replacement = "Cellular Debris")
  compartment <- str_replace_all(string = compartment, pattern = "Cyto", replacement = "Cytosolic Fraction")
  # compartment = compartment_group
  treatment = str_match(string = bnames,pattern = "([CR])_")[,2]
  treatment <- str_replace_all(string = treatment, pattern = "C", replacement = "Control")
  treatment <- str_replace_all(string = treatment, pattern = "R", replacement = "Rifampicin")
  # treatment = treatment_group
  # bnames = c("C.CW","C.Cyto","C.CD","R.CW","R.Cyto","R.CD")
  # Line below is still hard-coded. :(
  temp2 = data.frame(Group = bnames,Intensity=means,emin=emin,emax=emax,SD=SDs,Compartment=compartment,Treatment=treatment)#Int_One=intens1,Int_Two=intens2,Int_Three=intens3,
  return(temp2)
}

prePlot_new = function(data,i){
  intensity = as.numeric(t(data[i,int_cols]))
  group = grouping_group
  compartment = compartment_group
  treatment = treatment_group
  temp = data.frame(Group=group,Intensity=intensity,Compartment=compartment,Treatment=treatment)
  return(temp)
}

# prePlot_new_Alt = function(data,i,int_string="Raw") {
#   int_string = toupper(int_string)
#   if (int_string == "RAW") {
#     cols = 110:127
#     else if (int_string == "LFQ") {
#       cols = 147:164
#     }
#     else {cat("Please choose LFQ or RAW for Intensity.")}
#   }
#   Intensity = as.numeric(t(data[i,cols]))
#   Group = rep(x = c("C.CW","C.Cyto","C.CD","R.CW","R.Cyto","R.CD"),each=3)
#   Time = rep(c(1,2,3),each=3)
#   Treatment = rep(c("C","R"),each=9)
#   temp = data.frame(Group=Group,Intensity=Intensity,Time=Time,Treatment=Treatment)
#   return(temp)
# }

prePlot_Alt = function(data,i,int_string="Raw"){
  int_string = toupper(int_string)
  if (int_string == "RAW"){
    cols = 110:127
  } else if (int_string == "LFQ"){
      cols = 147:164
  } else {
    cat("Please choose LFQ or RAW for Intensity.")
  }
  # cols = 1:18
  # cols = grep(pattern = "LFQ",x = colnames(data)) # Change "Intensity" to "LFQ"
  # cols = grep(pattern = "Intensity\\.",x = colnames(data)) # Change "Intensity" to "LFQ"
  temp = data[i,cols]
  # temp = temp[,-1] #remove
  blocks = list(as.numeric(temp[1:3]),as.numeric(temp[4:6]),as.numeric(temp[7:9]),as.numeric(temp[10:12]),as.numeric(temp[13:15]),as.numeric(temp[16:18]))
  means = as.numeric(lapply(blocks,mean))
  means1= means[1:3]
  means2= means[4:6]
  ratios= log2(means2/means1)
  SDs = as.numeric(lapply(blocks,sd))
  emin= means-SDs
  emax= means+SDs
  bnames = c("C.CW","C.Cyto","C.CD","R.CW","R.Cyto","R.CD")
  temp2 = data.frame(Group = bnames,Intensity=means,emin=emin,emax=emax,SD=SDs,Time=rep(c(1,2,3),times=2),Treatment=rep(c("Control","Rifampicin"),each=3))
  return(temp2)
}

capitalise = function(name){
  return(paste0(toupper(substr(name, 1, 1)), substr(name, 2, nchar(name))))
}

# Trying to make a plot of expression ratios for special MarR
# ratioPlot_Alt2 = function(data,i){
#   # cols = 1:18
#   # cols = grep(pattern = "LFQ",x = colnames(data)) # Change "Intensity" to "LFQ"
#   cols = grep(pattern = "Intensity\\.",x = colnames(data)) # Change "Intensity" to "LFQ"
#   temp = data[i,cols]
#   # temp = temp[,-1] #remove
#   blocks = list(as.numeric(temp[1:3]),as.numeric(temp[4:6]),as.numeric(temp[7:9]),as.numeric(temp[10:12]),as.numeric(temp[13:15]),as.numeric(temp[16:18]))
#   means = as.numeric(lapply(blocks,mean))
#   means1= means[1:3]
#   means2= means[4:6]
#   ratios= log2(means2/means1) 
#   bnames= c("CW","Cyto","CD")
#   temp2 = data.frame(Group = bnames, Ratios=ratios, Time=c(1,2,3),Treatment="log2(R/C)")
#   
#   qplot(x = Time, y = Ratios, data = temp2, geom = c("point","path"))
#   
#   B = ggplot(data = temp2,aes(x = Time, y=Intensity,fill=Treatment)) + geom_smooth(stat = "identity") + ggtitle(paste(geneName)) + scale_fill_manual(values=c("#F8766D")) + theme(text = element_text(size = 25), title = element_text(face="bold",size=40), legend.title = element_text(size=25), legend.key.size = unit(x = 2,units = "cm")) + xlab("Time Points")
#   return(B)
# }


#   A = ggplot() + geom_bar(aes(x = bnames, y = means))
#   tempLong = gather(data = temp,key = experiment,value = intensity,1:18)[]
#   A = ggplot(data = tempLong,aes(x=experiment,y=intensity))
# }

patt="MYCS2[ ]+(.*) OS=.*GN=(.*) PE=.*"

# Following section works perfectly, but takes a long time depending on how many proteins you are plotting...
####################################

mce_and_pkng <- c("A0QQK3","A0QQV1","A0R4N1","A0R4N5","A0R4N6","A0R4N8","A0R6G4")

orig4 <- c("A0R6G4","A0R4N8","A0R4N6","A0QQK3")
new4 <- c("A0R624","A0R043","A0QPX4","A0R1I5") #"A0QPM7" Gluconate Transporter

proteins_in_Fig5 <- c(orig4,new4)


output_ePlots = checkMake(output,"Expression_Plots/")
# output_ABC <- checkMake(output_ePlots,"ABC_Transporters")

data = smallProt %>% filter(Majority.protein.IDs %in% proteins_in_Fig5)
# data = export2Perseus[,colIndex] %>% filter(Majority.protein.IDs=="A0QU51")

###########################################################
compartment_hash <- hash("CW"="Cell Wall","Cyto"="Cytosolic Fraction", "CD"="Cellular Debris")  #,key=c("CW","Cyto","CD"), values=c("Cell Wall","Cytosolic Fraction","Cellular Debris"))
#Remember to check which prePlot function is being used (use prePlot_Alt for proteins not making it through selction criteria)
tot = nrow(data)
plot_optional_text <- "03_straightLine_stars"
tempScale = 1
isPDF <- F
for (i in 1:tot) {
  protID = data$Majority.protein.IDs[i]
#   pr.row = grep(protID,proteome$Entry)
#   gene = str_split(string = proteome$Gene.names[pr.row],pattern = " ")[[1]][1]
#   protName = proteome$Protein.names[pr.row]
  protFasta = data$Fasta.headers[i]
  mat = str_match(string = protFasta,pattern = patt)
  protName = mat[,2]
  geneName = mat[,3]
  shortFasta = paste(
    paste("Description:", mat[,2], sep=" "),
    #paste("Gene Name:", mat[,3], sep=" "),
    sep="\n"
  )
  
  # B = ratioPlot_Alt2(data,i) # Use only when plotting ratios from raw intensity values
  # A = prePlot_Alt(data,i)
  # A = prePlot(data,i)
  
  
  A = prePlot_new(data,i)
  B = prePlot(data,i)
  # A = prePlot_new_Alt(data,i,"LFQ")   # For data not passing thresholds. "LFQ" or "Raw" for LFQ or raw intensity use.
  # B = prePlot_Alt(data,i)             # For data not passing thresholds. "LFQ" or "Raw" for LFQ or raw intensity use.
  
  # comPlot = joicomPlot
  # B = ggplot(data = A,aes(x = as.factor(Time), y=Intensity, fill=Treatment)) + geom_bar(position = "dodge",stat = "identity") + geom_errorbar(aes(ymin=emin,ymax=emax), width=0.2,position=position_dodge(0.9)) + ggtitle(paste(geneName)) + scale_fill_manual(values=c("#00BFC4","#F8766D")) + theme(text = element_text(size = 35/tempScale), title = element_text(face="bold",size=60/tempScale), legend.title = element_text(size=45/max(tempScale-1,1)), legend.key.size = unit(x = 2/tempScale, units="cm")) + xlab("Time Points") #+ facet_wrap(~Time)
  # B = ggplot(data = A,aes(x = as.factor(Time), fill=Treatment)) + geom_bar(data = A, position = "dodge",stat = "identity") + geom_errorbar(aes(ymin=emin,ymax=emax), width=0.2,position=position_dodge(0.9)) + ggtitle(paste(geneName)) + scale_fill_manual(values=c("#00BFC4","#F8766D")) + xlab("Time Points") #+ facet_wrap(~Time)
  # + theme(text = element_text(size = 35/tempScale), title = element_text(face="bold",size=60/tempScale), legend.title = element_text(size=45/max(tempScale-1,1)), legend.key.size = unit(x = 2/tempScale, units="cm"))
  
  dodge_width = 0.9
  D = ggplot(data = A, aes(x = as.factor(Compartment), y = Intensity, fill=Treatment)) + 
    geom_bar(data = B, position = "dodge",stat = "identity") + 
    geom_bar(data = B, position = "dodge",stat = "identity", colour="black", size = 1, show.legend=FALSE) +
    geom_point(position = position_dodge((width = dodge_width+0.1)),aes(),size = 1,show.legend=F) + 
    theme_classic() + 
    ggtitle(paste(capitalise(geneName))) + xlab("Cellular Fraction") + ylab("iBAQ Intensity") + labs(caption = shortFasta) + 
    scale_x_discrete(expand = c(0,0.6)) + 
    theme(plot.margin = unit(c(1,1,2,1.5),"cm"), 
          axis.title.y=element_text(vjust=2.5), 
          axis.title.x=element_text(vjust=-2.5), 
          legend.key.size = unit(1, "cm"), 
          legend.title = element_text("Treatment\n"), 
          legend.key = element_rect(colour="black", size = 1), 
          panel.grid = element_blank(), text = element_text(size = 16)) + 
    scale_fill_grey(start = 0.5, end = 0.9) + 
    geom_errorbar(data = B, aes(ymin=emin,ymax=emax), width=0.2, size = 0.75, position=position_dodge(0.9))  #+ geom_crossbar(data = B, aes(y=Intensity,ymax=Intensity,ymin=Intensity, colour=Treatment), width=dodge_width, position = position_dodge(width = dodge_width),show_guide=T)#+ axis.ticks.length #theme(legend.key = element_rect(fill = Treatment))# + facet_wrap(~Time) + 
  # D
  
  E <- D
  signif_label.df.carry <- data.frame()
  signif_line.df.carry <- data.frame()
  
  allIntMax <- max(A$Intensity)
  comp_short <- c("CW","CD","Cyto")
  for (comp in comp_short) {
    pos <- grep(comp, comp_short)
    comp_long <- compartment_hash[[comp]]
    compIndex <- grep(pattern = as.character(protID),x = filter(data_all_pValC,Compartment==comp)$Majority.protein.IDs)
    if (length(compIndex==1)) {
      comp_pVal <- round(filter(data_all_pValC,Compartment==comp)[[compIndex,"pVals"]],digits = 3)
      comp_stars <- filter(data_all_pValC,Compartment==comp)[[compIndex,"Stars"]]
      compIntMax <- max(filter(A,Compartment==comp_long)$Intensity)
      
      halfwidth <- 0.1
      y_line_adjust <- allIntMax*0.02
      y_label_adjust <- y_line_adjust*1.5
      
      x <- c(-halfwidth,halfwidth) + pos
      y <- rep(compIntMax + y_line_adjust,2)
      # r <- 0.15
      # t <- seq(0, 180, by = 1) * pi / 180
      # x <- r * cos(t) + pos
      # y <- r*allIntMax*0.25 * sin(t) + compIntMax + allIntMax*0.02
      signif_line.df <- data.frame(Compartment = x, Intensity = y, Treatment=rep("Rifampicin",length(x)))
      
      # signif_line.df <- data.frame(Compartment=c(pos-0.3,pos+0.3), Treatment=c("Control","Rifampicin"), Intensity=rep((compIntMax + allIntMax*0.05),2))
      signif_label.df <- data.frame(Compartment=c(comp_long), Treatment=c("Control"), Intensity=c(compIntMax + y_label_adjust), pVal=comp_pVal, Stars=comp_stars)
      
      # signif_line.df.carry <- rbind(signif_line.df.carry,signif_line.df)
      signif_label.df.carry <- rbind(signif_label.df.carry,signif_label.df)
      E <- E + geom_line(data = signif_line.df, aes(x = Compartment, Intensity), lty=1)
    }
  }
  
  
  E <- E +
    # geom_line(data = signif_line.df.carry, aes(x = Compartment, Intensity)) +
    geom_text(data = signif_label.df.carry, aes(as.factor(Compartment), Intensity, label=Stars))
  
  E
  # print(B)
  # B = ggplot(data = A,aes(x=as.factor(Time),weight=Intensity,fill=Treatment)) + geom_bar(position = "dodge") + facet_wrap(~Treatment) + ggtitle(paste(protID,geneName,sep="/")) + geom_errorbar(aes(x=Intensity,ymin=emin,ymax=emax))
  # C = arrangeGrob(B,sub=textGrob(as.character(shortFasta), x = 0, hjust = -0.1, vjust=0.1, gp = gpar(fontface = "italic", fontsize = 18)))
  # ggsave(filename = paste(output,geneName,".jpg",sep=""),plot = C,scale = Scale)
  plot_type <- NULL
  if (isPDF){plot_type <- cairo_pdf}
  ggsave(filename = paste(output_ePlots,plot_optional_text,protID,"_",geneName,".tiff",sep=""),width = plot_width, height = plot_height, units = plot_units, dpi = plot_resolution, plot = E,scale = tempScale,limitsize = F, device=plot_type)
  print(paste(100*i/tot,"% complete"))
}
#dev.off()
##################

# ```

# Use STRINGdb to perform GO term enrichment and produce STRING graphs with up/down regulation (can I show when a protein is both up and down? Should I? Perhaps plot interactions of proteins which are shared separately?)

# ```{r,echo=FALSE}
library(STRINGdb)

# Bear in mind this order may break, might need to define string.db first before this function can be defined,
my_plot_network = function (string_ids, payload_id = NULL, string.db=string.db, required_score = NULL, add_link = TRUE, add_summary = TRUE, file=NULL) {
  "\nDescription:\n  Plots an image of the STRING network with the given proteins.\n\nInput parameters:\n  \"string_ids\"        a vector of STRING identifiers\n  \"payload_id\"        an identifier of payload data on the STRING server (see method post_payload for additional informations)\n  \"score_threshold\"   a threshold on the score that overrides the default score_threshold, that we use only for the picture\n  \"add_link\"          parameter to specify whether you want to generate and add a short link to the relative page in STRING. \n                      As default this option is active but we suggest to deactivate it in case one is generating many images (e.g. in a loop). \n                      Deactivating this option avoids to generate and store a lot of short-urls on our server.\n  \"add_summary\"       parameter to specify whether you want to add a summary text to the picture. This summary includes a p-value and the number of proteins/interactions.\n\nAuthor(s):\n   Andrea Franceschini\n"
  if (is.null(required_score)) 
    required_score = string.db$score_threshold
  img = string.db$get_png(string_ids = string_ids, network_flavor = "confidence", file = NULL, payload_id = payload_id, required_score = required_score)
  if (!is.null(img)) {
    # jpeg(filename = file)
    plot(1:(dim(img)[2]), type = "n", xaxt = "n", yaxt = "n", 
         xlab = "", ylab = "", ylim = c(1, dim(img)[1]), xlim = c(1, 
                                                                  (dim(img)[2])), asp = 1)
    if (add_summary) 
      mtext(string.db$get_summary(string_ids), cex = 0.7)
    if (add_link) 
      mtext(string.db$get_link(string_ids, payload_id = payload_id, required_score = required_score), cex = 0.7, side = 1)
    rasterImage(img, 1, 1, (dim(img)[2]), (dim(img)[1]))
    # dev.off()
  }
}

myStringPlots = function(data,string.db,output,timeChar){
  payload = string.db$post_payload(stringIds = data$STRING_id,colors = data$color) # New object carying the "Halo" info for colouring the graph.
  # Code to enable saving of plots...
  # my_plot_network(string_ids = data$STRING_id,payload_id = payload)
#  filename = paste(output,timeChar,"_STRING.jpg",sep="")
#  my_plot_network(string_ids = data$STRING_id,payload_id = payload,required_score = NULL,add_link = T,add_summary = F,string.db = string.db,file = filename)
  # dev.off()
  # Save plot
  # Code to enable saving of plots...
  # string.db$plot_ppi_enrichment(string_ids = data$STRING_id, quiet=TRUE )
  # Save plot
  return(payload)
}

myCollect = function(lol,num){ # List of Lists, number of Item to collect
  lapply(X = lol,FUN = function(x){x[[num]]})
}

# A good idea but doesn't seem to work so I no longer define the function but rather just let it run.
# initiateString = function(proteome,speciesID=246196){
# cat(paste("Species ID is:",speciesID,"\n"))

# Set things up
string.db <- STRINGdb$new( version="10", species=246196,score_threshold=0, input_directory="" ) # Create a new instance of the STRINGdb (R6) class that behaves more like objects in Ruby than what you are likely familiar with in R
identifiers_proteome = string.db$map(my_data_frame = proteinGroups_clean_norm_iBAQ_cut_ComBat,my_data_frame_id_col_names = "Protein.IDs",removeUnmappedRows = T)
identifiers_CW_bg = string.db$map(my_data_frame = data_CW_added,my_data_frame_id_col_names = "Protein.IDs",removeUnmappedRows = T)
identifiers_Cyto_bg = string.db$map(my_data_frame = data_Cyto_added,my_data_frame_id_col_names = "Protein.IDs",removeUnmappedRows = T)
identifiers_CD_bg = string.db$map(my_data_frame = data_CD_added,my_data_frame_id_col_names = "Protein.IDs",removeUnmappedRows = T)
# identifiers_proteome = string.db$map(my_data_frame = proteome,my_data_frame_id_col_names = "Entry",removeUnmappedRows = T)

# string.db$set_background(background_vector = identifiers_proteome$STRING_id) # Make sure the proteome I used is the background for enrichment p-value calculations later.

# cat(paste("\nFinished setting up STRING instance.\n"))
# Sys.sleep(5)
# return(string.db)
# }

myStringFunction = function(data,background_identifiers=identifiers_proteome,string.db=string.db,output,timeChar){
  string.db$set_background(background_vector = background_identifiers$STRING_id) 
  # Make sure those proteins entering the t-test comparison are used as the background for enrichment p-value calculations.
  
  data = string.db$map(my_data_frame = data,my_data_frame_id_col_names = "Majority.protein.IDs",removeUnmappedRows = T) # Map data to internal STRING_IDs
  data = string.db$add_proteins_description(data) # Add columns with annotations (similar, but not identical with, Uniprot Fasta descriptions)
  data = string.db$add_diff_exp_color(screen = data,logFcColStr = "logFC") # Add colour by FC
  
  payload = myStringPlots(data = data,string.db = string.db,output = output,timeChar = timeChar) # For when wanting to generate interaction and PPI Enrichment plots
  print(paste("Payload: ",as.character(payload)))
  
  GO_bp = string.db$get_enrichment(string_ids = data$STRING_id,category = "Process",methodMT = "fdr")
  GO_fun = string.db$get_enrichment(string_ids = data$STRING_id,category = "Function",methodMT = "fdr")
  GO_com = string.db$get_enrichment(string_ids = data$STRING_id,category = "Component",methodMT = "fdr")
  KEGG = string.db$get_enrichment(string_ids = data$STRING_id,category = "KEGG",methodMT = "fdr")
  
  outSTRING = list(data=data,GO_bp=GO_bp,GO_com=GO_com,GO_fun=GO_fun,KEGG=KEGG)
  return(outSTRING)
}

output_GO = checkMake(output,"GO_Enrichment/")

# # string.db = initiateString(proteome = proteome,speciesID = 246196)
CW = myStringFunction(data = data_CW_pValC, background_identifiers = identifiers_CW_bg, string.db = string.db,output_GO,"CW")
Cyto = myStringFunction(data = data_Cyto_pValC, background_identifiers = identifiers_Cyto_bg, string.db = string.db,output_GO,"Cyto")
CD = myStringFunction(data = data_CD_pValC, background_identifiers = identifiers_CD_bg, string.db = string.db,output_GO,"CD")
lol = list(CW=CW,Cyto=Cyto,CD=CD)

get_payload <- function(data,id_col_name = "Majority.protein.IDs",fold_change_col = "logFC"){
  data = string.db$map(my_data_frame = data,my_data_frame_id_col_names = id_col_name,removeUnmappedRows = T) # Map data to internal STRING_IDs
  data = string.db$add_proteins_description(data) # Add columns with annotations (similar, but not identical with, Uniprot Fasta descriptions)
  data = string.db$add_diff_exp_color(screen = data,logFcColStr = "logFC") # Add colour by FC
  payload = string.db$post_payload(stringIds = data$STRING_id,colors = data$color)
  return(payload)
}

# Make a new dataframe with the present/absent proteins added.
if (doPresAbs) {
  
  data_CW_pValC_pa = rbind(
    select(.data = data_CW_pValC,
           sapply(c("Majority.protein.IDs","logFC","pVals"),function(x){
             grep(pattern = x,x = colnames(data_CW_pValC))
           }) #c("Majority.protein.IDs","logFC","pVals")),
    ),
    data_CW_pa
  )
  data_Cyto_pValC_pa = rbind(
    select(.data = data_Cyto_pValC,
           sapply(c("Majority.protein.IDs","logFC","pVals"),function(x){
             grep(pattern = x,x = colnames(data_Cyto_pValC))
           }) #c("Majority.protein.IDs","logFC","pVals")),
    ),
    data_Cyto_pa
  )
  data_CD_pValC_pa = rbind(
    select(.data = data_CD_pValC,
           sapply(c("Majority.protein.IDs","logFC","pVals"),function(x){
             grep(pattern = x,x = colnames(data_CD_pValC))
           }) #c("Majority.protein.IDs","logFC","pVals")),
    ),
    data_CD_pa
  )

## Do this with pres/abs data included ##
CW_pa = myStringFunction(data = data_CW_pValC_pa, string.db = string.db,output_GO,"CW_pa")
Cyto_pa = myStringFunction(data = data_Cyto_pValC_pa,string.db = string.db,output_GO,"Cyto_pa")
CD_pa = myStringFunction(data = data_CD_pValC_pa,string.db = string.db,output_GO,"CD_pa")
lol_pa = list(CW=CW_pa,Cyto=Cyto_pa,CD=CD_pa)

payload_CW = get_payload(data_CW_pValC_pa)
payload_Cyto = get_payload(data_Cyto_pValC_pa)
payload_CD = get_payload(data_CD_pValC_pa)
}

payload_CW = get_payload(data_CW_pValC)
payload_Cyto = get_payload(data_Cyto_pValC)
payload_CD = get_payload(data_CD_pValC)

## Do this with full proteome as background
# CW_fp = myStringFunction(data = data_CW_pValC,string.db = string.db,output,"CW")
# Cyto_fp = myStringFunction(data = data_Cyto_pValC,string.db = string.db,output,"Cyto")
# CD_fp = myStringFunction(data = data_CD_pValC,string.db = string.db,output,"CD")
# lol_fp = list(CW=CW_fp,Cyto=Cyto_fp,CD=CD_fp)

titles = c("Biological Processes","Cellular Compartments","Molecular Functions","KEGG Pathways")
colours = c("blue","brown","red","yellow")

############# Define GO Term Levels
# 
# library(GO.db)
# 
# getAllBPChildren <- function(goids){
#   ans <- unique(unlist(mget(goids, GOBPCHILDREN), use.names=FALSE))
#   ans <- ans[!is.na(ans)]
# }
# 
# level1_BP_terms <- getAllBPChildren("GO:0008150")     # 23 terms
# level2_BP_terms <- getAllBPChildren(level1_BP_terms)  # 256 terms
# level3_BP_terms <- getAllBPChildren(level2_BP_terms)  # 3059 terms
# level4_BP_terms <- getAllBPChildren(level3_BP_terms)  # 9135 terms
# level5_BP_terms <- getAllBPChildren(level4_BP_terms)  # 15023 terms
# 
# # library(org.Hs.eg.db)
# # level5_genes <- mget(intersect(level5_BP_terms,
# #                                keys(org.Hs.egGO2EG)),
#                      org.Hs.egGO2EG)
#############
# addProtNum <- function(holder=CW, levelx_terms=level1_BP_terms, GO_type="bp"){
#   GO_type = paste0("GO_",as.character(GO_type))
#   data = holder[[GO_type]]
#   data.frame(
#     Num_Proteins = apply(X = levelx_terms[1],FUN = function(GO_term){grep(pattern = as.character(GO_term),x = data$term_id)})
#                          )
# }
# level1_BP_terms$proteins <- apply(X = levelgrep
# CW$GO_bp["GO level"] <- grep()
names = c("CW","Cyto","CD")
# names = c("background")

if (doPresAbs) {
  lol = lol
  suffix_linker = "_pa_"
} else {
  lol = lol_pa
  suffix_linker = "_"
}

for (l in 1:length(lol)){
  GOs = lol[[paste(names[l])]][-1]
  for (i in 1:length(GOs)){
    data = GOs[[i]] %>%
      arrange((desc(pvalue))) %>%
      filter(pvalue_fdr<=0.05) %>%
      filter(proteins<100)
    p = max(nrow(data)-14,1) # Max number of items less one (i.e. 15 items to display)
    data = data[p:nrow(data),]
    data$term_ordered = factor((str_wrap(data$term_description,width = 30)),levels = (str_wrap(data$term_description,width = 30)))
    # data$order = as.numeric(1:length(data))
    p = ggplot(data = data,mapping = aes(x=term_ordered,y=as.integer(hits))) + geom_bar(stat = "identity",fill=colours[i],colour=colours[i]) + coord_flip() + theme(legend.position = "none",text = element_text(size = 30)) + xlab("Term Descriptions") + ylab("Number of Dysregulated Proteins") + ggtitle(titles[i])
    ggsave(paste(output_GO,"",names[l],"_",titles[i],suffix_linker,"100.jpg",sep = ""),plot = p,scale = Scale)
    #+ guide_legend(postition=NULL)#+ scale_x_reverse()#+ scale_x_discrete(labels = data$term_description)
  }
}

# S_data = myCollect(lol,1)
# GO_bp = myCollect(lol,2)
# GO_com = myCollect(lol,3)
# GO_fun = myCollect(lol,4)
# KEGG = myCollect(lol,5)
# 
# dysregulatedProteins.genes = lapply(X = S_data,FUN = function(x){x$preferred_name})
# dysregulatedProteins.genes.unlist = unique(unlist(dysregulatedProteins.genes))
# 
# ####
# # Regulons
# devR_Reg = read.table(file = "C:/Users/Alexander/Desktop/RifAnalysisTemp/DevR_Regulon.tsv",header = F,sep = "\t")
# colnames(devR_Reg) <- c("id","locus_id","gene")
# devR_Reg = mutate(devR_Reg, uniprot = sapply(X = devR_Reg$locus_id,FUN = function(x){proteome$Entry[grep(pattern = x,x = proteome$Gene.names)]}))
# 
# ####
# # Set up the data for use with pathview
# alls = full_join(
#   x = CW$data[,c("Majority.protein.IDs","logFC")],
#   y = full_join(Cyto$data[,c("Majority.protein.IDs","logFC")],CD$data[,c("Majority.protein.IDs","logFC")],by="Majority.protein.IDs")
# )
# all2 = alls[,2:4]
# colnames(all2) <- c("CW","Cyto","CD")
# ##
# # library(UniProt.ws)
# # up <- UniProt.ws(taxId=246196)
# # select(up, uniprots, "ENTREZ_GENE")
# rownames(all2) <- alls$Majority.protein.IDs
# all2[is.na(all2)] <- 0
# 
# # Regulons
# devR_Reg[which(devR_Reg$uniprot %in% alls$Majority.protein.IDs),]
# devR_Reg[which(devR_Reg$uniprot %in% S_data$CD$Majority.protein.IDs),]
# 
# # Select KEGG pathways of interest.
# paths = lapply(X = KEGG,FUN = function(x){filter(x,(pvalue<0.05 & proteins<50))$term_id})
# paths = c("00860")
# pathsUnique = unique(unlist(paths))
# 
# ###
# stringLists = lapply(X = S_data,FUN = function(x){x$STRING_id})
# pathways = lapply(X = KEGG,FUN = function(x){x$term_id[1:10]})
# 
# eh.bp = string.db$enrichment_heatmap(genesVectors = stringLists,vectorNames = c("CW","Cyto","CD"),title = "Biological Processes",enrichmentType = "Process",pvalue_threshold = pVal.CutOff,cexRow = 0.5,cexCol = 0.01)
# eh.com  = string.db$enrichment_heatmap(genesVectors = stringLists,vectorNames = c("CW","Cyto","CD"),title = "Cellular Component",enrichmentType = "Component",pvalue_threshold = pVal.CutOff,limit = 20)
# ###
# 
# library(KEGGREST)
# library(pathview)
# for (i in 1:length(pathsUnique)){
#   pv.out <- pathview(gene.data = all2, pathway.id = pathsUnique[i], gene.idtype = "kegg", species = "msm", out.suffix = "same_T-match_T", keys.align = "y", kegg.native = T, match.data = T, multi.state = T, same.layer = F)
# }



# #################################################################

# # get_STRING_species(version="10", species_name="smegmatis") # Use this when looking for your organisms identifier or when seeing if your organism is listed
# string.db <- STRINGdb$new( version="10", species=246196,score_threshold=0, input_directory="" ) # Create a new instance of the STRINGdb (R6) class that behaves more like objects in Ruby than what you are likely familiar with in R
# # STRINGdb$methods() # Use this to see which methods are available for your R6 class object.
# # STRINGdb$help("benchmark_ppi_pathway_view") # Use this to get help with those methods.
# identifiers_proteome = string.db$map(my_data_frame = proteome,my_data_frame_id_col_names = "Entry",removeUnmappedRows = T)
# string.db$set_background(background_vector = identifiers_proteome$STRING_id)
# identifiers_CW = string.db$map(my_data_frame = data_CW_pValC,my_data_frame_id_col_names = "Majority.protein.IDs",removeUnmappedRows = T) # Add a column that translates an identifying column in your data to an internal STRING ID.
# # head(identifiers_CW)
# identifiers_CW = string.db$add_proteins_description(identifiers_CW) # Add columns with annotations (similar, but not identical with, Uniprot Fasta descriptions)
# 
# # string.db$plot_network(string_ids = identifiers_CW$STRING_id,add_summary = T,add_link = T) # The extra options specified are the defaults anyway.
# 
# # identifiers_CW %>%
# #   mutate(FC_Colour = ifelse(logFC>0,yes = "Yellow",no = "Red")) # Not used right now.
# 
# identifiers_CW = string.db$add_diff_exp_color(screen = identifiers_CW,logFcColStr = "logFC")
# payload_id_CW = string.db$post_payload(stringIds = identifiers_CW$STRING_id,colors = identifiers_CW$color)
# string.db$plot_network(string_ids = identifiers_CW$STRING_id,payload_id = payload_id_CW)
# 
# string.db$plot_ppi_enrichment(string_ids = identifiers_CW$STRING_id, quiet=TRUE )
# 
# # Gene Ontology and KEGG Pathway Enrichment
# GO_bp_CW = string.db$get_enrichment(string_ids = identifiers_CW$STRING_id,category = "Process",methodMT = "fdr")
# KEGG_CW = string.db$get_enrichment(string_ids = identifiers_CW$STRING_id,category = "KEGG",methodMT = "fdr")
# ```



# Colour code each protein by up/down regulation
# In progress...
# ```{r echo=FALSE}
# kegg = export2Perseus[,c("Majority.protein.IDs","")]
# ```


# Extra Bits of code
# ```{r,echo=FALSE}
# ggplotColours <- function(n=6, h=c(0, 360) +15){
#   if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
#     hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
# }
# 
# ggplotColours(2)
# # [1] "#F8766D" "#00BFC4"
# 
# ggplotColours(3)
# # [1] "#F8766D" "#00BA38" "#619CFF"

# library(colorspace)
# pal = choose_palette()

# ```
# test_background <- myStringFunction(data = test_data_background,string.db = string.db,output = output,timeChar = "background")

addFasta <- function(df,ref_df,id_col="Majority.protein.IDs",ref_id_col="Majority.protein.IDs",ref_fasta_col="Fasta.headers"){ #Use Uniprot accession to match Fastas and cbind
  ind <- sapply(df[,id_col],function(x){grep(x,ref_df[,ref_id_col])})
  fastas <- as.character(ref_df[,ref_fasta_col][ind])
  new_df <- cbind(df,Fasta.headers=fastas)
  return(new_df)
}

dysregulatedProteins_pVal <- cbind(
  addFasta(ref_df = proteinGroups_clean,
      df = rbind(data_CD_pValC_pa,data_CW_pValC_pa,data_Cyto_pValC_pa)),
  Compartment = c(rep("CD",nrow(data_CD_pValC_pa)),rep("CW",nrow(data_CW_pValC_pa)),rep("Cyto",nrow(data_Cyto_pValC_pa)))
)
  
dys_cols <- c("Majority.protein.IDs","Fasta.headers","pVals","logFC")
dysregulatedProteins_FC   <- rbind(data_CD_FC[,dys_cols],data_CW_FC[,dys_cols],data_Cyto_FC[,dys_cols])
dysregulatedProteins_FC[,"Compartment"] <- c(rep("CD",nrow(data_CD_FC)),rep("CW",nrow(data_CW_FC)),rep("Cyto",nrow(data_Cyto_FC)))

write.table(x = dysregulatedProteins_FC,file = paste0(output,"dysregulatedProteins_FC.tsv"),quote = T,sep = "\t",row.names = F)
write.table(x = dysregulatedProteins_pVal,file = paste0(output,"dysregulatedProteins_pVal.tsv"),quote = T,sep = "\t",row.names = F)


transporters_ABC <- dysregulatedProteins_pVal %>% filter(Compartment=="CW") %>% filter(grepl(pattern="ABC|[Pp]ermease|[Tt]ransporter",Fasta.headers))
# Rif_vs_CW <- read.table(file = "G:/2017/TB/Clemens/Paper_with_Alex/Alex/Rif_vs_CW.csv",header = T,sep = ",",quote = "\"",as.is = T)
# newlyDysregulated <- setdiff(
#   union(
#     union(CW$data$Majority.protein.IDs,Cyto$data$Majority.protein.IDs),
#     CD$data$Majority.protein.IDs),
#   union(
#     union(Rif_vs_CW$T1,Rif_vs_CW$T2),
#     Rif_vs_CW$T3)
# )
# write.table(x = newlyDysregulated,file = "C:/Users/user/Desktop/newlyDysregulated.txt",quote = T,sep = "\t",row.names = F)
