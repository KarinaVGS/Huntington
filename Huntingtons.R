#************************************************************
# Project:  		Huntington's Disease (wild type and abnormal)
# Course:			Introducción a Bioinformática
# Professor:  		Claudia Rangel
# Midterm: 			October 2, 2019
# Data 2 groups: 
#					B) WiTy (WT_1, WT_2, WT_3, WT_4)
#					C) polyQ (polyQ_1, polyQ_2, polyQ_3, polyQ_4)
# Total of 8 microarrays
# ***********************************************************

# --------------------------------------------------------------
# R-packages: oligo, affy, gplots
# --------------------------------------------------------------
library(oligo)


# --------------------------------------------------------------
# Data: 8 ".CEL files" 
# --------------------------------------------------------------

setwd("C:/Users/gonza/Documents/School/Bioinformatics/PF_CEL")
dir.celfiles=getwd()
myfiles=list.celfiles(dir.celfiles, full.names = TRUE)
length(myfiles)

# --------------------------------------------------------------
# Read in data
# --------------------------------------------------------------
rawdata=read.celfiles(myfiles)
sN = sampleNames(rawdata)

#------------------------------------------------------------------
# a) reorder the data diles accordingly
# b) label the cel files properly
# c) create colors vector
# d) create index to arrange samples per group
#------------------------------------------------------------------

indx.samples = c(5:8,1:4)
rawdata=rawdata[, indx.samples]
data.labels=c(paste("WildType",1:4, sep=""), paste("polyQ", 1:4, sep=""))
sampleNames(rawdata)=data.labels
sample.colors=c(rep("chartreuse3",4), rep("darkviolet", 4))


#------------------------------------------------------------------
# QC plots
#------------------------------------------------------------------
# These plots should be generated using raw data but also on normalized
# data to see the normalization effects
# Core = transcription

oligo::boxplot(rawdata,target="core",col=sample.colors, main="Raw Data Distribution", cex=0.8)
oligo::hist(rawdata,target="core",col=sample.colors, main="Raw data Density curves", cex=0.8)


# --------------------------------------------------------------
# Normalization & Background correction
# -------------------------------------------------------------- 
# Background correction and Normalization using RMA, read the 
# help for RMA (from oligo package) to set all parameters up
# NOTE: To perform normalization within groups (on replicates)
# it must be done prior to RMA and therefore your variable won´t be "data"
# as it shows in the code. It will now be "data.norm" slide 44 lecture 6

# normalizing each group
data.norm = rawdata
qnorm= normalize(rawdata[,1:4], method="quantile")
pm(data.norm)[,1:4] = pm(qnorm)
qnorm= normalize(rawdata[,5:8], method="quantile")
pm(data.norm)[,5:8] = pm(qnorm)
# plot to see the effect
oligo::boxplot(data.norm,target="core", col= sample.colors, main = "Quantile-normalized by group", cex.axis=.8, las=2)
oligo::hist(data.norm,target="core", col= sample.colors, main = "Quantile-normalized by group", cex=.8, lwd=2)


# normalizing all samples 
data.RMA=oligo::rma(data.norm,target="core")
exp.matrix=exprs(data.RMA)
# plot to see the effect
oligo::boxplot(exp.matrix, target="core", col=sample.colors, main = "RMA-Normalized Data", las =2, ylab = "Expression levels")

# The plotDensity function only works for the affy package 
# so, once affy is loaded, it overwrites some of the oligo 
# functions. Hence, if you need to run previous code lines,
# you will need to re-load oligo.
library(affy)
plotDensity(exp.matrix, col=sample.colors, main = "RMA-Normalized Data", xlab = "Expression Levels")
save(data.RMA, sample.colors, file="DataNorm.RData")

# --------------------------------------------------------------
# Classification: Differential expression
# a) Install and load all required libraries
# b) Create the design matrix
# c) Create the contrasts matrix
# d) Fit the linear models (one per gene)
# e) Plots: volcanoplot and heat map
# f) Obtain the gene expression data for further analyses
# -------------------------------------------------------------- 
library(limma)
library(gtools)
library(gplots)
library(GO.db)
library(annotate)

# Needs mouse library
library(mouse4302.db)
library(mogene10stprobeset.db)

annDB = "mogene10stprobeset.db"

NoTranscripts = 22632 

design = matrix(rep(0,2*8), nrow=8)
colnames(design) = c("WildType", "polyQ")
design[1:4,1] = 1
design[5:8,2] = 1

contrasts.list = c("polyQ vs. WildType");
cont.matrix = makeContrasts("polyQ - WildType", levels=design)


fit = lmFit(exp.matrix,design)
fitC = contrasts.fit(fit, cont.matrix)
fitCB = eBayes(fitC)
TT = topTable(fit=fitCB, coef=1, adjust = "fdr", sort.by = "logFC", number=NoTranscripts, genelist=fit$genes)



# Using B-statistic (B) and log fold-change (M)
# ----------------------------------------------
B = 3
M = 1

# Identify differentially expressed genes according 
# to B and M (variable "selected")
selected = TT[(TT$B>B & abs(TT$logFC)>M),]
# selected = TT[(TT$B>B),]

# Label the rows by searching the gene symbols
probenames = rownames(selected)
entID = getEG(probenames,data = annDB)
sym = getSYMBOL(probenames,data = annDB)
diff.expressed = data.frame(sym, selected)
indx = !is.na(diff.expressed$sym)
diff.expressed = diff.expressed[indx,]


# Extract that list of differentially expressed genes 
write.csv(diff.expressed, file="Bstat_Diff_Expressed.csv")	
htmlpage(genelist = list(entID[indx]), filename = "DiffExp.HTML", title=paste("Differentially Expressed Genes: WildType vs. polyQ"), othernames = diff.expressed, table.head = c("Entrez ID", colnames(diff.expressed)), table.center = TRUE)	


# Volcano Plot (x0 and x1 are the min and max values for log fold-change plun epsilon)
# (y0 and y1 are the min and max for B). This is a color plot	
x0 = -6.5
x1 = 6.5
y0 = -8
y1 = 20
volcanoplot(fitCB, col="blue", ylim=c(y0,y1), xlim=c(x0,x1), coef=1, main="Differential Expression: WildType vs. polyQ", cex.lab=1.3, ylab="B")
par(new=T)
abline(v=-M, col="brown", ylab="", xlab="")
par(new=T)
abline(v=M, col="brown", ylab="", xlab="")
par(new=T)
abline(h=B, col="black", ylab="", xlab="")
par(new=T)
ind1 = abs(fitCB$coef)>M
ind2 = fitCB$lods >B
ind3 = (fitCB$coef>M & fitCB$lods>B)
ind4 = (fitCB$coef< -M & fitCB$lods>B)
x = as.matrix(fitCB$coef[ind1])
y = as.matrix(fitCB$lods[ind1])
plot(x, y, col="magenta",ylim=c(y0,y1), xlim=c(x0,x1),main="", pch = "*",xlab="", ylab="",cex.lab=1.2)
x = as.matrix(fitCB$coef[ind2])
y = as.matrix(fitCB$lods[ind2])
par(new=T)
plot(x, y, col="orange",  ylim=c(y0,y1), xlim=c(x0,x1), main="", pch = 20,xlab="", ylab="",cex.lab=1.2)
x = as.matrix(fitCB$coef[ind3])
y = as.matrix(fitCB$lods[ind3])
par(new=T)
plot(x, y, col="red",  ylim=c(y0,y1), xlim=c(x0,x1), main="", pch = 20,xlab="", ylab="",cex.lab=1.2)
x = as.matrix(fitCB$coef[ind4])
y = as.matrix(fitCB$lods[ind4])
par(new=T)
plot(x, y, col="darkgreen", ylim=c(y0,y1), xlim=c(x0,x1), main="", pch = 20,xlab="", ylab="",cex.lab=1.2)



# -------------------------------------------------------------#
# (1) Hierarchical Clustering
# -------------------------------------------------------------#
# There will be as many heatmaps as contrast of interest, in this example 3:
# LungCa vs BrCa, LungCa vs TNBrCa, and TNBrCa vs BrCa. Make sure labels and 
# titles are set accordingly.
# If you want to generate this plot in a pdf file remove "#" from lines 196 and 199 below
data.clus = exp.matrix[rownames(diff.expressed),]
genelabels = diff.expressed$sym	
#pdf(file="Bstat_Heatmap.pdf", height=16, width=8)
par(oma = c(3,1,3,4),mar=c(12,5,2,2)+0.1)
ind.hmap = heatmap.2(data.clus, col=greenred(75), scale="row", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=.55, cexCol=.7, main = "LungCa vs BrCa", labRow = genelabels, ColSideColors = sample.colors)

