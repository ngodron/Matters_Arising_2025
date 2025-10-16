library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggbeeswarm)

ASS=commandArgs(trailingOnly=TRUE)[1]
COVfile=commandArgs(trailingOnly=TRUE)[2]
MGEfile=commandArgs(trailingOnly=TRUE)[3]
KKfile=commandArgs(trailingOnly=TRUE)[4]
OUTfile=commandArgs(trailingOnly=TRUE)[5]

minlim=0.1

# Parse the COVERAGE output table
COVtable=read.table(COVfile,sep="\t",header=FALSE, col.names=c("contig.id","startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq"))
COVtable[which(COVtable$meandepth == 0),"meandepth"]=0.1

# Parse the MGE detection output table
MGEtable=read.table(MGEfile,sep="\t",header=FALSE, fill=TRUE)[,c(1,2)]
colnames(MGEtable)<-c("IS","contig.id")

# dirty modification to derepplicate contig ids.
MGEtable$IS=MGEtable$contig.id
MGEtable=unique(MGEtable)

# add coverage and size info
MGEtable=cbind(MGEtable,COVtable[grep(paste(MGEtable$contig.id, collapse="|"),COVtable$contig.id),c(3,7)])


# Parse the Kraken output table, and select contigs > 500 bp
KKtmp=read.table(KKfile,sep="_",header=FALSE, fill=TRUE)
KKtable=KKtmp[,c(2,3,4)]
colnames(KKtable) = c("prefix","contig.id","Taxon")
KKtable$contig.id=paste0(KKtable$prefix,"_",KKtable$contig.id)
KKtable=cbind(KKtable[order(KKtable$contig.id),],COVtable[order(COVtable$contig.id),c(3,7)])


# Extract taxonomic names (up to genus)
ALLgenera=sub(pattern=" .*",replacement="",x=KKtable$Taxon,perl=TRUE)
DEREPgenera=rle(ALLgenera[order(ALLgenera)])


#Enterobacetriaceae taxa
Enterotaxa="Enterobacter|Escherichia|Salmonella|Klebsiella|Proteus|Citrobacter"
#For expected Actinomycetes
Exptaxon="Amycolatopsis|Saccharomonospora|Streptomyces" 
#For exprected Bacillota
Exptaxon="Bacill|Clostridi|Mediterraneibacter|Enterocloster|Thomasclavelia|Staphylo|Streptocco" 

# Select the contigs matching the expected taxon
Exptaxonc = KKtable[grep(Exptaxon,KKtable$Taxon),c(4,5)]

# Remove the expected taxon
DEREPgenera$lengths<-DEREPgenera$lengths[-grep(Exptaxon,DEREPgenera$values)]
DEREPgenera$values<-DEREPgenera$values[-grep(Exptaxon,DEREPgenera$values)]

# Find specifically Enterobactriaceae or Escherichia contigs
Enterotaxac = KKtable[grep(Enterotaxa,KKtable$Taxon),c(4,5)]

# Remove Enterobacteriaceae taxa
DEREPgenera$lengths<-DEREPgenera$lengths[-grep(Enterotaxa,DEREPgenera$values)]
DEREPgenera$values<-DEREPgenera$values[-grep(Enterotaxa,DEREPgenera$values)]

#Select all other taxa contigs
Othertaxac=KKtable[grep(paste0(DEREPgenera$values,collapse="|"),KKtable$Taxon),c(4,5)]


Expnb=nrow(Exptaxonc)
Entnb=nrow(Enterotaxac)
ISnb=nrow(MGEtable)

data <- data.frame(
 all   = c(rep(ASS,Expnb+Entnb+ISnb) ),
 contig  = c( rep("Expected taxon",Expnb), rep("Enterobacteriaceae",Entnb), rep("IS1-carrying",ISnb) ),
 value = c(Exptaxonc$meandepth,Enterotaxac$meandepth,MGEtable$meandepth)
)


# Compute the max coverage value for ARG or MGE-containing contig, to scale the graph at the next 100th value  
#maxylim=max(KKtable$meandepth)
#maxylim=(1+as.integer(maxylim/100))*100
maxylim=5000


png(file=OUTfile,width=680,height=850)
par(mar=c(1,1,1,1))

data %>%
  ggplot( aes(x=all, y=value, color=contig)) +
  scale_y_log10(limits=c(0.09,maxylim)) +
  geom_beeswarm(method = "center", size=c(rep(6,Expnb+Entnb),rep(10,ISnb)), pch=c(rep(1,Expnb+Entnb),rep(16,ISnb)), stroke=2.5 ) +
  theme(text=element_text(size=40)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800","magenta")) +
  xlab("")

dev.off()


