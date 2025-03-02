library(dplyr)

ASS=commandArgs(trailingOnly=TRUE)[1]
ARGfile=paste0("02.ASSEMBLIES/",ASS,"/",ASS,".ARGS.coverages")
MGEfile=paste0("02.ASSEMBLIES/",ASS,"/",ASS,".MGES.coverages")
KKfile=paste0("02.ASSEMBLIES/",ASS,"/",ASS,".kk_slow.out")
OUTfile=paste0("03.COVERAGEPLOTS/",ASS,"_Contigs_coverage.png")

# Parse the ARG detection output table
ARGtable=read.table(ARGfile,sep="_",header=FALSE, fill=TRUE, col.names=c("contig size","contig coverage"))
#colnames(ARGtable) = c("contig size","contig coverage")

# Parse the MGE detection output table
MGEtable=read.table(MGEfile,sep="_",header=FALSE, fill=TRUE, col.names=c("contig size","contig coverage"))

# Find contigs containing both ARG and MGE.
AMtable=inner_join(ARGtable,MGEtable)

# Parse the Kraken output table, and select contigs > 150 bp
KKtmp=read.table(KKfile,sep="_",header=FALSE, fill=TRUE)
KKtable=KKtmp[KKtmp[,5]>150,c(5,7,8)]
colnames(KKtable) = c("contig size (>150 bp)","contig coverage","Taxon")


# Extract taxonomic names (up to genus)
ALLgenera=sub(pattern=" .*",replacement="",x=KKtable[,3],perl=TRUE)
DEREPgenera=rle(ALLgenera[order(ALLgenera)])

# Find the most matched taxon and select the corresponding contigs
Best1stmatch=DEREPgenera$values[which.max(DEREPgenera$lengths)]
Best1stc = KKtable[grep(Best1stmatch,KKtable[,3]),c(1,2)]

# Remove the most matched taxon
DEREPgenera$values<-DEREPgenera$values[-(which.max(DEREPgenera$lengths))]
DEREPgenera$lengths<-DEREPgenera$lengths[-(which.max(DEREPgenera$lengths))]

# Find the second most matched taxon and select the corresponding contigs
Best2ndmatch=DEREPgenera$values[which.max(DEREPgenera$lengths)]
Best2ndc = KKtable[grep(Best2ndmatch,KKtable[,3]),c(1,2)]

# Remove the second most matched taxon
DEREPgenera$values<-DEREPgenera$values[-(which.max(DEREPgenera$lengths))]
DEREPgenera$lengths<-DEREPgenera$lengths[-(which.max(DEREPgenera$lengths))]

# Find the third most matched taxon and select the corresponding contigs
Best3rdmatch=DEREPgenera$values[which.max(DEREPgenera$lengths)]
Best3rdc = KKtable[grep(Best3rdmatch,KKtable[,3]),c(1,2)]

# Compute the max coverage value for contigs, to scale the graph at the next 100th value  
maxylim=max(ARGtable[,2],MGEtable[,2],KKtable[,2])
maxylim=(1+as.integer(maxylim/100))*100

# create the figure
png(file=OUTfile,width=1300,height=700)
par(mar=c(5,4.5,4,2), cex=2, cex.lab=1.2)

# Plot the contigs
#------------------------------------------------------------------------
colone="grey"
coltwo="#cea30d" # Staphylococcal specific color

# Plot contigs assigned to none of the three most matched taxon (others)
plot(x=KKtable[,1],y=KKtable[,2], ylim=c(0.5,maxylim), main=ASS, log="xy", pch=20, cex=1.5, col="grey", yaxp=c(1,200,1), xlab=colnames(KKtable)[1], ylab=colnames(KKtable)[2])

# Generic
#points(Best1stc, pch=21, cex=1.5, bg=white)
#points(Best2ndc, pch=20, cex=1.5, col=black)
#points(Best3rdc, pch='+', cex=2, col="black")

# for ERR212931
#points(Best2ndc, pch=21, cex=1.5, bg=colone)
#points(Best1stc, pch=21, cex=1.5, bg=coltwo)

# for SRR1522630
points(Best2ndc, pch=21, cex=1.5, bg=coltwo)
points(Best1stc, pch=21, cex=1.5, bg=colone)

points(MGEtable, pch=21, cex=1.5, bg="orange")
points(ARGtable, pch=21, cex=1.5, bg="cyan1")
points(AMtable,cex=1.6, pch=21, bg="red")

# Print legend
#------------------------------------------------------------------------
# Generic
#legend(x="right",inset=-0.1,xpd=TRUE, legend=c("MGE-containing","ARG-containing","MGE+ARG-containing",Best1stmatch,Best2ndmatch,Best3rdmatch,'Other'),pch=c(21,21,21,21,21,43,20),pt.cex=1.1,pt.bg=c("orange","cyan1","red","mediumseagreen","khaki1","black","grey"),col=c("black","black","black","black","black","black","grey"))

# for ERR212931
#legend(x="bottomright",inset=c(-0.15,0), ncol=3, bty="n", xpd=TRUE, legend=c("ARG-containing","MGE+ARG-containing",Best2ndmatch,Best1stmatch,'Other'),pch=c(21,21,21,21,20),pt.cex=1.1,pt.bg=c("cyan1","red",coltwo,colone,"grey"),col=c("black","black","black","black","grey"))

# for SRR1522630
legend(x="bottomright",inset=c(-0.15,0), ncol=3, bty="n", xpd=TRUE, legend=c("ARG-containing",Best2ndmatch,Best1stmatch,'Other'),pch=c(21,21,21,20),pt.cex=1.1,pt.bg=c("cyan1",coltwo,colone,"grey"),col=c("black","black","black","grey"))

dev.off()
