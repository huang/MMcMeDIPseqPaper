#!/usr/bin/Rscript

# script runs the various QC steps for the medusa pipeline, in addition it creates wig files - utilises the MEDIPS bioconductor package.
args <- commandArgs();
sample<-args[length(args)-5]
species<-args[length(args)-4]
refName<-args[length(args)-3]
windowSize<-as.numeric(args[length(args)-2])
strands<-as.numeric(args[length(args)-1])
path2output<-args[length(args)]


library(MEDIPS)
if (species == "Mouse"){
	bsGenome<-paste("BSgenome.Mmusculus.UCSC.",refName,sep="")
	library(bsGenome,character.only=T)
	genome <- bsGenome
} else if (species == "Human"){
	bsGenome<-paste("BSgenome.Hsapiens.UCSC.",refName,sep="")
	library(bsGenome,character.only=T)
	genome <- bsGenome
} else if (species == "Dog"){
	bsGenome<-paste("BSgenome.Cfamiliaris.UCSC.",refName,sep="")
	library(bsGenome,character.only=T)
	genome <- bsGenome
} else if (species == "Chimp"){
	bsGenome<-paste("BSgenome.Ptroglodytes.UCSC.",refName,sep="")
	library(bsGenome,character.only=T)
	genome <- bsGenome
} else if (species == "Macaca"){
	bsGenome<-paste("BSgenome.Mmulatta.NCBI.",refName,sep="")
	library(bsGenome,character.only=T)
	genome <- bsGenome
} else {
        stop("Species must be Mouse, Human, Chimp, Macaca or Dog. Currently using mouse mm9, human hg19, chimp panTro2, macaca mmul1 and dog canFam2")
}

setwd(path2output)
sampleWig.set<-MEDIPS.createSet(BSgenome=genome,file=paste(path2output,"/",sample,".bed",sep=""),window_size=windowSize)
#sampleWig.set<-MEDIPS.createSet(BSgenome=genome,file="MMC_2pos.bed",window_size=100)
MEDIPS.exportWIG(Set=sampleWig.set,file=paste(sample,"_rpkm.wig",sep=""),format="rpkm")
#MEDIPS.exportWIG(Set=sampleWig.set,file=paste("MMC_2pos","_rpkm.wig",sep=""),format="rpkm")
if(strands == 0) {
#function (file = NULL, BSgenome = NULL, nit = 10, nrit = 1, empty_bins = TRUE, 
#	      rank = FALSE, extend = 0, shift = 0, window_size = 500, uniq = 0.001, 
#	          chr.select = NULL, paired = F, isSecondaryAlignment = FALSE, 
#	          simpleCigar = TRUE) 
sr<-MEDIPS.saturation(file=paste(path2output,"/",sample,".bed",sep=""),BSgenome=genome,uniq=1e-3,extend=0,shift=0,window_size=windowSize,nit=10,rank=T)
#sr<-MEDIPS.saturation(file="MMC_2pos.bed",BSgenome=genome,uniq=1e-3,extend=0,shift=0,window_size=100,nit=10,rank=T)
png(file=paste(sample,"_saturation.png",sep=""),type="cairo")
MEDIPS.plotSaturation(sr)
dev.off()

#cr<-q(file="MMC_2pos.bed",BSgenome=genome,uniq=1e-3,extend=0,shift=0) --> What is q() function?
#function (file = NULL, BSgenome = NULL, pattern = "CG", extend = 0, 
#    shift = 0, uniq = 0.001, chr.select = NULL, paired = F, isSecondaryAlignment = FALSE, 
#    simpleCigar = TRUE) 
cr<-MEDIPS.seqCoverage(file=paste(path2output,"/",sample,".bed",sep=""),pattern="CG",BSgenome=genome,uniq=1e-3,extend=0,shift=0)
png(file=paste(sample,"_coveragePie.png",sep=""),type="cairo")
#MEDIPS.plotSeqCoverage(seqCoverageObj=cr, main="Sequence pattern coverage", type="pie", cov.level = c(0,1,2,3,4,5))
MEDIPS.plotSeqCoverage(cr,type="pie",cov.level=c(0,1,5,10,20,50))
dev.off()
#function (file = NULL, BSgenome = NULL, extend = 0, shift = 0, uniq = 0.001, chr.select = NULL, paired = F) 
er<-MEDIPS.CpGenrich(file=paste(path2output,"/",sample,".bed",sep=""),BSgenome=genome,uniq=1e-3,extend=0,shift=0)
write.table(er,paste(sample,"_enrichment.txt",sep=""),sep="\t",quote=F)
}



#BiocManager::install(version = "3.10")
# 
# library(MEDIPS)
# library("BSgenome.Mmusculus.UCSC.mm10",character.only=T)
# genome <- "BSgenome.Mmusculus.UCSC.mm10"
# setwd("/home/jhuang/DATA/Data_Nicole_MeDIP/alg")
# sample="MMC_2neg"
# 
# sampleWig.set<-MEDIPS.createSet(BSgenome=genome,file=paste(sample,".bed",sep=""),window_size=100)
# MEDIPS.exportWIG(Set=sampleWig.set,file=paste(sample,"_rpkm.wig",sep=""),format="rpkm")
# 
# sr<-MEDIPS.saturation(file=paste(sample,".bed",sep=""),BSgenome=genome,uniq=1e-3,extend=0,shift=0,window_size=windowSize,nit=10,rank=T)
# png(file=paste(sample,"_saturation.png",sep=""),type="cairo")
# MEDIPS.plotSaturation(sr)
# dev.off()
# 
# cr=MEDIPS.seqCoverage(file=paste(sample,".bed",sep=""), pattern="CG", BSgenome=genome, uniq=1e-3,extend=0,shift=0)
# png(file=paste(sample,"_coveragePie.png",sep=""),type="cairo")
# MEDIPS.plotSeqCoverage(cr,type="pie",cov.level=c(0,1,5,10,20,50))
# dev.off()
# 
# er<-MEDIPS.CpGenrich(file=paste(sample,".bed",sep=""),BSgenome=genome,uniq=1e-3,extend=0,shift=0)
# write.table(er,paste(sample,"_enrichment.txt",sep=""),sep="\t",quote=F)






# pkgs <- c(
#     "MEDIPS", "S4Vectors", "IRanges", "GenomicRanges", "DelayedArray",
#     "XVector", "GenomicAlignments", "ShortRead",
#     "VariantAnnotation", "AnnotationHub", "GGtools",
#     "ggbio", "ChemmineR", "InteractionSet", "flowCore",
#     "GenomicTuples", "CNEr", "MultiAssayExperiment",
#     "genomeIntervals", "TFBSTools", "IWTomics", "spliceSites",
#     "podkat", "kebabs", "matter", "dada2",
#     "ClusterSignificance", "gespeR", "HiTC", "tigre", "soGGi"
# )
# update <- intersect(rownames(installed.packages()), pkgs)
# #biocLite(update, type="source", ask=FALSE)
# BiocManager::install(update, type="source", ask=FALSE)
# 
# library(MEDIPSData)
# library("BSgenome.Hsapiens.UCSC.hg19")
# bam.file.hESCs.Rep1.MeDIP = system.file("extdata", "hESCs.MeDIP.Rep1.chr22.bam", package="MEDIPSData")
# #er=MEDIPS.CpGenrich(file=bam.file.hESCs.Rep1.MeDIP, BSgenome="BSgenome.Hsapiens.UCSC.hg19", chr.select="chr22", 
# er=MEDIPS.CpGenrich(file=bam.file.hESCs.Rep1.MeDIP, BSgenome="BSgenome.Hsapiens.UCSC.hg19", chr.select="chr22",extend="150")

