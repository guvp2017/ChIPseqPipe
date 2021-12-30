########################
#  Plot BedCovCompare  #
#  Berkley Gryder      #
#  Feb. 2019           #
########################

##1. Import parameters. Load data. 
  #to test the script, uncomment and use these parameters
    #setwd("\\\\ads.case.edu/rc/SOM_GENE_BEG33/ChIP_seq/hg38/projects/OsteoMet/")
    #file.cov_scaled = "MG63_H3K27ac_MG63.3_H3K27ac_p-7_BCC/MG63_H3K27ac_MG63.3_H3K27ac.RPM.bed"
    #file.bedoverlap = "MG63_H3K27ac_MG63.3_H3K27ac_p-7_BCC/MG63_H3K27ac_MG63.3_H3K27ac.UpMet.TSS.bed"
    #name.overlap = "UpMet.TSS"
    #dir.output = "MG63_H3K27ac_MG63.3_H3K27ac_p-7_BCC"
    #name.sample1 = "MG63_H3K27ac"
    #name.sample2 = "MG63.3_H3K27ac"
    #"MG63_H3K27ac_MG63.3_H3K27ac_p-7_BCC/MG63_H3K27ac_MG63.3_H3K27ac.UpMet.TSS.sortinstitch.RPM.bed"

  #inputs from shell script runBedCovCompare.sh or runBedCovCompareSpiked.sh
  args = commandArgs(TRUE)
    file.cov_scaled = args[1]
    file.bedoverlap = args[2]
    name.overlap = args[3]
    dir.output = args[4]
    name.sample1 = args[5]
    name.sample2 = args[6]
  
  #define outputs
  output.file.prefix = paste(dir.output,"/",name.sample1,"_",name.sample2,sep="")
  output.file.noverlap = paste(output.file.prefix,"_no_",name.overlap,".MA.pdf",sep="")
  output.file.overlap  = paste(output.file.prefix,"_in_",name.overlap,".MA.pdf",sep="")
  output.file.bed  = paste(output.file.prefix,"_",name.overlap,".overlap.AB.delta.l2fc.bed",sep="")
  output.file.up = paste(output.file.prefix,".l2fc_UP.bed",sep="")
  output.file.dn = paste(output.file.prefix,".l2fc_DN.bed",sep="")
  output.file.bednoverlap  = paste(output.file.prefix,"_no_",name.overlap,".l2fc.bedgraph",sep="")
  output.file.bedoverlap  = paste(output.file.prefix,"_in_",name.overlap,".l2fc.bedgraph",sep="")
  
  bedCov <- read.table(file.cov_scaled, sep="\t", header=F) 
  bedOverlap <- read.table(file.bedoverlap, sep="\t", header=F) 
  
##2. Calculate delta log2 fold change
  
  library(ggplot2)
  ggplot(bedCov, aes(x=log2(bedCov[,4]+1),y=log2(bedCov[,5]+1)))+geom_point(alpha=0.1)
  ggplot(bedCov, aes(x=log2(bedCov[,4]+1)))+geom_histogram(bins = 50)
  bedCov$log2fc = log2(bedCov[,5]+1)-log2(bedCov[,4]+1); bedCov$log2fc = round(bedCov$log2fc, digits = 5)
  bedCov$rankL2FC = rank(bedCov$log2fc)
    bedCov.up=subset(bedCov,bedCov$log2fc>1)
    bedCov.dn=subset(bedCov,bedCov$log2fc<(-1))
    write.table(x = bedCov.up,file = output.file.up, col.names = F, row.names = F, quote = F,sep="\t")
    write.table(x = bedCov.dn,file = output.file.dn, col.names = F, row.names = F, quote = F,sep="\t")
  
  rankl2fcplot = ggplot(bedCov, aes(x=rankL2FC,y=log2fc))+geom_point()+geom_point(data=bedCov.up,color="red")+geom_point(data=bedCov.dn,color="blue")+
    geom_text(x=(max(bedCov$rank)-max(bedCov$rank)/10), y=max(bedCov$log2fc), label=paste(nrow(bedCov.up)," regions up, l2fc > 1"),vjust = "inward", hjust = "inward",color="red")+
    geom_text(x=(max(bedCov$rank)/10), y=min(bedCov$log2fc), label=paste(nrow(bedCov.dn)," regions down, l2fc < -1"),vjust = "inward", hjust = "inward",color="blue")+
    labs(y = paste("log2FC ",name.sample1," - ",name.sample2))
  
    output.file.rank = paste(output.file.prefix,".rankl2fc.pdf",sep="")
    ggsave(output.file.rank, plot = rankl2fcplot ,width = 5.5, height = 4)
  
##3. Split bedCov by overlapping genomic feature
  
  bedOverlap$V4[bedOverlap$V4>1] <- 1 #just in case of multiple overlapping features
  
  bedCov$overlap = bedOverlap[,4];bedCov$overlap <- gsub(0, paste("not in ",name.overlap,sep=""), bedCov$overlap);bedCov$overlap <- gsub(1, paste("overlaps ",name.overlap,sep=""), bedCov$overlap)
  bedCov$peaksize = bedCov[,3]-bedCov[,2]
  bedCov$delta = bedCov[,5]-bedCov[,4]
  
  bedCov$max = do.call(pmax, bedCov[4:5])
  bedCov.noverlap = subset(bedCov, bedCov$overlap == paste("not in ",name.overlap,sep=""))
  bedCov.overlap = subset(bedCov, bedCov$overlap == paste("overlaps ",name.overlap,sep=""))
  bedout = bedCov[,c("V1","V2","V3","overlap","V4","V5","delta","log2fc")]; write.table(x = bedout,file = output.file.bed, col.names = F, row.names = F, quote = F,sep="\t")
  bedout.noverlap = bedCov.noverlap[,c(1,2,3,6)]; write.table(x = bedout.noverlap,file = output.file.bednoverlap, col.names = F, row.names = F, quote = F,sep="\t")
  bedout.overlap  = bedCov.overlap[,c(1,2,3,6)]; write.table(x = bedout.overlap,file = output.file.bedoverlap, col.names = F, row.names = F, quote = F,sep="\t")
  
##4. Plot

  # sample stats
  library(ggplot2)
  peaktext= paste(dir.output,"\n",paste(nrow(bedCov.noverlap)," peaks not in ",name.overlap,sep=""),"\n",paste(nrow(bedCov.overlap)," peaks overlap ",name.overlap,sep=""),sep="")
  peaksizeplot = ggplot(bedCov, aes(x=as.character(overlap),y=log10(peaksize), fill=as.character(overlap),color=as.character(overlap)))+geom_boxplot(alpha=0.5)+
    theme_bw()+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+ scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    labs(subtitle = peaktext)  

  output.file.peaksize = paste(output.file.prefix,".",name.overlap,".peaksizes.pdf",sep="")
  ggsave(output.file.peaksize, plot = peaksizeplot ,width = 5.5, height = 4)
  
  # MA plots split into overlapping and non-overlapping groups
  #install.packages("LSD")
  library(LSD)
  
  log2bedCovlimit = log2(max(bedCov$max+1))
  
  name.samples.noverlap = paste("Peaks not in ",name.overlap,sep="")
  pdf(output.file.noverlap) 
  comparisonplot(log2(bedCov.noverlap$max+1),bedCov.noverlap$log2fc, colpal="ylgnbu", main = name.samples.noverlap, add.density = TRUE, ... = abline(h = 0.5), xlim = c(2,log2bedCovlimit),ylim = c(-5,5))
  dev.off()
  
  name.samples.overlap = paste("Peaks found in ",name.overlap,sep="")
  pdf(output.file.overlap) 
  comparisonplot(log2(bedCov.overlap$max+1),bedCov.overlap$log2fc, colpal="ylgnbu", main = name.samples.overlap, add.density = TRUE, ... = abline(h = 0.5), xlim = c(2,log2bedCovlimit),ylim = c(-5,5))
  dev.off()
  
  #delta - captures size as well as change.
  #comparisonplot(bedCov$V4,bedCov$delta, colpal="ylgnbu", main = dir.output, add.density = TRUE, ... = abline(h = 0.5))#, xlim = c(2,log2bedCovlimit),ylim = c(-5,5))
  
  
  #box plots for direct comparisons
  ttestpVal = round(t.test(bedCov.noverlap$log2fc, bedCov.overlap$log2fc, paired = F)$p.value, 4)
  l2fcboxplot = ggplot(bedCov, aes(x=overlap,y=log2fc))+geom_violin(aes(fill=overlap),alpha=0.5)+geom_boxplot(aes(fill=overlap),width = 0.5, outlier.shape = NA)+theme_bw()+geom_hline(yintercept = 0)+
    annotate("text",x=1.5,y=2,label=paste("Welch t-test, p = ",ttestpVal)) +annotate("segment", x = 1.1, xend = 1.9, y = 1.8, yend = 1.8)+ggtitle(paste(name.sample1,"v.",name.sample2,sep=" "))+ scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  output.file.l2fcboxplot = paste(output.file.prefix,".",name.overlap,".l2fcboxplot.pdf",sep="")
  ggsave(output.file.l2fcboxplot, plot = l2fcboxplot ,width = 5, height = 5)
  
  ttestpVal.sample1Cov = round(t.test(bedCov.noverlap$V4, bedCov.overlap$V4, paired = F)$p.value, 4)
  logCovboxplot = ggplot(bedCov, aes(x=overlap,y=log2(V4)))+geom_violin(aes(fill=overlap),alpha=0.5)+geom_boxplot(aes(fill=overlap),width = 0.3, outlier.shape = NA)+theme_bw()+coord_flip()+ylab(paste("log2(",name.sample1," RRPM)",sep = ""))+
    annotate("text",x=1.5,y=4,label=paste("Welch t-test, p = ",ttestpVal.sample1Cov)) +annotate("segment", x = 1.1, xend = 1.9, y = 1.8, yend = 1.8)+ggtitle(paste(name.sample1,"v.",name.sample2,sep=" "))+ scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  output.file.logCovboxplot = paste(output.file.prefix,".",name.overlap,".logCovboxplot.pdf",sep="")
  ggsave(output.file.logCovboxplot, plot = logCovboxplot ,width = 7, height = 4)
  
  #box plots for direct comparisons
  ttestpVal.delta = round(t.test(bedCov.noverlap$delta, bedCov.overlap$delta, paired = F)$p.value, 4)
  deltaboxplot = ggplot(bedCov, aes(x=overlap,y=delta))+geom_violin(aes(fill=overlap),alpha=0.5)+geom_boxplot(aes(fill=overlap),width = 0.5, outlier.shape = NA)+theme_bw()+geom_hline(yintercept = 0)+
    annotate("text",x=1.5,y=100,label=paste("Welch t-test, p = ",ttestpVal.delta)) +annotate("segment", x = 1.1, xend = 1.9, y = 1.8, yend = 1.8)+ggtitle(paste(name.sample1,"v.",name.sample2,sep=" "))+ scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  output.file.deltaboxplot = paste(output.file.prefix,".",name.overlap,".deltaboxplot.pdf",sep="")
  ggsave(output.file.deltaboxplot, plot = deltaboxplot ,width = 5, height = 5)
  
  #library(grid)
  #library(gridExtra)
  #grid.arrange(arrangeGrob(peaksizeplot,l2fcboxplot,logCovboxplot,ncol=2))
  
