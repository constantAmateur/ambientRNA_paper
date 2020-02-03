##############
#************#
#* Preamble *#
#************#
##############
setwd('~/trueHome/Projects/SoupX')

#TODO
#Need a convincingish biological example for tumour data.
#Colour scheme for expression in 1B
#Referencing
#Details of how figures were made in methods

#Supplementary Figures:
# Optimum number of UMIs for calculating BG - Done
# Lobster for Drop-seq data - Done
# Show that soup doesn't vary strongly between cells. - Drop
# Showing correlation between effective contamination fraction and inferred one. - Done
# Selection of contamination gene plot. - ?
# Example showing that HB genes work well. - Done

#############
# Libraries #
#############

library(ggplot2)
library(cowplot)
library(reshape2)
library(beeswarm)
#library(cowplot)
library(wesanderson)
#library(reshape2)
library(ComplexHeatmap)
#library(BiocParallel)
library(Seurat)
#library(Matrix)
#library(dplyr)
#library(ggrepel)
#library(DropletUtils)
library(SoupX)
#colSums = Matrix::colSums
#rowSums = Matrix::rowSums
library(Matrix)
import('Code/plotParts.R',all=TRUE)

##########
# Params #
##########

#Where to store the results
plotDir='Results/Virgin'
#Recalculate everything from scratch even if we have the data?
forceRecalc=FALSE
#Manifest with information about Kidney Samples
fp_cMani = 'Data/Tumour/Table5.tsv'
#Path for species-mix experiment
fp_idealHg = 'Data/Mixture/hg19'
fp_idealMm = 'Data/Mixture/mm10'
fp_mixDropSeq = 'Data/DropEstData/dropseq_thousand.rds'
fp_mixInDropHg = 'Data/DropEstData/indrop_mixture_human.rds'
fp_mixInDropMm = 'Data/DropEstData/indrop_mixture_mouse.rds'
#The PBMC data
fp_pbmc = 'Data/PBMC_4k'
#Regex for soup specific genes
rhoEstGenes = '^HB[BGMQDAZE][12]?_'
#Define populations to use to estimate rho.
markers = list(HK = '^(EIF[0-9]|RPL[0-9]|RPS[0-9]|RPN1|POLR[0-9]|SNX[0-9]|HSP[AB][0-9]|H1FX|H2AF[VXYZ]|PRKA|NDUF[ABCSV]|ATP[0-9]|PSM[ABCDEFG][0-9]|UBA[0-9]|UBE[0-9]|USP[0-9]|TXN|SLC[0-9])',
               MT = '^MT-',
               HLA = '^HLA-',
               HEM = '^HB[BGMQDAZE][12]?_'
               )
#Colour palette
colPal = as.character(wes_palette('Darjeeling1',n=5))
colPal = c('#CC3333','#6699CC',colPal[-(1:2)])
#Set1
colPalLarge = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999')
#Set3
#colPalLarge = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9')
#Manual
colPalLarge = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666')
#Diverging colour palette
divColPal = rev(c('#8c510a','#d8b365','#f6e8c3','#f5f5f5','#c7eae5','#5ab4ac','#01665e'))
divColPal = rev(c('#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4'))
rasterRes = 300
#Cut-off for when a cluster is considered to map to its cleaned partner 1-1
matchFrac = 0.8
#How many cores to use
nCores = 16

#############
# Functions #
#############
#' Run standard analysis
standard10X = function(dat,nPCs=50,res=1.0){
  srat = CreateSeuratObject(dat)
  srat = NormalizeData(srat)
  srat = ScaleData(srat)
  srat = FindVariableFeatures(srat)
  srat = RunPCA(srat)
  srat = RunUMAP(srat,dims=seq(nPCs))
  srat = RunTSNE(srat,dims=seq(nPCs))
  srat = FindNeighbors(srat,dims=seq(nPCs))
  srat = FindClusters(srat,res=res)
  return(srat)
}
#' Normal plots
#'
#' The idea being that this will save a pdf and png version.  If splitRaster is set to TRUE, it will also save a two-part figure where part of the plot is rasterised and saved as a png and the other bits are saved as a pdf.
#' All this works by assuming that the doPlot function takes one parameter that can be 'All', 'Text', or 'Raster'
makePlots = function(baseNom,doPlot,splitRaster=FALSE,...){
  #Normal plots
  pdf(paste0(baseNom,'.pdf'),...,useDingbats=FALSE)
    doPlot()
  dev.off()
  png(paste0(baseNom,'.png'),...,res=rasterRes,units='in')
    doPlot()
  dev.off()
  #Split raster and text
  if(splitRaster){
    pdf(paste0(baseNom,'_nonRaster.pdf'),...,useDingbats=FALSE)
      doPlot('Text')
    dev.off()
    png(paste0(baseNom,'_raster.png'),...,res=rasterRes,units='in')
      doPlot('Raster')
    dev.off()
  }
}


#############
# Schematic #
#############
#Use PBMC data
scPBMC = load10X(fp_pbmc)
##########################
# Soup profile estimation
#Elbow plot
x = scPBMC$nDropUMI
o = order(-x)
x = x[o]
pdf(file.path('Results','Schematic_fitBG1.pdf'))
  plot(seq_along(x),x,
       log='xy',
       type='l',
       frame.plot=FALSE,
       ylab='nUMIs',
       xlab='Droplets')
  lims = range(which(x<10))
  polygon(c(lims[1],lims[1],lims[2],lims[2]),
          c(1,max(x),max(x),1),
          col='#11111133')
  text(lims[1]/2,mean(x),labels=sprintf('%d droplets with %d counts',sum(x<10),sum(x[x<10])))
dev.off()
#Background barplot
x = scPBMC$soupProfile
x = x[order(-x$est),]
#Top and bottom 5
xx = x[c(1:5,tail(which(x$est>0),5)),]
#Percentage, one sig-fig
xx$perc = sprintf('%.01g',xx$est*100)
pdf(file.path('Results','Schematic_fitBG2.pdf'))
  y = barplot(xx$est,names.arg=rownames(xx),las=2)
  text(y[,1],max(xx$est)*1.1,xx$perc,xpd=TRUE,srt=90)
dev.off()
###########################
# Estimating contamination
scPBMC$metaData = cbind(scPBMC$metaData,PBMC_DR[rownames(scPBMC$metaData),])
colnames(scPBMC$metaData) = gsub('Cluster','cluster',colnames(scPBMC$metaData))
igGenes = c('IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHD','IGHE','IGHM',
            'IGLC1','IGLC2','IGLC3','IGLC4','IGLC5','IGLC6','IGLC7',
            'IGKC')
useToEst = estimateNonExpressingCells(scPBMC,list(IG=igGenes),setNames(PBMC_DR$Cluster,rownames(PBMC_DR)))
gg = plotMarkerMap(scPBMC,igGenes,PBMC_DR,useToEst=useToEst)
x = plotMarkerMap(scPBMC,igGenes,PBMC_DR)$data
labs = x[x$qVals<0.05,]
labs = labs[match(unique(labs$Cluster),labs$Cluster),]
mids = aggregate(cbind(RD1,RD2) ~ Cluster,data=x,FUN=mean)
gg1 = ggplot(x,aes(RD1,RD2,colour=qVals<0.05)) +
  #Plot the non-sig ones first
  geom_point(data=x[x$qVals>=0.05,],size=0.5) +
  geom_point(data=x[x$qVals<0.05,],size=0.5) +
  #Label cluster mids
  geom_label(data=mids,inherit.aes=FALSE,aes(RD1,RD2,label=paste0('Clust',Cluster))) +
  #Label one of the bad ones randomly
  geom_label(data=labs,aes(label=Cluster),alpha=0.4)+
  geom_density2d(inherit.aes=FALSE,aes(RD1,RD2,group=Cluster),colour='black') +
  scale_colour_manual(values=c(`FALSE`='#888888',`TRUE`='#cc0000')) +
  guides(colour=FALSE) +
  theme_void()
pdf(file.path('Results','Schematic_rho1.pdf'),useDingbats=FALSE)
  plot(gg1)
dev.off()
tmp = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')
tmp = c('#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#000000')
tmp = c('#f7f7f7','#cccccc','#969696','#525252')
x$igGenes = colSums(scPBMC$toc[igGenes,])[rownames(x)]
x$useToEst = useToEst[rownames(x),1]
x[!useToEst[,1],'igGenes'] = 0
doPlot = function(e){
  plot(x$RD1[x$useToEst],x$RD2[x$useToEst],
       type='n',
       frame.plot=FALSE,
       xlab='',
       ylab='',
       axes=FALSE)
  points(x$RD1[x$useToEst],x$RD2[x$useToEst],
         pch=19,
         #col='black',
         col=colAlpha(tmp[x$igGenes[x$useToEst]+1],1.0))
}
gg2 = ggplot(x,aes(RD1,RD2,colour=as.character(igGenes))) + 
  geom_point(data=x[!x$useToEst,],colour='grey',alpha=0.5) +
  geom_point(data=x[x$useToEst,]) +
  ggtitle(sprintf('%d IG counts, %g expected soup',
                  sum(scPBMC$toc[igGenes,useToEst[,1]]),
                  sum(scPBMC$toc[,useToEst[,1]])*sum(scPBMC$soupProfile[igGenes,'est'])
                  )) +
  scale_colour_manual(values=tmp[seq(7)]) +
  guides(colour=FALSE) +
  theme_void()
pdf(file.path('Results','Schematic_rho2.pdf'),useDingbats=FALSE)
  #plot(gg2)
  doPlot()
dev.off()
####################
# Correcting counts
scPBMC = calculateContaminationFraction(scPBMC,useToEst=useToEst,nonExpressedGeneList=list(IG=igGenes))
out = adjustCounts(scPBMC,setNames(PBMC_DR$Cluster,rownames(PBMC_DR)),verbose=2,nCores=nCores)
#Now get the genes to display and total for cluster 6
toDisp = c('LYZ','IGKC','S100A8','S100A9','CD74')
tgtClust = '6'
x = rowSums((scPBMC$toc-out)[,scPBMC$metaData$cluster==tgtClust])
#Round to integers and fake y-distribution
cnts = round(x[toDisp])
#Add one for total
cnts = c(cnts,round(1000))
df = data.frame(gene = rep(seq_along(c(toDisp,'Total')),cnts)+rnorm(sum(cnts),sd=0.15),
                yDist = 0.5+rnorm(sum(cnts),sd=0.15))
df$gene = as.character(df$gene)
pdf(file.path(plotDir,'Schematic_adjust1.pdf'),height=4)
  plot(df$gene,df$yDist,
       ylim=c(-1,2),
       pch=19,
       cex=0.2)
dev.off()
#Make a matrix with a random 10 cells for each gene
spCnts = (scPBMC$toc-out)[toDisp,scPBMC$metaData$cluster==tgtClust]
x = t(as.matrix(round(spCnts[,sample(ncol(spCnts),10)])))
y = t(as.matrix(scPBMC$toc[colnames(x),rownames(x)]))
rownames(x) = paste0('Cell',c(1:5,(ncol(spCnts)-4):ncol(spCnts)))





###########################
# Species mix experiments #
###########################
#Load each of the different technologies
scMixes = list()
##########
# DropSeq
ideal = readMM('Data/DropSeqMixture/SRR1748411/quants_mat.mtx.gz')
ideal = t(ideal)
ideal@x = round(ideal@x)
rownames(ideal) = read.table('Data/DropSeqMixture/SRR1748411/quants_mat_cols.txt')[,1]
colnames(ideal) = read.table('Data/DropSeqMixture/SRR1748411/quants_mat_rows.txt')[,1]
hgGenes= grep('ENSG',rownames(ideal),value=TRUE)
mmGenes = grep('ENSM',rownames(ideal),value=TRUE)
hgCnts = colSums(ideal[hgGenes,])
mmCnts = colSums(ideal[mmGenes,])
#Designate cells, excluding doublets
cells = colnames(ideal)[colSums(ideal) >= 5000 & colSums(ideal>0) >=300 & pmin(mmCnts,hgCnts)/colSums(ideal)<0.2 ]
scMix = SoupChannel(ideal,ideal[,cells],dataType='DropSeq',keepDroplets=TRUE,soupRange=c(10,100))
scMix$hgCnts = colSums(scMix$tod[hgGenes,])
scMix$mmCnts = colSums(scMix$tod[mmGenes,])
isHuman = scMix$hgCnts[colnames(scMix$toc)] > scMix$mmCnts[colnames(scMix$toc)]
isDoublet = pmin(scMix$hgCnts[colnames(scMix$toc)], scMix$mmCnts[colnames(scMix$toc)])>=1000
scMix$isHuman = isHuman
scMix$isDoublet = isDoublet
scMix$hgGenes = hgGenes
scMix$mmGenes = mmGenes
scMix$seurat = standard10X(scMix$toc) 
scMix = setClusters(scMix,scMix$seurat@active.ident)
scMixes[['DropSeq']]=scMix
######
# 10X
hg = Read10X(fp_idealHg)
mm = Read10X(fp_idealMm)
#Combine
ideal = rbind(hg,mm)
#Use cut-off to select cells
cells = colnames(ideal)[colSums(ideal) >= 1000 & colSums(ideal>0) >= 500]
#Now construct the soupchannels
scMix = SoupChannel(ideal,ideal[,cells],dataType='10X',keepDroplets=TRUE)
hgGenes = grep('^hg19_',rownames(scMix$tod),value=TRUE)
mmGenes = grep('^mm10_',rownames(scMix$tod),value=TRUE)
#Store useful bits of information
scMix$hgCnts = colSums(scMix$tod[hgGenes,])
scMix$mmCnts = colSums(scMix$tod[mmGenes,])
#Calculate the contamination
isHuman = scMix$hgCnts[colnames(scMix$toc)] > scMix$mmCnts[colnames(scMix$toc)]
isDoublet = pmin(scMix$hgCnts[colnames(scMix$toc)], scMix$mmCnts[colnames(scMix$toc)])>=1000
scMix$isHuman = isHuman
scMix$isDoublet = isDoublet
scMix$hgGenes = hgGenes
scMix$mmGenes = mmGenes
scMix$seurat = standard10X(scMix$toc) 
scMix = setClusters(scMix,scMix$seurat@active.ident)
scMixes[['10X']]=scMix
########################
# Calculate rho for all
for(tNom in names(scMixes)){
  scMix = scMixes[[tNom]]
  message(sprintf("Processing for technology %s",tNom))
  #Define the genes to use for estimation
  geneSets = list(hg=scMix$hgGenes,
                  mm=scMix$mmGenes)
  #We could do this manually, but the automated thing gets it 100% right
  useToEst = estimateNonExpressingCells(scMix,nonExpressedGeneList=geneSets,clusters=FALSE)
  #Calculate a global rho
  scMix$globalRho = calculateContaminationFraction(scMix,nonExpressedGeneList=geneSets,useToEst=useToEst)$metaData$rho[1]
  #And now the fine-grained estimation
  scMix = calculateContaminationFraction(scMix,nonExpressedGeneList=geneSets,useToEst=useToEst,verbose=TRUE,cellSpecificEstimates=TRUE,stanParams= list(chains = 1, warmup = 8000, iter = 48000, cores = 1))
  scMixes[[tNom]] = scMix
}
################
# Lobster plot
#Make one version of this for all techs
for(tNom in names(scMixes)){
  scMix = scMixes[[tNom]]
  df = data.frame(hg=scMix$hgCnts,mm=scMix$mmCnts,isCell=colnames(scMix$tod) %in% colnames(scMix$toc))
  df$isDoublet = FALSE
  df[rownames(scMix$metaData),'isDoublet'] = scMix$isDoublet[rownames(scMix$metaData)]
  #df = data.frame(hg=hg_umiCounts,mm=mm_umiCounts,qval=qval)
  df$M = log10(0.5*(df$hg+df$mm))
  df$A = log10(df$hg/df$mm)
  ext = ceiling(1.01*max(abs(df$A[is.finite(df$A)])))
  df$A[df$A == -Inf] = -ext
  df$A[df$A == Inf] = ext
  #Drop the doublets
  df = df[!df$isDoublet,]
  #Estimate global rho for two populations
  rho_hg = log10(mean(scMix$metaData$rho[!scMix$isDoublet & scMix$isHuman])/2)
  rho_mm = log10(mean(scMix$metaData$rho[!scMix$isDoublet & !scMix$isHuman])/2)
  #Make the plot
  doPlot = function(mode=c('All','Raster','Text')){
    mode = match.arg(mode)
    if(mode=='Raster'){
      plot(df$M,df$A,
           type='n',
           xaxt='n',
           yaxt='n',
           frame.plot=FALSE,
           xlab='',
           ylab=''
           )
    }else{
      plot(df$M,df$A,
           type='n',
           frame.plot=FALSE,
           xlab='log10(Average UMIs)',
           ylab='log10(Human/Mouse)'
           )
      legend('topright',
             legend=c('Yes','No'),
             pch=19,
             col=c('black',grey(0.7)),
             title='Droplet contains cell?',
             bty='n',
             xpd=TRUE)
    }
    if(mode!='Text'){
      points(df$M,df$A,
           cex=0.1,
           pch=19,
           col=ifelse(df$isCell,'black',grey(0.7))
           )
    }
  }
  makePlots(file.path(plotDir,sprintf('lobster_%s',tNom)),doPlot,width=3.5,height=3.5,pointsize=8,splitRaster=TRUE)
  ###gg_lobster = ggplot(df,aes(M,A,colour=isCell)) +
  ###  geom_point(size=1.0) +
  ###  #geom_hline(yintercept=c(-rho_hg,rho_mm)) +
  ###  labs(colour='Is a cell?') +
  ###  xlab('log10(Average UMIs)') +
  ###  ylab('log10(Human/Mouse)') +
  ###  scale_colour_manual(values=colPal[1:2]) +
  ###  guides(colour=FALSE) + 
  ###  theme_bw(base_size=18)
  ###pdf(file.path(plotDir,sprintf('lobster_%s.pdf',tNom)))
  ###plot(gg_lobster)
  ###dev.off()
  ###png(file.path(plotDir,sprintf('lobster_%s.png',tNom)),width=7,height=7,units='in',res=rasterRes)
  ###plot(gg_lobster)
  ###dev.off()
}
########################
# Soup/Cell correlation 
#Get just the human counts in human cells and mouse counts in mouse cells
for(tNom in names(scMixes)){
  scMix = scMixes[[tNom]]
  spFrac = scMix$soupProfile$counts
  clFrac = rowSums(scMix$toc)
  #Sub-sample to match
  rat = sum(spFrac)/sum(clFrac)
  clFrac = rbinom(length(clFrac),clFrac,rat)
  doPlot = function(mode=c('All','Raster','Text')){
    mode = match.arg(mode)
    x=spFrac
    y=clFrac
    plot(x,
         y,
         log='xy',
         type='n',
         las=1,
         xlab='Counts in background',
         ylab='Counts in cells',
         frame.plot=FALSE
         )
    if(mode!='Text'){
      abline(0,1,col='black')
      points(x,
             y,
             cex=0.1,
             pch=19,
             #col=ifelse(df$species=='human',colPal[1],colPal[2])
             col='black'
             )
    }
    cc=cor(x,y)
    text(1,1000,bquote(R^2 == .(format(cc,digits=2))),adj=c(0,0.5))
    #text(-5,-2.3,label=labs[2])
  }
  makePlots(file.path(plotDir,sprintf('idealCorrelation_%s',tNom)),doPlot,width=3.5,height=3.5,pointsize=8,splitRaster=FALSE)
  ###dat = scMix$toc
  ###isHuman = scMix$hgCnts[colnames(dat)] > scMix$mmCnts[colnames(dat)]
  ###isDoublet = scMix$isDoublet[colnames(dat)]
  ###hgSoup = scMix$soupProfile[scMix$hgGenes,'est']
  ###mmSoup = scMix$soupProfile[scMix$mmGenes,'est']
  ####dat = ideal[,idealCells]
  ####isHuman = hgCnts[idealCells]>mmCnts[idealCells]
  ####isDoublet = pmin(hgCnts[idealCells],mmCnts[idealCells])>=1000
  ###hgCell = (dat[scMix$hgGenes,isHuman & !isDoublet])
  ###hgCell = Matrix::rowSums(hgCell)/sum(hgCell)
  ####hgSoup = idealSoup$est[grep('hg',rownames(idealSoup))]
  ###hgSoup = hgSoup/sum(hgSoup)
  ###mmCell = (dat[scMix$mmGenes,!isHuman & !isDoublet])
  ###mmCell = Matrix::rowSums(mmCell)/sum(mmCell)
  ####mmSoup = idealSoup$est[grep('mm',rownames(idealSoup))]
  ###mmSoup = mmSoup/sum(mmSoup)
  ###df = data.frame(cellFrac = c(hgCell,mmCell),
  ###                soupFrac = c(hgSoup,mmSoup),
  ###                species = rep(c('human','mouse'),c(length(hgCell),length(mmCell)))
  ###                )
  ####Drop the top expressed genes as they will skew correlation and likely are biological 
  ###w = which(df$cellFrac < quantile(df$cellFrac,0.99) & df$soupFrac < quantile(df$soupFrac,0.99))
  ###tmp = df[w,]
  ###hgCor = cor(tmp$cellFrac[tmp$species=='human'],tmp$soupFrac[tmp$species=='human'])
  ###mmCor = cor(tmp$cellFrac[tmp$species=='mouse'],tmp$soupFrac[tmp$species=='mouse'])
  ###labs = c(sprintf('R^2 == %.02f (Human)',hgCor),
  ###         sprintf('R^2 == %.02f (Mouse)',mmCor))
  ####Alterantive total correlation
  ###totCor = cor(tmp$cellFrac,tmp$soupFrac)
  ###labs = sprintf('R^2 = %.02f',totCor)
  ####Just a straight correlation plot
  ###doPlot = function(textOnly=FALSE){
  ###  x=log10(df$cellFrac)
  ###  y=log10(df$soupFrac)
  ###  plot(x,
  ###       y,
  ###       type='n',
  ###       las=1,
  ###       xlab='Average expression in cells',
  ###       ylab='Average expression in background',
  ###       frame.plot=FALSE
  ###       )
  ###  if(!textOnly){
  ###    abline(0,1,col='black')
  ###    points(x,
  ###           y,
  ###           cex=0.1,
  ###           pch=19,
  ###           #col=ifelse(df$species=='human',colPal[1],colPal[2])
  ###           col='black'
  ###           )
  ###  }
  ###  text(-5,-2,labs[1])
  ###  #text(-5,-2.3,label=labs[2])
  ###}
  ###pdf(file.path(plotDir,sprintf('idealCorrelation_%s_textOnly.pdf',tNom)),width=3.5,height=3.5,pointsize=8,useDingbats=FALSE)
  ###  doPlot(TRUE)
  ###dev.off()
  ###pdf(file.path(plotDir,sprintf('idealCorrelation_%s.pdf',tNom)),width=3.5,height=3.5,pointsize=8,useDingbats=FALSE)
  ###  doPlot()
  ###dev.off()
  ###png(file.path(plotDir,sprintf('idealCorrelation_%s.png',tNom)),width=3.5,height=3.5,units='in',res=rasterRes,pointsize=8)
  ###  doPlot()
  ###dev.off()
  ###gg_idealCor = ggplot(df,aes(log10(cellFrac),log10(soupFrac))) +
  ###  geom_point(aes(colour=species),size=1.0,alpha=1) +
  ###  geom_abline(intercept=0,slope=1) +
  ###  annotate('text',x=-5,y=-2,label= labs[1],parse=TRUE) +
  ###  annotate('text',x=-5,y=-2.3,label= labs[2],parse=TRUE) +
  ###  xlab('log10(cell expression)') +
  ###  ylab('log10(soup expression)') +
  ###  scale_colour_manual(values=colPal[1:2]) + 
  ###  guides(colour=FALSE) +
  ###  theme_bw(base_size=18)
  ####The MA plot instead 
  ####gg_idealCor = ggplot(df,aes(log10(0.5*(cellFrac+soupFrac)),log10(cellFrac/soupFrac))) +
  ####  geom_point(aes(colour=species)) +
  ####  geom_hline(yintercept=0) +
  ####  annotate('text',x=-2,y=-2,label= labs[1],parse=TRUE) +
  ####  annotate('text',x=-2,y=-2.3,label= labs[2],parse=TRUE) +
  ####  xlab('log10(Average expression)') +
  ####  ylab('log10(cell/soup)') +
  ####  xlim(-6,NA) +
  ####  guides(colour=FALSE) +
  ####  theme_bw(base_size=18)
  ####  #ggtitle(sprintf('Gene correlation between cell and soup\n R=%g,%g for human,mouse',cor(hgCell,hgSoup),cor(mmCell,mmSoup)))
  ###pdf(file.path(plotDir,sprintf('idealCorrelation_%s.pdf',tNom)))
  ###plot(gg_idealCor)
  ###dev.off()
  ###png(file.path(plotDir,sprintf('idealCorrelation_%s.png',tNom)),width=7,height=7,units='in',res=rasterRes)
  ###plot(gg_idealCor)
  ###dev.off()
}
####################
# Cell specific rho
#Make one plot that contains everything
df = lapply(names(scMixes),function(e) {
              tmp = scMixes[[e]]$metaData
              tmp$isDoublet = scMixes[[e]]$isDoublet[rownames(tmp)]
              tmp$isHuman = scMixes[[e]]$isHuman[rownames(tmp)]
              tmp$techName = e
              tmp})
df = do.call(rbind,df)
df = df[!df$isDoublet,]
doPlot = function(mode=c('All','Raster','Text')){
  mode = match.arg(mode)
  #Transparency code
  alpha=0.4
  hexAlpha = toupper(as.hexmode(round(alpha*255)))
  #Define layout
  zones=matrix(c(1,2), ncol=2, byrow=TRUE)
  layout(zones, widths=c(1/10,9/10), heights=c(1))
  #Make density plots by tech
  s = lapply(split(df$rho,df$techName),density)
  maxDen = 100
  maxRho = 0.2
  tmp = density(df$rho)
  par(mar=c(5,0,4,0))
  plot(0,0,type='n',
       xlim=c(-maxDen,0),
       ylim=c(0,maxRho),
       frame.plot=FALSE,
       xaxt='n',
       yaxt='n',
       xlab='Frequency',
       ylab='')
  for(i in seq_along(s)){
    lines(-s[[i]]$y*maxDen/max(s[[i]]$y),s[[i]]$x,col=colPal[i])
    polygon(c(0,-s[[i]]$y*maxDen/max(s[[i]]$y),0),c(0,s[[i]]$x,max(s[[i]]$x)),col=paste0(colPal[i],hexAlpha))
  }
  #And the plots of rho
  par(mar=c(5,4,4,2))
  plot(0,0,type='n',
       xlim=range(log10(df$nUMIs)),
       ylim=c(0,maxRho),
       las=1,
       frame.plot=FALSE,
       xlab='log10(#UMIs)',
       ylab='Contamination Fraction')
  tNoms = names(s)
  for(i in seq_along(tNoms)){
    tNom=tNoms[i]
    w = which(df$techName==tNom)
    x = log10(df$nUMIs)[w]
    y = df$rho[w]
    points(x,y,
           pch=19,
           col=paste0(colPal[i],hexAlpha),
           cex=0.3)
    sm = smooth.spline(x,y,nknots=5)
    lines(sm$x,sm$y,col=colPal[i])
  }
}
makePlots(file.path(plotDir,'rhoMixTrend'),doPlot,width=4.5,height=3.5,pointsize=8,splitRaster=FALSE)
#gg = ggplot(df,aes(x=log10(nUMIs),rho,colour=techName)) + 
#  geom_point(size=0.5) +
#  geom_smooth() +
#  scale_colour_manual(values=colPal[seq_len(length(scMixes))]) + 
#  xlab('log10(#UMIs)') +
#  ylab('Contamination Fraction')
##Add marginal distribution on left
#yPlot = ggplot(df,aes(rho,fill=techName)) + 
#  geom_density(alpha=0.5) +
#  coord_flip() +
#  scale_y_reverse() +
#  scale_fill_manual(values=colPal[seq_len(length(scMixes))]) + 
#  guides(colour=FALSE,fill=FALSE) +
#  theme_nothing()
#ggp = plot_grid(yPlot,gg,align='v',axis='r',nrow=1,ncol=2,rel_widths=c(1,4))
###pdf(file.path(plotDir,'rhoMixTrend.pdf'),width=4.5,height=3.5,useDingbats=FALSE,pointsize=8)
###  doPlot()
###dev.off()
###png(file.path(plotDir,'rhoMixTrend.png'),width=4.5,height=3.5,units='in',res=rasterRes,pointsize=8)
###  doPlot()
###dev.off()
#####################
# Species retention 
#Get the corrected counts
df = list()
for(tNom in names(scMixes)){
  scMix = scMixes[[tNom]]
  #Use global estimate
  scMix = setContaminationFraction(scMix,scMix$globalRho)
  out = adjustCounts(scMix,verbose=2,nCores=nCores)
  tmp = data.frame(mFracNew = colSums(out[scMix$mmGenes,])/colSums(out),
                   mFracOld = colSums(scMix$toc[scMix$mmGenes,])/scMix$metaData[colnames(out),'nUMIs'],
                   hFracNew = colSums(out[scMix$hgGenes,])/colSums(out),
                   hFracOld = colSums(scMix$toc[scMix$hgGenes,])/scMix$metaData[colnames(out),'nUMIs'],
                   isHuman = scMix$isHuman[colnames(out)],
                   isDoublet = scMix$isDoublet[colnames(out)])
  #tmp$isHuman = tmp$hFracOld>0.5
  #tmp$isDoublet = FALSE
  tmp$mRatio = tmp$mFracNew/tmp$mFracOld
  tmp$hRatio = tmp$hFracNew/tmp$hFracOld
  tmp = melt(tmp[!tmp$isDoublet,],id.vars='isHuman',measure.vars=c('mRatio','hRatio'))
  tmp$crossSpecies = xor(tmp$isHuman,grepl('h',tmp$variable))
  tmp$techName = tNom
  df[[tNom]] = tmp
}
df = do.call(rbind,df)
doPlot = function(mode=c('All','Raster','Text')){
  mode = match.arg(mode)
  #Block left 0-2 (with buffer)
    #Block inner-left 0.2-1 (with buffer)
      #0.3 - 0.9
    #Block inner-right 1-1.8 (with buffer)
      #1.1 - 1.7
  #Block right 2-4 (with buffer)
    #
  plot(0,0,
       xlim = c(0,4),
       ylim = c(0,1.2),
       las=1,
       type='n',
       xaxt='n', 
       xlab='',
       ylab='Fractional Change',
       frame.plot=FALSE
       )
  axis(1,at=c(1,3),labels=c('True\nCounts','Contaminating\nCounts'),padj=0.5)
  abline(h=1)
  #Within species
  #10X
  sets = list(list(cs=FALSE,tn='10X',col=colPal[1],mid=0.6),
              list(cs=FALSE,tn='DropSeq',col=colPal[2],mid=1.4),
              list(cs=TRUE,tn='10X',col=colPal[1],mid=2.6),
              list(cs=TRUE,tn='DropSeq',col=colPal[2],mid=3.4))
  for(set in sets){
    x = df$value[df$crossSpecies==set$cs & df$techName==set$tn]
    points(set$mid+ runif(length(x),-.25,.25),
           x,
           col=colAlpha(set$col,0.1),
           pch=19,
           cex=0.1
    )
    boxplot(x,
            col=NA,
            at=set$mid,
            boxwex=1.0,
            border=set$col,
            outline=FALSE,
            add=TRUE,
            xaxt='n', 
            yaxt='n', 
            bty='n',
            frame.plot=FALSE
    )
  }
  legend('topright',
         legend=c('10X','DropSeq'),
         fill=colPal,
         bty='n',
         title='Experiment',
         xpd=TRUE)
}
makePlots(file.path(plotDir,'speciesMixPerformance'),doPlot,width=3,height=4,pointsize=8,splitRaster=FALSE)
#gg = ggplot(df,aes(crossSpecies,value,colour=techName)) +
#  geom_hline(yintercept=1.0,colour='black') +
#  geom_point(position=position_jitterdodge(jitter.width=0.25),size=0.1,alpha=0.4) +
#  geom_boxplot(outlier.colour=NA,alpha=0.5) + 
#  xlab('Cross species') +
#  ylab('Fraction retained') +
#  scale_colour_manual(values=rep(colPal[seq_len(length(scMixes))],2)) +
#  theme_classic() +
#  guides(colour=FALSE)
#pdf(file.path(plotDir,'speciesMixPerformance.pdf'),width=2.5,height=3.5,pointsize=8,useDingbats=FALSE)
#  doPlot()
#dev.off()
#png(file.path(plotDir,'speciesMixPerformance.png'),width=2.5,height=3.5,units='in',res=rasterRes,pointsize=8)
#  doPlot()
#dev.off()
###########################################
# Compare effective rho with estimated one
df = list()
for(tNom in names(scMixes)){
  scMix = scMixes[[tNom]]
  trueRho = scMix$metaData$rho
  out = adjustCounts(scMix,verbose=2,nCores=nCores)
  effRho = colSums(scMix$toc-out)/colSums(scMix$toc)
  df[[tNom]] = data.frame(nUMIs = scMix$metaData$nUMIs,
                          trueRho = trueRho,
                          effRho = effRho,
                          tNom = tNom)
}
df = do.call(rbind,df)
doPlot = function(mode=c('All','Raster','Text')){
  tmp = c(df$trueRho,df$effRho)
  plot(0,0,
       type='n',
       frame.plot=FALSE,
       xlim=range(tmp),
       ylim=range(tmp),
       xlab='True contamination',
       ylab='Effective contamination',
       log='xy')
  abline(0,1,col='black')
  points(df$trueRho,df$effRho,
         pch=19,
         cex=0.3,
         col=ifelse(df$tNom=='DropSeq',colPal[1],colPal[2])
         )
}
makePlots(file.path(plotDir,'speciesMixEffectiveRho'),doPlot,width=7,height=7,pointsize=8,splitRaster=FALSE)
###xx = log10(scMix$metaData$nUMIs)
###plot(xx,trueRho,pch=19,cex=0.1)
###points(xx,effRho,pch=19,cex=0.1,col='black')
###segments(xx,trueRho,xx,effRho,col=ifelse(trueRho>effRho,'red','green'))
###abline(h=scMix$globalRho)
###plot(xx,effRho,pch=19,cex=0.1)
####points(xx,effRho,pch=19,cex=0.1,col='black')
###segments(xx,rep(scMix$globalRho,length(xx)),xx,effRho,col=ifelse(effRho>scMix$globalRho,'green','red'))
###abline(h=scMix$globalRho)
###################################################
#### Comparison of different methods of correction
###noClustFit = adjustCounts(sc,FALSE,'multinomial',verbose=2)
###subFit = adjustCounts(sc,clusters,'subtraction',verbose=2)
###multFit = adjustCounts(sc,clusters,'multinomial',verbose=2)
###tt = adjustCounts(sc,FALSE,'multinomial',verbose=2)
###pFit = adjustCounts(sc,clusters,'soupOnly',verbose=2)
####Check the effective rho from cluster expansion to that from direct estimation
###directRho = sc$metaData$rho
###effectiveRho = colSums(sc$toc-multFit)/colSums(sc$toc)
####Compare the estimation from multinomial to subtraction
###a = as.matrix(sc$toc-subFit)
###b = as.matrix(sc$toc-multFit)
###diff = sapply(seq(ncol(subFit)),function(e) cor(a[,e],b[,e]))
####Fraction of cross-speciesness
###df = list()
###tgtFits = list(subFit=subFit,multFit=multFit,pFit=pFit)
###for(tNom in names(tgtFits)){
###  tgtFit = tgtFits[[tNom]]
###  tmp = data.frame(mFracNew = colSums(tgtFit[sc$mmGenes,])/colSums(tgtFit),
###                   mFracOld = colSums(sc$toc[sc$mmGenes,])/colSums(sc$toc),
###                   hFracNew = colSums(tgtFit[sc$hgGenes,])/colSums(tgtFit),
###                   hFracOld = colSums(sc$toc[sc$hgGenes,])/colSums(sc$toc),
###                   isHuman = sc$isHuman[colnames(tgtFit)],
###                   isDoublet = sc$isDoublet[colnames(tgtFit)],
###                   cellID = colnames(tgtFit)
###  )
###  #tmp$isHuman = tmp$hFracOld>0.5
###  #tmp$isDoublet = FALSE
###  tmp$mRatio = tmp$mFracNew/tmp$mFracOld
###  tmp$hRatio = tmp$hFracNew/tmp$hFracOld
###  tmp = melt(tmp[!tmp$isDoublet,],id.vars=c('cellID','isHuman'),measure.vars=c('mRatio','hRatio'))
###  tmp$crossSpecies = xor(tmp$isHuman,grepl('h',tmp$variable))
###  tmp$techName = tNom
###  df[[tNom]] = tmp
###}
###df = do.call(rbind,df)
###df$clusterID = bu_clusters[df$cellID]
###df$nUMIs = cut(sc$metaData[df$cellID,'nUMIs'],10)
###gg = ggplot(df,aes(crossSpecies,value,colour=techName)) +
###  geom_point(aes(shape=isHuman),position=position_jitterdodge(jitter.width=0.1)) +
###  geom_boxplot(aes(shape=isHuman),outlier.colour=NA,alpha=0.5) #+ 
###  #scale_colour_manual(values=rep(colPal[seq_len(length(tgtFits))],2))
############################################
#### Optimal droplets for soup calculation
####Define "true background".  Cross-species counts in cells
###scMix = scMixes[['10X']]
###hgBG = rowSums(scMix$toc[scMix$mmGenes,scMix$isHuman & !scMix$isDoublet & scMix$metaData$nUMIs > 1000])
###hgBG = hgBG/sum(hgBG)
###mmBG = rowSums(scMix$toc[scMix$hgGenes,!scMix$isHuman & !scMix$isDoublet & scMix$metaData$nUMIs > 1000])
###mmBG = mmBG/sum(mmBG)
###out = list()
###for(i in seq(100)){
###  w = which(scMix$nDropUMIs==i)
###  hgSoup = rowSums(scMix$tod[scMix$mmGenes,w])
###  hgSoup = hgSoup/sum(hgSoup)
###  mmSoup = rowSums(scMix$tod[scMix$hgGenes,w])
###  mmSoup = mmSoup/sum(mmSoup)
###  out[[i]] = data.frame(nUMIs=i,
###                        cor = c(cor(hgSoup,hgBG),cor(mmSoup,mmBG),cor(c(hgSoup,mmSoup),c(hgBG,mmBG))),
###                        species = c('human','mouse','both')
###                        )
###}
###out = do.call(rbind,out)
####Keep only the both version
###out = out[out$species=='both',]
###doPlot = function(e){
###  plot(out$nUMIs,out$cor,
###       type='n',
###       xlab='Droplet UMI count',
###       ylab='Correlation with true contamination',
###       las=1,
###       frame.plot=FALSE)
###  points(out$nUMIs,out$cor,
###         pch=19,
###         col='black',
###         cex=0.7)
###}
###makePlots(file.path(plotDir,'speciesMixOptimalRho'),doPlot,width=7,height=7,pointsize=8,splitRaster=FALSE)




#############
# PBMC data #
#############
#Demonstration of what genes to use
scPBMC = load10X(fp_pbmc)
data(PBMC_DR)
scPBMC$metaData = cbind(scPBMC$metaData,PBMC_DR)
scPBMC$DR=c('RD1','RD2')
colnames(scPBMC$metaData) = gsub('Cluster','clusters',colnames(scPBMC$metaData))
###Test the mixed model version
##library(reshape2)
##df = melt(as.matrix(scPBMC$toc))
##df$nUMIs = scPBMC$metaData[df$Var2,'nUMIs']
##df$soup = scPBMC$soupProfile[df$Var1,'est']
##df$exp = df$nUMIs*df$soup
###Drop those rows that have 0 in soup, they contribute nothing to the variance in likelihood with rho
##df = df[df$soup>0,]
###Keep only the top soup genes
##tgts = rownames(scPBMC$soupProfile)[order(-scPBMC$soupProfile$est)]
##tgts = tgts[seq(200)]
###Sub sample
##ddf = df[df$Var2 %in% sample(colnames(scPBMC$toc),100),]
##ddf = ddf[ddf$Var1 %in% tgts,]
##fit = glmer(value ~ 0 + (exp | Var2),data=ddf,family=poisson(link='identity'),start=1,verbose=1)
#####################
# Gene/Cells to use
#Now show which ones we'll use for the two sets
igGenes = c('IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHD','IGHE','IGHM',
            'IGLC1','IGLC2','IGLC3','IGLC4','IGLC5','IGLC6','IGLC7',
            'IGKC')
geneSets = list(LYZ='LYZ',IG=igGenes)
#First plot the clustering
df = scPBMC$metaData
doPlot = function(mode=c('All','Raster','Text')){
  mode = match.arg(mode)
  plot(df$RD1,df$RD2,
       type='n',
       frame.plot=FALSE,
       xlab='tSNE1',
       ylab='tSNE2'
       )
  cols = colPalLarge[seq_along(unique(scPBMC$metaData$clusters)) %% length(colPalLarge) +1]
  cols='black'
  addDensityContours(df$RD1,df$RD2,scPBMC$metaData$clusters,nSplits=20,splitsToUse=2,col=cols)
  points(df$RD1,df$RD2,
         cex=0.3,
         pch=19,
         col=colPalLarge[as.numeric(scPBMC$metaData$clusters) %% length(colPalLarge) +1],
         lwd=0.1,
         #col=grey(0.6)
         )
  mids = aggregate(cbind(RD1,RD2) ~ clusters,data=df,FUN=mean)
  text(mids$RD1,mids$RD2,
       labels=mids$clusters,
       cex=1.5
       )
}
makePlots(file.path(plotDir,'usePBMC_clusters'),doPlot,width=3.5,height=3.5,pointsize=8)
#pdf(file.path(plotDir,'usePBMC_clusters.pdf'))
#doPlot()
#dev.off()
#png(file.path(plotDir,'usePBMC_clusters.png'),width=7,height=7,units='in',res=rasterRes)
#doPlot()
#dev.off()
#Now the usage
for(tgt in names(geneSets)){
  useToEst = estimateNonExpressingCells(scPBMC,geneSets[tgt])
  gg = plotMarkerMap(scPBMC,geneSets[[tgt]],useToEst=useToEst)
  df = gg$data
  df$cluster = scPBMC$metaData[rownames(df),'clusters']
  doPlot = function(mode=c('All','Raster','Text')){
    mode = match.arg(mode)
    plot(df$RD1,df$RD2,
         type='n',
         frame.plot=FALSE,
         xlab='tSNE1',
         ylab='tSNE2'
         )
    #Work out which are the approved clusters
    appClust = unique(df$cluster[df$useToEst])
    clust = factor(df$cluster)
    cols = ifelse(levels(clust) %in% appClust,'black','black')
    lwds = ifelse(levels(clust) %in% appClust,2,1)
    #shadeCols = ifelse(levels(clust) %in% appClust,grey(0.8,0.7),NA)
    shadeCols = ifelse(levels(clust) %in% appClust,NA,NA)
    #shadeCols = ifelse(levels(clust) %in% appClust,colAlpha('darkgreen',0.2),NA)
    addDensityContours(df$RD1,df$RD2,df$cluster,nSplits=20,splitsToUse=2,col=cols,lwd=lwds,shadeCol=shadeCols)
    #Add the ones with NAs as dots
    w = is.na(df$logRatio)
    points(df$RD1[w],df$RD2[w],cex=0.1,pch=19,col='black')
    #Work out colour gradients
    cols = circlize::colorRamp2(seq(-2,2,length.out=length(divColPal)),divColPal)
    points(df$RD1[!w],df$RD2[!w],
           cex=0.7,
           #lwd=0.1,
           pch=19,
           #col = 'black',
           col=cols(df$logRatio[!w]),
           #col=grey(0.6)
           )
    addColorBar('topleft',col=divColPal,breaks=c(-2,2),title='logRatio')
  }
  makePlots(file.path(plotDir,sprintf('usePBMC_%s',tgt)),doPlot,width=3.5,height=3.5,pointsize=8)
}
##############
# Annotation
cMap = c('0'='MNP',
         '1'='CD8 T-Cell',
         '2'='CD8 T-Cell',
         '3'='CD4 T-Cell',
         '4'='B-Cell',
         '5'='CD4 T-Cell',
         '6'='NK',
         '7'='B-Cell',
         '8'='NK',
         '9'='MNP',
         '10'='MNP',
         '11'='?')
PBMC_DR$Annotation = factor(cMap[as.character(PBMC_DR$Cluster)])
mids = lapply(split(PBMC_DR[,1:2],PBMC_DR$Annotation),apply,2,mean)
mids = cbind(as.data.frame(do.call(rbind,mids)),Annotation=names(mids))
mids[1,1:2]= mids[1,1:2]+5
#Don't label the ?
mids = mids[mids$Annotation!='?',]
doPlot = function(mode=c('All','Raster','Text')){
  mode = match.arg(mode)
  df = PBMC_DR
  colnames(df) = gsub('Cluster','clusters',colnames(df))
  plot(df$RD1,df$RD2,
       type='n',
       frame.plot=FALSE,
       xlab='tSNE1',
       ylab='tSNE2'
       )
  #Work out which are the approved clusters
  clust = factor(df$clusters)
  shadeColMap = setNames(colAlpha(colPalLarge[seq_along(unique(df$Annotation))],0.4),levels(factor(df$Annotation)))
  shadeCols = as.numeric(factor(df$Annotation[match(levels(clust),df$clusters)]))
  shadeCols = colAlpha(colPalLarge[shadeCols],0.4)
  shadeCols = shadeColMap[df$Annotation[match(levels(clust),df$clusters)]]
  #shadeCols = ifelse(levels(clust) %in% appClust,colAlpha('darkgreen',0.2),NA)
  addDensityContours(df$RD1,df$RD2,df$clusters,nSplits=20,splitsToUse=2,shadeCol=shadeCols)
  points(df$RD1,df$RD2,
         cex=0.2,
         pch=19,
         col = 'black',
         bg = colAlpha(grey(0.8),0.8)
         )
  #Precise placement
  points(-25,20,bg=shadeColMap['MNP'],lwd=1,cex=2,pch=22)
  text(-25,20,'MNP',adj=-0.2)
  points(-10,-30,bg=shadeColMap['NK'],lwd=1,cex=2,pch=22)
  text(-10,-30,'NK',adj=-0.2)
  points(25,-15,bg=shadeColMap['B-Cell'],lwd=1,cex=2,pch=22)
  text(25,-15,'B-Cell',xpd=NA,adj=-0.2)
  points(0,30,bg=shadeColMap['CD4 T-Cell'],lwd=1,cex=2,pch=22)
  text(0,30,'CD4 T-Cell',adj=-0.2)
  points(20,30,bg=shadeColMap['CD8 T-Cell'],lwd=1,cex=2,pch=22)
  text(20,30,'CD8 T-Cell',adj=-0.2)
  #Legend
  #legend(x='topleft',
  #       legend = levels(factor(df$Annotation)),
  #       #legend=c('B-Cell','CD4 T-Cell','CD8 T-Cell','MNP','NK'),
  #       fill = colAlpha(colPalLarge[seq_along(levels(factor(df$Annotation)))],0.4),
  #       title = 'Cell Type',
  #       bty='n',
  #       xpd=NA)
}
makePlots(file.path(plotDir,'PBMC_annotation'),doPlot,splitRaster=FALSE,width=3.5,height=3.5,pointsize=8)
#pdf(file.path(plotDir,'PBMC_annotation.pdf'))
#doPlot()
#dev.off()
#png(file.path(plotDir,'PBMC_annotation.png'),width=7,height=7,units='in',res=rasterRes)
#doPlot()
#dev.off()
#######################
# Show DE of key genes
useToEst = estimateNonExpressingCells(scPBMC,geneSets['IG'])
scPBMC = calculateContaminationFraction(scPBMC,geneSets['IG'],useToEst)
out = adjustCounts(scPBMC,verbose=2,nCores=nCores)
bad = standard10X(scPBMC$toc,res=1.0)
good = standard10X(out,res=1.0)
#Get the reference annotation in each cluster
bad@meta.data$Annotation = PBMC_DR[rownames(bad@meta.data),'Annotation']
good@meta.data$Annotation = PBMC_DR[rownames(good@meta.data),'Annotation']
#Get the consensu annotation for each cluster
tt = table(bad@active.ident,bad@meta.data$Annotation)
bad@meta.data$clusterAnn = setNames(colnames(tt)[apply(tt,1,which.max)],rownames(tt))[as.character(bad@active.ident)]
tt = table(good@active.ident,good@meta.data$Annotation)
good@meta.data$clusterAnn = setNames(colnames(tt)[apply(tt,1,which.max)],rownames(tt))[as.character(good@active.ident)]
#Get shared cluster IDs
clust = data.frame(barcode = names(good@active.ident),
                   goodClust = as.character(good@active.ident),
                   badClust = as.character(bad@active.ident[names(good@active.ident)]),
                   stringsAsFactors=TRUE)
tt = aggregate(goodClust ~ badClust,FUN=table,data=clust)
#Number moved
rowSums(tt$goodClust)-apply(tt$goodClust,1,max)
#Make a heatmap of the move, normalised to fraction of source cell
nMap = tt$goodClust[order(as.numeric(as.character(tt$badClust))),]
rownames(nMap) = paste0('dirtyCluster',as.character(tt$badClust)[order(as.numeric(as.character(tt$badClust)))])
nMap = nMap[,order(as.numeric(colnames(nMap)))]
colnames(nMap) = paste0('cleanCluster',colnames(nMap))
#Save a count version
clMap = nMap
#Normalise it
nMap = nMap/rowSums(nMap)
#Make the heatmap
col = c('#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#000000')
col = circlize::colorRamp2(seq(0,1,length.out=length(col)),col)
hm = Heatmap(nMap[seq(nrow(nMap),1),],
             col = col,
             name = 'Intersection',
             cluster_rows=FALSE,
             cluster_columns=FALSE,
             row_title='Uncorrected Clusters',
             row_title_side = 'left',
             row_names_side = 'left',
             column_title='Corrected Clusters',
             column_title_side='bottom',
             row_labels = rev(gsub('.*Cluster','',rownames(nMap))),
             column_labels = gsub('.*Cluster','',colnames(nMap)),
             row_title_gp = gpar(fontsize=8),
             row_names_gp = gpar(fontsize=8),
             column_title_gp = gpar(fontsize=8),
             column_names_gp = gpar(fontsize=8),
             heatmap_legend_param = list(title_gp = gpar(fontsize=8),
                                         labels_gp = gpar(fontsize=8))
             )
doPlot = function(e){
  print(hm)
}
makePlots(file.path(plotDir,'PBMC_clusterChange'),doPlot,splitRaster=FALSE,width=3.5,height=3.5,pointsize=8)
#pdf(file.path(plotDir,'PBMC_clusterChange.pdf'))
#print(hm)
#dev.off()
#png(file.path(plotDir,'PBMC_clusterChange.png'),width=7,height=7,units='in',res=rasterRes)
#print(hm)
#dev.off()
##############################################################################
# Show that on average cleaning up soup increases the logFC between clusters.
#Full cross-mapping list
cMap = apply(clMap,1,function(e) colnames(clMap)[e>0])
#Use the good clustering the default one in both maps
good@meta.data$cMap = good@active.ident
bad@meta.data$cMap = good@meta.data[rownames(bad@meta.data),'seurat_clusters']
bad = SetIdent(bad,rownames(bad@meta.data),bad@meta.data$cMap)
####This will set everything correctly, except those that are split
###cMap = setNames(colnames(tt$goodClust)[apply(tt$goodClust,1,which.max)],as.character(tt$badClust))
###bad@meta.data$cMap = cMap[as.character(bad@active.ident)]
###good@meta.data$cMap = as.character(good@active.ident)
####Set this to the active ident
###good = SetIdent(good,rownames(good@meta.data),good@meta.data$cMap)
###bad = SetIdent(bad,rownames(bad@meta.data),bad@meta.data$cMap)
#Get markers and compare.  This is slow, so save after it's done
tgtFile=file.path(plotDir,'PBMC_clusterMarkers.RDS')
if(!file.exists(tgtFile)){
  gMarks = FindAllMarkers(good)
  bMarks = FindAllMarkers(bad)
  saveRDS(list(gMarks=gMarks,bMarks=bMarks),tgtFile)
}else{
  tmp = readRDS(tgtFile)
  bMarks = tmp$bMarks
  gMarks = tmp$gMarks
}
#Get complete set of old/new
mMarks = bMarks
m = match(with(gMarks,interaction(cluster,gene)),with(bMarks,interaction(cluster,gene)))
mMarks = rbind(mMarks,gMarks[is.na(m),])
mMarks = mMarks[,c('gene','cluster')]
mMarks$oldFC = NA
mMarks$newFC = NA
mMarks$inOld = with(mMarks,interaction(cluster,gene)) %in% with(bMarks,interaction(cluster,gene))
mMarks$inNew = with(mMarks,interaction(cluster,gene)) %in% with(gMarks,interaction(cluster,gene))
#Get the old logFC from stractch
for(cl in unique(mMarks$cluster)){
  w = which(mMarks$cluster==cl)
  a = apply(bad@assays$RNA@data[mMarks[w,'gene'],bad@active.ident==cl],1,function(e) log(mean(expm1(e))+1))
  b = apply(bad@assays$RNA@data[mMarks[w,'gene'],bad@active.ident!=cl],1,function(e) log(mean(expm1(e))+1))
  mMarks$oldFC[w] = a-b
  a = apply(good@assays$RNA@data[mMarks[w,'gene'],good@active.ident==cl],1,function(e) log(mean(expm1(e))+1))
  b = apply(good@assays$RNA@data[mMarks[w,'gene'],good@active.ident!=cl],1,function(e) log(mean(expm1(e))+1))
  mMarks$newFC[w] = a-b
}
tmp = mMarks[!is.na(mMarks$oldFC),]
tmp = tmp[sign(tmp$newFC)==sign(tmp$oldFC),]
tmp = tmp[abs(log(tmp$newFC/tmp$oldFC))<1,]
x = tmp$oldFC
y = tmp$newFC/tmp$oldFC
doPlot = function(e){
  layout(matrix(c(1,2),ncol=2),widths=c(10,1))
  par(mar=c(5.1,4.1,4.1,0))
  plot(x,y,
       las=1,
       type='n',
       xlab='Uncorrected Fold Change',
       ylab='CorrectedFC/UncorrectedFC',
       frame.plot=FALSE)
  col = c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b')
  hb = addHexbin(x,y,colSet=c(grey(0.95),grey(0)),xbins=100,logCounts=TRUE,maxCnts=1000)
  abline(h=1,col='black',lty=3)
  #Add special ones as points
  #w = which(xor(tmp$inOld,tmp$inNew))
  #points(x[w],y[w],
  #     col=ifelse(tmp$inOld[w],'red','green'),
  #     pch=19,
  #     cex=0.2)
  addColorBar('topright',col=hb$col,breaks=hb$breaks,ticks=0:3,title='log10(counts)')
  #Add right plot giving density
  par(mar=c(5.1,0,4.1,0))
  den = density(tmp$newFC/tmp$oldFC)
  plot(den$y,den$x,type='n',
       frame.plot=FALSE,
       xaxt='n',
       yaxt='n',
       xlab='Frequency',
       ylab='',
  xpd=NA)
  lines(den$y,den$x,col='black')
  polygon(c(0,den$y,0),c(0,den$x,0),col=grey(0.8,alpha=0.6))
  abline(h=1,col='black',lty=3)
}
makePlots(file.path(plotDir,'PBMC_foldChange'),doPlot,splitRaster=FALSE,width=5,height=3.5,pointsize=8)
####Main plot area
###gg = ggplot(tmp,aes(oldFC,(newFC/oldFC))) +
###  stat_binhex(aes(fill=log10(..count..)),bins=100) +
###  geom_hline(yintercept=1,colour='red') +
###  geom_point(data=tmp[xor(tmp$inOld,tmp$inNew),],aes(colour=ifelse(inOld,'Uncorrected','Corrected')),size=0.1,alpha=0.6) +
###  scale_y_continuous(position = "right") +
###  xlab('Uncorrected Fold Change') +
###  ylab('CorrectedFC/UncorrectedFC') +
###  guides(colour = guide_legend(override.aes = list(size=3))) +
###  labs(colour='Exclusive to')
####Add marginal densities
###xMarg = ggplot(tmp,aes(oldFC)) +
###  geom_density(fill='black',alpha=0.3) + 
###  geom_density(data=tmp[xor(tmp$inOld,tmp$inNew),],aes(fill=ifelse(inOld,'Uncorrected','Corrected')),alpha=0.5) +
###  theme_void() +
###  guides(fill=FALSE)
###yMarg = ggplot(tmp,aes(newFC/oldFC)) +
###  geom_density(fill='black',alpha=0.3) + 
###  geom_density(data=tmp[xor(tmp$inOld,tmp$inNew),],aes(fill=ifelse(inOld,'Uncorrected','Corrected')),alpha=0.5) +
###  geom_vline(xintercept=1,colour='red') +
###  theme_void() +
###  coord_flip() +
###  scale_y_reverse() +
###  guides(fill=FALSE)
#### Align histograms with scatterplot
###aligned_x_hist = align_plots(xMarg, gg, align = "v")[[1]]
###aligned_y_hist = align_plots(yMarg, gg, align = "h")[[1]]
###pg = plot_grid(
###    #NULL
###    #, aligned_x_hist
###    aligned_y_hist
###    , gg
###    , ncol = 2
###    , nrow = 1
###    #, rel_heights = c(0.2, 1)
###    , rel_widths = c(0.2, 1)
###  )
###pdf(file.path(plotDir,'PBMC_foldChange.pdf'),width=10)
###plot(pg)
###dev.off()
###png(file.path(plotDir,'PBMC_foldChange.png'),width=10,height=7,units='in',res=rasterRes)
###plot(pg)
###dev.off()
#gg = ggplot(mMarks,aes(oldFC,newFC/oldFC)) +
#  stat_binhex(aes(fill=log10(..count..)),bins=50) +
#  geom_hline(yintercept=1,colour='red') +
#  geom_point(data=mMarks[xor(mMarks$inOld,mMarks$inNew),],aes(colour=ifelse(inOld,'Uncorrected','Corrected')),size=0.2) +
#  xlab('Uncorrected Fold Change') +
#  ylab('CorrectedFC/UncorrectedFC') +
#  labs(colour='Exclusive to')
#pdf(file.path(plotDir,'PBMC_foldChange.pdf'))
#plot(gg)
#dev.off()
#png(file.path(plotDir,'PBMC_foldChange.png'),width=7,height=7,units='in',res=rasterRes)
#plot(hm)
#dev.off()
###################################
# Compare LYZ in MNP versus others
tgtGenes='LYZ'
tgtTypes='MNP'
dePlot = function(tgtGenes,tgtTypes){
  colSet = c('black',grey(0.6))
  expGood = split(colSums(good@assays$RNA@data[tgtGenes,,drop=FALSE]),good@meta.data$clusterAnn %in% tgtTypes)
  expBad = split(colSums(bad@assays$RNA@data[tgtGenes,,drop=FALSE]),bad@meta.data$clusterAnn %in% tgtTypes)
  #Plot it
  df = data.frame(src=rep(c('Corrected','Uncorrected'),each=length(bad@active.ident)),
                  exp = c(unlist(expGood),unlist(expBad)),
                  type = ifelse(grepl('^TRUE.',names(c(unlist(expGood),unlist(expBad)))),tgtTypes[1],'Other'))
  #Block left 0-2 (with buffer)
    #Block inner-left 0.2-1 (with buffer)
      #0.3 - 0.9
    #Block inner-right 1-1.8 (with buffer)
      #1.1 - 1.7
  #Block right 2-4 (with buffer)
    #
  plot(0,0,
       xlim = c(0,4),
       ylim = c(0,max(df$exp)),
       las=1,
       type='n',
       xaxt='n', 
       xlab='',
       ylab=sprintf('log(%s Expression)',tgtGenes),
       frame.plot=FALSE
       )
  axis(1,at=c(1,3),labels=c('Uncorrected','Corrected'))
  #Define the plots
  #Within species
  #10X
  sets = list(list(cs='Uncorrected',tn=tgtTypes,col=colSet[1],mid=0.6),
              list(cs='Uncorrected',tn='Other',col=colSet[2],mid=1.4),
              list(cs='Corrected',tn=tgtTypes,col=colSet[1],mid=2.6),
              list(cs='Corrected',tn='Other',col=colSet[2],mid=3.4))
  for(set in sets){
    x = df$exp[df$src == set$cs & df$type==set$tn]
    points(set$mid+ runif(length(x),-.25,.25),
           x,
           col=colAlpha(set$col,0.1),
           pch=19,
           cex=0.1
    )
    boxplot(x,
            col=NA,
            at=set$mid,
            boxwex=1.0,
            border=set$col,
            outline=FALSE,
            add=TRUE,
            xaxt='n', 
            yaxt='n', 
            bty='n',
            frame.plot=FALSE
    )
  }
  #legend('topright',
  #       legend=c(tgtTypes,'Other'),
  #       fill=colPal,
  #       bty='n',
  #       title='CellGroup',
  #       xpd=NA)
}
doPlot = function(e) dePlot('LYZ','MNP')
makePlots(file.path(plotDir,'PBMC_exampleA'),doPlot,splitRaster=FALSE,width=3.5,height=3.5,pointsize=8)
##The bar plot
#gg = dePlot('LYZ','MNP')
#gg = gg + ggtitle("Change in fold-change in LYZ")
#pdf(file.path(plotDir,'PBMC_example1A.pdf'))
#plot(gg)
#dev.off()
#png(file.path(plotDir,'PBMC_example1A.png'),width=7,height=7,units='in',res=rasterRes)
#plot(gg)
#dev.off()
#And the umap
gg = plotChangeMap(scPBMC,good@assays$RNA@counts,'LYZ')
df = gg$data
doPlot = function(e){
  plot(df$RD1,df$RD2,
       type='n',
       frame.plot=FALSE,
       xlab='tSNE1',
       ylab='tSNE2')
  col = c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704')
  #col = c(grey(0.95),'black')
  cc = circlize::colorRamp2(seq(0,1,length.out=length(col)),col)#,transparency = 0.2)
  #Do as in the earlier panel
  w = is.na(df$relChange)
  points(df$RD1[w],df$RD2[w],cex=0.1,pch=19,col='black')
  points(df$RD1[!w],df$RD2[!w],
         col=cc(df$relChange[!w]),
         #col='black',
         pch=19,
         #lwd=0.01,
         cex=0.7)
  addColorBar('topleft',col=col,breaks=c(1,0),ticks=5)
}
makePlots(file.path(plotDir,'PBMC_exampleB'),doPlot,splitRaster=FALSE,width=3.5,height=3.5,pointsize=8)
###gg = gg + ggtitle("Change in expression of LYZ")
###pdf(file.path(plotDir,'PBMC_example1B.pdf'))
###plot(gg)
###dev.off()
###png(file.path(plotDir,'PBMC_example1B.png'),width=7,height=7,units='in',res=rasterRes)
###plot(gg)
###dev.off()
#####################
# Genes to consider
gg = plotMarkerDistribution(scPBMC) +
  scale_size(range=c(0.5,2))
doPlot = function(e) plot(gg)
makePlots(file.path(plotDir,'guessPBMCgenes'),doPlot,width=7,height=4,pointsize=8)
#pdf(file.path(plotDir,'guessPBMCgenes.pdf'))
#plot(gg)
#dev.off()
#png(file.path(plotDir,'guessPBMCgenes.png'),width=7,height=7,units='in',res=rasterRes)
#plot(gg)
#dev.off()





###############
# Tumour Data #
###############
cMani = read.table(fp_cMani,sep='\t',header=TRUE)
cMani$Path = file.path('Data/KidneyTumour',cMani$Label)
#Keep just the tumour/NR channels
cMani = cMani[cMani$TissueDiseaseState!='Normal',]
#Exclude RCC3, it's full of failure
cMani = cMani[cMani$Experiment!='RCC3',]
#Fix path
cMani$Path = gsub('KidneyTumour','Tumour',cMani$Path)
#Define markers
geneList = list(
  igGenes = c('IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHD','IGHE','IGHM',
              'IGLC1','IGLC2','IGLC3','IGLC4','IGLC5','IGLC6','IGLC7',
              'IGKC'),
  hgGenes = c('HBA1','HBA2','HBB','HBD','HBE1','HBG1','HBG2','HBM','HBQ1','HBZ'),
  mastGenes = c('TPSB2','TPSAB1'),
  colGenes = c('COL1A1','COL1A2','COL3A1'),
  s100Genes = c('S100A8','S100A9')
  )
#Load each channel
scTum = list()
#Need something else for i=10,11,12
tgtSet='hgGenes'
for(i in seq(nrow(cMani))){
  #Only one sample doesn't really play nice with hgGenes
  if(i %in% c(10,11,12)){
    tgtSet='igGenes'
  }else{
    tgtSet='hgGenes'
  }
  message(cMani$Label[i])
  scTum[[i]] = load10X(cMani$Path[i])
  ute = estimateNonExpressingCells(scTum[[i]],nonExpressedGeneList = geneList[tgtSet])
  if(sum(ute)==0)
    ute = estimateNonExpressingCells(scTum[[i]],nonExpressedGeneList = geneList[tgtSet],clusters=FALSE,maximumContamination=0.6,FDR=0.05)
  plotMarkerMap(scTum[[i]],geneSet=geneList[[tgtSet]],useToEst=ute)
  scTum[[i]] = calculateContaminationFraction(scTum[[i]],nonExpressedGeneList=geneList[tgtSet],useToEst=ute)
}
#Make a spot the difference
dirty = lapply(scTum,function(e) e$toc)
cleaned = lapply(scTum,adjustCounts,verbose=2,nCores=nCores)
for(i in seq(nrow(cMani))){
  colnames(cleaned[[i]]) = paste0(cMani$Label[i],'___',colnames(cleaned[[i]]))
  colnames(dirty[[i]]) = paste0(cMani$Label[i],'___',colnames(dirty[[i]]))
}
dirty = do.call(cbind,dirty)
cleaned = do.call(cbind,cleaned)
dirty = standard10X(dirty,res=1.0)
cleaned = standard10X(cleaned,res=1.0)
saveRDS(list(dirty,cleaned),file.path(plotDir,'Tumour_seurat_init.RDS'))
###################
# Basic source map
#Annotation for tSNE.  Extremely quick and rough.
#mrks = quickMarkers(dirty@assays$RNA@counts,setNames(dirty@meta.data$seurat_clusters,rownames(dirty@meta.data)),N=Inf)
ann = c(
        '0' = 'T',
        '1' = 'NK',
        '2' = 'T',
        '3' = 'MNP',
        '4' = 'MNP',
        '5' = 'T',
        '6' = 'MNP',
        '7' = 'MNP',
        '8' = 'ccRCC',
        '9' = 'Endothelium',
        '10' = 'T',
        '11' = 'T',
        '12' = 'ccRCC',
        '13' = 'NK',
        '14' = 'NK',
        '15' = 'ccRCC',
        '16' = 'MNP',
        '17' = 'MNP',
        '18' = 'NK',
        '19' = 'T',
        '20' = 'T',
        '21' = 'RBC',
        '22' = 'Fibroblast',
        '23' = 'Wilms',
        '24' = 'Endothelium',
        '25' = 'MNP',
        '26' = 'pRCC',
        '27' = 'ccRCC',
        '28' = 'RBC',
        '29' = 'B',
        '30' = 'Fibroblast',
        '31' = 'ccRCC',
        '32' = 'B',
        '33' = 'ccRCC',
        '34' = 'B',
        '35' = 'MNP',
        '36' = 'ccRCC'
        )
#ann = c(
#        '0' = 'T',
#        '1' = 'T',
#        '2' = 'MNP',
#        '3' = 'MNP',
#        '4' = 'MNP',
#        '5' = 'MNP',
#        '6' = 'NK',
#        '7' = 'Endothelium',
#        '8' = 'T',
#        '9' = 'T',
#        '10' = 'NK',
#        '11' = 'ccRCC',
#        '12' = 'NK',
#        '13' = 'ccRCC',
#        '14' = 'NK',
#        '15' = 'ccRCC',
#        '16' = 'MNP',
#        '17' = 'T',
#        '18' = 'T',
#        '19' = 'Fibroblast',
#        '20' = 'RBC',
#        '21' = 'MNP',
#        '22' = 'Wilms',
#        '23' = 'ccRCC',
#        '24' = 'MNP',
#        '25' = 'RBC',
#        '26' = 'pRCC',
#        '27' = 'Endothelium',
#        '28' = 'ccRCC',
#        '29' = 'B',
#        '30' = 'B',
#        '31' = 'Fibroblast',
#        '32' = 'ccRCC',
#        '33' = 'ccRCC',
#        '34' = 'ccRCC',
#        '35' = '?',
#        '36' = '?'
#        )
#df = cbind(dirty@meta.data,dirty@reductions$tsne@cell.embeddings)
#gg = ggplot(df,aes(tSNE_1,tSNE_2)) +
#  geom_point(aes(colour=orig.ident),size=0.1) +
#  scale_colour_manual(values=colPalLarge[seq(7)]) +
#  guides(colour = guide_legend(override.aes = list(size=1))) +
#  labs(colour='Donor') +
#  xlab('tSNE1') +
#  ylab('tSNE2') +
#  theme_classic()
#Order tweaked to make colours sensible
colPalTum = c('T'='#B3DE69',
              RBC='#FB8072',
              NK="#8DD3C7",
              Endothelium="#FCCDE5",
              ccRCC="#FDB462",
              pRCC="#BC80BD",
              MNP='#FFFFB3',
              Wilms='#BEBADA',
              Fibroblast="#D9D9D9",
              B='#80B1D3')
doPlot = function(mode=c('All','Raster','Text')){
  mode = match.arg(mode)
  df = cbind(dirty@meta.data,dirty@reductions$tsne@cell.embeddings)
  df$cluster = df$seurat_clusters
  df$Annotation = ann[as.character(df$seurat_clusters)]
  plot(df$tSNE_1,df$tSNE_2,
       type='n',
       frame.plot=FALSE,
       xlab='tSNE1',
       ylab='tSNE2'
       )
  #Work out which are the approved clusters
  wKeep = df$Annotation!='?'
  clust = factor(df$cluster[wKeep])
  #shadeColMap = setNames(colAlpha(colPalTum[seq_along(unique(df$Annotation[wKeep]))],0.4),unique(df$Annotation[wKeep]))
  shadeColMap = setNames(colAlpha(colPalTum,0.4),names(colPalTum))
  shadeCols = shadeColMap[ann[levels(clust)]]
  addDensityContours(df$tSNE_1[wKeep],df$tSNE_2[wKeep],clust,shadeCol=shadeCols,useKS=TRUE)
  points(df$tSNE_1,df$tSNE_2,
         cex=0.1,
         pch=19,
         col = grey(0.2,0.2)
         )
  #Precise placement
  annPos = list(ccRCC=c(-50,-30),
                pRCC=c(-38,-42),
                Wilms=c(-48,40),
                NK=c(-22,38),
                Endothelium=c(10,50),
                B=c(-30,45),
                'T'=c(35,25),
                MNP=c(35,-35),
                RBC=c(20,-50),
                Fibroblast=c(-30,-50))
  for(nom in names(annPos)){
    points(annPos[[nom]][1],annPos[[nom]][2],
           bg = shadeColMap[nom],
           lwd = 1,
           cex = 2,
           pch =22)
    text(annPos[[nom]][1],annPos[[nom]][2],nom,
         adj = ifelse(nchar(nom)<3,ifelse(nchar(nom)<2,-0.8,-0.5),-0.2))
  }
}
makePlots(file.path(plotDir,'Tumour_ann'),doPlot,splitRaster=FALSE,width=4.5,height=4.5,pointsize=8)
#############################
# Compare cluster membership
#Get shared cluster IDs
clust = data.frame(barcode = rownames(cleaned@meta.data),
                   goodClust = as.character(cleaned@meta.data$seurat_clusters),
                   badClust = as.character(dirty@meta.data[rownames(cleaned@meta.data),'seurat_clusters']),
                   stringsAsFactors=TRUE)
tt = aggregate(goodClust ~ badClust,FUN=table,data=clust)
#Number moved
rowSums(tt$goodClust)-apply(tt$goodClust,1,max)
#Make a heatmap of the move, normalised to fraction of source cell
nMap = tt$goodClust[order(as.numeric(as.character(tt$badClust))),]
rownames(nMap) = paste0('dirtyCluster',as.character(tt$badClust)[order(as.numeric(as.character(tt$badClust)))])
nMap = nMap[,order(as.numeric(colnames(nMap)))]
colnames(nMap) = paste0('cleanCluster',colnames(nMap))
#Save a count version
clMap = nMap
#Normalise it
nMap = nMap/rowSums(nMap)
#Make the heatmap
col = c('#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#000000')
col = circlize::colorRamp2(seq(0,1,length.out=length(col)),col)
hm = Heatmap(nMap[seq(nrow(nMap),1),],
             col = col,
             name = 'Intersection',
             cluster_rows=FALSE,
             cluster_columns=FALSE,
             row_title='Uncorrected Clusters',
             row_title_side = 'left',
             row_names_side = 'left',
             column_title='Corrected Clusters',
             column_title_side='bottom',
             row_labels = rev(gsub('.*Cluster','',rownames(nMap))),
             column_labels = gsub('.*Cluster','',colnames(nMap)),
             row_title_gp = gpar(fontsize=8),
             row_names_gp = gpar(fontsize=8),
             column_title_gp = gpar(fontsize=8),
             column_names_gp = gpar(fontsize=8),
             heatmap_legend_param = list(title_gp = gpar(fontsize=8),
                                         labels_gp = gpar(fontsize=8))
             )
doPlot = function(e) {print(hm)}
makePlots(file.path(plotDir,'Tumour_clusterChange'),doPlot,splitRaster=FALSE,width=4.5,height=4.0,pointsize=8)
##############################################################################
# Show that on average cleaning up soup increases the logFC between clusters.
#Full cross-mapping list
cMap = apply(clMap,1,function(e) colnames(clMap)[e>0])
#Use the good clustering the default one in both maps
cleaned@meta.data$cMap = cleaned@meta.data$seurat_clusters
dirty@meta.data$cMap = cleaned@meta.data[rownames(dirty@meta.data),'seurat_clusters']
dirty = SetIdent(dirty,rownames(dirty@meta.data),dirty@meta.data$cMap)
####Split them into those that have 1-1 correspondent and those that don't
###badToGood = data.frame(badClust = rownames(nMap),
###                       goodClust = colnames(nMap)[apply(nMap,1,which.max)],
###                       mapFrac = apply(nMap,1,max))
###compMap = list(c0=c('d0'),
###               c1=c('d2','d5'),
###               c2=c('d3'),
###               c3=c('d4','d35'),
###               c4=c('d7'),
###               c5=c('d13'),
###               c6=c('d9'),
###               c7=c('d6'),
###               c8=c('d10'),
###               c9=c('d11'),
###               d1=c('c10','c14'),
###               c11=c('d12'),
###               c12=c('d18','d19'),
###               c13=c('d14'),
###               d8=c('c15','c24'),
###               c16=c('d15'),
###               c17=c('d22','d30'),
###               c18=c('d17'),
###               c19=c('d20'),
###               c20=c('d25'),
###               c21=c('d16'),
###               c22=c('d28'),
###               c23=c('d23'),
###               c25=c('d21'),
###               c26=c('d24'),
###               c27=c('d26'),
###               c28=c('d27'),
###               c29=c('d32','d34'),
###               c30=c('d29'),
###               c31=c('d31'),
###               c32=c('d33'),
###               c33=c(),
###               c34=c('d36')
###               )
####Make each of these into a connonical cluster
###cleanToGold = paste0('c',levels(cleaned@active.ident))
###m = match(cleanToGold,names(compMap))
###m[is.na(m)] = rep(seq_along(compMap),lengths(compMap))[match(cleanToGold[is.na(m)],unlist(compMap))]
###cleanToGold = setNames(paste0('g',m),cleanToGold)
###dirtyToGold = paste0('d',levels(dirty@active.ident))
###m = match(dirtyToGold,names(compMap))
###m[is.na(m)] = rep(seq_along(compMap),lengths(compMap))[match(dirtyToGold[is.na(m)],unlist(compMap))]
###dirtyToGold = setNames(paste0('g',m),dirtyToGold)
###good = SetIdent(cleaned,rownames(cleaned@meta.data),cleanToGold[paste0('c',as.character(cleaned@meta.data$seurat_clusters))])
###bad = SetIdent(dirty,rownames(dirty@meta.data),dirtyToGold[paste0('d',as.character(dirty@meta.data$seurat_clusters))])
#Get markers and compare
tgtFile=file.path(plotDir,'Tumour_clusterMarkers.RDS')
if(!file.exists(tgtFile)){
  gMarks = FindAllMarkers(cleaned)
  bMarks = FindAllMarkers(dirty)
  saveRDS(list(gMarks=gMarks,bMarks=bMarks),tgtFile)
}else{
  tmp = readRDS(tgtFile)
  bMarks = tmp$bMarks
  gMarks = tmp$gMarks
}
#Get complete set of old/new
mMarks = bMarks
m = match(with(gMarks,interaction(cluster,gene)),with(bMarks,interaction(cluster,gene)))
mMarks = rbind(mMarks,gMarks[is.na(m),])
mMarks = mMarks[,c('gene','cluster')]
mMarks$oldFC = NA
mMarks$newFC = NA
mMarks$inOld = with(mMarks,interaction(cluster,gene)) %in% with(bMarks,interaction(cluster,gene))
mMarks$inNew = with(mMarks,interaction(cluster,gene)) %in% with(gMarks,interaction(cluster,gene))
#Get the logFC from stractch
for(cl in unique(mMarks$cluster)){
  w = which(mMarks$cluster==cl)
  a = apply(dirty@assays$RNA@data[mMarks[w,'gene'],dirty@active.ident==cl],1,function(e) log(mean(expm1(e))+1))
  b = apply(dirty@assays$RNA@data[mMarks[w,'gene'],dirty@active.ident!=cl],1,function(e) log(mean(expm1(e))+1))
  mMarks$oldFC[w] = a-b
  a = apply(cleaned@assays$RNA@data[mMarks[w,'gene'],cleaned@active.ident==cl],1,function(e) log(mean(expm1(e))+1))
  b = apply(cleaned@assays$RNA@data[mMarks[w,'gene'],cleaned@active.ident!=cl],1,function(e) log(mean(expm1(e))+1))
  mMarks$newFC[w] = a-b
}
tmp = mMarks[!is.na(mMarks$oldFC),]
tmp = tmp[sign(tmp$newFC)==sign(tmp$oldFC),]
tmp = tmp[abs(log(tmp$newFC/tmp$oldFC))<1,]
x = tmp$oldFC
y = tmp$newFC/tmp$oldFC
doPlot = function(e){
  layout(matrix(c(1,2),ncol=2),widths=c(10,1))
  par(mar=c(5.1,4.1,4.1,0))
  plot(x,y,
       las=1,
       type='n',
       xlab='Uncorrected Fold Change',
       ylab='CorrectedFC/UncorrectedFC',
       frame.plot=FALSE)
  col = c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b')
  hb = addHexbin(x,y,colSet=c(grey(0.95),grey(0)),xbins=100,logCounts=TRUE,maxCnts=2000)
  abline(h=1,col='black',lty=3)
  #Add special ones as points
  #w = which(xor(tmp$inOld,tmp$inNew))
  #points(x[w],y[w],
  #     col=ifelse(tmp$inOld[w],'red','green'),
  #     pch=19,
  #     cex=0.2)
  addColorBar('topright',col=hb$col,breaks=hb$breaks,ticks=0:3,title='log10(counts)')
  #Add right plot giving density
  par(mar=c(5.1,0,4.1,0))
  den = density(tmp$newFC/tmp$oldFC)
  plot(den$y,den$x,type='n',
       frame.plot=FALSE,
       xaxt='n',
       yaxt='n',
       xlab='Frequency',
       ylab='',
  xpd=NA)
  lines(den$y,den$x,col='black')
  polygon(c(0,den$y,0),c(0,den$x,0),col=grey(0.8,alpha=0.6))
  abline(h=1,col='black',lty=3)
}
makePlots(file.path(plotDir,'Tumour_foldChange'),doPlot,splitRaster=FALSE,width=5,height=3.5,pointsize=8)
########################
# MarkerMap for Example
tgts=c('SAA1','COL1A1','SAA2','HBB','S100A9')
for(tgt in tgts){
  df = list()
  for(i in seq_along(scTum)){
    tmp = cleaned@assays$RNA@counts[,grep(cMani$Label[i],colnames(cleaned))]
    colnames(tmp) = gsub('.*___','',colnames(tmp))
    #Get gene in old and new
    old = scTum[[i]]$toc[tgt,]
    new = tmp[tgt,]
    relChange = (old-new)/old
    df[[i]] = plotChangeMap(scTum[[i]],tmp,tgt,logData=FALSE)
    df[[i]]$data$relChange = relChange
    #x = df[[i]]$data[!is.na(relChange),]
  }
  df = lapply(seq_along(df),function(i) {e=df[[i]]$data;rownames(e) = paste0(cMani$Label[i],'___',rownames(e));e})
  df = do.call(rbind,df)
  #Update UMAP coordinates
  m = match(rownames(df),colnames(dirty))
  df$RD1 = dirty@reductions$tsne@cell.embeddings[m,1]
  df$RD2 = dirty@reductions$tsne@cell.embeddings[m,2]
  #Order
  df = df[order(!is.na(df$relChange)),]
  doPlot = function(e){
    plot(df$RD1,df$RD2,
         type='n',
         frame.plot=FALSE,
         xlab='tSNE1',
         ylab='tSNE2')
    col = c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704')
    #col = c(grey(0.95),'black')
    cc = circlize::colorRamp2(seq(0,1,length.out=length(col)),col)#,transparency = 0.2)
    #Do as in the earlier panel
    w = is.na(df$relChange)
    points(df$RD1[w],df$RD2[w],cex=0.01,pch=19,col=grey(0.8))
    points(df$RD1[!w],df$RD2[!w],
           col=colAlpha(cc(df$relChange[!w]),1.0),
           #col='black',
           pch=19,
           #lwd=0.1,
           cex=0.3)
    #addColorBar('topleft',col=col,breaks=c(1,0),ticks=5)
  }
  makePlots(file.path(plotDir,sprintf('Tumour_change_%s',tgt)),doPlot,splitRaster=FALSE,width=4.5,height=4.5,pointsize=8)
  #Truncate to prevent -Inf -> NA coloured dots
  #df$relChange[which(df$relChange< -2)] = -2
  #df$relChange[which(df$relChange>0)] = 0
  #df$relChange = 10**df$relChange
  #Make plot
  ###makePlot = function(){
  ###  heatCols = c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704')
  ###  cc = circlize::colorRamp2(seq(min(df$relChange,na.rm=TRUE),max(df$relChange,na.rm=TRUE),length.out=length(heatCols)),heatCols)
  ###  df = df[order(!is.na(df$relChange),df$relChange),]
  ###  par(mar=c(5,4,4,2+2)+0.1)
  ###  plot(df$RD1,df$RD2,
  ###       frame.plot=FALSE,
  ###       col = ifelse(is.na(df$relChange),'grey',cc(df$relChange)),
  ###       xlab='ReducedDim1',
  ###       ylab='ReducedDim2',
  ###       pch=19,
  ###       cex=0.1)
  ###  nBlocks = 101
  ###  nTicks = 5 
  ###  xx = seq(0,-2,length.out=nBlocks)
  ###  labs = rep(NA,length(xx))
  ###  w = round(seq(1,nBlocks,length.out=nTicks))
  ###  labs[w] = as.character(round(xx[w],2))
  ###  tmp = legend(x=max(df$RD1),
  ###               y=mean(range(df$RD2)),
  ###               legend=labs, #Tick markers
  ###               fill=cc(xx), #Fill blocks
  ###               border=NA, #With no border
  ###               y.intersp=0.1, #Bunch up blocks to make continuousish
  ###               adj = c(0.4,0.5),#Move text closer to bar
  ###               cex=0.6, #Size of text
  ###               xjust=0, #To right of max
  ###               yjust=0.5, #Center around y location
  ###               bty='n', #No border
  ###               xpd=TRUE #Stick outside plot area
  ###  )
  ###  #Add the title
  ###  text(tmp$rect$left+0.5*tmp$rect$w, #In the middle of the bar in x-direction
  ###        tmp$rect$top,
  ###        adj = c(0.5,-1), #Move up a bit
  ###        labels='log10(RelChange)',
  ###        xpd=TRUE #Don't clip
  ###  )
  ###}
  ###gg = ggplot(df,aes(RD1,RD2)) +
  ###  geom_point(aes(col=10**relChange),size=0.1) +
  ###  xlab('ReducedDim1') +
  ###  ylab('ReducedDim2') +
  ###  labs(colour='log10(relChange)') + 
  ###  ggtitle(sprintf('Change in expression of %s',tgt))
  ###gg = gg + scale_colour_gradientn(colours=c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704'))
  ###pdf(file.path(plotDir,sprintf('Tumour_example%s.pdf',tgt)),width=10)
  ###plot(gg)
  ####makePlot()
  ###dev.off()
  ###png(file.path(plotDir,sprintf('Tumour_example%s.png',tgt)),width=10,height=7,units='in',res=rasterRes)
  ###plot(gg)
  ####makePlot()
  ###dev.off()
}
#################################################################
# Calculate the shannon entropy before and after changes locally
#How many random bootstraps
nBoots = 250
#How many cells from each batch
nCells = 100
nNgbs = 100
#Find neighbours
tgtObjs = list('uncorrected'=dirty,'corrected'=cleaned)
tgtObjs = lapply(tgtObjs,FindNeighbors,k.param=nNgbs)
cellIDs = rownames(tgtObjs[[1]]@meta.data)
batch = gsub('_[0-9]$','',gsub('___.*','',rownames(tgtObjs[[1]]@meta.data)))
out = matrix(0,nrow=nBoots,ncol=2)
pb = txtProgressBar(max=nBoots)
for(ii in seq(nBoots)){
  setTxtProgressBar(pb,ii)
  #Pick random cells
  tgtCells = unlist(lapply(split(cellIDs,batch),sample,nCells))
  jj=0
  for(srat in list(dirty,cleaned)){
    jj = jj + 1
    #Get the neighbours
    cNgbs = as.matrix(srat@graphs$RNA_nn[tgtCells,]>0)
    ent = 0
    kk=0
    for(tgtCell in tgtCells){
      kk = kk +1
      #cNgbs = which(srat@graphs$RNA_nn[tgtCell,]>0)
      #Get batch frequencies
      lCnts = table(batch[which(cNgbs[kk,])])
      lFreq = lCnts/sum(lCnts)
      ent = ent - sum(lFreq *log(lFreq))
    }
    out[ii,jj] = ent
  }
}
close(pb)
##Make a nice beehive plot of the results
doPlot = function(e){
  #plot(0,0,
  #     type='n',
  #     xaxt='n',
  #     xlim = c(0.42,2.58),
  #     ylim = c(min(out)-12,max(out)+12),
  #     xlab='',
  #     ylab='Batch Entropy',
  #     frame.plot=FALSE)
  #axis(1,at=c(1,2),labels=c('Uncorrected','Corrected'))
  #beeswarm(list(Uncorrected=out[,1],Corrected=out[,2]),
  #         cex=0.2,
  #         pch=19,
  #         col=grey(c(0,0.6)),
  #         spacing=3,
  #         add=FALSE)
  #Stop making things more complicated,  just do the same as the others
  plot(0,0,
       xlim = c(0.5,2.5),
       ylim = range(out),
       las=1,
       type='n',
       xaxt='n', 
       xlab='',
       ylab='Batch Entropy',
       frame.plot=FALSE
       )
  axis(1,at=c(1,2),labels=c('Uncorrected','Corrected'))
  #Define the plots
  #Within species
  #10X
  sets = list(list(cs='Uncorrected',tn=1,col=grey(0),mid=1.0),
              list(cs='Uncorrected',tn=2,col=grey(0),mid=2.0))
  for(set in sets){
    x = out[,set$tn]
    points(set$mid+ runif(length(x),-.25,.25),
           x,
           col=colAlpha(set$col,0.5),
           pch=19,
           cex=0.1
    )
    boxplot(x,
            col=NA,
            at=set$mid,
            boxwex=1.0,
            border=set$col,
            outline=FALSE,
            add=TRUE,
            xaxt='n', 
            yaxt='n', 
            bty='n',
            frame.plot=FALSE
    )
  }
}
makePlots(file.path(plotDir,'Tumour_entropy'),doPlot,splitRaster=FALSE,width=3,height=4.0,pointsize=8)
##########################################
# Example showing that HB genes work well
i=16
ute = estimateNonExpressingCells(scTum[[i]],nonExpressedGeneList = geneList[tgtSet])
gg = plotMarkerMap(scTum[[i]],geneSet=geneList[[tgtSet]],useToEst=ute)
#df = scTum[[i]]$metaData
##First plot the clustering
#makePlot = function(){
#  plot(df$tSNE1,df$tSNE2,
#       frame.plot=FALSE,
#       col = colPalLarge[as.numeric(df$clusters) %% length(colPalLarge) +1],
#       pch=19,
#       cex=0.3,
#       xlab='tSNE1',
#       ylab='tSNE2'
#       )
#  mids = aggregate(cbind(tSNE1,tSNE2) ~ clusters,data=df,FUN=mean)
#  text(mids$tSNE1,mids$tSNE2,
#       labels=mids$clusters,
#       cex=1.5
#       )
#}
doPlot = function(e) plot(gg)
makePlots(file.path(plotDir,'Tumour_egHB'),doPlot,splitRaster=FALSE,width=7,height=7,pointsize=8)


#########
# Other #
#########
#####################################
# Calculate correlation coef for all
getCor = function(sc){
  spFrac = sc$soupProfile$counts
  clFrac = rowSums(sc$toc)
  #Drop the extreme genes
  w = which(spFrac>quantile(spFrac,0.99) | clFrac > quantile(clFrac,0.99))
  spFrac = spFrac[-w]
  clFrac = clFrac[-w]
  #Sub-sample to match
  rat = sum(spFrac)/sum(clFrac)
  clFrac = rbinom(length(clFrac),clFrac,rat)
  cc = cor(spFrac,clFrac)
  return(cc)
}
#Apply it to everything
ccs = c(sapply(scMixes,getCor),getCor(scPBMC),sapply(scTum,getCor))


#' Automatically calculate the contamination fraction
#'
#' The idea of this method is that genes that are highly expressed in the soup and be marker genes for some population can be used to estimate the background.  Marker genes are identified using the tfidf method (see \code{\link{quickMarkers}}).  The contamination fraction is then calculated at the cluster level for each of these genes and clusters are then aggressively pruned to remove those that give implausible estimates.  A final estimate is based on the average across these estimates, under the assumption that biased estimates will not cluster, while true ones will (around the true value).
#'
#' Note that the interpretation of the tf-idf cut-off in this context is that a cut-off t implies that a marker gene has the property that geneFreqGlobal < exp(-t/geneFreqInClust).
#'
#' @param sc The SoupChannel object.
#' @param nMarks How many marker genes to use. These must also pass the tfidf cut-off and be highly expressed in the soup, so there may end up being fewer than this number.
#' @param tfidfMin Minimum value of tfidf to accept for a marker gene.
#' @param soupMin Minimum value of expression in the background to accept for marker.
#' @param rhoMax
#' @param FDR
#' @param clustPerGene
#' @param maxClustFrac
#' @param minCnts
autoCont = function(sc,nMarks = 50,tfidfMin=0.8,soupMin=1e-4,rhoMax=0.5,FDR=0.1,clustPerGene=1,maxClustFrac=0.5,minCnts=1){
  #First collapse by cluster
  s = split(rownames(sc$metaData),sc$metaData$clusters)
  tmp = do.call(cbind,lapply(s,function(e) rowSums(sc$toc[,e,drop=FALSE])))
  ssc = sc 
  ssc$toc = tmp
  ssc$metaData = data.frame(nUMIs = colSums(tmp),row.names=colnames(tmp))
  ###################
  # Get best markers
  #Get the top N soup Genes
  tgts = ssc$soupProfile[order(ssc$soupProfile$est,decreasing=TRUE),]
  tgts = rownames(tgts)[tgts$est>soupMin]
  #Refine this to the best markers we can manage
  mrks = quickMarkers(sc$toc,sc$metaData$clusters,N=Inf)
  #Keep only the ones that are high in soup
  mrks = mrks[mrks$gene %in% tgts,]
  #And only the most specific entry for each gene
  mrks = mrks[order(mrks$gene,-mrks$tfidf),]
  mrks = mrks[!duplicated(mrks$gene),]
  #Order by tfidif maxness
  mrks = mrks[order(-mrks$tfidf),]
  #Apply tf-idf cut-off
  mrks = mrks[mrks$tfidf > tfidfMin,]
  tgts = head(mrks$gene,nMarks)
  ############################
  # Get estimates in clusters
  #Get which ones we'd use and where with canonical method
  tmp = as.list(tgts)
  names(tmp) = tgts
  ute = estimateNonExpressingCells(sc,tmp,maximumContamination=rhoMax,FDR=FDR)
  m = rownames(sc$metaData)[match(rownames(ssc$metaData),sc$metaData$clusters)]
  ute = t(ute[m,])
  colnames(ute) = rownames(ssc$metaData)
  #Now calculate the observed and expected counts for each cluster for 
  expCnts = outer(ssc$soupProfile$est,ssc$metaData$nUMIs)
  rownames(expCnts) = rownames(ssc$soupProfile)
  colnames(expCnts) = rownames(ssc$metaData)
  expCnts = expCnts[tgts,]
  #And the observed ones
  obsCnts = ssc$toc[tgts,]
  #Filter out the shite
  #Get the p-value for this being less than 1
  pp = ppois(obsCnts,expCnts*rhoMax,lower.tail=TRUE)
  qq = p.adjust(pp,method='BH')
  qq = matrix(qq,nrow=nrow(pp),ncol=ncol(pp),dimnames=dimnames(pp))
  #Get the cluster level ratio
  rhos = obsCnts/expCnts
  #Index in range
  rhoIdx = t(apply(rhos,1,function(e) order(order(e))))
  #Make a data.frame with everything
  dd = data.frame(gene = rep(rownames(ute),ncol(ute)),
                  passNonExp = as.vector(ute),
                  rho = as.vector(rhos),
                  rhoIdx = as.vector(rhoIdx),
                  obsCnt = as.vector(obsCnts),
                  expCnt = as.vector(expCnts),
                  qVal = as.vector(qq)
                  )
  dd$pass = dd$obsCnt >= minCnts & 
    dd$qVal < FDR & 
    dd$rhoIdx <= min(clustPerGene,floor(ncol(rhoIdx)*maxClustFrac)) & 
    dd$passNonExp
  #I think the best way to do this is based on the density.
  tmp = density(dd$rho[dd$pass])
  dd$rhoEst = tmp$x[which.max(tmp$y)]
  return(dd)
}
