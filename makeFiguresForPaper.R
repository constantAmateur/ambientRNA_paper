#' Calculates the results for SoupX paper.
#####################
# Working directory #
#####################
setwd('~/Projects/SoupX')

#############
# Libraries #
#############
library(ggplot2)
library(cowplot)
library(reshape2)
library(beeswarm)
library(wesanderson)
library(ComplexHeatmap)
library(Seurat)
library(SoupX)
library(Matrix)
library(splatter)
import('Code/plotParts.R',all=TRUE)

##########
# Params #
##########
#Where to store the results
plotDir='Results/'
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

#############
# Functions #
#############
#' Run standard analysis
standard10X = function(dat,nPCs=50,res=1.0,verbose=FALSE){
  srat = CreateSeuratObject(dat)
  srat = NormalizeData(srat,verbose=verbose)
  srat = ScaleData(srat,verbose=verbose)
  srat = FindVariableFeatures(srat,verbose=verbose)
  srat = RunPCA(srat,verbose=verbose)
  srat = RunUMAP(srat,dims=seq(nPCs),verbose=verbose)
  srat = RunTSNE(srat,dims=seq(nPCs),verbose=verbose)
  srat = FindNeighbors(srat,dims=seq(nPCs),verbose=verbose)
  srat = FindClusters(srat,res=res,verbose=verbose)
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
PBMC_DR = PBMC_metaData
##########################
# Soup profile estimation
#Elbow plot
x = scPBMC$nDropUMI
o = order(-x)
x = x[o]
pdf(file.path(plotDir,'Schematic_fitBG1.pdf'))
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
pdf(file.path(plotDir,'Schematic_fitBG2.pdf'))
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
pdf(file.path(plotDir,'Schematic_rho1.pdf'),useDingbats=FALSE)
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
pdf(file.path(plotDir,'Schematic_rho2.pdf'),useDingbats=FALSE)
  #plot(gg2)
  doPlot()
dev.off()
##################
# Automated genes
#Get marker genes that meet criteria
scPBMC$metaData$clusters = scPBMC$metaData$cluster
mrks = quickMarkers(scPBMC$toc, scPBMC$metaData$clusters, N = Inf)
mrks = mrks[order(mrks$gene, -mrks$tfidf), ]
mrks = mrks[!duplicated(mrks$gene), ]
mrks = mrks[order(-mrks$tfidf), ]
mrks = mrks[mrks$tfidf>1,]
#Show top 3 per cluster
forPlot = lapply(split(mrks,mrks$cluster),head,n=3)
for(cl in names(forPlot)){
  message(sprintf("Genes for cluster %s",cl))
  cat(paste(forPlot[[cl]]$gene,collapse='\n'))
  cat('\n')
}
tmp = autoEstCont(scPBMC)
rhoProbes = seq(0, 1, 0.001)
#Make the histogram
pdf(file.path(plotDir,'Schematic_autoRho.pdf'),useDingbats=FALSE)
  hist(tmp$fit$dd$rhoEst[tmp$fit$dd$rhoEst<1],
       n=100,
       xlab='Contamination fraction estimate',
       ylab='Frequency',
       main=''
       )
  lines(rhoProbes,tmp$fit$posterior*10)
  abline(v=tmp$metaData$rho[1],col='red')
dev.off()
####################
# Correcting counts
scPBMC = calculateContaminationFraction(scPBMC,useToEst=useToEst,nonExpressedGeneList=list(IG=igGenes))
out = adjustCounts(scPBMC,setNames(PBMC_DR$Cluster,rownames(PBMC_DR)),verbose=2)
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
scMix = SoupChannel(ideal,ideal[,cells],dataType='DropSeq',calcSoupProfile=FALSE)
scMix = estimateSoup(scMix,soupRange=c(10,100),keepDroplets=TRUE)
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
scMix = SoupChannel(ideal,ideal[,cells],dataType='10X',calcSoupProfile=FALSE)
scMix = estimateSoup(scMix,keepDroplets=TRUE)
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
  scMix = autoEstCont(scMix)
  #Define the genes to use for estimation
  geneSets = list(hg=scMix$hgGenes,
                  mm=scMix$mmGenes)
  #We could do this manually, but the automated thing gets it 100% right
  useToEst = estimateNonExpressingCells(scMix,nonExpressedGeneList=geneSets)
  #Calculate a global rho
  scMix$globalRho = calculateContaminationFraction(scMix,nonExpressedGeneList=geneSets,useToEst=useToEst)$metaData$rho[1]
  #And now the fine-grained estimation.  We have enough power we don't need to use STAN.
  #Human cells
  tmp = colSums(scMix$toc[geneSets$mm,scMix$isHuman])/(colSums(scMix$toc[,scMix$isHuman])*sum(scMix$soupProfile[geneSets$mm,'est']))
  scMix$metaData[scMix$isHuman,'rho']=tmp
  #Mouse cells
  tmp = colSums(scMix$toc[geneSets$hg,!scMix$isHuman])/(colSums(scMix$toc[,!scMix$isHuman])*sum(scMix$soupProfile[geneSets$hg,'est']))
  scMix$metaData[!scMix$isHuman,'rho']=tmp
  #scMix = calculateContaminationFraction(scMix,nonExpressedGeneList=geneSets,useToEst=useToEst,verbose=TRUE,cellSpecificEstimates=TRUE,stanParams= list(chains = 1, warmup = 8000, iter = 48000, cores = 1),forceAccept=TRUE)
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
#####################
# Species retention 
#Get the corrected counts
df = list()
for(tNom in names(scMixes)){
  scMix = scMixes[[tNom]]
  #Use global estimate
  scMix = setContaminationFraction(scMix,scMix$globalRho)
  out = adjustCounts(scMix,verbose=2)
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
###########################################
# Compare effective rho with estimated one
df = list()
for(tNom in names(scMixes)){
  scMix = scMixes[[tNom]]
  trueRho = scMix$metaData$rho
  out = adjustCounts(scMix,verbose=2)
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



#############
# PBMC data #
#############
#Demonstration of what genes to use
scPBMC = load10X(fp_pbmc)
PBMC_DR = PBMC_metaData
scPBMC$metaData = cbind(scPBMC$metaData,PBMC_DR)
scPBMC$DR=c('RD1','RD2')
colnames(scPBMC$metaData) = gsub('Cluster','clusters',colnames(scPBMC$metaData))
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
#######################
# Show DE of key genes
useToEst = estimateNonExpressingCells(scPBMC,geneSets['IG'])
scPBMC = autoEstCont(scPBMC)
scPBMC = calculateContaminationFraction(scPBMC,geneSets['IG'],useToEst)
out = adjustCounts(scPBMC,verbose=2)
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
##############################################################################
# Show that on average cleaning up soup increases the logFC between clusters.
#Full cross-mapping list
cMap = apply(clMap,1,function(e) colnames(clMap)[e>0])
#Use the good clustering the default one in both maps
good@meta.data$cMap = good@active.ident
bad@meta.data$cMap = good@meta.data[rownames(bad@meta.data),'seurat_clusters']
bad = SetIdent(bad,rownames(bad@meta.data),bad@meta.data$cMap)
#Get markers and compare.  This is slow, so save after it's done
tgtFile=file.path(plotDir,'PBMC_clusterMarkers.RDS')
if(!file.exists(tgtFile)){
  gMarks = FindAllMarkers(good)
  bMarks = FindAllMarkers(bad)
  #saveRDS(list(gMarks=gMarks,bMarks=bMarks),tgtFile)
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
#####################
# Genes to consider
gg = plotMarkerDistribution(scPBMC) +
  scale_size(range=c(0.5,2))
doPlot = function(e) plot(gg)
makePlots(file.path(plotDir,'guessPBMCgenes'),doPlot,width=7,height=4,pointsize=8)




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
  tmp = tryCatch(autoEstCont(scTum[[i]]),
                 error=function(e) NULL,
                 warning=function(e) NULL)
  if(is.null(tmp)){
    tmp = tryCatch(autoEstCont(scTum[[i]],tfidfMin=0.8),
                   error=function(e) NULL,
                   warning = function(e) NULL)
    if(is.null(tmp)){
      tmp =  tryCatch(autoEstCont(scTum[[i]],tfidfMin=0.6),
                     error=function(e) NULL,
                     warning = function(e) NULL) 
      if(is.null(tmp)){
        tmp=scTum[[i]]
        tmp$fit=list(rhoEst=NA)
        warning("Can't auto-calculate contamination.")
      }
    }
  }
  tmp$autoEst$rhoEst=tmp$fit$rhoEst
  scTum[[i]]=tmp
  scTum[[i]] = calculateContaminationFraction(scTum[[i]],nonExpressedGeneList=geneList[tgtSet],useToEst=ute)
  message(sprintf('Auto %g Marker %g',scTum[[i]]$autoEst$rhoEst,scTum[[i]]$metaData$rho[1]))
}
#Make a comparison plot
auto = sapply(scTum,function(e) e$autoEst$rhoEst)
man = sapply(scTum,function(e) e$metaData$rho[1])
pdf(file.path(plotDir,'autoComparison.pdf'),width=7,height=7)
  plot(auto,man,
       ylab='Contamination fraction (manual)',
       xlab='Contamination fraction (auto)',
       frame.plot=TRUE,
       pch=19,
       cex=0.5
       )
  abline(0,1,lty=2)
  cc = cor(auto[!is.na(auto)],man[!is.na(auto)])
  text(x=0.03,y=0.45,bquote(R^2 == .(format(cc,digits=2))))
dev.off()
#Make a spot the difference
dirty = lapply(scTum,function(e) e$toc)
cleaned = lapply(scTum,adjustCounts,verbose=2)
for(i in seq(nrow(cMani))){
  colnames(cleaned[[i]]) = paste0(cMani$Label[i],'___',colnames(cleaned[[i]]))
  colnames(dirty[[i]]) = paste0(cMani$Label[i],'___',colnames(dirty[[i]]))
}
dirty = do.call(cbind,dirty)
cleaned = do.call(cbind,cleaned)
dirty = standard10X(dirty,res=1.0)
cleaned = standard10X(cleaned,res=1.0)
#saveRDS(list(dirty,cleaned),file.path(plotDir,'Tumour_seurat_init.RDS'))
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
#Get markers and compare
tgtFile=file.path(plotDir,'Tumour_clusterMarkers.RDS')
if(!file.exists(tgtFile)){
  gMarks = FindAllMarkers(cleaned)
  bMarks = FindAllMarkers(dirty)
  #saveRDS(list(gMarks=gMarks,bMarks=bMarks),tgtFile)
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
doPlot = function(e) plot(gg)
makePlots(file.path(plotDir,'Tumour_egHB'),doPlot,splitRaster=FALSE,width=7,height=7,pointsize=8)






#########
# Liver #
#########
#Load meta-data
mDat = read.table('Data/FetalLiverMetaData.tsv',sep='\t',header=TRUE)
#Samples to load
srcs = list.files('Data/FetalLiver')
srcs = setNames(file.path('Data/FetalLiver',srcs),srcs)
scLiver = list()
out = list()
for(i in seq_along(srcs)){
#for(i in seq_along(dd$Path)){
  message(sprintf("Loading data for 10X channel %s from %s",names(srcs)[i], srcs[i]))
  sc = load10X(srcs[i])
  if(ncol(sc$metaData)==1){
    message(sprintf("No clustering information found, performing basic clustering on channel %d",i))
    tmp = standard10X(sc$toc)
    sc = setClusters(sc,tmp@meta.data$seurat_clusters)
  }
  #These ones have too little data to infer contamination fraction well.  Default to the average of 10%
  if(i %in% c(1,2,3,4)){
    sc = setContaminationFraction(sc,0.1)
  }else{
    sc = autoEstCont(sc)
  }
  tmp = adjustCounts(sc)
  colnames(tmp) = paste0(names(srcs)[i],'_',colnames(tmp))
  scLiver[[i]] = sc
  out[[i]] = tmp
}
dirty = do.call(cbind,lapply(scLiver,function(e) e$toc))
colnames(dirty) = paste0(rep(names(srcs),sapply(scLiver,function(e) ncol(e$toc))),'_',colnames(dirty))
clean = do.call(cbind,out)
#Keep only those in mDat
dirty = dirty[,mDat$cellID]
clean = clean[,mDat$cellID]
#Define ordering of types
cTypes = c('Early Erythroid','Mid  Erythroid','Late Erythroid','EI macrophage','VCAM1+ EI macrophage','Hepatocyte','Fibroblast','Endothelial cell','HSC/MPP','pre-B cell','Pre pro B cell','pro-B cell','B cell','MEMP','Mast cell','Megakaryocyte','Neutrophil-myeloid progenitor','Monocyte-DC precursor','pDC precursor','Monocyte','DC1','DC2','Mono-Mac','Kupffer Cell','Mono-NK','ILC progenitor','Early lymphoid/T lymphocyte','NK')
############
# UMAP plot
pdf(file.path(plotDir,'liverUMAP.pdf'),width=7,height=7)
gg = Seurat:::SingleDimPlot(mDat[,c('UMAP1','UMAP2','cellType')],c(1,2),col.by = 'cellType',label=TRUE,repel=TRUE)+guides(colour=FALSE)
plot(gg)
dev.off()
#############################
# Expression comparison plot
doPlot = function(mode=c('All','Raster','Text'),tgtGns='HBB',geneLab=tgtGns,ymax=7){
  mode = match.arg(mode)
  dd = mDat
  dd$HBB = log(1+1e4*colSums(clean[tgtGns,,drop=FALSE])/colSums(clean))
  dd$HBB_old = log(1+1e4*colSums(dirty[tgtGns,,drop=FALSE])/colSums(dirty))
  dd2 = data.frame(src=rep(c('Corrected','Uncorrected'),each=nrow(dd)),
                   exp = c(dd$HBB,dd$HBB_old),
                   type = rep(dd$cellType,times=2))
  #Define ordering
  #cTypes = unique(dd$cell.labels)
  par(mar=c(8,5,0,0))
  plot(0,
       xlim=c(0,2*length(cTypes)+1),
       ylim=c(0,ymax),
       las=1,
       type='n',
       xaxt='n',
       yaxt=ifelse(mode=='Raster','n','s'),
       xlab='',
       ylab=ifelse(mode=='Raster','',sprintf('log(%s expression)',geneLab)),
       frame.plot=FALSE
       )
  if(mode!='Raster')
    axis(1,at=2*seq_along(cTypes)-1,label=cTypes,las=2)
  #Add boxplots one at at time
  if(mode!='Text'){
    for(i in seq_along(cTypes)){
      mid = i*2-1
      w = which(dd2$type==cTypes[i])
      y = dd2$exp[w]
      xBase = mid+0.3*ifelse(dd2$src[w]=='Corrected',-1,1)
      w2 = dd2$src[w]=='Corrected'
      xx = max(min(1,200/(length(w2))),0.005)
      #message(sprintf("Cell type %s has %d types and alpha of %g",cTypes[i],length(w2),xx))
      points(xBase+runif(length(xBase),-.25,.25),
             y,
             #col=plotrix::smoothColors('black',alpha=255*max(min(1,.1/log(length(w2))),0.01)),
             col=plotrix::smoothColors('black',alpha=255*xx),
             pch=19,
             cex=0.1
             )
      boxplot(y[w2],
              col=NA,
              at=xBase[w2][1],
              boxwex=1.0,
              border='black',
              outline=FALSE,
              add=TRUE,
              xaxt='n',
              yaxt='n',
              bty='n',
              frame.plot=FALSE
              )
      boxplot(y[!w2],
              col=NA,
              at=xBase[!w2][1],
              boxwex=1.0,
              border='black',
              outline=FALSE,
              add=TRUE,
              xaxt='n',
              yaxt='n',
              bty='n',
              frame.plot=FALSE
              )
    }
  }
}
tmpFun = function(mode='All') {return(doPlot(mode,'HBB'))}
makePlots(file.path(plotDir,'liverComparison'),tmpFun,splitRaster=TRUE,width=9,height=7)
tmpFun = function(mode='All') {return(doPlot(mode,'GYPA',ymax=5))}
makePlots(file.path(plotDir,'liverComparisonGYPA'),tmpFun,splitRaster=TRUE,width=9,height=7)


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
  if(rat<1){
    clFrac = rbinom(length(clFrac),clFrac,rat)
  }else{
    spFrac = rbinom(length(spFrac),spFrac,1/rat)
  }
  cc = cor(spFrac,clFrac)
  return(cc)
}
#Apply it to everything
ccs = c(sapply(scMixes,getCor),getCor(scPBMC),sapply(scTum,getCor),sapply(scLiver,getCor))






