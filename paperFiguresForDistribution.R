##############
#************#
#* Preamble *#
#************#
##############

#################
# Standard Init #
#################

#For reproducability
set.seed(123)
#NOTE: This is the essential part of my usual .Rprofile
options(menu.graphics=FALSE)
options(stringsAsFactors=FALSE)


#############
# Libraries #
#############

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

##########
# Params #
##########

#Working directory
workDir = '~/trueHome/scratch/SoupX/'
#Where to store the results
plotDir='Output'
#Recalculate everything from scratch even if we have the data?
forceRecalc=FALSE
#Manifest with information about Kidney Samples
fp_cMani = 'KidneyTumour/Table5.tsv'
#Path for ideal experiment
fp_idealHg = 'Mixture/hg19'
fp_idealMm = 'Mixture/mm10'
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
rasterRes = 300

###########################
#*************************#
#* Load and process data *#
#*************************#
###########################

####################
# Load shared data #
####################

setwd(workDir)
#Move data to directory
cMani = read.table(fp_cMani,sep='\t',header=TRUE)
cMani$Path = file.path('KidneyTumour',cMani$Label)
#Keep just the tumour/NR channels
cMani = cMani[cMani$TissueDiseaseState!='Normal',]
#Exclude RCC3, it's full of failure
cMani = cMani[cMani$Experiment!='RCC3',]
#Load the symbol table
symbs = read.table(file.path(cMani$Path[1],'raw_gene_bc_matrices','GRCh38','genes.tsv'),sep='\t',header=FALSE)
symbs$merged = paste0(symbs$V2,'_',symbs$V1)
#Convert markers to gene names
rhoMarkers = lapply(markers,grep,symbs$merged,value=TRUE)
#Simplify sample names
N1 = as.integer(gsub('.*_ldc_([0-9])_([0-9])','\\1',cMani$Label))
N2 = as.integer(gsub('.*_ldc_([0-9])_([0-9])','\\2',cMani$Label))
#This looks a bit magical, but just constructs a unique "replicate number" for each
simpNames = setNames(paste0(gsub('_.*','',cMani$Label),ifelse(grepl('_R_',cMani$Label),'Rest_','_'),(N1-1)*2+N2),cMani$Label)
superSimpNames = setNames(gsub('Rest','',paste0(gsub('_.*','',simpNames),ifelse(grepl('RCC._[34]',simpNames),'_2',''))),names(simpNames))


######################
# Process with SoupX #
######################

#Create or load the ideal data
if(forceRecalc | !file.exists(file.path(plotDir,'ideal_scl.RDS'))){
  #Have to manually construct the scl object
  hg = Read10X(fp_idealHg)
  mm = Read10X(fp_idealMm)
  #Combine
  ideal = rbind(hg,mm)
  #Use cut-off to select cells
  cells = colnames(ideal)[colSums(ideal) >= 1000 & colSums(ideal>0) >= 500]
  #Now construct the soupchannels
  sc = SoupChannel(ideal,ideal[,cells],'mixture',dataType='10X',keepDroplets=TRUE)
  #And construct the scl object
  sclIdeal = SoupChannelList(list(sc))
  #Add a global table of droplets
  sclIdeal$tod = do.call(cbind,lapply(sclIdeal$channels,function(e) e$tod))
  sclIdeal$hgCnts = colSums(sclIdeal$tod[grep('^hg19_',rownames(sclIdeal$tod)),])
  sclIdeal$mmCnts = colSums(sclIdeal$tod[grep('^mm10_',rownames(sclIdeal$tod)),])
  #Calculate the contamination
  isHuman = sclIdeal$hgCnts[colnames(sclIdeal$toc)] > sclIdeal$mmCnts[colnames(sclIdeal$toc)]
  isDoublet = pmin(sclIdeal$hgCnts[colnames(sclIdeal$toc)], sclIdeal$mmCnts[colnames(sclIdeal$toc)])>=1000
  sclIdeal$isHuman = isHuman
  sclIdeal$isDoublet = isDoublet
  #The manual way
  sclIdeal$channels$mixture$candidateNonExpressedGenes = data.frame(row.names=rownames(sclIdeal$tod),
                                                                    mu=rep(NA,nrow(sclIdeal$tod)),
                                                                    phi=NA,
                                                                    isCandidate=TRUE)
  #Create masks
  masks = data.frame(human=grepl('^mm10_',rownames(sclIdeal$tod)),
                     mouse=grepl('^hg19_',rownames(sclIdeal$tod)),
                     row.names=rownames(sclIdeal$tod))
  masks = as.matrix(masks[,ifelse(isHuman,'human','mouse')])
  colnames(masks) = colnames(sclIdeal$toc)
  sclIdeal$channels$mixture$unexpressedMatrix = masks
  sclIdeal = findUnexpressedGenesByCell(sclIdeal,knownUnexpressedMatrix = masks)
  sclIdeal = calculateContaminationFraction(sclIdeal)
  sclIdeal = interpolateCellContamination(sclIdeal)
  #Store the gold standard
  sclIdeal$channels$mixture$metaData$rhosGold = sclIdeal$channels$mixture$metaData$rhos
  sclIdeal$channels$mixture$rhoGroupedGold = sclIdeal$channels$mixture$rhoGrouped
  sclIdeal$channels$mixture$unexpressedMatrixGold = masks
  #Process it the standard way and using manual intervention
  sclIdeal = findCandidateSoupGenes(sclIdeal)
  sclIdeal = findUnexpressedGenesByCell(sclIdeal,maxContamination=0.2)
  sclIdeal = calculateContaminationFraction(sclIdeal)
  sclIdeal = interpolateCellContamination(sclIdeal)
  #The manual way
  #sclIdeal = calculateContaminationFraction(sclIdeal,'mixture',
  #                                          list(mm=grep('^mm10_',rownames(sclIdeal$toc),value=TRUE),
  #                                               hg=grep('^hg19_',rownames(sclIdeal$toc),value=TRUE)),
  #                                          t(data.frame(hg=isHuman & !isDoublet,
  #                                                     mm=!isHuman & !isDoublet)))
  ##Infer the cell specific contamination
  #sclIdeal = interpolateCellContamination(sclIdeal,'mixture')
  #Clean it up
  sclIdeal = strainCells(sclIdeal)
  sclIdeal = adjustCounts(sclIdeal)
  saveRDS(sclIdeal,file.path(plotDir,'ideal_scl.RDS'))
}else{
  sclIdeal = readRDS(file.path(plotDir,'ideal_scl.RDS'))
}
#Create or load the scl object
if(!forceRecalc & file.exists(file.path(plotDir,'scl.RDS'))){
  sclReal= readRDS(file.path(plotDir,'scl.RDS'))
}else{
  sclReal = load10X(cMani$Path)
  #Work out which genes to use in each channel
  igGenes = c('IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHD','IGHE','IGHM',
              'IGLC1','IGLC2','IGLC3','IGLC4','IGLC5','IGLC6','IGLC7',
              'IGKC')
  hgGenes = c('HBA1','HBA2','HBB','HBD','HBE1','HBG1','HBG2','HBM','HBQ1','HBZ')
  mastGenes = c('TPSB2','TPSAB1')
  colGenes = c('COL1A1','COL1A2','COL3A1') 
  s100Genes = c('S100A8','S100A9')
  #The modern way
  #Manually specify which genes to consider
  sclReal = findCandidateSoupGenes(sclReal,includeGenes = c(igGenes,hgGenes,mastGenes,colGenes,s100Genes),nSoupGenes=0)
  sclReal = findUnexpressedGenesByCell(sclReal,maxContamination=1.0,pThr=0.01)
  sclReal = calculateContaminationFraction(sclReal)
  sclReal = interpolateCellContamination(sclReal)
  #Do additional tweaks on the unexpressed gene list for gold standard
  estGenes = list('Channel1'=list(HG=hgGenes,MAST=mastGenes),
                  'Channel2'=list(HG=hgGenes,MAST=mastGenes),
                  'Channel3'=list(HG=hgGenes),
                  'Channel4'=list(HG=hgGenes),
                  'Channel5'=list(HG=hgGenes),
                  'Channel6'=list(HG=hgGenes),
                  'Channel7'=list(HG=hgGenes),
                  'Channel8'=list(IG=igGenes),
                  'Channel9'=list(HG=hgGenes),
                  'Channel10'=list(IG=igGenes),
                  'Channel11'=list(IG=igGenes),
                  'Channel12'=list(IG=igGenes),
                  'Channel13'=list(HG=hgGenes,MAST=mastGenes),#IGFBP7 looks the only reasonable guess
                  'Channel14'=list(HG=hgGenes,MAST=mastGenes),#as above
                  'Channel15'=list(IG=igGenes,HG=hgGenes),
                  'Channel16'=list(HG=hgGenes),
                  'Channel17'=list(HG=hgGenes),
                  'Channel18'=list(HG=hgGenes),
                  'Channel19'=list(HG=hgGenes,IG=igGenes),
                  'Channel20'=list(HG=hgGenes),
                  'Channel21'=list(HG=hgGenes),
                  'Channel22'=list(HG=hgGenes),
                  'Channel23'=list(HG=hgGenes,MAST=mastGenes),
                  'Channel24'=list(HG=hgGenes,MAST=mastGenes))
  for(channel in names(estGenes)){
    tmp = sclReal$channels[[channel]]$unexpressedMatrix
    for(i in which(!(rownames(tmp) %in% unlist(estGenes[[channel]]))))
      tmp[i,]=FALSE
    sclReal$unexpressedMatrix= tmp
  }
  #Save the gold standard
  for(i in seq_along(sclReal$channels)){
    sclReal$channels[[i]]$metaData$rhosGold = sclReal$channels[[i]]$metaData$rhos
    sclReal$channels[[i]]$rhoGroupedGold = sclReal$channels[[i]]$rhoGrouped
    sclReal$channels[[i]]$unexpressedMatrixGold = sclReal$channels[[i]]$unexpressedMatrix
  }
  #Repeat using the standard gene selection method.  Calibrate aggressively to be sure we correct everything
  sclReal = findCandidateSoupGenes(sclReal)
  sclReal = findUnexpressedGenesByCell(sclReal,maxContamination=1.0,pThr=0.01)
  sclReal = calculateContaminationFraction(sclReal)
  sclReal = interpolateCellContamination(sclReal,interpolationMethod='fixed')
  #scl = inferNonExpressedGenes(scl)
  #pdf(file.path(plotDir,'inferredNonExpressed.pdf'),width=14)
  #for(nom in names(scl$channels)){
  #  plot(plotMarkerDistribution(scl,nom)+ggtitle(nom))
  #}
  #dev.off()
  ##These are decided by looking at plotMarkerDistribution(scl,'ChannelX')
  #estGenes = list('Channel1'=list(HG=hgGenes,MAST=mastGenes),
  #                'Channel2'=list(HG=hgGenes,MAST=mastGenes),
  #                'Channel3'=list(HG=hgGenes),
  #                'Channel4'=list(HG=hgGenes),
  #                'Channel5'=list(HG=hgGenes),
  #                'Channel6'=list(HG=hgGenes),
  #                'Channel7'=list(HG=hgGenes),
  #                'Channel8'=list(IG=igGenes),
  #                'Channel9'=list(HG=hgGenes),
  #                'Channel10'=list(IG=igGenes),
  #                'Channel11'=list(IG=igGenes),
  #                'Channel12'=list(IG=igGenes),
  #                'Channel13'=list(HG=hgGenes,MAST=mastGenes),#IGFBP7 looks the only reasonable guess
  #                'Channel14'=list(HG=hgGenes,MAST=mastGenes),#as above
  #                'Channel15'=list(IG=igGenes,HG=hgGenes),
  #                'Channel16'=list(HG=hgGenes),
  #                'Channel17'=list(HG=hgGenes),
  #                'Channel18'=list(HG=hgGenes),
  #                'Channel19'=list(HG=hgGenes,IG=igGenes),
  #                'Channel20'=list(HG=hgGenes),
  #                'Channel21'=list(HG=hgGenes),
  #                'Channel22'=list(HG=hgGenes),
  #                'Channel23'=list(HG=hgGenes,MAST=mastGenes),
  #                'Channel24'=list(HG=hgGenes,MAST=mastGenes))
  ##Now calculate the channel contamination
  #pdf(file.path(plotDir,'contaminationEstimates.pdf'),width=14)
  #for(nom in names(scl$channels)){
  #  scl = calculateContaminationFraction(scl,nom,estGenes[[nom]])
  #  plot(plotChannelContamination(scl,nom)+ggtitle(nom))
  #}
  #dev.off()
  ##Decide which ones to just use global estimate
  #useGlobal = rep(FALSE,length(scl$channels))
  #useGlobal[c(1,7,8,9,10,11,12,19,20)] = TRUE
  #names(useGlobal) = names(scl$channels)
  ##Interpolate final rhos
  #for(nom in names(scl$channels)){
  #  scl = interpolateCellContamination(scl,nom,useGlobal=useGlobal[nom])
  #}
  ##Now do both corrections
  sclReal = strainCells(sclReal)
  sclReal = adjustCounts(sclReal,pCut=0.1)
  saveRDS(sclReal,file.path(plotDir,'scl.RDS'))
}
####Make the Seurat object and process it
###fast_tsne_path = '/home/ubuntu/FIt-SNE/bin/fast_tsne'
###clusterData = function(toc,mDat=NULL,nPCs=100){
###  srat = CreateSeuratObject(toc)
###  srat = NormalizeData(srat)
###  if(!is.null(mDat))
###    srat = AddMetaData(srat,mDat)
###  srat = FindVariableGenes(srat,do.plot=FALSE)
###  srat = ScaleData(srat,genes.use=srat@var.genes)
###  srat = RunPCA(srat,pcs.compute=nPCs,do.print=FALSE)
###  srat = RunUMAP(srat,reduction.use='pca',dims.use=seq(nPCs))
###  srat = RunTSNE(srat,reduction.use='pca',dims.use=seq(nPCs),tsne.method='FIt-SNE',nthreads=4,reduction.name='FItSNE',reduction.key='FItSNE_',fast_tsne_path=fast_tsne_path)
###  srat = FindClusters(srat,reduction.type='pca',dims.use=seq(nPCs),resolution=1,print.output=FALSE)
###  srat = AddMetaData(srat,data.frame(srat@dr$umap@cell.embeddings))
###  srat = AddMetaData(srat,data.frame(srat@dr$FItSNE@cell.embeddings))
###  srat@misc$markers = quickMarkers(srat@data,srats[[i]]@meta.data$res.1,N=Inf)
###  #Calculate fraction of expression in each of the gene sets
###  for(nom in names(geneSets)){
###    m = match(geneSets[[nom]],gsub('___.*','',rownames(srat@data)))
###    m = m[!is.na(m)]
###    srat@meta.data[,nom] = colSums(srat@raw.data[rownames(srat@data)[m],])/colSums(srat@raw.data)
###  }
###  srat@scale.data=NULL
###  return(srat)
###}
scl = sclReal
nPCs = 30
perplexity= 30
if(forceRecalc | !file.exists(file.path(plotDir,'sratAdjusted.RDS'))){
  mDat = cMani[as.integer(gsub('Channel([0-9]+)___.*','\\1',colnames(scl$atoc))),]
  rownames(mDat) = colnames(scl$atoc)
  srat = CreateSeuratObject(scl$atoc)
  srat = NormalizeData(srat)
  mDat = cMani[as.integer(gsub('Channel([0-9]+)___.*','\\1',colnames(srat@data))),]
  rownames(mDat) = colnames(srat@data)
  srat = AddMetaData(srat,mDat)
  srat = FindVariableGenes(srat,do.plot=FALSE)
  srat = ScaleData(srat)
  srat = RunPCA(srat,pcs.compute=nPCs,do.print=FALSE,maxit=1000)
  srat = RunTSNE(srat,dims.use=seq(nPCs),perplexity=perplexity)
  srat = FindClusters(srat,dims.use=seq(nPCs),resolution=1,print.output=0,save.SNN=TRUE,plot.SNN=FALSE)
  srat@scale.data=NULL
  saveRDS(srat,file.path(plotDir,'sratAdjusted.RDS'))
}
#Make the Seurat object for the unedited  
if(forceRecalc | !file.exists(file.path(plotDir,'sratSoggy.RDS'))){
  srat = CreateSeuratObject(scl$toc)
  srat = NormalizeData(srat)
  mDat = cMani[as.integer(gsub('Channel([0-9]+)___.*','\\1',colnames(srat@data))),]
  rownames(mDat) = colnames(srat@data)
  srat = AddMetaData(srat,mDat)
  srat = FindVariableGenes(srat,do.plot=FALSE)
  srat = ScaleData(srat)
  srat = RunPCA(srat,pcs.compute=nPCs,do.print=FALSE,maxit=1000)
  srat = RunTSNE(srat,dims.use=seq(nPCs),perplexity=perplexity)
  srat = FindClusters(srat,dims.use=seq(nPCs),resolution=1,print.output=0,save.SNN=TRUE,plot.SNN=FALSE)
  srat@scale.data=NULL
  saveRDS(srat,file.path(plotDir,'sratSoggy.RDS'))
}
#Finally, the Seurat object for the strained
if(forceRecalc | !file.exists(file.path(plotDir,'sratStrained.RDS'))){
  srat = createCleanedSeurat(scl)
  mDat = cMani[as.integer(gsub('Channel([0-9]+)___.*','\\1',colnames(srat@data))),]
  rownames(mDat) = colnames(srat@data)
  srat = AddMetaData(srat,mDat)
  srat = FindVariableGenes(srat,do.plot=FALSE)
  srat = ScaleData(srat)
  srat = RunPCA(srat,pcs.compute=nPCs,do.print=FALSE,maxit=1000)
  srat = RunTSNE(srat,dims.use=seq(nPCs),perplexity=perplexity)
  srat = FindClusters(srat,dims.use=seq(nPCs),resolution=1,print.output=0,save.SNN=TRUE,plot.SNN=FALSE)
  srat@scale.data=NULL
  saveRDS(srat,file.path(plotDir,'sratStrained.RDS'))
}
#Load the soggy and adjusted objects
soggy = readRDS(file.path(plotDir,'sratSoggy.RDS'))
adjusted = readRDS(file.path(plotDir,'sratAdjusted.RDS'))
#Add the annotation for soggy
annMap = c('0'='T',
          '1'='CD8+ T',
          '2'='MNP',
          '3'='MNP',
          '4'='MNP',
          '5'='CD8+ T',
          '6'='T',
          '7'='ccRCC',
          '8'='ccRCC',
          '9'='Endothelium',
          '10'='ccRCC',
          '11'='CD8+ T',
          '12'='NK',
          '13'='Dead',
          '14'='NK',
          '15'='RBC',
          '16'='?',
          '17'='ccRCC',
          '18'='T',
          '19'='?',
          '20'='MAST',
          '21'='NK',
          '22'='Wilms',
          '23'='Fibroblasts',
          '24'='Endothelium',
          '25'='pRCC',
          '26'='?',
          '27'='Proliferating T',
          '28'='B',
          '29'='Fibroblasts',
          '30'='Proliferating MNP',
          '31'='Activated B',
          '32'='?')

######################
#********************#
#* Generate figures *#
#********************#
######################

######################
# Basic cluster plot #
######################
pdf(file.path(plotDir,'tsne_clusters.pdf'))
TSNEPlot(soggy,do.label=TRUE,no.legend=TRUE)
dev.off()
soggy@meta.data$annotation = annMap[as.character(soggy@meta.data$res.1)]
pdf(file.path(plotDir,'tsne_annotation.pdf'))
TSNEPlot(soggy,group.by='annotation',do.label=TRUE,no.legend=TRUE)
dev.off()


##################################
# Gene selection validation plot #
##################################
#Work out what the fraction of 
isHuman = sclIdeal$hgCnts[colnames(sclIdeal$toc)] > sclIdeal$mmCnts[colnames(sclIdeal$toc)]
isDoublet = pmin(sclIdeal$hgCnts[colnames(sclIdeal$toc)], sclIdeal$mmCnts[colnames(sclIdeal$toc)])>=1000
#Get the fraction of markers that are human and mouse
markFrac =  apply(sclIdeal$channels$mixture$unexpressedMatrix,2,function(e) c(sum(grepl('^hg19_',rownames(sclIdeal$channels$mixture$unexpressedMatrix))[e]),sum(grepl('mm10_',rownames(sclIdeal$channels$mixture$unexpressedMatrix))[e])))
rownames(markFrac) = c('human','mouse')
#Get fraction correct
trueFrac = markFrac['human',]/colSums(markFrac)
trueFrac[isHuman] = 1 -trueFrac[isHuman]
#Drop doublets
trueFrac = trueFrac[!isDoublet]
nUMIs = sclIdeal$metaData[names(trueFrac),'nUMIs']
#plot(log10(nUMIs),trueFrac)
#Now look at the difference in inferred contamination
df = sclIdeal$channels$mixture$rhoGroupedGold
df$estCoal = sclIdeal$channels$mixture$rhoGrouped$est
df$ratio = df$estCoal/df$est
#Exclude doublets
df = df[!(df$nUMIs %in% sclIdeal$metaData[names(isDoublet)[isDoublet],'nUMIs']),]
#Plot the log ratio of estimated contamination
gg_validation = ggplot(df,aes(log10(nUMIs),log10(ratio))) +
  geom_hex(bins=100) +
  geom_hline(yintercept=0,colour='red')
#And the actual values
df = melt(df[,c('est','estCoal','nUMIs')],id.vars=c('nUMIs'))
gg_validation = ggplot(df,aes(log10(nUMIs),value,colour=variable)) +
  geom_point(alpha=1/1)
#Get the counts of the genes used to inspect
inspectCounts = sclIdeal$toc[rownames(sclIdeal$channels$mixture$unexpressedMatrix),]
#Get the ratios with limits in a sensible way
getRatioWithLimits = function(scl,useCells=NULL,cellGroupings=NULL,tgtSoupCntsPerGroup=100){
  if(is.null(useCells))
    useCells=colnames(scl$toc)
  #Get the estimation matricies for each condition
  useToEstCoal = as.matrix(scl$unexpressedMatrix[,useCells])
  useToEstGold = as.matrix(scl$unexpressedMatrixGold[,useCells])
  #Determine a globally appropriate cell grouping if not supplied
  if(is.null(cellGroupings)){
    o = order(scl$metaData[useCells,'nUMIs'])
    #The expected counts per cell
    expCntsGold = colSums(useToEstGold*scl$soupProfile[rownames(useToEstGold),'est'])*scl$metaData[useCells,'nUMIs']
    expCntsCoal = colSums(useToEstCoal*scl$soupProfile[rownames(useToEstCoal),'est'])*scl$metaData[useCells,'nUMIs']
    expCnts = pmin(expCntsGold,expCntsCoal)
    #Split groups so that the expected soup count is respected
    cellGroupings = cut_interval(cumsum(expCnts[o]),length=tgtSoupCntsPerGroup)
    #Reverse ordering
    cellGroupings = cellGroupings[match(seq_along(o),o)]
    #Drop empty levels
    cellGroupings = factor(as.character(cellGroupings))
  }
  #Get total UMIs in each grouping (and global)
  nUMIs = sapply(split(scl$metaData[useCells,'nUMIs'],cellGroupings),sum)
  nUMIs['Global'] = sum(scl$metaData[useCells,'nUMIs'])
  #Estimation
  #Gold
  obsSoupCntsGold = SoupX:::estRateLims(colSums(scl$toc[rownames(useToEstGold),useCells]*useToEstGold),scl$metaData[useCells,'nUMIs'])*scl$metaData[useCells,'nUMIs']
  expSoupCntsGold = colSums(useToEstGold*scl$soupProfile[rownames(useToEstGold),'est'],na.rm=TRUE)*scl$metaData[useCells,'nUMIs']
  globRhoGold = SoupX:::estGrouped(obsSoupCntsGold$est,expSoupCntsGold,scl$metaData[useCells,'nUMIs'],rep('Global',length(useCells)))
  groupRhosGold = SoupX:::estGrouped(obsSoupCntsGold$est,expSoupCntsGold,scl$metaData[useCells,'nUMIs'],cellGroupings)
  rhoGroupedGold = rbind(globRhoGold,groupRhosGold)
  #Coal
  obsSoupCntsCoal = SoupX:::estRateLims(colSums(scl$toc[rownames(useToEstCoal),useCells]*useToEstCoal),scl$metaData[useCells,'nUMIs'])*scl$metaData[useCells,'nUMIs']
  expSoupCntsCoal = colSums(useToEstCoal*scl$soupProfile[rownames(useToEstCoal),'est'],na.rm=TRUE)*scl$metaData[useCells,'nUMIs']
  globRhoCoal = SoupX:::estGrouped(obsSoupCntsCoal$est,expSoupCntsCoal,scl$metaData[useCells,'nUMIs'],rep('Global',length(useCells)))
  groupRhosCoal = SoupX:::estGrouped(obsSoupCntsCoal$est,expSoupCntsCoal,scl$metaData[useCells,'nUMIs'],cellGroupings)
  rhoGroupedCoal = rbind(globRhoCoal,groupRhosCoal)
  #Compare across groups
  rhoRatio = list()
  for(i in seq(nrow(rhoGroupedGold))){
    nn = nUMIs[rownames(rhoGroupedGold)[i]]
    #Sample rho in gold
    rhoGold = rbinom(1000,nn,rhoGroupedGold[i,'obsSoupCnts']/nn)/rhoGroupedGold[i,'expSoupCnts']
    rhoCoal = rbinom(1000,nn,rhoGroupedCoal[i,'obsSoupCnts']/nn)/rhoGroupedCoal[i,'expSoupCnts']
    rhoRatio[[rownames(rhoGroupedGold)[i]]] = rhoCoal/rhoGold
  }
  #Get median and confidence
  df = sapply(lapply(rhoRatio,log2),quantile,c(0.025,0.5,0.975),na.rm=TRUE)
  df = data.frame(est=df[2,],
                  lower=df[1,],
                  upper=df[3,],
                  rhoGold = rhoGroupedGold$est,
                  lowerGold = rhoGroupedGold$lower,
                  upperGold = rhoGroupedGold$upper,
                  rhoTest = rhoGroupedCoal$est,
                  lowerTest = rhoGroupedCoal$lower,
                  upperTest = rhoGroupedCoal$upper,
                  row.names=rownames(rhoGroupedGold),
                  nUMIs=rhoGroupedGold$nUMIs)
  gg = ggplot(df,aes(log10(nUMIs),est)) +
    geom_errorbar(aes(ymin=lower,ymax=upper),size=0.1) +
    geom_hex(bins=100) +
    geom_hline(yintercept=0,col='red',linetype=2)
  return(gg)
}
plots = lapply(sclReal$channels,getRatioWithLimits)
#Make a plot where the two are shown on an absolute scale
df = lapply(plots,function(e) e$data)
df = do.call(rbind,df)
df$Channel = gsub('\\..*','',rownames(df))
ddf = melt(df[,c('rhoGold','rhoTest','nUMIs','Channel')],id.vars=c('Channel','nUMIs'),value.name='rho')
tmp = melt(df[,c('lowerGold','lowerTest','nUMIs','Channel')],id.vars=c('Channel','nUMIs'),value.name='rho')
ddf$lower = tmp$rho
tmp = melt(df[,c('upperGold','upperTest','nUMIs','Channel')],id.vars=c('Channel','nUMIs'),value.name='rho')
ddf$upper = tmp$rho
ddf$variable = ifelse(ddf$variable=='rhoGold',
                      'Manual',
                      'Automatic'
                      )
##Add the ideal to the real
#tmp = sclIdeal$channels$mixture$rhoGroupedGold
#m = match(tmp$nUMIs,sclIdeal$metaData[names(sclIdeal$isDoublet),'nUMIs'])
#tmp = tmp[which(!sclIdeal$isDoublet[m]),]
#ddf = rbind(ddf,data.frame(Channel='Discovery',nUMIs=tmp$nUMIs,variable='Manual',rho=tmp$est,lower=tmp$lower,upper=tmp$upper))
#tmp = sclIdeal$channels$mixture$rhoGrouped
#m = match(tmp$nUMIs,sclIdeal$metaData[names(sclIdeal$isDoublet),'nUMIs'])
#tmp = tmp[which(!sclIdeal$isDoublet[m]),]
#ddf = rbind(ddf,data.frame(Channel='Discovery',nUMIs=tmp$nUMIs,variable='Automatic',rho=tmp$est,lower=tmp$lower,upper=tmp$upper))
gg_sel_all = ggplot(ddf,aes(log10(nUMIs),rho,colour=variable)) +
  geom_point(size=0.1) +
  geom_errorbar(aes(ymin=lower,ymax=upper)) + 
  geom_smooth() + 
  ylim(0,1) +
  ylab('Contamination fraction') +
  xlab('log10(nUMIs)') + 
  labs(colour='Gene Selection Method') +
  facet_wrap(~Channel)
pdf(file.path(plotDir,'AllSelectionPlots.pdf'),width=14,height=14)
plot(gg_sel_all)
dev.off()
png(file.path(plotDir,'AllSelectionPlots.png'),width=14,height=14,units = 'in',res=rasterRes)
plot(gg_sel_all)
dev.off()
#Same but just for the ideal
tmp = sclIdeal$channels$mixture$rhoGroupedGold
m = match(tmp$nUMIs,sclIdeal$metaData[names(sclIdeal$isDoublet),'nUMIs'])
tmp = tmp[which(!sclIdeal$isDoublet[m]),]
ddf = data.frame(Channel='Discovery',nUMIs=tmp$nUMIs,variable='Manual',rho=tmp$est,lower=tmp$lower,upper=tmp$upper)
tmp = sclIdeal$channels$mixture$rhoGrouped
m = match(tmp$nUMIs,sclIdeal$metaData[names(sclIdeal$isDoublet),'nUMIs'])
tmp = tmp[which(!sclIdeal$isDoublet[m]),]
ddf = rbind(ddf,data.frame(Channel='Discovery',nUMIs=tmp$nUMIs,variable='Automatic',rho=tmp$est,lower=tmp$lower,upper=tmp$upper))
gg_sel_all = ggplot(ddf,aes(log10(nUMIs),rho,colour=variable)) +
  geom_point(size=0.1) +
  geom_errorbar(aes(ymin=lower,ymax=upper)) + 
  geom_smooth() + 
  ylim(0,1) +
  ylab('Contamination fraction') +
  xlab('log10(nUMIs)') + 
  labs(colour='Gene Selection Method')
pdf(file.path(plotDir,'IdealSelectionPlots.pdf'),width=14,height=14)
plot(gg_sel_all)
dev.off()
png(file.path(plotDir,'IdealSelectionPlots.png'),width=14,height=14,units='in',res=rasterRes)
plot(gg_sel_all)
dev.off()

################
# Lobster plot #
################
df = data.frame(hg=sclIdeal$hgCnts,mm=sclIdeal$mmCnts,isCell=colnames(sclIdeal$tod) %in% colnames(sclIdeal$toc))
#df = data.frame(hg=hg_umiCounts,mm=mm_umiCounts,qval=qval)
df$M = log10(0.5*(df$hg+df$mm))
df$A = log10(df$hg/df$mm)
ext = ceiling(1.01*max(abs(df$A[is.finite(df$A)])))
df$A[df$A == -Inf] = -ext
df$A[df$A == Inf] = ext
#Pick out the doublets
df$isDoublet = pmin(df$hg,df$mm)>=1000
#Drop the doublets
df = df[!df$isDoublet,]
#Estimate global rho for two populations
w = which(df$isCell & !df$isDoublet & df$hg>df$mm)
rho_hg = (sum(df[w,'mm'])/sum(df[w,c('hg','mm')]))/sum(sclIdeal$soupMatrix[grep('^mm10_',rownames(sclIdeal$soupMatrix))])
w = which(df$isCell & !df$isDoublet & df$hg<df$mm)
rho_mm = (sum(df[w,'hg'])/sum(df[w,c('hg','mm')]))/sum(sclIdeal$soupMatrix[grep('^hg19_',rownames(sclIdeal$soupMatrix))])
#Make the plot
gg_lobster = ggplot(df,aes(M,A,colour=isCell)) +
  geom_point(size=1.0) +
  #geom_hline(yintercept=c(-log10(rho_hg/2),log10(rho_mm/2))) +
  labs(colour='Is a cell?') +
  xlab('log10(Average UMIs)') +
  ylab('log10(Human UMIs/Mouse UMIs)') +
  scale_colour_manual(values=colPal[1:2]) +
  guides(colour=FALSE) + 
  theme_bw(base_size=18)
pdf(file.path(plotDir,'lobster.pdf'))
plot(gg_lobster)
dev.off()
png(file.path(plotDir,'lobster.png'),width=7,height=7,units='in',res=rasterRes)
plot(gg_lobster)
dev.off()


##################################
# Soup/Cell correlation in ideal #
##################################
#Get just the human counts in human cells and mouse counts in mouse cells
dat = sclIdeal$toc
isHuman = sclIdeal$hgCnts[colnames(dat)] > sclIdeal$mmCnts[colnames(dat)]
isDoublet = pmin(sclIdeal$hgCnts[colnames(dat)], sclIdeal$mmCnts[colnames(dat)])>=1000
hgSoup = sclIdeal$soupMatrix[grep('^hg19_',rownames(sclIdeal$soupMatrix))]
mmSoup = sclIdeal$soupMatrix[grep('^mm10_',rownames(sclIdeal$soupMatrix))]
#dat = ideal[,idealCells]
#isHuman = hgCnts[idealCells]>mmCnts[idealCells]
#isDoublet = pmin(hgCnts[idealCells],mmCnts[idealCells])>=1000
hgCell = (dat[grep('^hg19_',rownames(dat)),isHuman & !isDoublet])
hgCell = Matrix::rowSums(hgCell)/sum(hgCell)
#hgSoup = idealSoup$est[grep('hg',rownames(idealSoup))]
hgSoup = hgSoup/sum(hgSoup)
mmCell = (dat[grep('^mm10_',rownames(dat)),!isHuman & !isDoublet])
mmCell = Matrix::rowSums(mmCell)/sum(mmCell)
#mmSoup = idealSoup$est[grep('mm',rownames(idealSoup))]
mmSoup = mmSoup/sum(mmSoup)
df = data.frame(cellFrac = c(hgCell,mmCell),
                soupFrac = c(hgSoup,mmSoup),
                species = rep(c('human','mouse'),c(length(hgCell),length(mmCell)))
                )
labs = c(sprintf('R^2 == %.02f (Human)',cor(hgCell,hgSoup)),
         sprintf('R^2 == %.02f (Mouse)',cor(mmCell,mmSoup)))
#Just a straight correlation plot
gg_idealCor = ggplot(df,aes(log10(cellFrac),log10(soupFrac))) +
  geom_point(aes(colour=species),size=1.0,alpha=1) +
  geom_abline(intercept=0,slope=1) +
  annotate('text',x=-6,y=-2,label= labs[1],parse=TRUE) +
  annotate('text',x=-6,y=-2.3,label= labs[2],parse=TRUE) +
  xlab('log10(cell expression)') +
  ylab('log10(soup expression)') +
  scale_colour_manual(values=colPal[1:2]) + 
  guides(colour=FALSE) +
  theme_bw(base_size=18)
#The MA plot instead 
#gg_idealCor = ggplot(df,aes(log10(0.5*(cellFrac+soupFrac)),log10(cellFrac/soupFrac))) +
#  geom_point(aes(colour=species)) +
#  geom_hline(yintercept=0) +
#  annotate('text',x=-2,y=-2,label= labs[1],parse=TRUE) +
#  annotate('text',x=-2,y=-2.3,label= labs[2],parse=TRUE) +
#  xlab('log10(Average expression)') +
#  ylab('log10(cell/soup)') +
#  xlim(-6,NA) +
#  guides(colour=FALSE) +
#  theme_bw(base_size=18)
#  #ggtitle(sprintf('Gene correlation between cell and soup\n R=%g,%g for human,mouse',cor(hgCell,hgSoup),cor(mmCell,mmSoup)))
pdf(file.path(plotDir,'idealCorrelation.pdf'))
plot(gg_idealCor)
dev.off()
png(file.path(plotDir,'idealCorrelation.png'),width=7,height=7,units='in',res=rasterRes)
plot(gg_idealCor)
dev.off()


#####################################
# Correlation for non-ideal samples #
#####################################
df = list()
for(i in seq_along(scl$channels)){
  #Get the soup expression profile
  df[[i]] = data.frame(Cells = rowSums(scl$channels[[i]]$toc)/sum(scl$channels[[i]]$metaData[colnames(scl$channels[[i]]$toc),'nUMIs']),
                       Soup = scl$soupMatrix[,i],
                       Sample = cMani$Label[i],
                       ShortName = simpNames[cMani$Label[i]]
                       )
  #Get the correlation co-efficient
  df[[i]]$Label = sprintf('%s (%.02f)',df[[i]]$ShortName,cor(df[[i]]$Cells,df[[i]]$Soup))
}
df = do.call(rbind,df)
#Create the plot
gg_cor = ggplot(df,aes(log10(Cells),log10(Soup))) +
  geom_abline(intercept=0,slope=1) +
  geom_point(size=0.5,alpha=1/2) +
  facet_wrap(~Label)
pdf(file.path(plotDir,'allCorrelations.pdf'),width=14,height=14)
plot(gg_cor)
dev.off()
png(file.path(plotDir,'allCorrelations.png'),width=14,height=14,units='in',res=rasterRes)
plot(gg_cor)
dev.off()


#########################################
# Optimal droplets for soup calculation #
#########################################
#Define "true background"
dat = sclIdeal$toc
isHuman = sclIdeal$hgCnts[colnames(dat)] > sclIdeal$mmCnts[colnames(dat)]
isDoublet = pmin(sclIdeal$hgCnts[colnames(dat)], sclIdeal$mmCnts[colnames(dat)])>=1000
hgBG = rowSums(dat[grep('mm10_',rownames(dat)),isHuman & !isDoublet & sclIdeal$metaData[colnames(sclIdeal$toc),'nUMIs'] > 1000])
hgBG = hgBG/sum(hgBG)
mmBG = rowSums(dat[grep('hg19_',rownames(dat)),!isHuman & !isDoublet & sclIdeal$metaData[colnames(sclIdeal$toc),'nUMIs'] > 1000])
mmBG = mmBG/sum(mmBG)
#Now calculate soup in each bin from 0 to 100
out=list()
for(i in seq(100)){
  w = which(sclIdeal$metaData$nUMIs==i)
  hgSoup = rowSums(sclIdeal$tod[grep('^mm10_',rownames(dat)),w])
  hgSoup = hgSoup/sum(hgSoup)
  mmSoup = rowSums(sclIdeal$tod[grep('^hg19_',rownames(dat)),w])
  mmSoup = mmSoup/sum(mmSoup)
  out[[i]] = data.frame(nUMIs=i,
                        cor = c(cor(hgSoup,hgBG),cor(mmSoup,mmBG)),
                        species = c('human','mouse')
                        )
}
out = do.call(rbind,out)
gg = ggplot(out,aes(nUMIs,cor,colour=species)) +
  geom_point() +
  scale_colour_manual(values=colPal[1:2]) +
  xlab('Droplet UMI count') +
  ylab('Correlation with true contamination')
pdf(file.path(plotDir,'OptimalSoupDroplets.pdf'))
plot(gg)
dev.off()
png(file.path(plotDir,'OptimalSoupDroplets.png'),width=7,height=7,units='in',res=rasterRes)
plot(gg)
dev.off()



#################################
# Soup comparison between cells #
#################################
dat = sclIdeal$toc
isHuman = sclIdeal$hgCnts[colnames(dat)] > sclIdeal$mmCnts[colnames(dat)]
isDoublet = pmin(sclIdeal$hgCnts[colnames(dat)], sclIdeal$mmCnts[colnames(dat)])>=1000
hgSoup = sclIdeal$soupMatrix[grep('^hg19_',rownames(sclIdeal$soupMatrix))]
mmSoup = sclIdeal$soupMatrix[grep('^mm10_',rownames(sclIdeal$soupMatrix))]
hgSoup = hgSoup/sum(hgSoup)
mmSoup = mmSoup/sum(mmSoup)
#Get correlation to soup in bins of cells
hgCellSoups = dat[grep('mm10_',rownames(dat)),isHuman & !isDoublet]
nCnts = Matrix::colSums(hgCellSoups)
o = order(nCnts)
hgCellSoups = hgCellSoups[,o]
nCnts = nCnts[o]
#Split into bins of roughly 1000 reads.
w = cut_interval(cumsum(nCnts),length=1000)
hgSoupCor = sapply(levels(w),function(e) cor(mmSoup,Matrix::rowSums(hgCellSoups[,w==e,drop=FALSE])/sum(hgCellSoups[,w==e])))
hgSoupCnts = sapply(split(nCnts,w),sum)
df = data.frame(Counts=hgSoupCnts,
                Correlation = hgSoupCor,
                Species='Human')
#Repeat for mouse
mmCellSoups = dat[grep('hg19_',rownames(dat)),!isHuman & !isDoublet]
nCnts = Matrix::colSums(mmCellSoups)
o = order(nCnts)
mmCellSoups = mmCellSoups[,o]
nCnts = nCnts[o]
#Split into bins of roughly 1000 reads.
w = cut_interval(cumsum(nCnts),length=1000)
mmSoupCor = sapply(levels(w),function(e) cor(hgSoup,Matrix::rowSums(mmCellSoups[,w==e,drop=FALSE])/sum(mmCellSoups[,w==e])))
mmSoupCnts = sapply(split(nCnts,w),sum)
df = rbind(df,
           data.frame(Counts=mmSoupCnts,
                Correlation = mmSoupCor,
                Species='Mouse')
           )
gg_cellsoup = ggplot(df,aes(Counts,Correlation)) +
  geom_point(aes(colour=Species)) +
  scale_colour_manual(values=colPal[1:2]) +
  xlab('Cell counts') +
  ylab('Correlation with soup') +
  ylim(0,1)
pdf(file.path(plotDir,'cellSoupComparison.pdf'))
plot(gg_cellsoup)
dev.off()
png(file.path(plotDir,'cellSoupComparison.png'),width=7,height=7,units='in',res=rasterRes)
plot(gg_cellsoup)
dev.off()



############################
# Spot the difference tSNE #
############################
# Make a spot the difference tSNE
df = soggy@meta.data
df$tSNE1 = soggy@dr$tsne@cell.embeddings[rownames(df),1]
df$tSNE2 = soggy@dr$tsne@cell.embeddings[rownames(df),2]
df$soup = 'Contaminated'
alt = adjusted@meta.data
alt$annotation = df$annotation
alt$tSNE1 = adjusted@dr$tsne@cell.embeddings[rownames(df),1]
alt$tSNE2 = adjusted@dr$tsne@cell.embeddings[rownames(df),2]
alt$soup = 'Cleaned'
df = rbind(df,alt[,colnames(df)])
#Adjust labels to emphasise different experiments and biological reps
df$labs = superSimpNames[df$Label]
#Set order of panels
df$soup = factor(df$soup,levels=c('Contaminated','Cleaned'))
gg_std = ggplot(df,aes(tSNE1,tSNE2,colour=labs)) +
  geom_point(size=0.3,alpha=1/2) +
  scale_colour_manual(values=colPalLarge) +
  labs(colour='Sample') + 
  guides(colour = guide_legend(override.aes = list(size=3,alpha=1))) +
  facet_grid(~soup)
pdf(file.path(plotDir,'spotTheDifference.pdf'),width=14)
plot(gg_std)
dev.off()
png(file.path(plotDir,'spotTheDifference.png'),width=14,height=7,units='in',res=rasterRes)
plot(gg_std)
dev.off()


############################
# Batch comparison heatmap #
############################
mats = list()
for(srat in list(soggy,adjusted)){
  #Get all channels
  channels = unique(srat@meta.data$Label)
  #Make a cluster membership matrix
  membTab = sapply(lapply(split(srat@meta.data$Label,srat@meta.data$res.1),unique),ltable,levels=channels)
  out = list()
  for(channel in channels){
    ngbFrac = rowSums(membTab[,srat@meta.data$res.1[srat@meta.data$Label==channel]])
    #Clean up names
    names(ngbFrac) = simpNames[names(ngbFrac)]
    out[[simpNames[channel]]] = ngbFrac/sum(srat@meta.data$Label==channel)
  }
  mats[[length(mats)+1]] = do.call(rbind,out)
}
names(mats) = c('Contamintaed','Cleaned')
#Define groups to group by
groups = ifelse(grepl('Wilms',rownames(mats[[1]])),
                'Wilms',
                ifelse(grepl('pRCC',rownames(mats[[1]])),
                       'pRCC',
                       'ccRCC'
                       )
                )
groups = factor(groups,levels=c('Wilms','ccRCC','pRCC'))
htList = Heatmap(mats[[1]],
        #column_title='Contaminated',
        row_names_side='left',
        show_heatmap_legend=FALSE,
        row_split = groups,
        row_gap = unit(2,'mm'),
        column_gap = unit(2,'mm'),
        column_split = groups,
        show_row_names = FALSE,
        show_column_names = FALSE,
        cluster_row=FALSE,
        cluster_col=FALSE) +
  Heatmap(mats[[2]],
          name='Membership\nFraction',
          show_row_names=FALSE,
          #column_title='Cleaned',
          row_split = groups,
          column_split = groups,
          row_gap = unit(2,'mm'),
          column_gap = unit(2,'mm'),
          show_column_names = FALSE,
          cluster_row=FALSE,
          cluster_col=FALSE)
pdf(file.path(plotDir,'batchComparison.pdf'),width=14,height=7)
draw(htList,
     column_title='Target Sample',
     row_title='Source Sample',
     gap = unit(2,'cm'),
     column_title_side='bottom')
dev.off()
png(file.path(plotDir,'batchComparison.png'),width=14,height=7,units='in',res=rasterRes)
draw(htList,
     column_title='Target Sample',
     row_title='Source Sample',
     gap = unit(2,'cm'),
     column_title_side='bottom')
dev.off()



######################
# MT frac comparison #
######################
mtSoggy = Matrix::colSums(exp(soggy@data[grep('^MT-.*',rownames(soggy@data)),])-1)/1e4
mtClean = Matrix::colSums(exp(adjusted@data[grep('^MT-.*',rownames(adjusted@data)),])-1)/1e4
df = data.frame(clean=mtClean,soggy=mtSoggy)
df$M = log(df$clean/df$soggy)
df$A = df$soggy
gg_mtComp = ggplot(df[df$clean<0.2 | df$soggy<0.2,],aes(soggy,clean-soggy)) +
  geom_point(size=0.5,alpha=1/2) +
  xlab('Contaminated MT frac') +
  ylab('Clean MT - Contaminated MT')
pdf(file.path(plotDir,'MTChange.pdf'))
plot(gg_mtComp)
dev.off()
png(file.path(plotDir,'MTChange.png'),width=7,height=7,units='in',res=rasterRes)
plot(gg_mtComp)
dev.off()


##################################
# Ideal data marker gene changes #
##################################
isHuman = sclIdeal$hgCnts[colnames(sclIdeal$toc)] > sclIdeal$mmCnts[colnames(sclIdeal$toc)]
isDoublet = pmin(sclIdeal$hgCnts[colnames(sclIdeal$toc)], sclIdeal$mmCnts[colnames(sclIdeal$toc)])>=1000
#Get the fraction of contamination from mouse in human
w = which(isHuman & !isDoublet)
#Looking in human cells
fracKeptBadHg = colSums(sclIdeal$atoc[grep('^mm10_',rownames(sclIdeal$atoc)),w])/colSums(sclIdeal$toc[grep('^mm10_',rownames(sclIdeal$toc)),w])
fracKeptGoodHg = colSums(sclIdeal$atoc[grep('^hg19_',rownames(sclIdeal$atoc)),w])/colSums(sclIdeal$toc[grep('^hg19_',rownames(sclIdeal$toc)),w])
#Looking in mouse cells
w = which(!isHuman & !isDoublet)
fracKeptBadMm = colSums(sclIdeal$atoc[grep('^hg19_',rownames(sclIdeal$atoc)),w])/colSums(sclIdeal$toc[grep('^hg19_',rownames(sclIdeal$toc)),w])
fracKeptGoodMm = colSums(sclIdeal$atoc[grep('^mm10_',rownames(sclIdeal$atoc)),w])/colSums(sclIdeal$toc[grep('^mm10_',rownames(sclIdeal$toc)),w])
#Convert to data frame 
df = data.frame(species = rep(c('Human','Mouse'),2*c(length(fracKeptBadHg),length(fracKeptBadMm))),
                moleculeType = rep(c('Contamination','Cellular','Contamination','Cellular'),c(length(fracKeptBadHg),length(fracKeptGoodHg),length(fracKeptBadMm),length(fracKeptGoodMm))),
                fracRetained = c((fracKeptBadHg),(fracKeptGoodHg),(fracKeptBadMm),(fracKeptGoodMm))
                )
#Set some ordering
df$moleculeType = factor(df$moleculeType,levels=c('Contamination','Cellular'))
df$species = factor(df$species,levels=c('Human','Mouse'))
#Make the plot
gg = ggplot(df,aes(moleculeType,fracRetained)) +
  geom_boxplot(aes(fill=species)) +
  scale_fill_manual(values=colPal[1:2]) +
  #guides(fill=FALSE) +
  xlab('mRNA Type') +
  ylab('Fraction retained counts after correction')
pdf(file.path(plotDir,'IdealCorrectionComparison.pdf'),width=4)
plot(gg)
dev.off()
png(file.path(plotDir,'IdealCorrectionComparison.png'),width=4,height=7,units='in',res=rasterRes)
plot(gg)
dev.off()



##################################
# Real data marker gene changes  #
##################################
cl = quickMarkers(soggy@data,soggy@meta.data$res.1,N=10)
#Get the fractional change for the markers that we care about
s = split(rownames(soggy@meta.data),soggy@meta.data$res.1)
fracChange = lapply(s,function(e){
                    a = rowSums(soggy@data[cl$gene,e,drop=FALSE]>0)
                    b = rowSums(adjusted@data[cl$gene,e,drop=FALSE]>0)
                    b/a
                })
fracChange = do.call(cbind,fracChange)
#Get only the "diagonal" entries
ww = cbind(match(cl$gene,rownames(fracChange)),match(cl$cluster,colnames(fracChange)))
cl$fracAfter = cl$geneFrequency * fracChange[ww]
df = melt(cl[,c(1,2,3,ncol(cl))],id.vars=c('gene','cluster'))
colnames(df)[3] = 'condition'
df$condition = ifelse(df$condition=='geneFrequency','Contaminated','Cleaned')
df$condition = factor(df$condition,levels=c('Contaminated','Cleaned'))
#Colour by the fraction by which they go down
mark = paste0(df$gene,'___',df$cluster)
tmp = sapply(split(df,mark),function(e) {
             a = e$value[e$condition=='Contaminated']
             b = e$value[e$condition=='Cleaned']
             (a-b)/a})
df$fracChange = tmp[mark]
#Dot and line alternative to heatmap
gg = ggplot(df,aes(condition,value,colour=fracChange>0.5)) +
  geom_point() +
  geom_line(aes(group=gene,colour=fracChange>0.5)) +
  xlab('')+
  ylab('Expression Fraction') +
  ylim(0,1) +
  guides(colour=FALSE) +
  scale_colour_manual(values=c('TRUE'=colPal[1],'FALSE'='black')) +
  facet_wrap(~cluster)
#Add label for those that change a lot
adf = df[df$condition=='Cleaned',]
gg = gg + geom_text_repel(data=adf,aes(label=gene,x=as.numeric(condition)+0.05,colour=fracChange>0.5),hjust=0,direction='y',max.iter=10000,min.segment.length=Inf,force=1)
pdf(file.path(plotDir,'ClusterChanges.pdf'),width=28,height=28)
plot(gg)
dev.off()
png(file.path(plotDir,'ClusterChanges.png'),width=28,height=28,units='in',res=rasterRes)
plot(gg)
dev.off()
#############################
# Just the ones that change
keepClusts = unique(df$cluster[df$fracChange>0.5])
df = df[df$cluster %in% keepClusts,]
df$clustNom = annMap[df$cluster]
gg = ggplot(df,aes(condition,value,colour=fracChange>0.5)) +
  geom_point() +
  geom_line(aes(group=gene,colour=fracChange>0.5)) +
  xlab('')+
  ylab('Expression Fraction') +
  ylim(0,1) +
  guides(colour=FALSE) +
  scale_colour_manual(values=c('TRUE'=colPal[1],'FALSE'='black')) +
  facet_wrap(~clustNom,ncol=4)
#Add label for those that change a lot
adf = df[df$condition=='Cleaned' & df$fracChange>0.5,]
#Manually tweak position for the main figure
adf$value[adf$clustNom=='Dead' & adf$gene=='DEFA3']=0.05
adf$value[adf$clustNom=='Dead' & adf$gene=='HBA1']=0.08
adf$value[adf$clustNom=='Dead' & adf$gene=='HBD']=0.11
adf$value[adf$clustNom=='MNP' & adf$gene=='SAA1']=0.06
adf$value[adf$clustNom=='NK' & adf$gene=='COL3A1']=0.03
adf$value[adf$clustNom=='NK' & adf$gene=='CTGF']=0.06
adf$value[adf$clustNom=='NK' & adf$gene=='CPXM1']=0.09
adf$value[adf$clustNom=='NK' & adf$gene=='COL1A1']=0.12
adf$value[adf$clustNom=='NK' & adf$gene=='COL1A2']=0.15
adf$value[adf$clustNom=='NK' & adf$gene=='TPSAB1']=0.18
adf$value[adf$clustNom=='NK' & adf$gene=='MYL9']=0.21
adf$value[adf$clustNom=='NK' & adf$gene=='RNASE1']=0.24
gg = gg + geom_text(data=adf,aes(label=gene,x=as.numeric(condition)+0.05,colour=fracChange>0.5),hjust=0)
pdf(file.path(plotDir,'ClusterChangesNegative.pdf'),width=14,height=7)
plot(gg)
dev.off()
png(file.path(plotDir,'ClusterChangesNegative.png'),width=14,height=7,units='in',res=rasterRes)
plot(gg)
dev.off()


##############################################
# Contamination fraction estimate by channel #
##############################################
nUMIs = sclReal$metaData[sclReal$metaData$isCell,'nUMIs']
fac = cut(nUMIs,50)
cellGroupings = split(fac,sclReal$metaData$channelName[sclReal$metaData$isCell])
tmp = sclReal
for(nom in names(tmp$channels)){
  tmp$channels[[nom]]$candidateNonExpressedGenes$isCandidate = rownames(tmp$channels[[nom]]$candidateNonExpressedGenes) %in% rownames(tmp$channels[[nom]]$unexpressedMatrix)
}
tmp = calculateContaminationFraction(tmp,cellGroupings=cellGroupings)
#Get a big data.frame
df = lapply(tmp$channels,function(e) e$rhoGrouped)
df = cbind(do.call(rbind,df),Channel = rep(names(df),sapply(df,nrow)))
df$Sample = simpNames[cMani$Label[as.integer(gsub('Channel','',df$Channel))]]
#Add the ideal samples
tmp = sclIdeal$channels$mixture$rhoGrouped
tmp$Channel = 'mixture'
tmp$Sample = 'Ideal'
#Drop doublets
m = match(sclIdeal$metaData[names(sclIdeal$isDoublet)[sclIdeal$isDoublet],'nUMIs'],tmp$nUMIs)
tmp = tmp[-m,]
df = rbind(df,tmp)
#Make the plot
globs = df[grep('Global',rownames(df)),]
#Grouped
groups = df[-grep('Global',rownames(df)),]
gg = ggplot(groups,aes(log10(nUMIs),log10(est))) +
  geom_hline(data=globs,aes(yintercept = log10(est)),linetype=1,colour='red',size=0.1) +
  geom_hline(data=globs,aes(yintercept = log10(upper)),linetype=2,colour='red',size=0.1) +
  geom_hline(data=globs,aes(yintercept = log10(lower)),linetype=2,colour='red',size=0.1) +
  #geom_hex(bins=20) +
  geom_point(size=0.5) +
  geom_line(size=0.5) +
  geom_errorbar(aes(ymin=log10(lower),ymax=log10(upper)),size=0.5) +
 ylab('Contamination Fraction') +
  xlab('log10(nUMIs/cell') +
  ylim(-3,0) +
  facet_wrap(~Sample)
pdf(file.path(plotDir,'AllContaminationPlots.pdf'),width=14,height=14)
plot(gg)
dev.off()
png(file.path(plotDir,'AllContaminationPlots.png'),width=14,height=14,units='in',res=rasterRes)
plot(gg)
dev.off()


#######################
# Ideal Contamination #
#######################
#Get the per-cell contamination
sc = sclIdeal$channels$mixture
sc = calculateContaminationFraction(sc,cellGroupings=factor(seq(ncol(sc$toc))))
df = sc$rhoGrouped[-1,]
rownames(df) = colnames(sc$toc)
df$isHuman = sclIdeal$isHuman
df$isDoublet = sclIdeal$isDoublet
df$Species = ifelse(df$isHuman,'Human','Mouse')
df$nUMIs = sc$metaData[rownames(df),'nUMIs']
df = df[!df$isDoublet,]
gg_idealRho = ggplot(df,aes(log10(nUMIs),log10(est),colour=Species)) +
  geom_point() + 
  geom_line() + 
  #geom_errorbar(aes(ymin=log10(lower),ymax=log10(upper))) + 
  ylab('log10(Contamination Fraction)') + 
  xlab('log10(nUMIs)') +
  scale_colour_manual(values=colPal[1:2]) +
  guides(colour=FALSE) + 
  theme_bw(base_size=18)
pdf(file.path(plotDir,'IdealContaminationPlot.pdf'),width=7,height=7)
plot(gg_idealRho)
dev.off()
png(file.path(plotDir,'IdealContaminationPlot.png'),width=7,height=7,units='in',res=rasterRes)
plot(gg_idealRho)
dev.off()


####################################
# Infering non-expressed gene plot #
####################################
genes = c('HBA1','HBA2','HBB','TPSAB1','TPSB2')
gg = plotMarkerDistribution(scl,'Channel24',genes) +
  ggtitle('pRCC_2')
pdf(file.path(plotDir,'pRCC_bimodalGenes.pdf'),width=14,height=7)
plot(gg)
dev.off()
png(file.path(plotDir,'pRCC_bimodalGenes.png'),width=14,height=7,units='in',res=rasterRes)
plot(gg)
dev.off()


#######################################
# Make a couple of showing off plots. #
#######################################
DR = as.data.frame(soggy@dr$tsne@cell.embeddings)
colnames(DR) = c('RD1','RD2')
pdf(file.path(plotDir,'ChangeInDistribution.pdf'),width=14)
for(gene in c('LYZ','CD74','HLA-DRA','IL32','TRAC','CD3D','S100A9','S100A8','LTB','NKG7','GNLY','CD4','CD8A','HBB','TPSAB1')){
  plot(plotChangeMap(scl,gene,DR,includePanels=c('Uncorrected','CorrectedCounts'))+ggtitle(gene))
}
dev.off()
#A more reduced set, presented cleanly
ggs = list()
for(gene in c('HBB','HLA-DRA','S100A9','CD74','SAA1','SAA2','COL1A2','COL1A1','COL3A1')){
  ggs[[gene]] = plotChangeMap(scl,gene,DR,includePanels=c('Uncorrected','CorrectedCounts'))
  df = ggs[[gene]]$data
  df$expressed = ifelse(df$data,'yes','no')
  df$correction = factor(ifelse(df$correction=='Uncorrected','Contaminated','Cleaned'),
                         levels=c('Contaminated','Cleaned'))
  ggs[[gene]] = ggplot(df,aes(RD1,RD2)) + 
    geom_point(aes(colour=expressed),size=0.3,alpha=1/2) +
    scale_colour_manual(values=c('yes'='#e60000','no'='#808080')) +
    guides(colour=FALSE) +
    facet_grid(~correction) +
    xlab('tSNE1') +
    ylab('tSNE2') +
    ggtitle(gene)
}
pdf(file.path(plotDir,'ChangeInDistributionKey.pdf'),width=28,height=28)
plot_grid(ggs[[1]],ggs[[2]],ggs[[3]],ggs[[4]],nrow=2,ncol=2)
dev.off()
png(file.path(plotDir,'ChangeInDistributionKey.png'),width=28,height=28,units='in',res=rasterRes)
plot_grid(ggs[[1]],ggs[[2]],ggs[[3]],ggs[[4]],nrow=2,ncol=2)
dev.off()
pdf(file.path(plotDir,'SAA1_correction.pdf'),width=14,height=7)
plot(ggs[['SAA1']])
dev.off()
png(file.path(plotDir,'SAA1_correction.png'),width=14,height=7,units='in',res=rasterRes)
plot(ggs[['SAA1']])
dev.off()
pdf(file.path(plotDir,'SAA2_correction.pdf'),width=14,height=7)
plot(ggs[['SAA2']])
dev.off()
png(file.path(plotDir,'SAA2_correction.png'),width=14,height=7,units='in',res=rasterRes)
plot(ggs[['SAA2']])
dev.off()
pdf(file.path(plotDir,'COL1A2_correction.pdf'),width=14,height=7)
plot(ggs[['COL1A2']])
dev.off()
png(file.path(plotDir,'COL1A2_correction.png'),width=14,height=7,units='in',res=rasterRes)
plot(ggs[['COL1A2']])
dev.off()
pdf(file.path(plotDir,'COL1A1_correction.pdf'),width=14,height=7)
plot(ggs[['COL1A1']])
dev.off()
png(file.path(plotDir,'COL1A1_correction.png'),width=14,height=7,units='in',res=rasterRes)
plot(ggs[['COL1A1']])
dev.off()
pdf(file.path(plotDir,'COL3A1_correction.pdf'),width=14,height=7)
plot(ggs[['COL3A1']])
dev.off()
png(file.path(plotDir,'COL3A1_correction.png'),width=14,height=7,units='in',res=rasterRes)
plot(ggs[['COL3A1']])
dev.off()

