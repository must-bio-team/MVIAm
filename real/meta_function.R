require('GEOquery')
require('affy')
require('annotate')
require('org.Hs.eg.db')
require('foreach')
require('impute')
require('dplyr')
require('sva')
require('glmnet')
require('beeswarm')
require('reshape2')
require('ROCR')
require('cba')
require('pheatmap')
require('ggplot2')
require('gridExtra')
require('RColorBrewer')


fixCustomCdfGeneIds = function(geneIds) {
  return(sub('_at', '', geneIds))}


fixGeoSampleNames = function(sampleNames) {
  sampleNames = paste0(toupper(sampleNames), '_')
  regexResult = regexpr('^GSM[0-9]+[^0-9]', sampleNames)
  sampleNamesNew = mapply(function(sampleName, matchLength) substr(sampleName, 1, matchLength-1),
                          sampleNames, attr(regexResult, 'match.length'))
  return(sampleNamesNew)}


fixCelSampleNames = function(sampleNames) {
  sampleNamesNew = gsub('\\.cel$', '', sampleNames, ignore.case=TRUE)
  return(sampleNamesNew)}


getGeneProbeMappingAffy = function(mappingFilePath) {
  mapping = read.table(mappingFilePath, sep='\t', header=TRUE, stringsAsFactors=FALSE)
  mappingUnique = unique(mapping[,c('Probe.Set.Name', 'Affy.Probe.Set.Name')])
  rownames(mappingUnique) = NULL
  colnames(mappingUnique) = c('geneId', 'probeSet')
  return(mappingUnique)}


getGeneProbeMappingDirect = function(featureDf, geneColname, probeColname='ID') {
  mapping = featureDf[,c(probeColname, geneColname)]
  mapping = mapping[apply(mapping, MARGIN=1, function(x) all(!is.na(x) & x!='')),]
  mapping = data.frame(lapply(mapping, as.character), stringsAsFactors=FALSE)
  colnames(mapping) = c('probeSet', 'geneId')
  return(mapping)}


getGeneProbeMappingAnno = function(featureDf, dbName, interName) {
  mappingProbeIntermediate = featureDf[!is.na(featureDf[,interName]) & featureDf[,interName]!='', c('ID', interName)]
  colnames(mappingProbeIntermediate) = c('probeSet', 'geneInter')
  mapTmp1 = eval(parse(text=dbName))
  mapTmp2 = mappedkeys(mapTmp1)
  mapTmp3 = as.list(mapTmp1[mapTmp2])
  geneId = do.call(c, mapTmp3)
  geneInter = do.call(c, mapply(function(inter, len) rep_len(inter, len), names(mapTmp3), sapply(mapTmp3, length),
                                SIMPLIFY=FALSE))
  if (dbName=='org.Hs.egUNIGENE2EG') {
    geneInter = sub('Hs.', '', geneInter, fixed=TRUE)}
  mappingIdInter = data.frame(geneId, geneInter, stringsAsFactors=FALSE)
  mapping = merge(mappingIdInter, mappingProbeIntermediate, by='geneInter', sort=FALSE)}


calcExprsByGene = function(eset, mapping) {
  geneIds = unique(mapping[,'geneId'])
  exprsByGene = matrix(nrow=length(geneIds), ncol=ncol(eset), dimnames=list(geneIds, sampleNames(eset)))
  for (geneId in geneIds) {
    exprsTmp = exprs(eset)[mapping[mapping[,'geneId']==geneId, 'probeSet'],, drop=FALSE]
    if (nrow(exprsTmp)==1) {
      exprsByGene[geneId,] = exprsTmp
    } else {
      exprsByGene[geneId,] = rowMedians(t(exprsTmp), na.rm=TRUE)}}
  return(exprsByGene)}


getSupportedPlatforms = function() {
  return(c('GPL180', 'GPL885', 'GPL887', 'GPL962', 'GPL1053', 'GPL1291', 'GPL1293', 'GPL1390',
           'GPL1708', 'GPL5645', 'GPL6254', 'GPL6480', 'GPL6884', 'GPL6947', 'GPL7015'))}


getStudyData = function(parentFolderPath, studyName, studyDataType, platformInfo) {
  cat(sprintf('Loading study %s...\n', studyName))
  # load the data, convert to gene id, normalize and transform where necessary
  if (studyDataType %in% c('affy_geo', 'affy_custom')) {
    require(platformInfo, character.only=TRUE)
    cwd = setwd(file.path(parentFolderPath, studyName))
    eset = justRMA(cdfname=platformInfo)
    setwd(cwd)
    featureNames(eset) = fixCustomCdfGeneIds(featureNames(eset))
    if (studyDataType=='affy_geo') {
      sampleNames(eset) = fixGeoSampleNames(sampleNames(eset))
    } else {
      sampleNames(eset) = fixCelSampleNames(sampleNames(eset))}
    
  } else if (studyDataType=='affy_series_matrix') {	
    mapping = getGeneProbeMappingAffy(file.path(parentFolderPath, paste0(platformInfo, '_mapping.txt')))
    esetOrig = getGEO(filename=file.path(parentFolderPath, paste0(studyName, '_series_matrix.txt')))
    exprs(esetOrig)[exprs(esetOrig)<=0] = min(exprs(esetOrig)[exprs(esetOrig)>0])
    exprs(esetOrig) = log2(exprs(esetOrig))
    exprsByGene = calcExprsByGene(esetOrig, mapping)
    rownames(exprsByGene) = fixCustomCdfGeneIds(rownames(exprsByGene))
    colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))
    eset = ExpressionSet(assayData=exprsByGene, phenoData=phenoData(esetOrig), experimentData=experimentData(esetOrig))
    
  } else if (studyDataType=='series_matrix') {
    supportedPlatforms = getSupportedPlatforms()
    
    if (!(platformInfo %in% supportedPlatforms)) {
      warning(sprintf('Study %s not loaded, because platform %s is not currently supported.', studyName, platformInfo))
      return(NA)}
    
    esetOrig = getGEO(filename=file.path(parentFolderPath, paste0(studyName, '_series_matrix.txt')))
    if (is.list(esetOrig) && length(esetOrig)==1) {
      esetOrig = esetOrig[[1]]}
    
    featureDf = pData(featureData(esetOrig))
    idx = sapply(featureDf, is.factor)
    featureDf[idx] = lapply(featureDf[idx], as.character)
    if (platformInfo=='GPL180') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egSYMBOL2EG', interName='GENE_SYM')
    } else if (platformInfo=='GPL885') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname='GENE')
    } else if (platformInfo=='GPL887') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname='GENE')
    } else if (platformInfo=='GPL962') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egUNIGENE2EG', interName='UNIGENE')
    } else if (platformInfo=='GPL1053') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egSYMBOL2EG', interName='GENE')
    } else if (platformInfo=='GPL1291') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname='Entrez Gene ID')
    } else if (platformInfo=='GPL1293') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname='Entrez Gene ID')
    } else if (platformInfo=='GPL1390') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egREFSEQ2EG', interName='GB_ACC')
    } else if (platformInfo=='GPL1708') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname='GENE')
    } else if (platformInfo=='GPL5645') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egSYMBOL2EG', interName='Gene Name')
    } else if (platformInfo=='GPL6254') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egENSEMBL2EG', interName='ENSEMBL_GENE_ID')
    } else if (platformInfo=='GPL6480') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname='GENE')
    } else if (platformInfo=='GPL6884') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname='Entrez_Gene_ID')
    } else if (platformInfo=='GPL6947') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname='Entrez_Gene_ID')
    } else if (platformInfo=='GPL7015') {
      mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egREFSEQ2EG', interName='GB_LIST')
    } else {
      warning(sprintf('Study %s not loaded, because platform %s is not currently supported.', studyName, platformInfo))
      return(NA)}
    
    exprsByGene = calcExprsByGene(esetOrig, mapping)
    if (any(is.na(exprsByGene))) {
      warning(sprintf('Imputing missing expression values for study %s.', studyName))
      resultImputed = impute.knn(exprsByGene)
      exprsByGene = resultImputed$data}
    colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))
    eset = ExpressionSet(assayData=exprsByGene, phenoData=phenoData(esetOrig), experimentData=experimentData(esetOrig))
    
  } else if (studyDataType=='eset_rds') {
    esetOrig = readRDS(file.path(parentFolderPath, paste0(studyName, '.rds')))
    
    featureDf = pData(featureData(esetOrig))
    if (platformInfo=='ready') {
      return(esetOrig)
    } else if (platformInfo=='rosetta') {
      mapping = getGeneProbeMappingDirect(featureDf, geneColname='EntrezGene.ID', probeColname='probe')
    } else {
      warning(sprintf('Study %s not loaded, because platform %s is not currently supported.', studyName, platformInfo))
      return(NA)}
    
    exprsByGene = calcExprsByGene(esetOrig, mapping)
    if (any(is.na(exprsByGene))) {
      warning(sprintf('Imputing missing expression values for study %s.', studyName))
      resultImputed = impute.knn(exprsByGene)
      exprsByGene = resultImputed$data}
    eset = ExpressionSet(assayData=exprsByGene, phenoData=phenoData(esetOrig), experimentData=experimentData(esetOrig))
    
  } else {
    warning(sprintf('Study %s not loaded, because data type %s is not currently supported.', studyName, studyDataType))
    eset = NA}
  return(eset)}


getStudyDataList = function(parentFolderPath, studyMetadata) {
  esetList = foreach(studyName=rownames(studyMetadata)) %do% {
    if (any(is.na(studyMetadata[studyName,]))) {
      NA
    } else {
      getStudyData(parentFolderPath, studyName, studyMetadata[studyName, 'studyDataType'],
                   studyMetadata[studyName, 'platformInfo'])}}
  names(esetList) = rownames(studyMetadata)
  return(esetList[!is.na(esetList)])}


cleanStudyData = function(esetList, sampleMetadata) {
  # select relevant samples, convert to matrix
  ematList = foreach(studyName=names(esetList)) %do% {
    keepIdx = colnames(esetList[[studyName]]) %in% sampleMetadata[sampleMetadata[,'study']==studyName, 'sample']
    exprs(esetList[[studyName]])[,keepIdx]}
  names(ematList) = names(esetList)
  return(ematList)}

makeMatchSampleMapping = function(metadata, subStudyNames, matchColname) {
  metadataNow = metadata[metadata[,'study'] %in% subStudyNames, c('study', 'sample', matchColname)]
  metadataNow = metadataNow[order(metadataNow[,'study'], decreasing=is.unsorted(subStudyNames)), c('sample', matchColname)]
  headFunc = function(x) x[1]
  mappingDf = metadataNow %>% group_by_(.dots=list(matchColname)) %>% summarise_each(funs(headFunc)) %>%
    data.frame(check.names=FALSE)
  mapping = mappingDf[,'sample']
  names(mapping) = mappingDf[,matchColname]
  return(mapping)}

mergeStudyData = function(ematList, sampleMetadata, batchColname='study', covariateName=NA,
                          batchCorrection=TRUE, parPrior=TRUE) {
  # merge data and perform cross-study normalization
  geneIds = Reduce(intersect, lapply(ematList, function(x) rownames(x)))
  ematList2 = foreach(studyName=names(ematList)) %do% {ematNow = ematList[[studyName]][geneIds,]}
  if (batchCorrection) {
    # if both one-color and two-color data is present, ComBat can fail catastrophically, if data is not scaled beforehand
    ematListScaled = lapply(ematList2, function(emat) (emat - mean(emat)) / sd(emat))
    ematMerged = do.call(cbind, ematListScaled)
    if (is.na(covariateName)) {
      covariateInfo = model.matrix(~rep_len(1, ncol(ematMerged)))
    } else {
      covariateInfo = model.matrix(~sampleMetadata[colnames(ematMerged), covariateName])}
    
    if (length(unique(sampleMetadata[colnames(ematMerged), batchColname]))>1) {
      ### combat ###
      ematMergedNorm = ComBat(ematMerged, batch=sampleMetadata[colnames(ematMerged), batchColname],
                              mod=covariateInfo, par.prior=parPrior)
      ### ber ###
      #batch=as.factor(sampleMetadata[colnames(ematMerged), batchColname])
      #ematMergedNorm = ber(t(ematMerged),batch)
      #ematMergedNorm = t(ematMergedNorm)
    } else {
      ematMergedNorm = ematMerged}
    
    return(ematMergedNorm)
  } else {
    return(do.call(cbind, ematList2))}}