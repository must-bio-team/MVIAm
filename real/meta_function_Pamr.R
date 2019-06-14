require('pamr')


crossValidateMergedPam = function(ematMerged, sampleMetadata, foldid, className='class', ...) {
	y = sampleMetadata[colnames(ematMerged), className]
	x = t(scale(t(ematMerged), center=TRUE, scale=FALSE))
	foldsTmp = foldid[colnames(ematMerged)]
	folds = split(seq_along(foldsTmp), foldsTmp)
	pamData = list(x=x, y=y, geneid=rownames(ematMerged))
	pamTrain = pamr.train(pamData, ...)
	cvFit = pamr.cv(pamTrain, pamData, folds=folds)
	return(cvFit)}


predictValidationDataPam = function(ematList, studyMetadata, sampleMetadata, discoverySampleNames, classesTrain, threshold,
												batchColname='study', covariateName=NA, className='class', predType='posterior', ...) {
	discoveryStudyNames = studyMetadata[studyMetadata[,'discovery'], 'study']
	validationStudyNames = studyMetadata[studyMetadata[,'validation'], 'study']
	
	predsList = foreach(validationStudyName=validationStudyNames) %do% {
		idxValidation = sampleMetadata[,'study']==validationStudyName & (sampleMetadata[,className] %in% classesTrain)
		validationSampleNames = sampleMetadata[idxValidation, 'sample']
		
		ematListNow = ematList[c(discoveryStudyNames, validationStudyName)]
		ematMergedDiscVal = mergeStudyData(ematListNow, sampleMetadata, batchColname=batchColname, covariateName=covariateName)
		ematMergedDisc = ematMergedDiscVal[,discoverySampleNames]
		
		pamData = list(x=ematMergedDisc, y=sampleMetadata[discoverySampleNames, className], geneid=rownames(ematMergedDiscVal))
		pamTrain = pamr.train(pamData, threshold=threshold, ...)
		pamPreds = pamr.predict(pamTrain, newx=ematMergedDiscVal[,validationSampleNames], threshold=threshold, type=predType)}
	
	names(predsList) = validationStudyNames
	return(predsList)}


writeConfusionCrossValidationPam = function(cvFitPam, idx, sampleMetadata, discoverySampleNames, classesTrain,
														  className='class', metaAnalysisName='metaAnalysis') {
	predsClass = factor(cvFitPam$yhat[,idx], levels=classesTrain)
	trueClass = factor(sampleMetadata[discoverySampleNames, className], levels=classesTrain)
	confus = table(trueClass, predsClass)
	write.csv(confus, file=sprintf('%s_pam_cv_thresh_%.3g_confusion.csv', metaAnalysisName, cvFitPam$threshold[idx]), quote=FALSE)}


writeConfusionValidationPam = function(predsListPam, threshold, sampleMetadata, className='class', classLevels=NA,
													metaAnalysisName='metaAnalysis') {
	if (is.na(classLevels[1])) {
		classLevels = colnames(predsListPam[[1]])}
	
	predsProb = do.call(rbind, predsListPam)
	predsClass = colnames(predsProb)[apply(predsProb, MARGIN=1, function(x) which.max(x))]
	predsFactor = factor(predsClass, levels=classLevels)
	trueClasses = factor(sampleMetadata[rownames(predsProb), className], levels=classLevels)
	confus = table(trueClasses, predsFactor)
	write.csv(confus, file=sprintf('%s_pam_val_thresh_%.3g_confusion.csv', metaAnalysisName, threshold), quote=FALSE)}


writeConfusionValidationEachPam = function(predsListPam, threshold, sampleMetadata, className='class', classLevels=NA,
														 metaAnalysisName='metaAnalysis') {
	if (is.na(classLevels[1])) {
		classLevels = colnames(predsListPam[[1]])}
	
	for (validationStudyName in names(predsListPam)) {
		predsProb = predsListPam[[validationStudyName]]
		predsClass = colnames(predsProb)[apply(predsProb, MARGIN=1, function(x) which.max(x))]
		predsFactor = factor(predsClass, levels=classLevels)
		trueClasses = factor(sampleMetadata[rownames(predsProb), className], levels=classLevels)
		confus = table(trueClasses, predsFactor)
		write.csv(confus, file=sprintf('%s_pam_val_%s_thresh_%.3g_confusion.csv', metaAnalysisName, validationStudyName, threshold),
					 quote=FALSE)}}

