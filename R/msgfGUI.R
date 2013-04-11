# TODO: Add comment
# 
# Author: Thomas
###############################################################################
library(pepmaps)
library(RGtk2)
library(ggplot2)
library(scales)
library(plyr)
library(cairoDevice)
library(grid)
library(mzR)
library(XML)
library(Biostrings)


filelist <- rGtkDataFrame(data.frame(Datafiles=character(), stringsAsFactors=FALSE))
gSignalConnect(filelist, 'row-inserted', f=function(treemodel, treepath, treeiter){
			validButtons()
		})
gSignalConnect(filelist, 'row-deleted', f=function(treemodel, treepath){
			validButtons()
		})

settingsSetNames <- c('New', 'Placeholder 1', 'Placeholder 2', 'Placeholder 3')

rawResults <- rGtkDataFrame(data.frame(SpecFile=character(), SpecID=character(), ScanNum=numeric(), FragMethod=character(), Precursor=numeric(), IsotopeError=numeric(), PrecursorError=numeric(), Charge=numeric(), Peptide=character(), Protein=character(), DeNovoScore=numeric(), MSGFScore=numeric(), SpecEValue=numeric(), EValue=numeric(), QValue=numeric(), PepQValue=numeric(), realPeptide=character(), '.vis'=logical(), stringsAsFactors=FALSE, check.names=FALSE))
#gSignalConnect(rawResults, 'row-inserted', f=function(treemodel, ...){
#			resultBox$setSensitive(TRUE)
#			peptides$setFrame(data.frame(Peptides=unique(rawResults[which(rawResults[,'.vis']), 'realPeptide'], stringsAsFactors=FALSE)))
#			proteins$setFrame(data.frame(Proteins=unique(rawResults[which(rawResults[,'.vis']), 'Protein'], stringsAsFactors=FALSE)))
#		})
#gSignalConnect(rawResults, 'row-deleted', f=function(treemodel, ...){
#			if(nrow(rawResults[]) == 0){
#				nPepID$setText('')
#				nProtID$setText('')
#				protIntersect$setText('')
#				charge$setText(paste('', '\t(min)\n', '', '\t(median)\n', '', '\t(max)', sep=''))
#				dev.set(FDRplot$getData('Device'))
#				grid.newpage()
#				dev.set(spectrumPlot$getData('Device'))
#				grid.newpage()
#				resultBox$setSensitive(FALSE)
#			} else {
#				peptides$setFrame(data.frame(Peptides=unique(rawResults[which(rawResults[,'.vis']), 'realPeptide'], stringsAsFactors=FALSE)))
#				proteins$setFrame(data.frame(Proteins=unique(rawResults[which(rawResults[,'.vis']), 'Protein'], stringsAsFactors=FALSE)))
#				sampleSelect$setActive(0)
#				plotFDR(model=rawResults, sample='All', device=FDRplot$getData('Device'))
#			}
#		})
rawResultsFilter <- rawResults$filter()
rawResultsFilter$setVisibleColumn(ncol(rawResults)-1)



samples <- rGtkDataFrame(data.frame(Samples='All', '.notAll'=FALSE, loaded=FALSE, database=NA, rawfile=NA, resultfile=NA, check.names=FALSE, stringsAsFactors=FALSE))
gSignalConnect(samples, 'row-inserted', f=function(treemodel, treepath, treeiter, user.data){
			nPepID$setText(length(peptides[]))
			nProtID$setText(length(proteins[]))
			protIntersect$setText(calcProtIntersect())
			chargeSum <- summary(rawResults[rawResults[,'.vis'], 'Charge'])
			charge$setText(paste(chargeSum['Min.'], '\t(min)\n', chargeSum['Median'], '\t(median)\n', chargeSum['Max.'], '\t(max)', sep=''))
		})
samplesFilter <- samples$filter()
samplesFilter$setVisibleColumn(1)

peptides <- rGtkDataFrame(data.frame(Peptides=character(), stringsAsFactors=FALSE))
proteins <- rGtkDataFrame(data.frame(Proteins=character(), stringsAsFactors=FALSE))
tempPeptideList <- rGtkDataFrame(data.frame(Peptides=character(), stringsAsFactors=FALSE))
tempSampleList <- rGtkDataFrame(data.frame(Samples=character(), stringsAsFactors=FALSE))

collectSettings <- function(){
	ans <- list()
	ans$tolerance <- paste(sub(',', '.', toleranceValue$getText()), toleranceUnit$getActiveText(), sep='')
	ans$isotopeError <- c(isotopeErrorLow$getValue(), isotopeErrorHigh$getValue())
	ans$tda <- tdaToggle$getMode()
	ans$fragmentation <- fragmentation$getActiveText()
	ans$instrument <- instrument$getActiveText()
	ans$protease <- enzyme$getActiveText()
	ans$termini <- termini$getValue()
	if(modification$getText() != 'Modification file...' & modification$getText() != ''){
		if(file.exists(modification$getText())){
			ans$modification <- modification$getText()
		} else {
			modificationWarning <- gtkMessageDialog(parent=mainWindow, flags='modal', type='warnings', buttons='ok', 'Modification file ignored', show=FALSE)
			modificationWarning['title'] <- 'File problem'
			gSignalConnect(modificationWarning, 'response', f=function(dialog, response, data){
						dialog$destroy()
					})
			modificationWarning$show()
		}
	} else {}
	ans$lengthRange <- c(lengthLow$getValue(), lengthHigh$getValue())
	ans$chargeRange <- c(chargeLow$getValue(), chargeHigh$getValue())
	ans$matches <- match$getValue()
	ans
}
validButtons <- function(){
	fileCheck <- nrow(filelist[drop=FALSE]) == 0
	databaseCheck <- database$getText() == 'Database fasta file...' | database$getText() == ''
	if(!databaseCheck){
		databaseCheck <- !file.exists(database$getText())
	} else {}
	modificationCheck <- modification$getText() == 'Modification file...' | modification$getText() == ''
	if(!modificationCheck){
		modificationCheck <- !file.exists(modification$getText())
	} else {
		modificationCheck <- FALSE
	}
	
	if(fileCheck){
		fileRemoveButton['sensitive'] <- FALSE
	} else {
		fileRemoveButton['sensitive'] <- TRUE
	}
	
	if(any(c(fileCheck, databaseCheck, modificationCheck))){
		analyzeButton['sensitive'] <- FALSE
	} else {
		analyzeButton['sensitive'] <- TRUE
	}
}
rawFilter <- function(data, sample, FDR, decoy=TRUE){
	if(decoy){
		ifelse(data[,1] == sample | sample == 'All', ifelse(data[,15] <= FDR, ifelse(!grepl('XXX_', data[,10]), TRUE, FALSE), FALSE), FALSE)
	} else {
		ifelse(data[,1] == sample | sample == 'All', ifelse(data[,15] <= FDR, TRUE, FALSE), FALSE)
	}
}
sampleFilter <- function(){
	ans <- sampleSelect$getActiveIter()
	if(ans$retval){
		ans <- samples$getValue(ans$iter, 0)$value
	} else {
		ans <- 'All'
	}
	ans
}
rawUpdated <- function(){
	if(nrow(rawResults[]) == 0){
		resultNotebook$setCurrentPage(resultNotebook$pageNum(summaryTab))
		nPepID$setText('')
		nProtID$setText('')
		protIntersect$setText('')
		charge$setText(paste('', '\t(min)\n', '', '\t(median)\n', '', '\t(max)', sep=''))
		peptides$setFrame(data.frame(Peptides=character(), stringsAsFactors=FALSE))
		proteins$setFrame(data.frame(Proteins=character(), stringsAsFactors=FALSE))
		tempPeptideList$setFrame(data.frame(Peptides=character(), stringsAsFactors=FALSE))
		tempSampleList$setFrame(data.frame(Samples=character(), stringsAsFactors=FALSE))
		dev.set(FDRplot$getData('Device'))
		grid.newpage()
		dev.set(spectrumPlot$getData('Device'))
		grid.newpage()
		resultBox$setSensitive(FALSE)
	} else {
		resultBox$setSensitive(TRUE)
		peptides$setFrame(data.frame(Peptides=unique(rawResults[which(rawResults[,'.vis']), 'realPeptide'], stringsAsFactors=FALSE)))
		proteins$setFrame(data.frame(Proteins=unique(rawResults[which(rawResults[,'.vis']), 'Protein'], stringsAsFactors=FALSE)))
		sampleSelect$setActive(0)
		plotFDR(model=rawResults, sample='All', device=FDRplot$getData('Device'))
	}
}
filterUpdated <- function(FDRChange=TRUE, SampleChange=TRUE){
	if(length(samples[, 'Samples']) > 1){
		rawResults[,'.vis'] <- rawFilter(rawResults[], sample=sampleFilter(), FDR=FDR$getValue())
		peptides$setFrame(data.frame(Peptides=unique(rawResults[which(rawResults[,'.vis']), 'realPeptide'], stringsAsFactors=FALSE)))
		proteins$setFrame(data.frame(Proteins=unique(rawResults[which(rawResults[,'.vis']), 'Protein'], stringsAsFactors=FALSE)))
		#plotPepNumber(rawResults, sampleFilter(), FDR$getValue(), pepNumberPlot$getData('Device'))
		nPepID$setText(length(peptides[]))
		nProtID$setText(length(proteins[]))
		chargeSum <- summary(rawResults[rawResults[,'.vis'], 'Charge'])
		charge$setText(paste(chargeSum['Min.'], '\t(min)\n', chargeSum['Median'], '\t(median)\n', chargeSum['Max.'], '\t(max)', sep=''))
		if(SampleChange){
			plotFDR(rawResults, sampleFilter(), device=FDRplot$getData('Device'))
		} else {}
		if(FDRChange){
			#FDR specific actions here
		} else {}
	} else {}
}
plotFDR <- function(model, sample, device){
	if(sample == 'All'){
		df <- model[]
	} else {
		df <- model[model[,'SpecFile'] == sample, ]
	}	
	df$Score <- -log(df$SpecEValue)
	df$decoy <- grepl('XXX_', df$Protein)
	p <- ggplot(data=df, aes(x=Score, group=decoy, fill=decoy)) + theme_bw()
	p <- p + geom_density(alpha=I(0.5)) + scale_fill_manual('', breaks=c(FALSE, TRUE), labels=c('Database          ', 'Decoy'), values=c('#468966', '#B64926'))
	p <- p + scale_x_continuous(limits=c(0, max(df$Score)), expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
	p <- p + theme(panel.border=element_blank(), panel.grid=element_blank(), axis.line=element_line(size=0.5), legend.position='bottom')
	dev.set(device)
	print(p)
}
plotPepNumber <- function(model, sample, FDR, device){
	if(sample == 'All'){
		df <- model[model[,'QValue'] <= FDR,]
	} else {
		df <- model[model[,'SpecFile'] == sample & model[,'QValue'] <= FDR, ]
	}
	df <- df[!grepl('XXX_', df$Protein),]
	nPep <- daply(df, .(Protein), function(x) length(unique(x$realPeptide)))
	p <- ggplot(data=data.frame(nPep=nPep)) + geom_histogram(aes(x=nPep)) + coord_flip() + scale_x_reverse()
	p <- p + theme_bw() + theme(panel.border=element_blank(), panel.grid=element_blank(), axis.line=element_line(size=0.5))
	dev.set(device)
	print(p)
}
calcProtIntersect <- function(){
	df <- rawResults[]
	df <- df[!grepl('XXX_', df$Protein),]
	df <- df[df$QValue <= FDR$getValue(), ]
	names(df)[1] <- sub('#', '', names(df)[1])
	df <- dlply(df, .(SpecFile), function(x) unique(x$Protein))
	length(Reduce(intersect, df))
}
traceIon <- function(data, index, mz, ppm){
	indexRange <- c(1, length(data))
	mzWin <- c(mz-(mz/1000000)*ppm, mz+(mz/1000000)*ppm)
	scan <- peaks(data, index)
	mzIndex <- which(scan[,1] > mzWin[1] & scan[,1] < mzWin[2])
	if(length(mzIndex) == 0){
		intensity <- c()
		retention <- c()
	} else {
		if(length(mzIndex) != 1){
			mzIndex <- mzIndex[which.max(scan[mzIndex, 2])]
		} else {}
		intensity <- c(scan[mzIndex, 2])
		retention <- c(header(data, index)$retentionTime)
	}
	indexBack <- index
	indexForward <- index
	doBreak <- FALSE
	nextIter <- 1
	while(as.logical(nextIter)){
		indexBack <- indexBack-1
		while(header(data, indexBack)$msLevel == 2){
			indexBack <- indexBack-1
			if(indexBack <= indexRange[1]){
				doBreak <- TRUE
				break
			}
		}
		if(doBreak) break
		if(indexBack < 1) break
		scan <- peaks(data, indexBack)
		mzIndex <- which(scan[,1] > mzWin[1] & scan[,1] < mzWin[2])
		if(length(mzIndex) == 0){
			nextIter <- nextIter + 1
		} else {
			if(length(mzIndex) != 1){
				mzIndex <- mzIndex[which.max(scan[mzIndex, 2])]
			} else {}
			intensity <- c(scan[mzIndex, 2], intensity)
			retention <- c(header(data, indexBack)$retentionTime, retention)
		}
		if(nextIter > 2){
			nextIter <- 0
		}
		if(indexBack <= indexRange[1]){
			nextIter <- 0
		}
	}
	nextIter <- 1
	doBreak <- FALSE
	while(as.logical(nextIter)){
		indexForward <- indexForward+1
		while(header(data, indexForward)$msLevel == 2){
			indexForward <- indexForward+1
			if(indexForward >= indexRange[2]){
				doBreak <- TRUE
				break
			}
		}
		if(doBreak) break
		scan <- peaks(data, indexForward)
		mzIndex <- which(scan[,1] > mzWin[1] & scan[,1] < mzWin[2])
		if(length(mzIndex) == 0){
			nextIter <- nextIter + 1
		} else {
			if(length(mzIndex) != 1){
				mzIndex <- mzIndex[which.max(scan[mzIndex, 2])]
			} else {}
			intensity <- c(intensity, scan[mzIndex, 2])
			retention <- c(retention, header(data, indexForward)$retentionTime)
		}
		if(nextIter > 2){
			nextIter <- 0
		}
		if(indexForward >= indexRange[2]){
			nextIter <- 0
		}
	}
	data.frame(intensity, retention)
}
plotSpectrum <- function(data, scan, seq, ppm=40, device, ionTrace=TRUE, ions='aby', neutralLosses=TRUE){
	if(missing(device)){
		device <- dev.cur()
	}
	index <- which(header(data)$acquisitionNum == scan)
	scanInfo <- header(data, index)
	spec <- peaks(data, index)
	spec <- data.frame(mz=spec[,1], intensity=spec[,2], ymin=0)
	if(scanInfo$msLevel > 1){
		precursor <- which.min(abs(spec$mz-scanInfo$precursorMZ))
		if(abs(spec$mz[precursor]-scanInfo$precursorMZ) < (scanInfo$precursorMZ/1000000)*ppm){
			precursor <- spec[precursor, , drop=FALSE]
		} else {
			precursor <- data.frame()
		}
	} else {
		precursor <- data.frame()
	}
	
	if(!missing(seq)){
		p <- ggplot() + geom_linerange(data=spec, aes(x=mz, ymax=intensity, ymin=ymin), colour=I('grey'))
		if(nrow(precursor) != 0){
			p <- p + geom_point(data=precursor, aes(x=mz, y=intensity), size=I(3), colour=I('#B64926'))
		}
		ionlab <- data.frame(spec, ion=NA)
		ionlist <- pepmaps:::fragPattern(seq, ions=ions, neutralLosses=neutralLosses)
		for(i in 1:nrow(ionlist)){
			ind <- which(ionlab$mz < ionlist$mz[i]+(ionlist$mz[i]/1000000)*ppm & ionlab$mz > ionlist$mz[i]-(ionlist$mz[i]/1000000)*ppm)
			if(length(ind) == 1){
				ionlab$ion[ind] <- ionlist$ion[i]
			} else if(length(ind) > 1){
				ind <- ind[which.min(abs(ionlab$mz[ind]-ionlist$mz[i]))]
				ionlab$ion[ind] <- ionlist$ion[i]
			} else {}
		}
		ionlab <- ionlab[!is.na(ionlab$ion),]
		if(nrow(ionlab) != 0){
			p <- p + geom_linerange(data=ionlab, aes(x=mz, ymax=intensity, ymin=ymin))
			p <- p + geom_text(data=ionlab, aes(x=mz, y=intensity, label=ion), colour=I('#468966'), vjust=-0.5, size=4)
		} else {}
	} else {}
	p <- p + theme_bw() + theme(panel.border=element_blank(), panel.grid=element_blank(), axis.line=element_line(size=0.5), legend.position='bottom')
	p <- p + scale_y_continuous('Intensity', expand=c(0,0), limits=c(0, max(spec$intensity)*1.1))
	
	if(ionTrace){
		parentIndex <- which(header(data)$acquisitionNum == scanInfo$precursorScanNum)
		parentInfo <- header(data, parentIndex)
		intensity <- peaks(data, parentIndex)
		parentData <- data.frame(retention=parentInfo$retentionTime, intensity=intensity[which.min(abs(intensity[, 1]-scanInfo$precursorMZ)), 2])
		trace <- traceIon(data, parentIndex, scanInfo$precursorMZ, ppm)
		parentLines <- data.frame(x=c(min(trace$retention), parentData$retention), xend=c(parentData$retention, parentData$retention), y=c(parentData$intensity, min(trace$intensity)), yend=c(parentData$intensity, parentData$intensity))
		q <- ggplot(data=trace, aes(x=retention, y=intensity)) + geom_line()
		q <- q + geom_segment(data=parentLines, aes(x=x, xend=xend, y=y, yend=yend), linetype=I(2))
		q <- q + geom_point(data=parentData, size=I(3), colour=I('#B64926'))
		q <- q + theme_minimal() + theme(panel.grid=element_blank(), axis.ticks=element_blank(), plot.background=element_rect(colour=alpha('#F5F5F5', 0.75), fill=alpha('#F5F5F5', 0.75)), plot.margin=unit(c(10, 10, 0, 0), 'points'))
		q <- q + scale_y_continuous('', expand=c(0,0), breaks=parentData$intensity) + scale_x_continuous('', expand=c(0,0), breaks=parentData$retention)
	} else {}
	dev.set(device)
	print(p)
	if(ionTrace){
		vp <- viewport(width=0.4, height=0.3, x=0.98, y=0.98, just=c('right', 'top'))
		print(q, vp=vp)	
	}
}
collectRawPlotSettings <- function(){
	ans <- list()
	ans$ppm <- as.numeric(ppmPlotSetting$getText())
	ans$neutralLosses <- neutralPlotSetting$getActive()
	ans$ionTrace <- tracePlotSetting$getActive()
	ans$ions <- paste(if(aIonPlotSetting$getActive()) 'a', if(bIonPlotSetting$getActive()) 'b', if(cIonPlotSetting$getActive()) 'c', if(xIonPlotSetting$getActive()) 'x', if(yIonPlotSetting$getActive()) 'y', if(zIonPlotSetting$getActive()) 'z', sep='')
	ans
}
getRunInfoFromMZID <- function(file){
	parsedMZID <- xmlTreeParse(file, useInternalNodes=TRUE)
	rootMZID <- xmlRoot(parsedMZID)
	database <- xmlAttrs(rootMZID[['DataCollection']][['Inputs']][['SearchDatabase']])['location']
	if(!file.exists(database)){
		database <- file.path(dirname(file), basename(database))
		if(!file.exists(database)){
			database <- NA
		} else {}
	} else {}
	rawfile <- xmlAttrs(rootMZID[['DataCollection']][['Inputs']][['SpectraData']])['location']
	if(!file.exists(rawfile)){
		rawfile <- file.path(dirname(file), basename(rawfile))
		if(!file.exists(rawfile)){
			rawfile <- NA
		} else {}
	} else {}
	resultfile <- file
	analysisPar <- list()
	for(i in 1:length(xmlSApply(rootMZID[['AnalysisProtocolCollection']][['SpectrumIdentificationProtocol']][['AdditionalSearchParams']], xmlName))){
		if(xmlName(rootMZID[['AnalysisProtocolCollection']][['SpectrumIdentificationProtocol']][['AdditionalSearchParams']][[i]]) == 'userParam'){
			param <- xmlAttrs(rootMZID[['AnalysisProtocolCollection']][['SpectrumIdentificationProtocol']][['AdditionalSearchParams']][[i]])
			if(param['name'] == 'TargetDecoyApproach'){
				analysisPar$tda <- as.logical(param['value'])
			} else if(param['name'] == 'MinIsotopeError'){
				analysisPar$isotopeError[1] <- as.numeric(param['value'])
			} else if(param['name'] == 'MaxIsotopeError'){
				analysisPar$isotopeError[2] <- as.numeric(param['value'])
			} else if(param['name'] == 'FragmentMethod'){
				analysisPar$fragmentation <- param['value']
			} else if(param['name'] == 'Instrument'){
				analysisPar$instrument <- param['value']
			} else if(param['name'] == 'NumTolerableTermini'){
				analysisPar$termini <- as.numeric(param['value'])
			} else if(param['name'] == 'NumMatchesPerSpec'){
				analysisPar$matches <- as.numeric(param['value'])
			} else if(param['name'] == 'MinPepLength'){
				analysisPar$lengthRange[1] <- as.numeric(param['value'])
			} else if(param['name'] == 'MaxPepLength'){
				analysisPar$lengthRange[2] <- as.numeric(param['value'])
			} else if(param['name'] == 'MinCharge'){
				analysisPar$chargeRange[1] <- as.numeric(param['value'])
			} else if(param['name'] == 'MaxCharge'){
				analysisPar$chargeRange[2] <- as.numeric(param['value'])
			} else {}
		} else {}
	}
	analysisPar$protease <- xmlAttrs(rootMZID[['AnalysisProtocolCollection']][['SpectrumIdentificationProtocol']][['Enzymes']][['Enzyme']][['EnzymeName']][['cvParam']])['name']
	analysisPar$tolerance <- as.numeric(xmlAttrs(rootMZID[['AnalysisProtocolCollection']][['SpectrumIdentificationProtocol']][['ParentTolerance']][['cvParam']])['value'])
	tolUnit <- xmlAttrs(rootMZID[['AnalysisProtocolCollection']][['SpectrumIdentificationProtocol']][['ParentTolerance']][['cvParam']])['unitName']
	if(tolower(tolUnit) == 'parts per million'){
		analysisPar$tolerance <- paste(analysisPar$tolerance, 'ppm', sep='')
	} else if(tolower(tolUnit) == 'dalton'){
		analysisPar$tolerance <- paste(analysisPar$tolerance, 'Da', sep='')
	} else {
		analysisPar$tolerance <- as.character(analysisPar$tolerance)
	}
	free(parsedMZID)
	list(database=database, rawfile=rawfile, resultfile=resultfile, parameters=analysisPar)
}
updateSampleList <- function(add){
	newSample <- data.frame(Samples=basename(add$rawfile), '.notAll'=TRUE, loaded=FALSE, database=add$database, rawfile=add$rawfile, resultfile=add$resultfile, check.names=FALSE, stringsAsFactors=FALSE)
	if(any(samples[,'Samples'] %in% newSample$Samples)){
		
	} else {}
	samples$appendRows(newSample)
}
findProteinSequence <- function(name){
	rawfiles <- unique(as.character(rawResults[rawResults[,'Protein'] == name, 1]))
	databases <- unique(samples[samples[, 'Samples'] %in% rawfiles, 'database'])
	seq <- list()
	for(i in 1:length(databases)){
		seqs <- readAAStringSet(databases[i])
		ans <- list(sequence=as.character(seqs[which(sapply(strsplit(names(seqs), ' '), '[', 1) == name)]), files=samples[which(samples[, 'database'] == databases[i]), 'Samples'])
		seq[[i]] <- ans
	}
	seq
}
findProteinEvidence <- function(name, sample, FDR){
	if(missing(sample)){
		ind <- which(rawResults[, 'Protein'] == name & rawResults[, 'QValue'] <= FDR)
		ans <- rawResults[ind, 'realPeptide']
	} else {
		ind <- which(rawResults[, 'Protein'] == name & rawResults[, 'SpecFile'] == sample & rawResults[, 'QValue'] <= FDR)
		ans <- rawResults[ind, 'realPeptide']
	}
	split(ind, ans)
}
findProteinEvidenceLocation <- function(sequence, evidence, collate=TRUE){
	ans <- list()
	for(i in 1:length(evidence)){
		m <- regexpr(names(evidence)[i], sequence)
		ans[[i]] <- c(m, m+attr(m, 'match.length')-1)
	}
	if(collate){
		m <- sort(unique(unlist(lapply(ans, function(x) seq(x[1], x[2])))))
		d <- diff(m)
		d1 <- which(d != 1)
		ans <- list()
		if(length(d1) == 0){
			ans[[1]] <- c(m[1], m[length(m)])
		} else {
			ans[[1]] <- c(m[1], m[d1[1]])
			if(length(d1) > 1){
				for(i in 2:length(d1)){
					ans[[i]] <- c(m[d1[i-1]+1], m[d1[i]])
				}
			}
			ans[[length(ans)+1]] <- c(m[d1[length(d1)]+1], m[length(m)])
		}
	} else {}
	ans
}
formatProteinPango <- function(sequence, evidence, evidenceHighlight, breakPosition=40){
	markupLocation <- findProteinEvidenceLocation(sequence, evidence, collate=TRUE)
	if(!missing(evidenceHighlight)){
		highlightLocation <- findProteinEvidenceLocation(sequence, evidence[evidenceHighlight], collate=FALSE)
	}
	sequence <- strsplit(sequence, '')[[1]]
	ans <- c()
	if(markupLocation[[1]][1] != 1){
		ans <- paste(sequence[1:(markupLocation[[1]][1]-1)], collapse='')
	} else {
		ans <- c()
	}
	for(i in 1:length(markupLocation)){
		if(!missing(evidenceHighlight) && markupLocation[[i]][1] <= highlightLocation[[1]][1] && markupLocation[[i]][2] >= highlightLocation[[1]][2]){
			if(markupLocation[[i]][1] != highlightLocation[[1]][1]) ans <- paste(ans, '<span color=\'#B64926\'>', paste(sequence[markupLocation[[i]][1]:(highlightLocation[[1]][1]-1)], collapse=''), '</span>', sep='')
			ans <- paste(ans, '<span bgcolor=\'#468966\'>', paste(sequence[highlightLocation[[1]][1]:highlightLocation[[1]][2]], collapse=''), '</span>', sep='')
			if(markupLocation[[i]][2] != highlightLocation[[1]][2]) ans <- paste(ans, '<span color=\'#B64926\'>', paste(sequence[(highlightLocation[[1]][2]+1):markupLocation[[i]][2]], collapse=''), '</span>', sep='')
		} else {
			ans <- paste(ans, '<span color=\'#B64926\'>', paste(sequence[markupLocation[[i]][1]:markupLocation[[i]][2]], collapse=''), '</span>', sep='')
		}
		if(i != length(markupLocation)){
			ans <- paste(ans, paste(sequence[(markupLocation[[i]][2]+1):(markupLocation[[i+1]][1]-1)], collapse=''), sep='')
		}
	}
	if(markupLocation[[length(markupLocation)]][2] < length(sequence)){
		ans <- paste(ans, paste(sequence[(markupLocation[[length(markupLocation)]][2]+1):length(sequence)], collapse=''), sep='')
	}
	if(breakPosition > 0 & breakPosition < length(sequence)){
		seqSplitted <- regmatches(ans, gregexpr('<.*?>', ans, perl=T), invert=TRUE)
		if(length(which(seqSplitted[[1]] == '') != 0)) {
			addEmpty <- which(seqSplitted[[1]] == '')
			seqSplitted[[1]] <- seqSplitted[[1]][-which(seqSplitted[[1]] == '')]
		} else {
			addEmpty <- c()
		}
		splitPosition <- split(rep(1:ceiling(length(sequence)/breakPosition), each=breakPosition)[1:length(sequence)], rep(1:length(seqSplitted[[1]]), nchar(seqSplitted[[1]])))
		lineBreakInsert <- mapply(function(x, y) paste(unlist(sapply(split(x, y), function(x) c(x, '\n'))), collapse=''), strsplit(seqSplitted[[1]], ''), splitPosition)
		for(i in 1:(length(lineBreakInsert)-1)){
			if(splitPosition[[i]][length(splitPosition[[i]])] == splitPosition[[i+1]][1]){
				lineBreakInsert[[i]] <- sub('\n$', '', lineBreakInsert[[i]], perl=T)
			} else {}
		}
		lineBreakInsert[[length(lineBreakInsert)]] <- sub('\n$', '', lineBreakInsert[[length(lineBreakInsert)]], perl=T)
		if(length(addEmpty) > 0){
			for(i in addEmpty){
				lineBreakInsert <- append(lineBreakInsert, '', i-1)
			}
		}
		lineBreakInsert <- list(lineBreakInsert)
		regmatches(ans, gregexpr('<.*?>', ans, perl=T), invert=TRUE) <- lineBreakInsert
	}
	ans
}
blastPut <- function(seq, db, ...){
	### Based on blast.pdb
	urlput <- paste("http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?CMD=Put&DATABASE=", db,"&PROGRAM=blastp&CLIENT=web&QUERY=", seq, sep="")
	txt <- scan(urlput, what="raw", sep="\n", quiet=TRUE)
	rid <- sub("^.*RID = " ,"",txt[ grep("RID =",txt) ])
	rid
}
blastCheck <- function(rid){
	urlget <- paste("http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=",rid, sep="")
	txt <- scan(urlget, what="raw", sep="\n", quiet=TRUE)
	status <- sub("^.*Status=" ,"",txt[ grep("Status=",txt) ])
	status
}
blastGet <- function(rid, browser=TRUE, ...){
	urlget <- paste("http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=",rid, sep="")
	if(browser){
		browseURL(urlget)
	}
}
formatPeptideSequence <- function(seq){
	ans <- strsplit(seq, '.', fixed=T)[[1]]
	ans <- paste(sub('-', '|', ans), collapse='--')
	modPos <- c()
	mod <- c()
	while(1){
		mods <- regexpr('\\+([\\d]|,)+', ans, perl=TRUE)
		if(mods == -1){
			break
		} else {
			modPos <- c(modPos, mods-1)
			mod <- c(mod, regmatches(ans, mods))
			regmatches(ans, mods) <- ''
		}
	}
	top <- list(rep(' ', nchar(ans)), rep(' ', nchar(ans)))
	bottom <- list(rep(' ', nchar(ans)), rep(' ', nchar(ans)))
	if(length(modPos) > 0){
		mod <- sub('+', '', mod, fixed=TRUE)
		mod <- sub(',', '.', mod, fixed=TRUE)
		for(i in length(modPos):1){
			start <- 2
			while(1){
				if(start > length(top)){
					top <- c(top, list(rep(' ', nchar(ans))))
				} else {}
				if(start > length(bottom)){
					bottom <- c(bottom, list(rep(' ', nchar(ans))))
				} else {}
				if(all(top[[start]][modPos[i]:(modPos[i]+nchar(mod[i])-1)] == ' ', na.rm=TRUE)){
					top[[start]][modPos[i]:(modPos[i]+nchar(mod[i])-1)] <- strsplit(mod, '')[[1]]
					retract <- start-1
					while(retract){
						top[[retract]][modPos[i]] <- '|'
						retract <- retract-1
					}
					break
				} else if(all(bottom[[start]][modPos[i]:(modPos[i]+nchar(mod[i])-1)] == ' ', na.rm=TRUE)){
					bottom[[start]][modPos[i]:(modPos[i]+nchar(mod[i])-1)] <- strsplit(mod, '')[[1]]
					retract <- start-1
					while(retract){
						bottom[[retract]][modPos[i]] <- '|'
						retract <- retract-1
					}
					break
				} else {
					start <- start+1
				}
			}
		}
	} else {}
	top <- paste(sapply(top, paste, collapse=''), '\n', sep='')
	bottom <- paste(sapply(bottom, paste, collapse=''), '\n', sep='')
	ans <- c(rev(top), paste(ans, '\n', sep=''), bottom)
	ans
}
createPeptideView <- function(seq){
	VBox <- gtkVBox(FALSE, 5)
	VBox['border-width'] <- 5
	alternativeRawSeq <- unique(rawResults[rawResults[,'realPeptide'] == seq & rawResults[,'.vis'], 'Peptide'])
	for(i in alternativeRawSeq){
		if(i != alternativeRawSeq[1]){
			separator <- gtkHSeparatorNew()
			separator['height-request'] <- 30
			VBox$packStart(separator, expand=FALSE, fill=FALSE)
		}
		peptideOverview <- gtkLabel(paste(formatPeptideSequence(i), collapse=''))
			font_descr <- pangoFontDescriptionNew()
			font_descr$setFamily('courier')
			peptideOverview$modifyFont(font_descr)
			peptideOverviewScroll <- gtkScrolledWindow()
				peptideOverviewScroll$setPolicy(hscrollbar.policy=GtkPolicyType[2], vscrollbar.policy=GtkPolicyType[3])
				peptideOverviewScroll$addWithViewport(peptideOverview)
				peptideOverviewScroll$getChild()$setShadowType(GtkShadowType[1])
			VBox$packStart(peptideOverviewScroll, expand=FALSE, fill=FALSE)
			
		peptideViewTable <- gtkTable()
			peptideViewTable['border-width'] <- 5
			
			nDetectionLabel <- gtkLabel('Number of times detected:')
			nDetectionLabel['xalign'] <- 1
			peptideViewTable$attach(nDetectionLabel, left.attach=0,1, top.attach=0,1)
			nDetectionValue <- gtkLabel(sum(rawResults[,'Peptide'] == i & rawResults[,'.vis']))
			nDetectionValue['xalign'] <- 0
			peptideViewTable$attach(nDetectionValue, left.attach=1,2, top.attach=0,1)
			
			pepLengthLabel <- gtkLabel('Length (residues):')
			pepLengthLabel['xalign'] <- 1
			peptideViewTable$attach(pepLengthLabel, left.attach=0,1, top.attach=1,2)
			pepLengthValue <- gtkLabel(nchar(seq))
			pepLengthValue['xalign'] <- 0
			peptideViewTable$attach(pepLengthValue, left.attach=1,2, top.attach=1,2)
			
			pepMassLabel <- gtkLabel('Mass:')
			pepMassLabel['xalign'] <- 1
			peptideViewTable$attach(pepMassLabel, left.attach=0,1, top.attach=2,3)
			pepMassValue <- gtkLabel(paste(sprintf('%.2f', pepMass(seq)), 'Da'))
			pepMassValue['xalign'] <- 0
			peptideViewTable$attach(pepMassValue, left.attach=1,2, top.attach=2,3)
			
			pepPILabel <- gtkLabel('Isoelectric point (pI):')
			pepPILabel['xalign'] <- 1
			peptideViewTable$attach(pepPILabel, left.attach=0,1, top.attach=3,4)
			pepPIValue <- gtkLabel(sprintf('%.2f', peppI(seq)))
			pepPIValue['xalign'] <- 0
			peptideViewTable$attach(pepPIValue, left.attach=1,2, top.attach=3,4)
			
			hydroPepLabel <- gtkLabel('Hydrophobicity (Q Value):')
			hydroPepLabel['xalign'] <- 1
			peptideViewTable$attach(hydroPepLabel, left.attach=0,1, top.attach=4,5)
			hydroPepValue <- gtkLabel(sprintf('%.1f', Qval(seq)))
			hydroPepValue['xalign'] <- 0
			peptideViewTable$attach(hydroPepValue, left.attach=1,2, top.attach=4,5)
			
			peptideViewTable$setRowSpacings(15)
			peptideViewTable$setColSpacings(5)
			VBox$packStart(peptideViewTable, expand=FALSE, fill=FALSE)
		
		peptideDF <- rawResults[rawResults[,'Peptide'] == i & rawResults[,'.vis'], ]
		peptideViewModel <- gtkTreeStore('gchararray', 'gchararray', 'gchararray')
			by(peptideDF, peptideDF$Protein, function(x){
						parent_iter <- peptideViewModel$append()
						peptideViewModel$setValue(parent_iter$iter, column=0, paste('Protein: ', x$Protein[1], sep=''))
						by(x, x$SpecFile, function(y){
									parent2_iter <- peptideViewModel$append(parent=parent_iter$iter)
									peptideViewModel$setValue(parent2_iter$iter, column=0, paste('Sample: ', y$SpecFile[1], sep=''))
									by(y, y$Charge, function(z){
												parent3_iter <- peptideViewModel$append(parent=parent2_iter$iter)
												peptideViewModel$setValue(parent3_iter$iter, column=0, paste('Charge: ', z$Charge[1], sep=''))
												for(i in 1:nrow(z)){
													child_iter <- peptideViewModel$append(parent=parent3_iter$iter)
													peptideViewModel$setValue(child_iter$iter, column=1, value=z[i, 'ScanNum'])
													peptideViewModel$setValue(child_iter$iter, column=2, value=z[i, 'QValue'])
												}
											})
								})
					})
			peptideTreeView <- gtkTreeView()
				peptideTreeView$insertColumnWithAttributes(-1, '', gtkCellRendererText(), text=0)
				peptideTreeView$insertColumnWithAttributes(-1, 'Scan Number', gtkCellRendererText(), text=1)
				peptideTreeView$insertColumnWithAttributes(-1, 'FDR', gtkCellRendererText(), text=2)
				peptideTreeView$setModel(peptideViewModel)
				peptideTreeViewScroll <- gtkScrolledWindow()
				peptideTreeView['height-request'] <- 200
				peptideTreeViewScroll$add(peptideTreeView)
				peptideTreeSelection <- peptideTreeView$getSelection()
				gSignalConnect(peptideTreeView, 'row-activated', function(...){
							if(!with(peptideTreeSelection$getSelected(), model$iterParent(iter)$retval)){
								proteinName <- with(peptideTreeSelection$getSelected(), model$getValue(iter, 0)$value)
								proteinName <- sub('Protein: ', '', proteinName)
								proteinPath <- gtkTreePathNewFromString(which(proteins[] == proteinName)-1)
								resultNotebook$setCurrentPage(resultNotebook$pageNum(proteinTab))
								proteinSelect$selectPath(proteinPath)
								protein$scrollToCell(proteinPath, use.align=TRUE, row.align=0.5)
								gSignalEmit(protein, 'cursor-changed')
							} else if(!with(peptideTreeSelection$getSelected(), model$iterHasChild(iter))){
								scanNum <- with(peptideTreeSelection$getSelected(), model$getValue(iter, 1)$value)
								sample <- sub('Sample: ', '', with(peptideTreeSelection$getSelected(), with(model$iterParent(model$iterParent(iter)$iter), model$getValue(iter, 0)$value)))
								proteinName <- sub('Protein: ', '', with(peptideTreeSelection$getSelected(), with(model$iterParent(model$iterParent(model$iterParent(iter)$iter)$iter), model$getValue(iter, 0)$value)))
								ind <- which(rawResults[,'SpecFile'] == sample & rawResults[,'ScanNum'] == scanNum & rawResults[,'Protein'] == proteinName)
								if(length(ind) == 1){
									rawPath <- gtkTreePathNewFromString(ind-1)
									resultNotebook$setCurrentPage(resultNotebook$pageNum(rawTab))
									rawSelect$selectPath(rawPath)
									raw$scrollToCell(rawPath)
									raw$rowActivated(rawPath, gtkTreeViewColumn())
								}
							}
						})
			VBox$packStart(peptideTreeViewScroll, expand=FALSE, fill=FALSE)
	}
	viewportVBox <- gtkViewport()
	viewportVBox$add(VBox)
	viewportVBox$setShadowType(GtkShadowType[1])
	return(viewportVBox)
}






mainWindow <- gtkWindow(show=FALSE)
mainWindow$setTitle('MS GF+')
topGroup <- gtkVBox(FALSE, 5)
mainWindow$add(topGroup)

toolbar <- gtkToolbar()
	topGroup$packStart(toolbar, expand=FALSE)
	importAction <- gtkAction(name='import', label='Import', tooltip='Import mzIdentML files', stock.id='gtk-add')
	gSignalConnect(importAction, 'activate', f=function(...){
				file <- gtkFileChooserDialog(title='Select MS datafiles', parent=NULL, action='open', 'gtk-ok', GtkResponseType['ok'], 'gtk-cancel', GtkResponseType['cancel'], show=FALSE)
				file['select-multiple'] <- TRUE
				gSignalConnect(file, 'response', f=function(dialog, response, data){
							dialog$hide()
						})
				response <- file$Run()
				if(response == GtkResponseType["ok"]){
					filenames <- unlist(file$getFilenames())
					for(i in 1:length(filenames)){
						if(Sys.info()["sysname"] == 'Windows'){
							filename <- paste('\"', filenames[i], '\"', sep='')
						} else {
							filename <- gsub(' ', '\\ ', filenames[i], fixed=T)
						}
						callConv <- paste('java -Xmx', 10000, 'M -cp ', R.home(component='library/pepmaps/java/MSGFplus.jar'), ' edu.ucsd.msjava.ui.MzIDToTsv -i ', filename, ' -o ', paste(filename, '.tsv', sep=''), ' -unroll 1 -showDecoy 1', sep='')
						system(callConv)
						if(length(scan(paste(filenames[i], '.tsv', sep=''), skip=1, nlines=1, what='character', quiet=T)) != 0){
							results <- read.table(paste(filenames[i], '.tsv', sep=''), sep='\t', header=TRUE, comment.char='', stringsAsFactors=FALSE)
							names(results) <- sub('#', '', scan(paste(filenames[i], '.tsv', sep=''), nlines=1, what=character(), quiet=TRUE))
							unlink(paste(filenames[i], '.tsv', sep=''))
							sNum <- regexpr('scanId=[\\d]+', results$SpecID, perl=TRUE)
							sNum <- substr(results$SpecID, sNum, sNum+attr(sNum, 'match.length')-1)
							sNum <- sub('scanId=', '', sNum)
							results$ScanNum <- sNum
							peptide <- strsplit(as.character(results$Peptide), '\\.')
							peptide <- unlist(peptide)[seq(2, length(peptide)*3, by=3)]
							peptide <- gsub('\\+([\\d]|,)+', '', peptide, perl=TRUE)
							results$realPeptide <- peptide
							results$'.vis' <- rawFilter(results, sample=sampleFilter(), FDR=FDR$getValue())
							rawResults$appendRows(results)
							rawUpdated()
							sampleInfo <- getRunInfoFromMZID(filenames[i])
							updateSampleList(add=sampleInfo)
						} else {}
					}
				} else {}
				file$destroy()
			})
	importToolButton <- importAction$createToolItem()
	toolbar$add(importToolButton)
	removeAction <- gtkAction(name='remove', label='Remove', tooltip='Remove analysis results', stock.id='gtk-remove')
	gSignalConnect(removeAction, 'activate', f=function(...){
				removeDialog <- gtkDialogNewWithButtons(title='Remove sample results', parent=mainWindow, flags=1, 'Remove', GtkResponseType['ok'], 'Cancel', GtkResponseType['cancel'], show=FALSE)
				removeDialogContent <- removeDialog$getContentArea()
				removeDialogScroll <- gtkScrolledWindow()
				removeDialogContent$packStart(removeDialogScroll)
				removeDialogList <- gtkTreeView(samplesFilter)
				removeDialogList$insertColumnWithAttributes(position=-1, title='MS datafiles', cell=gtkCellRendererText(), text=0)
				removeDialogList$setHeadersVisible(FALSE)
				removeDialogList['height-request'] <- 100
				removeDialogSelect <- removeDialogList$getSelection()
				removeDialogSelect$setMode('multiple')
				removeDialogScroll$add(removeDialogList)
				gSignalConnect(removeDialog, 'response', f=function(dialog, response, ...){
							if(response == GtkResponseType['ok']){
								samplesToRemoveIndex <- sapply(removeDialogSelect$getSelectedRows()$retval, function(x) as.numeric(samplesFilter$convertPathToChildPath(x)$toString())+1)
								samplesToRemove <- samples[samplesToRemoveIndex, 1]
								allIndex <- which(samplesToRemove == 'All')
								if(length(allIndex) != 0){
									samplesToRemoveIndex <- samplesToRemoveIndex[-allIndex]
									samplesToRemove <- samplesToRemove[-allIndex]
								}
								rawRowsToRemove <- which(rawResults[,'SpecFile'] %in% samplesToRemove)
								samples$setFrame(data.frame(samples[-samplesToRemoveIndex, , drop=FALSE]))
								rawResults$setFrame(data.frame(rawResults[-rawRowsToRemove, , drop=FALSE]))
								rawUpdated()
							}
							dialog$Destroy()
						})
				removeDialog$show()
			})
	removeToolButton <- removeAction$createToolItem()
	toolbar$add(removeToolButton)
	exportAction <- gtkAction(name='export', label='Export', tooltip='Export results', stock.id='gtk-save')
	exportToolButton <- exportAction$createToolItem()
	toolbar$add(exportToolButton)
	toolbar$add(gtkSeparatorToolItem())
	sampleInfoAction <- gtkAction(name='sample.info', label='Samples', tooltip='View sample information', stock.id='gtk-info')
	sampleInfoToolButton <- sampleInfoAction$createToolItem()
	toolbar$add(sampleInfoToolButton)

mainGroup <- gtkHBox(FALSE, 5)
mainGroup$setBorderWidth(5)
topGroup$packStart(mainGroup)

controlGroup <- gtkVBox(FALSE, 5)
mainGroup$packStart(controlGroup, expand=FALSE, fill=FALSE)

importFrame <- gtkFrame('Data files to use')
controlGroup$packStart(importFrame)
	importBox <- gtkVBox(FALSE, 5)
	importBox$setBorderWidth(5)
	importFrame$add(importBox)
	
	datafileSelect <- gtkHBox(FALSE, 5)
	importBox$packStart(datafileSelect)
	
	filesScroll <- gtkScrolledWindow()
	datafileSelect$packStart(filesScroll)
	files <- gtkTreeView(filelist)
		files$insertColumnWithAttributes(position=-1, title='MS datafiles', cell=gtkCellRendererText(), text=0)
		files$setHeadersVisible(FALSE)
		files['height-request'] <- 100
		fileSelect <- files$getSelection()
		fileSelect$setMode('multiple')
	filesScroll$add(files)
	
	datafileSelectButtons <- gtkVBox(FALSE, 5)
	datafileSelect$packStart(datafileSelectButtons, expand=FALSE, fill=FALSE)
		fileAddButton <- gtkButton(label='Add')
		fileAddButton['width-request'] <- 65
		gSignalConnect(fileAddButton, 'clicked', f=function(widget, ...){
				file <- gtkFileChooserDialog(title='Select MS datafiles', parent=mainWindow, action='open', 'gtk-ok', GtkResponseType['ok'], 'gtk-cancel', GtkResponseType['cancel'], show=FALSE)
					file['select-multiple'] <- TRUE
					gSignalConnect(file, 'response', f=function(dialog, response, data){
								if(response == GtkResponseType['ok']){
									files <- as.character(dialog$getFilenames())
									filelist$appendRows(data.frame(Datafiles=files, stringsAsFactors=FALSE))
									setwd(dirname(files[1]))
								}
								dialog$destroy()
							})
					file$show()
				})
		fileRemoveButton <- gtkButton(label='Remove')
		fileRemoveButton['width-request'] <- 65
		fileRemoveButton['sensitive'] <- FALSE
		gSignalConnect(fileRemoveButton, 'clicked', f=function(widget, ...){
					selected <- sapply(fileSelect$getSelectedRows()$retval, function(x) as.numeric(x$toString())+1)
					filelist$setFrame(filelist[-selected, , drop=FALSE])
				})
		datafileSelectButtons$packStart(fileAddButton, expand=FALSE, fill=FALSE)
		datafileSelectButtons$packStart(fileRemoveButton, expand=FALSE, fill=FALSE)
		
	databaseSelect <- gtkHBox(FALSE, 5)
	importBox$packStart(databaseSelect, expand=FALSE, fill=FALSE)
	
	database <- gtkEntry()
	database['height-request'] <- 24
	database$setText('Database fasta file...')
	gSignalConnect(database, 'changed', f=function(widget, ...){
				validButtons()
			})
	databaseImport <- gtkButton(label='Browse')
	databaseImport['width-request'] <- 65
	gSignalConnect(databaseImport, 'clicked', f=function(widget, ...){
				file <- gtkFileChooserDialog(title='Select MS datafiles', parent=mainWindow, action='open', 'gtk-ok', GtkResponseType['ok'], 'gtk-cancel', GtkResponseType['cancel'], show=FALSE)
				gSignalConnect(file, 'response', f=function(dialog, response, data){
							if(response == GtkResponseType['ok']){
								file <- as.character(dialog$getFilenames())
								database$setText(file)
							}
							dialog$destroy()
						})
				file$show()
			})
	databaseSelect$packStart(database)
	databaseSelect$packStart(databaseImport, expand=FALSE, fill=FALSE)
	
settingsFrame <- gtkFrame('Settings')
controlGroup$packStart(settingsFrame, expand=FALSE, fill=FALSE)
	settingsBox <- gtkVBox(FALSE, 15)
	settingsBox$setBorderWidth(5)
	settingsFrame$add(settingsBox)
	
	settingsSet <- gtkHBox(FALSE, 5)
	settingsBox$packStart(settingsSet, expand=FALSE, fill=FALSE)
	settingsSetChooser <- gtkComboBoxNewText()
	sapply(settingsSetNames, settingsSetChooser$appendText)
	settingsSetSaveButton <- gtkButton(label='Save')
	settingsSetSaveButton['width-request'] <- 65
	settingsSet$packStart(settingsSetChooser)
	settingsSet$packStart(settingsSetSaveButton, expand=FALSE, fill=FALSE)
	
	settingsTable <- gtkTable(rows=11, columns=2, homogeneous=FALSE)
	settingsTable$setColSpacings(10)
	settingsTable$setRowSpacings(10)
	settingsBox$packStart(settingsTable)
		toleranceLabel <- gtkLabel('Tolerance:')
		toleranceLabel['xalign'] <- 1
		settingsTable$attach(toleranceLabel, left.attach=0,1, top.attach=0,1)
		toleranceBox <- gtkHBox(FALSE, 5)
			toleranceValue <- gtkEntry()
			toleranceValue['height-request'] <- 24
			toleranceValue['width-chars'] <- 5
			toleranceValue$setAlignment(1)
			toleranceValue$setText(20)
			toleranceUnit <- gtkComboBoxNewText()
			sapply(c('ppm', 'Da'), toleranceUnit$appendText)
			toleranceUnit$setActive(0)
			toleranceBox$packStart(toleranceValue, expand=FALSE, fill=FALSE)
			toleranceBox$packStart(toleranceUnit, expand=FALSE, fill=FALSE)
		settingsTable$attach(toleranceBox, left.attach=1,2, top.attach=0,1)
		
		isotopeErrorLabel <- gtkLabel('Isotope error range:')
		isotopeErrorLabel['xalign'] <- 1
		settingsTable$attach(isotopeErrorLabel, left.attach=0,1, top.attach=1,2)
		isotopeErrorBox <- gtkHBox(FALSE, 5)
			isotopeErrorLow <- gtkSpinButton(min=-5, max=10, step=1)
			isotopeErrorLow['height-request'] <- 24
			isotopeErrorLow['width-chars'] <- 3
			isotopeErrorLow$setAlignment(1)
			isotopeErrorLow$setValue(0)
			gSignalConnect(isotopeErrorLow, 'value-changed', f=function(widget, ...){
						min <- widget$getValue()
						if(min > isotopeErrorHigh$getValue()){
							isotopeErrorHigh$setRange(min, 10)
						}
					})
			isotopeErrorHigh <- gtkSpinButton(min=0, max=10, step=1)
			isotopeErrorHigh['height-request'] <- 24
			isotopeErrorHigh['width-chars'] <- 3
			isotopeErrorHigh$setAlignment(1)
			isotopeErrorHigh$setValue(1)
			isotopeErrorBox$packStart(isotopeErrorLow, expand=FALSE, fill=FALSE)
			isotopeErrorBox$packStart(gtkLabel(' - '), expand=FALSE, fill=FALSE)
			isotopeErrorBox$packStart(isotopeErrorHigh, expand=FALSE, fill=FALSE)
		settingsTable$attach(isotopeErrorBox, left.attach=1,2, top.attach=1,2)
		
		tdaLabel <- gtkLabel('Target-Decoy:')
		tdaLabel['xalign'] <- 1
		settingsTable$attach(tdaLabel, left.attach=0,1, top.attach=2,3)
		tdaToggle <- gtkCheckButton()
		tdaToggle$setActive(TRUE)
		settingsTable$attach(tdaToggle, left.attach=1,2, top.attach=2,3)
		
		fragmentationLabel <- gtkLabel('Fragmentation method:')
		fragmentationLabel['xalign'] <- 1
		settingsTable$attach(fragmentationLabel, left.attach=0,1, top.attach=3,4)
		fragmentation <- gtkComboBoxNewText()
		sapply(c('From spectrum', 'CID', 'ETD', 'HCD', 'Merge'), fragmentation$appendText)
		fragmentation$setActive(0)
		settingsTable$attach(fragmentation, left.attach=1,2, top.attach=3,4)
		
		instrumentLabel <- gtkLabel('Instrument:')
		instrumentLabel['xalign'] <- 1
		settingsTable$attach(instrumentLabel, left.attach=0,1, top.attach=4,5)
		instrument <- gtkComboBoxNewText()
		sapply(c('Low-res LCQ/LTQ', 'High-res LTQ', 'TOF'), instrument$appendText)
		instrument$setActive(0)
		settingsTable$attach(instrument, left.attach=1,2, top.attach=4,5)
		
		enzymeLabel <- gtkLabel('Enzyme:')
		enzymeLabel['xalign'] <- 1
		settingsTable$attach(enzymeLabel, left.attach=0,1, top.attach=5,6)
		enzyme <- gtkComboBoxNewText()
		sapply(c('Unspecific', 'Trypsin', 'Chymotrypsin', 'Lys-C', 'Lys-N', 'Glu-C', 'Arg-C', 'Asp-N', 'alphaLP', 'None'), enzyme$appendText)
		enzyme$setActive(1)
		settingsTable$attach(enzyme, left.attach=1,2, top.attach=5,6)
		
		terminiLabel <- gtkLabel('Tolerable termini:')
		terminiLabel['xalign'] <- 1
		settingsTable$attach(terminiLabel, left.attach=0,1, top.attach=6,7)
		terminiBox <- gtkHBox(FALSE, 5)
			termini <- gtkSpinButton(min=0, max=2, step=1)
			termini$setValue(2)
			termini['height-request'] <- 24
			termini['width-chars'] <- 3
			termini$setAlignment(1)
			terminiBox$packStart(termini, expand=FALSE, fill=FALSE)
		settingsTable$attach(terminiBox, left.attach=1,2, top.attach=6,7)
		
		modificationLabel <- gtkLabel('Custom modification file:')
		modificationLabel['xalign'] <- 1
		settingsTable$attach(modificationLabel, left.attach=0,1, top.attach=7,8)
		modificationBox <- gtkHBox(FALSE, 5)
			modification <- gtkEntry()
			modification['height-request'] <- 24
			modification['width-chars'] <- 15
			modification$setText('Modification file...')
			gSignalConnect(modification, 'changed', f=function(widget, ...){
						validButtons()
					})
			modificationImport <- gtkButton(label='Browse')
			modificationImport['width-request'] <- 65
			gSignalConnect(modificationImport, 'clicked', f=function(widget, ...){
						file <- gtkFileChooserDialog(title='Select MS datafiles', parent=mainWindow, action='open', 'gtk-ok', GtkResponseType['ok'], 'gtk-cancel', GtkResponseType['cancel'], show=FALSE)
						gSignalConnect(file, 'response', f=function(dialog, response, data){
									if(response == GtkResponseType['ok']){
										file <- as.character(dialog$getFilenames())
										modification$setText(file)
									}
									dialog$destroy()
								})
						file$show()
					})
			modificationBox$packStart(modification)
			modificationBox$packStart(modificationImport, expand=FALSE, fill=FALSE)
		settingsTable$attach(modificationBox, left.attach=1,2, top.attach=7,8)
		
		lengthLabel <- gtkLabel('Peptide length:')
		lengthLabel['xalign'] <- 1
		settingsTable$attach(lengthLabel, left.attach=0,1, top.attach=8,9)
		lengthBox <- gtkHBox(FALSE, 5)
			lengthLow <- gtkSpinButton(min=1, max=100, step=1)
			lengthLow['height-request'] <- 24
			lengthLow['width-chars'] <- 3
			lengthLow$setAlignment(1)
			lengthLow$setValue(6)
			gSignalConnect(lengthLow, 'value-changed', f=function(widget, ...){
						min <- widget$getValue()
						if(min > lengthHigh$getValue()){
							lengthHigh$setRange(min, 100)
						}
					})
			lengthHigh <- gtkSpinButton(min=6, max=100, step=1)
			lengthHigh['height-request'] <- 24
			lengthHigh['width-chars'] <- 3
			lengthHigh$setAlignment(1)
			lengthHigh$setValue(40)
			lengthBox$packStart(lengthLow, expand=FALSE, fill=FALSE)
			lengthBox$packStart(gtkLabel(' - '), expand=FALSE, fill=FALSE)
			lengthBox$packStart(lengthHigh, expand=FALSE, fill=FALSE)
		settingsTable$attach(lengthBox, left.attach=1,2, top.attach=8,9)
		
		chargeLabel <- gtkLabel('Peptide charge:')
		chargeLabel['xalign'] <- 1
		settingsTable$attach(chargeLabel, left.attach=0,1, top.attach=9,10)
		chargeBox <- gtkHBox(FALSE, 5)
			chargeLow <- gtkSpinButton(min=1, max=10, step=1)
			chargeLow['height-request'] <- 24
			chargeLow['width-chars'] <- 3
			chargeLow$setAlignment(1)
			chargeLow$setValue(2)
			gSignalConnect(chargeLow, 'value-changed', f=function(widget, ...){
						min <- widget$getValue()
						if(min > chargeHigh$getValue()){
							chargeHigh$setRange(min, 10)
						}
					})
			chargeHigh <- gtkSpinButton(min=2, max=10, step=1)
			chargeHigh['height-request'] <- 24
			chargeHigh['width-chars'] <- 3
			chargeHigh$setAlignment(1)
			chargeHigh$setValue(3)
			chargeBox$packStart(chargeLow, expand=FALSE, fill=FALSE)
			chargeBox$packStart(gtkLabel(' - '), expand=FALSE, fill=FALSE)
			chargeBox$packStart(chargeHigh, expand=FALSE, fill=FALSE)
		settingsTable$attach(chargeBox, left.attach=1,2, top.attach=9,10)
		
		matchLabel <- gtkLabel('Matches per spectrum:')
		matchLabel['xalign'] <- 1
		settingsTable$attach(matchLabel, left.attach=0,1, top.attach=10,11)
		matchBox <- gtkHBox(FALSE, 5)
			match <- gtkSpinButton(min=1, max=100, step=1)
			match['height-request'] <- 24
			match['width-chars'] <- 3
			match$setAlignment(1)
			matchBox$packStart(match, expand=FALSE, fill=FALSE)
		settingsTable$attach(matchBox, left.attach=1,2, top.attach=10,11)

analyzeBox <- gtkHBox(FALSE, 5)
controlGroup$packStart(analyzeBox, expand=FALSE, fill=FALSE)
analyzeProgress <- gtkProgressBar()
analyzeProgress['height-request'] <- 24
analyzeBox$packStart(analyzeProgress)
analyzeButton <- gtkButton(label='Analyze')
analyzeButton['width-request'] <- 65
analyzeButton['sensitive'] <- FALSE
gSignalConnect(analyzeButton, 'clicked', f=function(widget, ...){
			analyzeProgress$setFraction(0)
			analyzeButton['sensitive'] <- FALSE
			while(gtkEventsPending()){
				gtkMainIteration()
			}
			datafiles <- filelist[]
			settings <- collectSettings()
			settings$database <- database$getText()
			settings$saveMZid <- TRUE
			settings$showDecoy <- TRUE
			for(i in 1:length(datafiles)){
				settings$file <- datafiles[i]
				analyzeProgress$setText(basename(settings$file))
				while(gtkEventsPending()){
					gtkMainIteration()
				}
				results <- do.call('MSGFplus', settings)
				peptide <- strsplit(as.character(results$Peptide), '\\.')
				peptide <- unlist(peptide)[seq(2, length(peptide)*3, by=3)]
				peptide <- gsub('\\+([\\d]|,)+', '', peptide, perl=TRUE)
				results$realPeptide <- peptide
				results$'.vis' <- rawFilter(results, sample=sampleFilter(), FDR=FDR$getValue())
				rawResults$appendRows(results)
				rawUpdated()
				mzidFile <- paste(sapply(strsplit(datafiles[i],"\\."), function(x) paste(x[1:(length(x)-1)], collapse=".")), '.mzid', sep='')
				sampleInfo <- getRunInfoFromMZID(mzidFile)
				updateSampleList(add=sampleInfo)
				analyzeProgress$setFraction(i/length(datafiles))
				while(gtkEventsPending()){
					gtkMainIteration()
				}
			}
			analyzeProgress$setText('Done...')
			analyzeButton['sensitive'] <- TRUE
		})
analyzeBox$packStart(gtkHBox(), expand=TRUE, fill=TRUE)
analyzeBox$packStart(analyzeButton, expand=FALSE, fill=FALSE)

resultBox <- gtkVBox(FALSE, 5)
resultBox$setSensitive(FALSE)
mainGroup$packStart(resultBox)

filterFrame <- gtkFrame(label='Filters')
filterBox <- gtkHBox(FALSE, 5)
filterBox['border-width'] <- 5
filterFrame$add(filterBox)
resultBox$packStart(filterFrame, expand=FALSE, fill=TRUE)
	filterBox$packStart(gtkLabel('Sample:'), expand=FALSE, fill=FALSE)
	sampleSelect <- gtkComboBoxNewWithModel(samples)
	sampleRenderer <- gtkCellRendererText()
	sampleSelect$packStart(sampleRenderer)
	sampleSelect$addAttribute(sampleRenderer, 'text', 0)
	sampleSelect$setActive(0)
	gSignalConnect(sampleSelect, 'changed', f=function(widget, ...){
				filterUpdated(FDRChange=FALSE)
			})
	filterBox$packStart(sampleSelect)
	filterBox$packStart(gtkHBox(), expand=TRUE, fill=TRUE)
	filterBox$packStart(gtkLabel('FDR cutoff:'), expand=FALSE, fill=FALSE)
	FDR <- gtkSpinButton(min=0, max=1, step=0.01)
	FDR['width-chars'] <- 5
	FDR['height-request'] <- 24
	FDR$setValue(0.01)
	gSignalConnect(FDR, 'value-changed', f=function(widget, ...){
				filterUpdated(SampleChange=FALSE)
			})
	filterBox$packStart(FDR, expand=FALSE, fill=FALSE)

resultNotebook <- gtkNotebook()
resultNotebook['width-request'] <- 700
resultBox$packStart(resultNotebook)
	summaryTab <- gtkHBox(FALSE, 5)
	summaryTab['border-width'] <- 5
	summaryLeft <- gtkVBox(FALSE, 5)
	summaryTab$packStart(summaryLeft)
	resultNotebook$appendPage(summaryTab, gtkLabel('Result summary'))
	FDRplot <- gtkDrawingArea()
	summaryLeft$packStart(FDRplot)
	#pepNumberPlot <- gtkDrawingArea()
	#summaryLeft$packStart(pepNumberPlot)
	#summaryRight <- gtkVBox(FALSE, 5)
	#summaryRight['width-request'] <- 100
	#summaryTab$packStart(summaryRight)
	
	summaryFrame <- gtkFrame(label='Summary')
	summaryLeft$packStart(summaryFrame, expand=FALSE, fill=FALSE)
	summaryTable <- gtkTable(rows=2, columns=2, homogeneous=FALSE)
	summaryTable['border-width'] <- 5
	summaryTable$setColSpacings(10)
	summaryTable$setRowSpacings(10)
	
	summaryFrame$add(summaryTable)
		nPepIDLabel <- gtkLabel('Identified peptides:')
		nPepIDLabel['xalign'] <- 1
		summaryTable$attach(nPepIDLabel, left.attach=0,1, top.attach=0,1)
		nPepID <- gtkLabel('')
		nPepID['xalign'] <- 0
		summaryTable$attach(nPepID, left.attach=1,2, top.attach=0,1)
		
		nProtIDLabel <- gtkLabel('Identified proteins:')
		nProtIDLabel['xalign'] <- 1
		summaryTable$attach(nProtIDLabel, left.attach=0,1, top.attach=1,2)
		nProtID <- gtkLabel('')
		nProtID['xalign'] <- 0
		summaryTable$attach(nProtID, left.attach=1,2, top.attach=1,2)
		
		protIntersectLabel <- gtkLabel('Protein intersection:')
		protIntersectLabel['xalign'] <- 1
		summaryTable$attach(protIntersectLabel, left.attach=0,1, top.attach=2,3)
		protIntersect <- gtkLabel('')
		protIntersect['xalign'] <- 0
		summaryTable$attach(protIntersect, left.attach=1,2, top.attach=2,3)
		
		chargeLabel <- gtkLabel('Charge:')
		chargeLabel['xalign'] <- 1
		chargeLabel['yalign'] <- 0
		summaryTable$attach(chargeLabel, left.attach=0,1, top.attach=3,4)
		charge <- gtkLabel('\t(min)\n\t(median)\n\t(max)')
		charge['xalign'] <- 0
		summaryTable$attach(charge, left.attach=1,2, top.attach=3,4)
	
	peptideTab <- gtkHBox(FALSE, 5)
	peptideTab['border-width'] <- 5
	peptideScroll <- gtkScrolledWindow()
		peptideTab$packStart(peptideScroll)
		peptide <- gtkTreeView(peptides)
		peptide$insertColumnWithAttributes(position=-1, title='Peptides', cell=gtkCellRendererText(), text=0)
		peptideSelection <- peptide$getSelection()
		gSignalConnect(peptide, 'cursor-changed', function(...){
					peptideSelected <- with(peptideSelection$getSelected(), model$getValue(iter, 0)$value)
					pepInfo <- createPeptideView(peptideSelected)
					if(!is.null(peptideViewScroll$getChild())) peptideViewScroll$getChild()$destroy()
					peptideViewScroll$add(createPeptideView(peptideSelected))
				})
		peptideScroll$add(peptide)
	peptideViewFrame <- gtkFrame('Information')
	peptideTab$packStart(peptideViewFrame)
	peptideViewScroll <- gtkScrolledWindow()
	peptideViewScroll['border-width'] <- 5
	peptideViewScroll$setPolicy(hscrollbar.policy=GtkPolicyType[2], vscrollbar.policy=GtkPolicyType[1])
	peptideViewFrame$add(peptideViewScroll)
	resultNotebook$appendPage(peptideTab, gtkLabel('Peptides'))
	
	proteinTab <- gtkVBox(FALSE, 5)
	proteinTab['border-width'] <- 5
	proteinTopBox <- gtkHBox(FALSE, 5)
	proteinTab$packStart(proteinTopBox)
	proteinScroll <- gtkScrolledWindow()
		proteinTopBox$packStart(proteinScroll)
		protein <- gtkTreeView(proteins)
		protein$insertColumnWithAttributes(position=-1, title='Protein', cell=gtkCellRendererText(), text=0)
		proteinScroll$add(protein)
		proteinSelect <- protein$getSelection()
		gSignalConnect(protein, 'cursor-changed', f=function(...){
					proteinName <- with(proteinSelect$getSelected(), model$getValue(iter, 0)$value)
					proteinSequenceFromDB <- findProteinSequence(proteinName)
					if(sampleFilter() == 'All'){
						evidenceList <- findProteinEvidence(proteinName, FDR=FDR$getValue())
					} else {
						evidenceList <- findProteinEvidence(proteinName, sample=sampleFilter(), FDR=FDR$getValue())
					}
					evidenceLocation <- findProteinEvidenceLocation(proteinSequenceFromDB[[1]]$sequence, evidenceList, collate=FALSE)
					tempPeptideSelection$unselectAll()
					nPeptideValue$setText(length(evidenceList))
					coverageValue$setText(paste(sprintf('%.1f', 100*length(unique(unlist(lapply(evidenceLocation, function(x) x[1]:x[2]))))/nchar(proteinSequenceFromDB[[1]]$sequence)), '%'))
					protLengthValue$setText(nchar(proteinSequenceFromDB[[1]]$sequence))
					protMassValue$setText(paste(sprintf('%.2f', pepMass(proteinSequenceFromDB[[1]]$sequence)), 'Da'))
					protPIValue$setText(sprintf('%.2f', peppI(proteinSequenceFromDB[[1]]$sequence)))
					hydroValue$setText(sprintf('%.1f', Qval(proteinSequenceFromDB[[1]]$sequence)))
					tempPeptideList$setFrame(data.frame(Peptides=paste(sapply(evidenceLocation, '[', 1), '-', sapply(evidenceLocation, '[', 2), ': ', names(evidenceList), sep=''), stringsAsFactors=FALSE))
					tempSampleList$setFrame(data.frame(Samples=unique(rawResults[rawResults[, 'Protein'] == proteinName, 'SpecFile']), stringsAsFactors=FALSE))
					sequenceLabel$setMarkup(formatProteinPango(proteinSequenceFromDB[[1]]$sequence, evidenceList))
				})
	proteinTopBoxRight <- gtkVBox(FALSE, 5)
	proteinTopBox$packStart(proteinTopBoxRight)
	proteinPeptideScroll <- gtkScrolledWindow()
		proteinTopBoxRight$packStart(proteinPeptideScroll)
		tempPeptide <- gtkTreeView(tempPeptideList)
		tempPeptide$insertColumnWithAttributes(position=-1, title='Peptide', cell=gtkCellRendererText(), text=0)
		gSignalConnect(tempPeptide, 'cursor-changed', f=function(...){
					proteinName <- with(proteinSelect$getSelected(), model$getValue(iter, 0)$value)
					proteinSequenceFromDB <- findProteinSequence(proteinName)
					if(sampleFilter() == 'All'){
						evidenceList <- findProteinEvidence(proteinName, FDR=FDR$getValue())
					} else {
						evidenceList <- findProteinEvidence(proteinName, sample=sampleFilter(), FDR=FDR$getValue())
					}
					evidenceSelection <- as.numeric(with(tempPeptideSelection$getSelected(), model$getPath(iter))$toString())+1
					sequenceLabel$setMarkup(formatProteinPango(proteinSequenceFromDB[[1]]$sequence, evidenceList, evidenceSelection))
				})
		gSignalConnect(tempPeptide, 'row-activated', function(...){
					selectedSequence <- strsplit(with(tempPeptideSelection$getSelected(), model$getValue(iter, 0)$value), ' ')[[1]][2]
					peptidePath <- gtkTreePathNewFromString(which(peptides[] == selectedSequence)-1)
					resultNotebook$setCurrentPage(resultNotebook$pageNum(peptideTab))
					peptideSelection$selectPath(peptidePath)
					peptide$scrollToCell(peptidePath, use.align=TRUE, row.align=0.5)
					gSignalEmit(peptide, 'cursor-changed')
				})
		tempPeptideSelection <- tempPeptide$getSelection()
		proteinPeptideScroll$add(tempPeptide)
	proteinSampleScroll <- gtkScrolledWindow()
		proteinTopBoxRight$packStart(proteinSampleScroll)
		tempSamples <- gtkTreeView(tempSampleList)
		tempSamples$insertColumnWithAttributes(position=-1, title='Samples', cell=gtkCellRendererText(), text=0)
		proteinSampleScroll$add(tempSamples)
	proteinBottomBox <- gtkHBox(FALSE, 5)
	proteinTab$packStart(proteinBottomBox, expand=FALSE, fill=FALSE)
	proteinSequence <- gtkFrame('Sequence:')
		proteinSequenceScroll <- gtkScrolledWindow()
			proteinSequenceScroll['border-width'] <- 5
			proteinSequenceScroll['vscrollbar-policy'] <- GtkPolicyType[1]
			proteinSequenceScroll['hscrollbar-policy'] <- GtkPolicyType[3]
			proteinSequence$add(proteinSequenceScroll)
		sequenceLabel <- gtkLabel()
			sequenceLabel['selectable'] <- TRUE
			sequenceLabel$setLineWrap(TRUE)
			sequenceLabel$setLineWrapMode('PANGO_WRAP_WORD')
			font_descr <- pangoFontDescriptionNew()
			font_descr$setFamily('courier')
			sequenceLabel$modifyFont(font_descr)
			proteinSequenceScroll$addWithViewport(sequenceLabel)
			proteinSequenceScroll$getChild()$setShadowType(GtkShadowType[1])
		proteinBottomBox$packStart(proteinSequence, expand=FALSE, fill=FALSE)
	proteinInformationFrame <- gtkFrame('Information:')
		proteinInformationTable <- gtkTable(rows=8, columns=2, homogeneous=FALSE)
			proteinInformationTable['border-width'] <- 5
			
			nPeptideLabel <- gtkLabel('Number of detected peptides:')
			nPeptideLabel['xalign'] <- 1
			proteinInformationTable$attach(nPeptideLabel, left.attach=0,1, top.attach=0,1)
			nPeptideValue <- gtkLabel()
			nPeptideValue['xalign'] <- 0
			proteinInformationTable$attach(nPeptideValue, left.attach=1,2, top.attach=0,1)
			
			coverageLabel <- gtkLabel('Coverage:')
			coverageLabel['xalign'] <- 1
			proteinInformationTable$attach(coverageLabel, left.attach=0,1, top.attach=1,2)
			coverageValue <- gtkLabel()
			coverageValue['xalign'] <- 0
			proteinInformationTable$attach(coverageValue, left.attach=1,2, top.attach=1,2)
			
			protLengthLabel <- gtkLabel('Length (residues):')
			protLengthLabel['xalign'] <- 1
			proteinInformationTable$attach(protLengthLabel, left.attach=0,1, top.attach=2,3)
			protLengthValue <- gtkLabel()
			protLengthValue['xalign'] <- 0
			proteinInformationTable$attach(protLengthValue, left.attach=1,2, top.attach=2,3)
			
			protMassLabel <- gtkLabel('Mass:')
			protMassLabel['xalign'] <- 1
			proteinInformationTable$attach(protMassLabel, left.attach=0,1, top.attach=3,4)
			protMassValue <- gtkLabel()
			protMassValue['xalign'] <- 0
			proteinInformationTable$attach(protMassValue, left.attach=1,2, top.attach=3,4)
			
			protPILabel <- gtkLabel('Isoelectric point (pI):')
			protPILabel['xalign'] <- 1
			proteinInformationTable$attach(protPILabel, left.attach=0,1, top.attach=4,5)
			protPIValue <- gtkLabel()
			protPIValue['xalign'] <- 0
			proteinInformationTable$attach(protPIValue, left.attach=1,2, top.attach=4,5)
			
			hydroLabel <- gtkLabel('Hydrophobicity (Q Value):')
			hydroLabel['xalign'] <- 1
			proteinInformationTable$attach(hydroLabel, left.attach=0,1, top.attach=5,6)
			hydroValue <- gtkLabel()
			hydroValue['xalign'] <- 0
			proteinInformationTable$attach(hydroValue, left.attach=1,2, top.attach=5,6)
			
			blastButton <- gtkButton('BLAST')
			blastButton['width-request'] <- 80
			blastButton['height-request'] <- 40
			gSignalConnect(blastButton, 'clicked', f=function(...){
						sequence <- findProteinSequence(with(proteinSelect$getSelected(), model$getValue(iter, 0)$value))
						RID <- blastPut(sequence[[1]]$sequence, db='nr')
						while(blastCheck(RID) != 'READY'){
							Sys.sleep(5)
							if(blastCheck(RID) != 'WAITING'){
								stop('Unknown problem occured')
							}
						}
						blastGet(RID)
					})
			proteinInformationTable$attach(blastButton, left.attach=0,2, top.attach=6,8, xoptions='', yoptions='')
			
			proteinInformationTable$setRowSpacings(15)
			proteinInformationTable$setColSpacings(5)
			proteinInformationFrame$add(proteinInformationTable)
		proteinBottomBox$packStart(proteinInformationFrame)
	resultNotebook$appendPage(proteinTab, gtkLabel('Proteins'))
	
	rawTab <- gtkVBox(FALSE, 5)
	rawTab['border-width'] <- 5
	rawExpander <- gtkExpander(label='Spectrum plot')
		rawPlotAndSettingsBox <- gtkVBox(FALSE, 5)
		rawExpander$add(rawPlotAndSettingsBox)
		rawSettingsExpander <- gtkExpander(label='Settings')
		rawSettingsExpander['border-width'] <- 5
			rawSettingsBox <- gtkHBox(FALSE, 5)
			rawSettingsBox$packStart(gtkHBox(), expand=TRUE, fill=TRUE)
			rawSettingsTable <- gtkTable(rows=3, columns=7, homogeneous=FALSE)
			rawSettingsTable$setColSpacing(0, 5)
			rawSettingsTable$setColSpacing(1, 70)
			rawSettingsTable$setColSpacing(2, 5)
				ppmLabel <- gtkLabel('ppm: ')
				ppmLabel['xalign'] <- 1
				rawSettingsTable$attach(ppmLabel, left.attach=0,1, top.attach=0,1)
				ppmPlotSetting <- gtkEntry()
				ppmPlotSetting$setText('20')
				ppmPlotSetting['height-request'] <- 24
				ppmPlotSetting['width-request'] <- 40
				gSignalConnect(ppmPlotSetting, 'activate', f=function(...){
							selection <- rawSelect$getSelected()
							sample <- with(selection, model$getValue(iter, 0)$value)
							if(!samples[samples[,'Samples']==sample, 'loaded']){
								raw$setData(key='msConnection', data=openMSfile(samples[samples[,'Samples']==sample, 'msfile']))
								samples[samples[,'Samples']==sample, 'loaded'] <- TRUE
							} else {}
							scan <- with(selection, model$getValue(iter, 2)$value)
							seq <- with(selection, model$getValue(iter, 16)$value)
							rawExpander$setExpanded(TRUE)
							plotSettings <- collectRawPlotSettings()
							plotSpectrum(data=raw$getData('msConnection'), scan=scan, seq=seq, device=spectrumPlot$getData('Device'), ppm=plotSettings$ppm, ionTrace=plotSettings$ionTrace, ions=plotSettings$ions, neutralLosses=plotSettings$neutralLosses)
						})
				rawSettingsTable$attach(ppmPlotSetting, left.attach=1,2, top.attach=0,1)
				
				neutralLabel <- gtkLabel('Neutral losses: ')
				neutralLabel['xalign'] <- 1
				rawSettingsTable$attach(neutralLabel, left.attach=0,1, top.attach=1,2)
				neutralPlotSetting <- gtkCheckButton()
				neutralPlotSetting$setActive(TRUE)
				gSignalConnect(neutralPlotSetting, 'toggled', f=function(...){
							selection <- rawSelect$getSelected()
							sample <- with(selection, model$getValue(iter, 0)$value)
							if(!samples[samples[,'Samples']==sample, 'loaded']){
								raw$setData(key='msConnection', data=openMSfile(samples[samples[,'Samples']==sample, 'msfile']))
								samples[samples[,'Samples']==sample, 'loaded'] <- TRUE
							} else {}
							scan <- with(selection, model$getValue(iter, 2)$value)
							seq <- with(selection, model$getValue(iter, 16)$value)
							rawExpander$setExpanded(TRUE)
							plotSettings <- collectRawPlotSettings()
							plotSpectrum(data=raw$getData('msConnection'), scan=scan, seq=seq, device=spectrumPlot$getData('Device'), ppm=plotSettings$ppm, ionTrace=plotSettings$ionTrace, ions=plotSettings$ions, neutralLosses=plotSettings$neutralLosses)
						})
				rawSettingsTable$attach(neutralPlotSetting, left.attach=1,2, top.attach=1,2)
				
				traceLabel <- gtkLabel('Ion trace: ')
				traceLabel['xalign'] <- 1
				rawSettingsTable$attach(traceLabel, left.attach=0,1, top.attach=2,3)
				tracePlotSetting <- gtkCheckButton()
				tracePlotSetting$setActive(TRUE)
				gSignalConnect(tracePlotSetting, 'toggled', f=function(...){
							selection <- rawSelect$getSelected()
							sample <- with(selection, model$getValue(iter, 0)$value)
							if(!samples[samples[,'Samples']==sample, 'loaded']){
								raw$setData(key='msConnection', data=openMSfile(samples[samples[,'Samples']==sample, 'msfile']))
								samples[samples[,'Samples']==sample, 'loaded'] <- TRUE
							} else {}
							scan <- with(selection, model$getValue(iter, 2)$value)
							seq <- with(selection, model$getValue(iter, 16)$value)
							rawExpander$setExpanded(TRUE)
							plotSettings <- collectRawPlotSettings()
							plotSpectrum(data=raw$getData('msConnection'), scan=scan, seq=seq, device=spectrumPlot$getData('Device'), ppm=plotSettings$ppm, ionTrace=plotSettings$ionTrace, ions=plotSettings$ions, neutralLosses=plotSettings$neutralLosses)
						})
				rawSettingsTable$attach(tracePlotSetting, left.attach=1,2, top.attach=2,3)
				
				ionsLabel <- gtkLabel('Ions:')
				ionsLabel['xalign'] <- 1
				rawSettingsTable$attach(ionsLabel, left.attach=2,3, top.attach=0,1)
				
				aIonLabel <- gtkLabel('a:')
				aIonLabel['xalign'] <- 1
				rawSettingsTable$attach(aIonLabel, left.attach=3,4, top.attach=0,1)
				aIonPlotSetting <- gtkCheckButton()
				aIonPlotSetting$setActive(TRUE)
				gSignalConnect(aIonPlotSetting, 'toggled', f=function(...){
							selection <- rawSelect$getSelected()
							sample <- with(selection, model$getValue(iter, 0)$value)
							if(!samples[samples[,'Samples']==sample, 'loaded']){
								raw$setData(key='msConnection', data=openMSfile(samples[samples[,'Samples']==sample, 'msfile']))
								samples[samples[,'Samples']==sample, 'loaded'] <- TRUE
							} else {}
							scan <- with(selection, model$getValue(iter, 2)$value)
							seq <- with(selection, model$getValue(iter, 16)$value)
							rawExpander$setExpanded(TRUE)
							plotSettings <- collectRawPlotSettings()
							plotSpectrum(data=raw$getData('msConnection'), scan=scan, seq=seq, device=spectrumPlot$getData('Device'), ppm=plotSettings$ppm, ionTrace=plotSettings$ionTrace, ions=plotSettings$ions, neutralLosses=plotSettings$neutralLosses)
						})
				rawSettingsTable$attach(aIonPlotSetting, left.attach=4,5, top.attach=0,1)
				
				bIonLabel <- gtkLabel('b:')
				bIonLabel['xalign'] <- 1
				rawSettingsTable$attach(bIonLabel, left.attach=3,4, top.attach=1,2)
				bIonPlotSetting <- gtkCheckButton()
				bIonPlotSetting$setActive(TRUE)
				gSignalConnect(bIonPlotSetting, 'toggled', f=function(...){
							selection <- rawSelect$getSelected()
							sample <- with(selection, model$getValue(iter, 0)$value)
							if(!samples[samples[,'Samples']==sample, 'loaded']){
								raw$setData(key='msConnection', data=openMSfile(samples[samples[,'Samples']==sample, 'msfile']))
								samples[samples[,'Samples']==sample, 'loaded'] <- TRUE
							} else {}
							scan <- with(selection, model$getValue(iter, 2)$value)
							seq <- with(selection, model$getValue(iter, 16)$value)
							rawExpander$setExpanded(TRUE)
							plotSettings <- collectRawPlotSettings()
							plotSpectrum(data=raw$getData('msConnection'), scan=scan, seq=seq, device=spectrumPlot$getData('Device'), ppm=plotSettings$ppm, ionTrace=plotSettings$ionTrace, ions=plotSettings$ions, neutralLosses=plotSettings$neutralLosses)
						})
				rawSettingsTable$attach(bIonPlotSetting, left.attach=4,5, top.attach=1,2)
				
				cIonLabel <- gtkLabel('c:')
				cIonLabel['xalign'] <- 1
				rawSettingsTable$attach(cIonLabel, left.attach=3,4, top.attach=2,3)
				cIonPlotSetting <- gtkCheckButton()
				gSignalConnect(cIonPlotSetting, 'toggled', f=function(...){
							selection <- rawSelect$getSelected()
							sample <- with(selection, model$getValue(iter, 0)$value)
							if(!samples[samples[,'Samples']==sample, 'loaded']){
								raw$setData(key='msConnection', data=openMSfile(samples[samples[,'Samples']==sample, 'msfile']))
								samples[samples[,'Samples']==sample, 'loaded'] <- TRUE
							} else {}
							scan <- with(selection, model$getValue(iter, 2)$value)
							seq <- with(selection, model$getValue(iter, 16)$value)
							rawExpander$setExpanded(TRUE)
							plotSettings <- collectRawPlotSettings()
							plotSpectrum(data=raw$getData('msConnection'), scan=scan, seq=seq, device=spectrumPlot$getData('Device'), ppm=plotSettings$ppm, ionTrace=plotSettings$ionTrace, ions=plotSettings$ions, neutralLosses=plotSettings$neutralLosses)
						})
				rawSettingsTable$attach(cIonPlotSetting, left.attach=4,5, top.attach=2,3)
				
				xIonLabel <- gtkLabel('x:')
				xIonLabel['xalign'] <- 1
				rawSettingsTable$attach(xIonLabel, left.attach=5,6, top.attach=0,1)
				xIonPlotSetting <- gtkCheckButton()
				gSignalConnect(xIonPlotSetting, 'toggled', f=function(...){
							selection <- rawSelect$getSelected()
							sample <- with(selection, model$getValue(iter, 0)$value)
							if(!samples[samples[,'Samples']==sample, 'loaded']){
								raw$setData(key='msConnection', data=openMSfile(samples[samples[,'Samples']==sample, 'msfile']))
								samples[samples[,'Samples']==sample, 'loaded'] <- TRUE
							} else {}
							scan <- with(selection, model$getValue(iter, 2)$value)
							seq <- with(selection, model$getValue(iter, 16)$value)
							rawExpander$setExpanded(TRUE)
							plotSettings <- collectRawPlotSettings()
							plotSpectrum(data=raw$getData('msConnection'), scan=scan, seq=seq, device=spectrumPlot$getData('Device'), ppm=plotSettings$ppm, ionTrace=plotSettings$ionTrace, ions=plotSettings$ions, neutralLosses=plotSettings$neutralLosses)
						})
				rawSettingsTable$attach(xIonPlotSetting, left.attach=6,7, top.attach=0,1)
				
				yIonLabel <- gtkLabel('y:')
				yIonLabel['xalign'] <- 1
				rawSettingsTable$attach(yIonLabel, left.attach=5,6, top.attach=1,2)
				yIonPlotSetting <- gtkCheckButton()
				yIonPlotSetting$setActive(TRUE)
				gSignalConnect(yIonPlotSetting, 'toggled', f=function(...){
							selection <- rawSelect$getSelected()
							sample <- with(selection, model$getValue(iter, 0)$value)
							if(!samples[samples[,'Samples']==sample, 'loaded']){
								raw$setData(key='msConnection', data=openMSfile(samples[samples[,'Samples']==sample, 'msfile']))
								samples[samples[,'Samples']==sample, 'loaded'] <- TRUE
							} else {}
							scan <- with(selection, model$getValue(iter, 2)$value)
							seq <- with(selection, model$getValue(iter, 16)$value)
							rawExpander$setExpanded(TRUE)
							plotSettings <- collectRawPlotSettings()
							plotSpectrum(data=raw$getData('msConnection'), scan=scan, seq=seq, device=spectrumPlot$getData('Device'), ppm=plotSettings$ppm, ionTrace=plotSettings$ionTrace, ions=plotSettings$ions, neutralLosses=plotSettings$neutralLosses)
						})
				rawSettingsTable$attach(yIonPlotSetting, left.attach=6,7, top.attach=1,2)
				
				zIonLabel <- gtkLabel('z:')
				zIonLabel['xalign'] <- 1
				rawSettingsTable$attach(zIonLabel, left.attach=5,6, top.attach=2,3)
				zIonPlotSetting <- gtkCheckButton()
				gSignalConnect(zIonPlotSetting, 'toggled', f=function(...){
							selection <- rawSelect$getSelected()
							sample <- with(selection, model$getValue(iter, 0)$value)
							if(!samples[samples[,'Samples']==sample, 'loaded']){
								raw$setData(key='msConnection', data=openMSfile(samples[samples[,'Samples']==sample, 'msfile']))
								samples[samples[,'Samples']==sample, 'loaded'] <- TRUE
							} else {}
							scan <- with(selection, model$getValue(iter, 2)$value)
							seq <- with(selection, model$getValue(iter, 16)$value)
							rawExpander$setExpanded(TRUE)
							plotSettings <- collectRawPlotSettings()
							plotSpectrum(data=raw$getData('msConnection'), scan=scan, seq=seq, device=spectrumPlot$getData('Device'), ppm=plotSettings$ppm, ionTrace=plotSettings$ionTrace, ions=plotSettings$ions, neutralLosses=plotSettings$neutralLosses)
						})
				rawSettingsTable$attach(zIonPlotSetting, left.attach=6,7, top.attach=2,3)
			
			rawSettingsBox$packStart(rawSettingsTable)
			rawSettingsBox$packStart(gtkHBox(), expand=TRUE, fill=TRUE)
			rawSettingsExpander$add(rawSettingsBox)
		rawPlotAndSettingsBox$packStart(rawSettingsExpander)
		rawTab$packStart(rawExpander, expand=FALSE, fill=FALSE)
		rawPlotBox <- gtkVBox(FALSE, 5)
		rawPlotBox['height-request'] <- 400
		rawPlotAndSettingsBox$packStart(rawPlotBox)
		rawExpander$setExpanded(TRUE)
		spectrumPlot <- gtkDrawingArea()
		rawPlotBox$packStart(spectrumPlot)
	rawScroll <- gtkScrolledWindow()
		rawTab$packStart(rawScroll)
		raw <- gtkTreeView(rawResultsFilter)
			for(i in 1:(ncol(rawResults)-2)){
				raw$insertColumnWithAttributes(position=-1, title=dimnames(rawResults)[[2]][i], cell=gtkCellRendererText(), text=i-1)
			}
		rawSelect <- raw$getSelection()
		rawScroll$add(raw)
		gSignalConnect(raw, 'row-activated', f=function(treeview, path, column, user.data){
					model <- treeview$getModel()
					selection <- model$getIter(path)$iter
					sample <- treeview$getModel()$getValue(selection, 0)$value
					if(!user.data$samples[user.data$samples[,'Samples']==sample, 'loaded']){
						treeview$setData(key='msConnection', data=openMSfile(user.data$samples[user.data$samples[,'Samples']==sample, 'rawfile']))
						user.data$samples[user.data$samples[,'Samples']==sample, 'loaded'] <- TRUE
					} else {}
					scan <- model$getValue(selection, 2)$value
					seq <- model$getValue(selection, 16)$value
					user.data$expander$setExpanded(TRUE)
					plotSettings <- collectRawPlotSettings()
					plotSpectrum(data=treeview$getData('msConnection'), scan=scan, seq=seq, device=user.data$spectrumPlot$getData('Device'), ppm=plotSettings$ppm, ionTrace=plotSettings$ionTrace, ions=plotSettings$ions, neutralLosses=plotSettings$neutralLosses)
					treeview$scrollToCell(path)
				}, data=list(samples=samples, expander=rawExpander, spectrumPlot=spectrumPlot))		
		gSignalConnect(rawSelect, 'changed', f=function(select, user.data){
					if(rawExpander$getExpanded()){
						selection <- select$getSelected()
						sample <- with(selection, model$getValue(iter, 0)$value)
						if(!user.data$samples[user.data$samples[,'Samples']==sample, 'loaded']){
							user.data$raw$setData(key='msConnection', data=openMSfile(user.data$samples[user.data$samples[,'Samples']==sample, 'rawfile']))
							user.data$samples[user.data$samples[,'Samples']==sample, 'loaded'] <- TRUE
						} else {}
						scan <- with(selection, model$getValue(iter, 2)$value)
						seq <- with(selection, model$getValue(iter, 16)$value)
						plotSettings <- collectRawPlotSettings()
						plotSpectrum(data=user.data$raw$getData('msConnection'), scan=scan, seq=seq, device=user.data$spectrumPlot$getData('Device'), ppm=plotSettings$ppm, ionTrace=plotSettings$ionTrace, ions=plotSettings$ions, neutralLosses=plotSettings$neutralLosses)
					} else {}
				}, data=list(samples=samples, raw=raw, spectrumPlot=spectrumPlot))
	resultNotebook$appendPage(rawTab, gtkLabel('Raw data'))

statusbar <- gtkStatusbar()
topGroup$packStart(statusbar, expand=FALSE, fill=FALSE)
mainWindow$show()

asCairoDevice(FDRplot)
FDRplot$setData(key='Device', data=dev.cur())
grid.newpage()
#asCairoDevice(pepNumberPlot)
#pepNumberPlot$setData(key='Device', data=dev.cur())
resultNotebook$setCurrentPage(3)
asCairoDevice(spectrumPlot)
spectrumPlot$setData(key='Device', data=dev.cur())
grid.newpage()
resultNotebook$setCurrentPage(0)
rawExpander$setExpanded(FALSE)