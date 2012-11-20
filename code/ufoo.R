setcolclass.df <- function(df, colclass, factor.levels) {
	ww <- options()$warn
	options(warn=-1)
	toCls <- function(x, cls) { do.call(paste("as", cls, sep = "."), list(x)) }
	df <- replace(df, , Map(toCls, x=df, cls=colclass))
	options(warn=ww)
	iix <- FALSE
	if(!missing(factor.levels)) { iix <- colclass == "factor" & !is.null(factor.levels) }
	if(any(iix)) {
		for(i in which(iix)) { levels(df[[i]]) <- factor.levels[[i]] }
	}
	return(df)
}


## courtesy of Matthew McCall
celfileDateHour <- function(filename) {
	require(affyio)
	h <- affyio::read.celfile.header(filename, info="full")
	#ddate <- grep("/", strsplit(h$DatHeader, " ")[[1]], value=TRUE)
	#ddate <- strsplit(ddate, split="/")[[1]]
	#CC <- ifelse(substr(ddate[3],1,1)=="9", "19", "20")
	ddate <- strsplit(h$ScanDate, "T")[[1]]
	names(ddate) <- c("day", "hour")
	return(ddate)
}

celfileChip <- function(filename) {
	require(affyio)
	h <- affyio::read.celfile.header(filename, info="full")
	return(as.character(h$cdfName))
}



