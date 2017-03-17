
library(Biobase)
library(PharmacoGx)

GDSC <- downloadPSet("GDSC")
CCLE <- downloadPSet("CCLE")
common <- intersectPSet(
    list( 'CCLE'=CCLE,
         'GDSC'=GDSC ),
    intersectOn=c("cell.lines", "drugs"),
    strictIntersect=TRUE)

rawToDataFrame <- function( dataset ){
    rawSensData <- common[[dataset]]@sensitivity$raw
    names( dimnames(rawSensData) ) <- c("drugCell", "doseID", "doseData")
    spNames <- strsplit( dimnames( rawSensData )[[1]], "_" )
    allCells <- unique( sapply( spNames, "[[", 1 ) )
    allDrugs <- unique( sapply( spNames, "[[", 2 ) )
    allData <- expand.grid(
        dimnames(rawSensData)[["doseID"]],
        allCells,
        allDrugs, stringsAsFactors=FALSE )
    colnames(allData) <- c("doseID", "cellLine", "drug")
    concatName <- with( allData,
                   paste( cellLine, drug, sep="_" ) )
    allData$concentration <- NA
    allData$viability <- NA
    for( i in seq_len( nrow( allData ) ) ){
        x <- dimnames(rawSensData)[[1]] %in% concatName[i]
        y <- dimnames(rawSensData)[[2]] %in% allData$doseID[i]
        if( any( x ) & any(y) ){
            allData$viability[i] <- rawSensData[x,y,"Viability"] 
            allData$concentration[i] <- rawSensData[x,y,"Dose"]
        }
    }
    allData <- na.omit(allData)
    allData$doseID <- as.factor( allData$doseID )
    allData$cellLine <- as.factor( allData$cellLine )
    allData$drug <- as.factor( allData$drug )
    allData$concentration <- as.numeric( allData$concentration )
    allData$viability <- as.numeric( allData$viability )
    allData$study <- dataset
    allData
}

rawSensitivityDf<- rbind( rawToDataFrame("CCLE"), rawToDataFrame("GDSC") )
rawSensitivityDf <- rawSensitivityDf[,c("cellLine", "drug", "doseID", "concentration", "viability", "study")]

sdList <- lapply( c("CCLE", "GDSC"), function(dataset){
    summarizedData <- common[[dataset]]@sensitivity$profiles
    keepCols <- c("ic50_published", "auc_published")
    summarizedData <- summarizedData[,keepCols]
    spNames <- as.data.frame( do.call(rbind, strsplit( rownames( summarizedData ), "_" )) )
    colnames( spNames ) <- c("cellLine", "drug")
    summarizedData <- cbind( spNames, summarizedData )
    colnames(summarizedData) <- gsub("_published", "", colnames( summarizedData ))
    rownames( summarizedData ) <- NULL
    summarizedData
} )
names(sdList) <- c("CCLE", "GDSC")

library(plyr)
joinedSumData <- join( sdList[[1]], sdList[[2]], by=c("cellLine", "drug") )
colnames(joinedSumData) <- c("cellLine", "drug", "ic50_CCLE", "auc_CCLE", "ic50_GDSC", "auc_GDSC")

write.table( joinedSumData, sep=",",
            quote=FALSE, col.names=TRUE,
            row.names=FALSE, file="summarizedPharmacoData.csv" )
write.table( rawSensitivityDf, sep=",",
            quote=FALSE, col.names=TRUE,
            row.names=FALSE, file="rawPharmacoData.csv" )

#library(dplyr)
#library(ggplot2)
#ggplot( dplyr:::filter( joinedSumData, drug == "Sorafenib" ), aes(-log10(ic50_GDSC), -log10(ic50_CCLE)) ) +
#    geom_point()
#dev.off()
