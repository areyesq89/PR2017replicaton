library(Biobase)
library(PharmacoGx)

GDSC <- downloadPSet("GDSC")
CCLE <- downloadPSet("CCLE")
common <- intersectPSet(
    list( 'CCLE'=CCLE,
         'GDSC'=GDSC ),
    intersectOn=c("cell.lines", "drugs"),
    strictIntersect=TRUE)

rawSensData <- common[["CCLE"]]@sensitivity$raw
str(dimnames(rawSensData))
names( dimnames(rawSensData) ) <- c("drugCell", "doseID", "doseData")
spNames <- strsplit( dimnames( rawSensData )[[1]], "_" )
allCells <- sapply( spNames, "[[", 1 )
allDrugs <- sapply( spNames, "[[", 2 )

allData <- expand.grid(
    dimnames(rawSensData)[["doseID"]],
    allCells,
    allDrugs, stringsAsFactors=FALSE )
colnames(allData) <- c("doseID", "cellLine", "drug")

concatName <- with( allData,
                   paste( cellLine, drug, sep="_" ) )

allData$concentration <- NA
allData$viability <- NA

for( i in seq_len( nrow( allData ) )[1:20] ){
    allData$viability[i] <- rawSensData[ concatName[i],
            allData$doseID[i],
            "Viability"] 
    allData$concentration[i] <- rawSensData[ concatName[i],
            allData$doseID[i],
            "Dose"] 
}
