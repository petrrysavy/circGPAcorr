# here we download which mirnas are in which GO terms

library(miRBaseConverter)

source("mirnaGO.R")

# here we have all data we can get in mirMiRNABase programatically
miRNATab=getMiRNATable(version="v22",species="hsa")
wait <- TRUE

for(rowI in 1:nrow(miRNATab)) {
    map <- loadMapMiRNAGoTerms()
    row <- miRNATab[rowI,]
    if(row$Mature1_Acc %in% map$Mature1_Acc) { next }
    print(paste("querying", row$Mature1_Acc))
    goTerms <- miRNAGO(row$Mature1_Acc)
    if(is.null(goTerms) | length(goTerms) == 0) { next }
    newRows <- cbind(row$Mature1, goTerms)
    newRows <- cbind(row$Mature1_Acc, newRows)
    colnames(newRows) <- c("Mature1_Acc","Mature1","goTerm")
    map <- rbind(map, newRows)
    write.csv(as.data.frame(map), file='mirna-go-terms.csv')
    if(wait) { Sys.sleep(50) }
}
