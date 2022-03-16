# for http request (was preinstalled on RCI)
library("httr")
# for regular expressions
library("stringr")

miRNAGO <- function(miRNABaseAccession) {
    print(paste("Downloading miRNA GO terms annotations for ", miRNABaseAccession))
    url <- paste("http://mirbase.org/cgi-bin/mature.pl?mature_acc=", miRNABaseAccession, sep="")
    r <- GET(url)
    if (r$status_code != 200) {
        stop("problem connecting to the miRNABAse database")
    }
    rnaCentral <- unique(str_extract_all(content(r, "text"), "URS[0-9a-fA-F]+_[0-9a-fA-F]+")[[1]])
    if(length(rnaCentral) != 1) {
        print(paste("there is no RNA-central identifier for", miRNABaseAccession))
        return(NULL)
    }
    print(paste("Accessing QuickGO for ", rnaCentral))
    url <- paste("https://www.ebi.ac.uk/QuickGO/services/annotation/search?geneProductId=",rnaCentral,"&limit=100&page=1", sep="")
    r <- GET(url)
    if (r$status_code != 200) {
        stop("problem connecting to the QuickGO database")
    }
    goTerms <- unique(str_extract_all(content(r, "text"), "GO:\\d{2,15}")[[1]])
    return(goTerms)
}

loadMapMiRNAGoTerms <- function() {
    tryCatch({
        map <- read.csv('mirna-go-terms.csv')
        map <- subset(map, select = -c(X) )
        colnames(map) <- c("Mature1_Acc","Mature1","goTerm")
        return(map)
    }, error = function(e) {
        data.frame(Mature1_Acc=character(0),Mature1=character(0),goTerm=character(0))
    })
}
