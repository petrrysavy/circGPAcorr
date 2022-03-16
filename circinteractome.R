# this code is used to download interacting miRNAs from the CircInteractome database
# unfortunatelly, CircInteractome is the only database working now (CircNet is down)
# and it has only "User-Input web interface"

# for http request (was preinstalled on RCI)
library("httr")
# for regular expressions
library("stringr")

#r <- GET("https://circinteractome.nia.nih.gov/api/v2/mirnasearch?circular_rna_query=hsa_circ_0007848&mirna_query=&submit=miRNA+Target+Search")
#content(r, "text")
# hsa-miR-xxxx
# regex hsa-(mir|miR|miRNA|mirna|microRNA|microrna|let|lin)-[[:alnum:]\\-]{1,10}

circInteractome <- function(circRNA) {
    print(paste("Downloading circRNA interactions for", circRNA))
    url <- paste("https://circinteractome.nia.nih.gov/api/v2/mirnasearch?circular_rna_query=", circRNA, "&mirna_query=&submit=miRNA+Target+Search", sep="")
    r <- GET(url)
    if (r$status_code != 200) {
        stop("problem connecting to the circinteractome database")
    }
    mirnas <- str_extract_all(content(r, "text"), "hsa-(mir|miR|miRNA|microRNA|let|lin)-[[:alnum:]\\-]{1,10}")[[1]]
    unique(mirnas)
}

circInteractomeList <- function(circRNAs) {
    result <- data.frame()
    
    for(circRNA in circRNAs) {
        miRNAs <- circInteractome(circRNA)
        newRows <- cbind(circRNA, miRNAs)
        result <- rbind(result, newRows)
    }
    
    result
}

loadMap <- function() {
    tryCatch({
        map <- read.csv('circrna-mirna-interactions.csv')
        subset(map, select = -c(X) )
    }, error = function(e) {
        data.frame(circRNA=character(0),miRNAs=character(0))
    })
}

circInteractomeListCached <- function(circRNAs, wait=TRUE) {
    for(circRNA in circRNAs) {
        map <- loadMap()
        if(circRNA %in% map$circRNA) { next }
        print(paste("querying", circRNA))
        miRNAs <- circInteractome(circRNA)
        newRows <- cbind(circRNA, miRNAs)
        map <- rbind(map, newRows)
        write.csv(as.data.frame(map), file='circrna-mirna-interactions.csv')
        if(wait) { Sys.sleep(50) }
    }
    
    map <- loadMap()
    return(map[map[,1] %in% circRNAs,])
}
