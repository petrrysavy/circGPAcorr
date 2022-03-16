library(multiMiR)

source('buildGraph.R')

getSingle <- function(mirna) {
    df <- get_multimir(mirna = mirna, table="validated", summary = TRUE)@data
    df <- df[!duplicated(df[, c("mature_mirna_id", "target_symbol")]), ]
    df <- df[,c("mature_mirna_id", "target_symbol")]
    df <- df[df[,"target_symbol"] %in% all_gene_list,]
    return(df)
}

df <- getSingle('hsa-let-7a-5p')

for(mirna in all_mirna_list[all_mirna_list != "hsa-let-7a-5p"]) {
    print(paste("Querying mirna ", mirna))
    tryCatch({
        df <- rbind(df, getSingle(mirna))
        write.csv(df, file='mirna-mrna-interactions.csv')
    }, error = function(e) {
        print(e)
        print("unable to get results")
    })
}

#     validated1 <- get_multimir(org     = "hsa",
#                          mirna   = "hsa-let-7a-5p",
#                          target  = all_gene_list,
#                          table   = "validated",
#                          summary = TRUE,
#                          predicted.cutoff.type = "p",
#                          predicted.cutoff      = 10,
#                          use.tibble = TRUE)
#     saveRDS(validated, file="validatedRaw.RDS")
#validated <- select(all_interactions, keytype = "type", keys = "validated", columns = columns(all_interactions))
#     validated <- validated[!duplicated(validated[, c("mature_mirna_id", "target_entrez")]), ]
