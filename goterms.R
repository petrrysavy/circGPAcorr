library(msigdbr)

if (file.exists("goterms.RData")) {
    load("goterms.RData")
} else {
    goTerms <- unlist(Term(GOTERM))
    goTermNames <- names(goTerms)

    #c5map <- read.table('goidmap.tsv', header=TRUE)
    m_df = msigdbr(species = "Homo sapiens") %>% dplyr::filter(gs_cat %in% c("C5"))
    m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name) # transform them into a form suitable for fgsea
    longnames <- paste("GO_", str_replace_all(toupper(goTerms), ' ', '_'), sep="")
    indices <- longnames %in% names(m_list)

    goTerms <- goTerms[indices]
    goTermNames <- goTermNames[indices]

    missing <- !names(m_list) %in% longnames
    missing <- names(m_list)[missing]
    goTermNames <- c(goTermNames, missing)

    missing <- setNames(as.list(missing), missing)
    goTerms <- c(goTerms, unlist(missing))
    #save(goTerms, goTermNames, longmanes, m_list, file = "goterms.RData")
}
