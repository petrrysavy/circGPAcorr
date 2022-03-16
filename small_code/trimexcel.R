library(openxlsx)
library(stringr)

for(file in list.files(".", pattern="hsa_circ_\\d+-short\\.xlsx")) {
    print(file)
    res <- read.xlsx(file)
    res <- res[res$fdr < 0.05,]
    write.xlsx(res, str_replace(file, "short", "trimmed"), overwrite=TRUE)
}
