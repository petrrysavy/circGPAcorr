source('annotate.R')

circs <- c("hsa_circ_0000005", "hsa_circ_0002816", "hsa_circ_0000006", "hsa_circ_0001897", "hsa_circ_0004624", "hsa_circ_0024604", "hsa_circ_0000228", "hsa_circ_0001540", "hsa_circ_0044708", "hsa_circ_0003583")

for(circ in circs) {
    muCirc <- buildMuCirc(circ)
    ImmuCirc <- Im %*% muCirc
    print(circ)
    print(sum(ImmuCirc))
}
