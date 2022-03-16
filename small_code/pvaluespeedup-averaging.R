library(stringr)

circs <- c("hsa_circ_0000005", "hsa_circ_0002816", "hsa_circ_0000006", "hsa_circ_0001897", "hsa_circ_0024604",  "hsa_circ_0000228", "hsa_circ_0001540", "hsa_circ_0044708", "hsa_circ_0003583")
sum <- 0
count <- 0
min <- 1e10
max <- 0

for(circ in circs) {
    circ <- str_replace_all(circ, "_", "-")
    print(circ)

    speedup <- read.csv(paste(circ, "-pvaluespeedup.dat", sep=""))
    sum <- sum + sum(speedup)
    min <- min(min, min(speedup))
    max <- max(max, max(speedup))
    count <- count + nrow(speedup)
}

print(sum/count)
# 2762.48
print(min)
print(max)
#15.88527
#7524.333

# 2847.669
# 15.88527
# 7524.333
