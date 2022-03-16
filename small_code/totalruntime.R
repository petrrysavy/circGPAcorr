library(openxlsx)

time <- 0
i <- 0

for(file in list.files(".", pattern="hsa_circ_\\d+\\.xlsx")) {
    print(i)
    i <- i + 1
    res <- read.xlsx(file)
    time <- time + sum(as.numeric(res$time))
}

print("Total time needed to calculate everything:")
print(time)
# 52947092
