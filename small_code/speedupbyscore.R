library(openxlsx)
library(ggplot2)
library(tikzDevice)
library(stringr)

circs <- c("hsa_circ_0000005", "hsa_circ_0002816", "hsa_circ_0000006", "hsa_circ_0001897", "hsa_circ_0024604", "hsa_circ_0000228", "hsa_circ_0001540", "hsa_circ_0044708", "hsa_circ_0003583")
connect <- c(34043, 34113, 33470, 302, 325, 1771, 5037, 1837, 1234)

allDataL <-  as.data.frame(matrix(nrow=0, ncol=3))

for(circ in circs) {
    res <- read.xlsx(paste(circ, ".xlsx", sep=""))
    pval <- as.numeric(res$pvalue)
    bootT <- as.numeric(res$bootstrapTime)
    pvalT <- as.numeric(res$pvalueTime)
    goSize <- as.numeric(res$goSize)
    print("outliers")
    print(sum(bootT > 1000)) # outlier, RCI probably put the job away 
    bootTinlier <- bootT < 1000
    included <- (pval < 1.0 | is.na(pval)) & bootTinlier & goSize < 1000
    pvalT <- pvalT[included]
    gosq <- goSize[included] * goSize[included]
    
    #circTeX <- str_replace_all(circ, "_", "\\\\_")
    circTeX <- circ
    print(sum(is.na(pvalT)))
    allDataL <- rbind(allDataL, cbind(gosq, pvalT, circTeX))
}
colnames(allDataL) <- c("gosq", "runtime", "circRNA")

allDataL$gosq <- as.numeric(allDataL$gosq)
allDataL$runtime <- as.numeric(allDataL$runtime)

maxim <- max(allDataL$runtime)
maxim2 <- max(allDataL$gosq)
#tikz(file = "runtimebyescore.tex", width = 15, height = 15, standAlone=TRUE)
plot <- ggplot(allDataL, aes(x=gosq, y=runtime, color=circRNA)) +
        geom_point(size = 0.5) +
        theme_bw() +
        theme(legend.position = "bottom") +
        #$\\mathbb{E}(s(c, g)) * \\| \\mathbf{g}^m \\|_1$
        labs( x = "Number of annotated mRNAs squared", y = "p-value calculation time (sec)") +
        scale_y_continuous(breaks=floor(seq(from=0, to=maxim, length.out=11))) +
        #scale_x_continuous(breaks=floor(seq(from=0, to=maxim2, length.out=11)))
        scale_x_continuous(breaks=floor(seq(from=0, to=1000000, length.out=11))) +
        expand_limits(x = 0, y = 0)

ggsave("runtimebyescoregg2.pdf", plot=plot, width = 22, height = 20, units = "cm")
ggsave("runtimebyescoregg2.png", plot=plot, type="cairo", width = 22, height = 20, dpi = 150, units = "cm")
#print(plot)
#dev.off()
