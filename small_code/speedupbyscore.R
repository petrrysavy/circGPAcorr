library(openxlsx)
library(ggplot2)
library(tikzDevice)
library(stringr)

circs <- c("hsa_circ_0000005", "hsa_circ_0002816", "hsa_circ_0000006", "hsa_circ_0001897", "hsa_circ_0024604", "hsa_circ_0000228", "hsa_circ_0001540", "hsa_circ_0044708", "hsa_circ_0003583")
connect <- c(34043, 34113, 33470, 302, 325, 1771, 5037, 1837, 1234)

allData <- as.data.frame(matrix(nrow=0, ncol=3))
allDataEs <- as.data.frame(matrix(nrow=0, ncol=3))
allDataC <- as.data.frame(matrix(nrow=0, ncol=3))
allDataL <-  as.data.frame(matrix(nrow=0, ncol=3))

for(circ in circs) {
    res <- read.xlsx(paste(circ, "-short.xlsx", sep=""))
    pval <- as.numeric(res$pvalue)
    pvalT <- as.numeric(res$pvalueTime)
    bootT <- as.numeric(res$bootstrapTime)
    score <- as.numeric(res$score)
    escore <- as.numeric(res$expectedScore)
    goSize <- as.numeric(res$goSize)
    print("outliers")
    print(sum(bootT > 300)) # outlier, RCI probably put the job away 
    bootTinlier <- bootT < 300
    included <- pval < 1.0 & bootTinlier
    bootT <- bootT[included]
    pvalT <- pvalT[included]
    score <- score[included]
    escore <- escore[included] * goSize[included]
    
    speedup <- bootT/pvalT
    
    #circTeX <- str_replace_all(circ, "_", "\\\\_")
    circTeX <- circ
    allData <- rbind(allData, cbind(score, speedup, circTeX))
    allDataEs <- rbind(allDataEs, cbind(escore, speedup, circTeX))
    allDataC <- rbind(allDataC, cbind(connect[which(circs == circ)], speedup, circTeX))
    allDataL <- rbind(allDataL, cbind(escore, pvalT, circTeX))
}
colnames(allData) <- c("score", "speedup", "circRNA")
colnames(allDataEs) <- c("expectedScore", "speedup", "circRNA")
colnames(allDataC) <- c("connect", "speedup", "circRNA")
colnames(allDataL) <- c("expectedScore", "runtime", "circRNA")

allData$score <- as.numeric(allData$score)
allData$speedup <- as.numeric(allData$speedup)
allDataEs$expectedScore <- as.numeric(allDataEs$expectedScore)
allDataEs$speedup <- as.numeric(allDataEs$speedup)
allDataC$connect <- as.numeric(allDataC$connect)
allDataC$speedup <- as.numeric(allDataC$speedup)
allDataL$expectedScore <- as.numeric(allDataL$expectedScore)
allDataL$runtime <- as.numeric(allDataL$runtime)

# maxim <- max(allData$score)
# tikz(file = "speedupbyscore.tex", width = 15, height = 15)
# plot <- ggplot(allData, aes(x=score, y=speedup, color=circRNA)) +
#         geom_point() +
#         theme(legend.position = "bottom") +
#         labs( x = "Statistics $s(c, g)$", y = "Relative speedup") +
#         scale_y_continuous(breaks=seq(from=0, to=8000, by=1000)) +
#         scale_x_continuous(breaks=floor(seq(from=0, to=maxim, length.out=11)))
# 
# ggsave("speedupbyscoregg.pdf", plot=plot)
# ggsave("speedupbyscoregg.png", plot=plot, type="cairo")
# print(plot)
# dev.off()
# 
# 
# maxim <- max(allDataEs$expectedScore)
# tikz(file = "speedupbyescore.tex", width = 15, height = 15)
# plot <- ggplot(allDataEs, aes(x=expectedScore, y=speedup, color=circRNA)) +
#         geom_point() +
#         theme(legend.position = "bottom") +
#         labs( x = "$E(s(c, g)) * number of annotated mRNAs$", y = "Relative speedup") +
#         scale_y_continuous(breaks=seq(from=0, to=8000, by=1000)) +
#         scale_x_continuous(breaks=floor(seq(from=0, to=maxim, length.out=11)))
# 
# ggsave("speedupbyescoregg.pdf", plot=plot)
# ggsave("speedupbyescoregg.png", plot=plot, type="cairo")
# print(plot)
# dev.off()
# 
# 
# tikz(file = "speedupbyconect.tex", width = 15, height = 15)
# plot <- ggplot(allDataC, aes(x=connect, y=speedup, color=circRNA)) +
#         geom_point() +
#         theme(legend.position = "bottom") +
#         labs( x = "Number of paths to a mRNA", y = "Relative speedup") +
#         scale_y_continuous(breaks=seq(from=0, to=8000, by=1000))
# 
# ggsave("speedupbyconnectgg.pdf", plot=plot)
# ggsave("speedupbyconnectgg.png", plot=plot, type="cairo")
# print(plot)
# dev.off()
# 
# 
# allDataEs$expectedScore <- 1/allDataEs$expectedScore
# maxim <- max(allDataEs$expectedScore)
# tikz(file = "speedupbyoescore.tex", width = 15, height = 15)
# plot <- ggplot(allDataEs, aes(x=expectedScore, y=speedup, color=circRNA)) +
#         geom_point() +
#         theme(legend.position = "bottom") +
#         labs( x = "$1/(E(s(c, g))*number of annotated mRNAs)$", y = "Relative speedup") +
#         scale_y_continuous(breaks=seq(from=0, to=8000, by=1000)) +
#         scale_x_continuous(breaks=floor(seq(from=0, to=maxim, length.out=11)))
# 
# ggsave("speedupbyoescoregg.pdf", plot=plot)
# ggsave("speedupbyoescoregg.png", plot=plot, type="cairo")
# print(plot)
# dev.off()



maxim <- max(allDataL$runtime)
maxim2 <- max(allDataL$expectedScore)
#tikz(file = "runtimebyescore.tex", width = 15, height = 15, standAlone=TRUE)
plot <- ggplot(allDataL, aes(x=expectedScore, y=runtime, color=circRNA)) +
        geom_point(size = 0.5) +
        theme_bw() +
        theme(legend.position = "bottom") +
        #$\\mathbb{E}(s(c, g)) * \\| \\mathbf{g}^m \\|_1$
        labs( x = "Expected score times the number of annotated mRNAs", y = "p-value calculation time (sec)") +
        scale_y_continuous(breaks=floor(seq(from=0, to=maxim, length.out=11))) +
        #scale_x_continuous(breaks=floor(seq(from=0, to=maxim2, length.out=11)))
        scale_x_continuous(breaks=floor(seq(from=0, to=130000, length.out=14))) +
        expand_limits(x = 0, y = 0)

ggsave("runtimebyescoregg.pdf", plot=plot, width = 22, height = 20, units = "cm")
ggsave("runtimebyescoregg.png", plot=plot, type="cairo", width = 22, height = 20, dpi = 150, units = "cm")
#print(plot)
#dev.off()
