library(openxlsx)
source('annotate.R')
source('goterms.R')

library(GGally)
library(network)
library(sna)

deltaReport <- function(circrna, file=paste(circrna, "-short.xlsx", sep=""), target_file=paste(circrna, "-delta-report.txt", sep=""), fdrcutoff=0.01,
                            plotsFolder="deltareport",
                            plotsFiltered1=paste(plotsFolder, "1", sep=""),
                            plotsFiltered2=paste(plotsFolder, "2", sep=""),
                            plotsFiltered5=paste(plotsFolder, "5", sep="")
                        ) {
    dir.create(plotsFolder)
    dir.create(plotsFiltered1)
    dir.create(plotsFiltered2)
    dir.create(plotsFiltered5)
    
    write(paste("Neighborhood of", circrna, "w. r. t. GO terms"), file=target_file, append=FALSE)
    write(paste("=================================================="), file=target_file, append=TRUE)
    
    res <- read.xlsx(file)
    goterms <- res$goTerms[res$fdr < fdrcutoff]
    muCirc <- buildMuCirc(circrna)
    
    circMi <- circ_mirna_map[which(circ_mirna_map$circRNA == circrna),]
    miM <- geneInteract[geneInteract[,1] %in% circMi$miRNAs,]
    
    for(goterm in goterms) {
        print(goterm)
        tryCatch({
            write(paste("\n\nGO term:", goterm, "\n===================================="), file=target_file, append=TRUE)
            
            goMu <- buildGOMu(goterm)
            goM <- buildGOM(goterm)
            
            ImmuCirc <- as.vector(Im %*% muCirc)
            deltaM <- ImmuCirc * goM
            deltaMu <- muCirc * (goMu + as.vector(t(Im) %*% goM))
            
            Mindices <- order(deltaM, decreasing=TRUE)[1:sum(deltaM != 0)]
            Muindices <- order(deltaMu, decreasing=TRUE)[1:sum(deltaMu != 0)]
            
            ms <- c()
            mus <- c()
            
            write(paste("----- miRNAs ------"), file=target_file, append=TRUE)
            for(mui in Muindices) {
                write(paste(all_mirna_list[mui], "( delta: ", deltaMu[mui], ")"), file=target_file, append=TRUE)
                mus <- c(mus, all_mirna_list[mui])
            }
            write(paste("----- mRNAs ------"), file=target_file, append=TRUE)
            for(mi in Mindices) {
                write(paste(all_gene_list[mi], "( delta: ", deltaM[mi], ")"), file=target_file, append=TRUE)
                ms <- c(ms, all_gene_list[mi])
            }
            
            # now plot the interaction network
            circMiL <- circMi[which(circMi$miRNAs %in% mus),]
            miML <- miM[which(miM$mature_mirna_id %in% mus & miM$target_symbol %in% ms),]
            
            colnames(circMiL) <- c("V1","V2")
            colnames(miML) <- c("V1","V2")
            e <- rbind(circMiL, miML)
            net <- network(e, directed = FALSE)

            x = data.frame(molecule = network.vertex.names(net))
            xcp = x
            xcp[x$molecule %in% all_mirna_list,] <- "miRNA"
            xcp[x$molecule == circrna,] <- "circRNA"
            xcp[x$molecule %in% all_gene_list,] <- "mRNA"
            net %v% "molecule" <- as.character(xcp$molecule)
            
            xsz = x
            for(i in 1:nrow(xcp)) {
                if(xcp[i,1] == "miRNA") xsz[i,1] <- deltaMu[which(all_mirna_list == xsz[i,1])]
                if(xcp[i,1] == "mRNA") xsz[i,1] <- deltaM[which(all_gene_list == xsz[i,1])]
                if(xcp[i,1] == "circRNA") xsz[i,1] <- max(c(deltaM, deltaMu)) + 1
            }
            colnames(xsz) <- c("delta")
            net %v% "delta" <- as.numeric(unlist(xsz))
            colfunc <- colorRampPalette(c("grey", "red"))
            net %v% "colr" <- colfunc(max(net %v% "delta"))[net %v% "delta"]
            
            pdf(paste(plotsFolder, "/", circrna, "-", goterm, ".pdf", sep=""))
            print(ggnet2(net, shape = "molecule", alpha = 0.75, size = "delta", edge.alpha = 0.5, label=TRUE, label.size = 2, color = "colr", size.legend = "Influence", shape.legend = "Molecule"))
            dev.off()
            
            delete.vertices(net, which(net %v% "delta" <= 1))
            pdf(paste(plotsFiltered1, "/", circrna, "-", goterm, ".pdf", sep=""))
            print(ggnet2(net, shape = "molecule", alpha = 0.75, size = "delta", edge.alpha = 0.5, label=TRUE, label.size = 2, color = "colr", size.legend = "Influence", shape.legend = "Molecule"))
            dev.off()
            
            delete.vertices(net, which(net %v% "delta" <= 2))
            pdf(paste(plotsFiltered2, "/", circrna, "-", goterm, ".pdf", sep=""))
            print(ggnet2(net, shape = "molecule", alpha = 0.75, size = "delta", edge.alpha = 0.5, label=TRUE, label.size = 2, color = "colr", size.legend = "Influence", shape.legend = "Molecule"))
            dev.off()
                    
            delete.vertices(net, which(net %v% "delta" <= 5))
            pdf(paste(plotsFiltered5, "/", circrna, "-", goterm, ".pdf", sep=""))
            print(ggnet2(net, shape = "molecule", alpha = 0.75, size = "delta", edge.alpha = 0.5, label=TRUE, label.size = 2, color = "colr", size.legend = "Influence", shape.legend = "Molecule"))
            dev.off()
        }, error=function(cond) {
            print("Error in delta report. Maybe an empty net?")
            print(cond)
            return(NULL)
        })
    }
}

args = commandArgs(trailingOnly=TRUE)
circ <- args[1]
deltaReport(circ)
