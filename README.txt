This repository contains source code and data needed to generate outputs
needed for the following paper:

Petr Ryšavý, Jiří Kléma, Michaela Dostálová Merkerová
circGPA: circRNA Functional Annotation Based on Probability-generating Functions
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04957-8

Outputs might by downloaded from: https://ida.fel.cvut.cz/~rysavy/circgpa/

=== Before we start ===
First, the following packages need to be installed in R. Open R terminal and install using the following commands:
install.packages("openxlsx")
install.packages("xlsx")
install.packages("readr")
install.packages("pracma")
install.packages("stringr")
install.packages("polynom")
install.packages("geometry")
install.packages("tictoc")
BiocManager::install("miRBaseConverter")
install.packages("GO.db")
install.packages("org.Hs.eg.db")
install.packages("httr")
install.packages("biomaRt")
install.packages("multiMiR")
install.packages("msigdbr")
install.packages("GGally")
install.packages("network")
install.packages("sna")
install.packages("DESeq2")
Next, download the source code with the associated graph. There is no need to build your own graph. A graph of interactions is included in the repository. To download the code, call in bash:
git clone https://github.com/petrrysavy/circGPAcorr.git

==== Run with the default graph ====
The repository you just downloaded comes with the graph of interactions used in the paper's experiments. To generate the output, go to bash and call
> ./circgpa-paper/run.sh hsa_circ_0000228
Some of the outputs that can be generated using this command are available at https://ida.fel.cvut.cz/~rysavy/circgpa/.

==== Run with a custom graph ====
If you want to run the circGPAcorr algorithm on your own interaction graph, you need to provide four inputs specifying the underlying graph (assuming a fixed ordering on miRNAs and mRNAs):
* Ammu - the adjacency matrix between the mRNAs and miRNAs. Row m, column mu is 1 iff mRNA m interacts with muRNA mu, 0 otherwise.
* Amuc - a binary vector, where 1 at position mu indicates that the miRNA mu interacts with the circRNA c of interest.
* gom - a binary vector of annotations of mRNAs - 1 if the mRNA m is annotated, 0 otherwise.
* gomu - a binary vector of annotations of miRNAs. Similar to gom.
===== Expression case =====
Additionally, two inputs are needed. In the case of the correlated version, you need to provide expression vectors for miRNAs and mRNAs. It is assumed that only downregulated expressions of miRNAs are used for upregulated circRNAs, and only upregulated expressions of mRNAs are used (and vice versa).
* mrnaExpressionVector - a vector of mRNA log-fold changes.
* mirnaExpressionVector - a vector of miRNA log-fold changes.
To run the code on the toy network from the paper, go to R in the circgpa-paper folder and call:
source('annotateExpressed.R')
gomu <- c(0,1,1);
gom <- c(1,1,1,0,0);
Amuc <- c(1,1,1);
Ammu <- matrix(c(0,1,1,1,1,1,0,0,1,1,0,0,1,1,0), nrow=5, ncol=3, byrow=TRUE);
mirnaExpressionVector <- c(3,2,0);
mrnaExpressionVector <- c(9,0,8,0,3);
annotateVectorized(Amuc , gomu, gom, Ammu, mirnaExpressionVector, mrnaExpressionVector)

===== Correlation case =====
Additionally, two inputs are needed. In the case of the correlated version, you need to provide expression vectors for miRNAs and mRNAs. It is assumed that only downregulated expressions of miRNAs are used for upregulated circRNAs, and only upregulated expressions of mRNAs are used (and vice versa).
* rhoMuCirc - matrix of weights of edges between the circRNA and miRNAs. Only edges showing the negative correlations should be considered.
* rhoMMuImHadammar - haddamar product of Im and matrix of weights of edges between all miRNAs and mRNAs. Only the edges showing negative correlations should be considered.
To run the code on the toy network from the paper, go to R in the circgpa-paper folder and call:
source('annotateExpressed.R')
gomu <- c(0,1,1);
gom <- c(1,1,1,0,0);
Amuc <- c(1,1,1);
Ammu <- matrix(c(0,1,1,1,1,1,0,0,1,1,0,0,1,1,0), nrow=5, ncol=3, byrow=TRUE);
rhoMuCirc <- c(0.7,0,0.3);
rhoMMuImHadammar <- matrix(c(0,0.1,0.3,0.7,0.9,0.2,0,0,0,0.2,0,0,0,0.6,0), nrow=5, ncol=3, byrow=TRUE);
annotateVectorized(Amuc , gomu, gom, Ammu, rhoMuCirc, rhoMMuImHadammar)

==== Known limitations ====
The algorithm requires a long double value with more exponent bits than common on regular desktop computers. The used long double type needs to accommodate (together with some room for operations with them) binomial coefficients up to the number of mRNAs over the size of the annotation term. See pvalue.cpp. 

==== Sources of the data for the default interaction graph ====
The included RData and CSV files contain a snapshot of interaction graph obtained
from the following databases:

[1] CircInteractome (https://circinteractome.nia.nih.gov/index.html)

Dudekula DB, Panda AC, Grammatikakis I, De S, Abdelmohsen K, and Gorospe M.
CircInteractome: A web tool for exploring circular RNAs and their interacting
proteins and microRNAs. RNA Biology, 2016, Jan 2;13(1):34-42


[2] multiMiR R package (http://multimir.org/)

Yuanbin Ru*, Katerina J. Kechris*, Boris Tabakoff, Paula Hoffman, Richard A. Radcliffe,
Russell Bowler, Spencer Mahaffey, Simona Rossi, George A. Calin, Lynne Bemis,
and Dan Theodorescu. (2014) The multiMiR R package and database: integration
of microRNA-target interactions along with their disease and drug associations.
Nucleic Acids Research, doi: 10.1093/nar/gku631.

[3] ENA Quick GO database (https://www.ebi.ac.uk/QuickGO/)

[4] MSigDB database C5 cathegory (http://www.gsea-msigdb.org/gsea/msigdb/)

Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette,
M. A., ... & Mesirov, J. P. (2005). Gene set enrichment analysis: a knowledge-based
approach for interpreting genome-wide expression profiles. Proceedings of the National
Academy of Sciences, 102(43), 15545-15550.
