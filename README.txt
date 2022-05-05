This repository contains source code and data needed to generate outputs
needed for the following paper:

# TODO the paper citation once published
Petr Ryšavý, Jiří Kléma, Michaela Dostálová Merkerová
circGPA: circRNA Functional Annotation Based on Probability-generating Functions

Outputs might by downloaded from: https://ida.fel.cvut.cz/~rysavy/circgpa/

To generate the results call in bash terminal
> run.sh hsa_circ_1234567
First, the following packages need to be installed in R:
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
