#!/bin/bash
echo $1
start=`date +%s`

ml R
Rscript run.R $1
Rscript writeTableForEnrichmentMap.R $1
Rscript delta_report.R $1

end=`date +%s`
runtime=$((end-start))
echo "Runtime of the program: $runtime s"
