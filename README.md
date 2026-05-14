# Differential abundance in functional screens

|Environment|Status|
|---|---|
|[BioC-release](https://bioconductor.org/packages/release/bioc/html/augere.screen.html)|[![Release OK](https://bioconductor.org/shields/build/release/bioc/augere.screen.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/augere.screen/)|
|[BioC-devel](https://bioconductor.org/packages/devel/bioc/html/augere.screen.html)|[![Devel OK](https://bioconductor.org/shields/build/devel/bioc/augere.screen.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/augere.screen/)|

Implements a pipeline function to generate parametrized Rmarkdown reports for differential abundance analyses of functional screen data.
The differential abundance analysis is done using `voom()` on the barcode-level counts, followed by consolidation to per-gene statistics for easier interpretation.
The report contains all of the R commands required to reproduce the analysis.
