# Vennerable

Vennerable provides Venn diagrams in R. 
It displays Venn and Euler diagrams for up to 9 different sets and using a variety of geometries. 
It allows the display of area-weighted Venn diagrams and allows fine graphical control over the result.

This package is not on CRAN and also requires some BioConductor packages. To install it, try something like
`source("https://bioconductor.org/biocLite.R");
biocLite(c("RBGL","graph"));
install.packages("devtools"); library(devtools);
install_github("Vennerable");
library(Vennerable");`

The function based documentations is sketchy. A better guide is the Venn vignette which you can see with
`library(Vennerable);vignette("Venn")`

