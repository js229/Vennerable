# Vennerable

Vennerable provides Venn diagrams in R. 
It displays Venn and Euler diagrams for up to 9 different sets and using a variety of geometries. 
It allows the display of area-weighted Venn diagrams and allows fine graphical control over the result.

This package needs a couple of BioConductor packages. Something like
`source("https://bioconductor.org/biocLite.R");
biocLite(c("RBGL","graph"))`
should get those.

This package may make it to CRAN one day, but it isn't there now so the easiest way to install it is with the `devtools` package:
`install.packages("devtools"); library(devtools);`. 

Finally you can actually install it with 
`install_github("Vennerable");
library(Vennerable);`

The function based documentation is sketchy. A better guide is the Venn vignette which you can see with
`vignette("Venn")`. 

