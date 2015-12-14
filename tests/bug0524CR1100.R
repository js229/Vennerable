library(Vennerable)
V4a <- Venn(SetNames=month.name[1:4],Weight=1:16)
CR4a <-  compute.CR(V4a)
VennGetFaceLabels(CR4a)
