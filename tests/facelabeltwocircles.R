library(Vennerable)
V3.big <- Venn(SetNames=LETTERS[1:3],Weight=2^(1:8))
Vmonth2.big <- V3.big[,c(1:2)]
plot(Vmonth2.big)
