library(Vennerable)
Weight<- c(0,0,0,0,0,0,0,1)
names(Weight) <- c("000","001","010","011","100","101","110","111")

V00000001 <- Venn(Weight=Weight,SetNames=LETTERS[1:3])

undebug(addSetToDrawing )
(compute.Venn(V00000001 ))

Weight <- c(0,0,0,1,0,1,0,1)
names(Weight) <- c("000","001","010","011","100","101","110","111")
V00000111 <- Venn(Weight=Weight,SetNames=LETTERS[1:3])
(compute.Venn(V00000111 ))

Weight <- c(0,1,0,0,0,0,0,1)
names(Weight) <- c("000","001","010","011","100","101","110","111")
V0001001<- Venn(Weight=Weight,SetNames=LETTERS[1:3])
plot(compute.Venn(V0001001))

Weight <- c(0,1,0,0,0,0,1,0)
names(Weight) <- c("000","001","010","011","100","101","110","111")
V0011000<- Venn(Weight=Weight,SetNames=LETTERS[1:3])
plot(compute.Venn(V0011000))

Weight <- c(0,1,1,0,0,1,0,0)
names(Weight) <- c("000","001","010","011","100","101","110","111")
V01100100<- Venn(Weight=Weight,SetNames=LETTERS[1:3])
try(compute.Venn(V01100100))
