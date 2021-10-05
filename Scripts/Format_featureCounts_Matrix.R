#Import featureCounts Count Matrix
Count_Matrix <- read.table("featureCounts_Count_Matrix.txt", header = T)

#Remove Chr, Start, End, Strand, and Length columns. N denotes the #final count column
Count_Matrix <- Count_Matrix[,c(1,7:N)]