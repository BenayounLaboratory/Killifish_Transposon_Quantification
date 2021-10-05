#Import Kallisto Count Matrix
Count_Matrix <- read.table("kallisto_count_matrix.txt", header = T, comment.char = "")

#Remove length and eff length columns. N denotes the final count #column
Count_Matrix <- Count_Matrix[,c(1,4:N)]