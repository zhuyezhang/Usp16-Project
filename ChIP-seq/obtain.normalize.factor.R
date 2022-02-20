random <- read.table("RANDOM.txt")


random2 <- random[,c(4,8,12,16,20,24,28,32)]

factors2 <-  c()

contain0 <- function(x){
  return(sum(x==0))
}

for(i in 1:dim(random2)[2]){
  random_select <- random2[,c(1,i)]
  ramdom_select2 <- random_select[apply(random_select,1,contain0) == 0,]
  factors2 <- c(factors2,sum(ramdom_select2[,2])/sum(ramdom_select2[,1]))
}

print(1/factors2)
