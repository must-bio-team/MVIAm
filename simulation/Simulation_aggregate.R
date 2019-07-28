## Setting path
WhereAmI <- "D:/path"

## Load package ##
library(sva)
library(ber)

## Load data ##
data1 <- read.table(paste0(WhereAmI,"Dataset1_matrix.txt"), sep = "\t", header = F)
data2 <- read.table(paste0(WhereAmI,"Dataset2_matrix.txt"), sep = "\t", header = F)
data3 <- read.table(paste0(WhereAmI,"Dataset3_matrix.txt"), sep = "\t", header = F)
sample <- read.table(paste0(WhereAmI,"Dataset_label.txt"), sep = "\t", header = F)

data1 <- read.table("Dataset1_matrix.txt",sep="\t",header=F)
data2 <- read.table("Dataset2_matrix.txt",sep="\t",header=F)
data3 <- read.table("Dataset3_matrix.txt",sep="\t",header=F)
sample <- read.table("Dataset_label.txt",sep = "\t",header = F)

dat<-t(rbind(data1,data2,data3))

# Using Combat_p to elimiate the batch effects
data_combat <- ComBat(dat, batch = sample[,2], mod = NULL, par.prior = TRUE, prior.plots = TRUE)
write.table(data_combat, file = paste0(WhereAmI,"x_meta1.txt"), sep = "\t")

# Using Combat_np to elimiate the batch effects
data_combat <- ComBat(dat, batch = sample[,2], mod = NULL, par.prior = FALSE, prior.plots = TRUE)
write.table(data_combat, file = paste0(WhereAmI,"x_meta2.txt"), sep = "\t")

# Using ber to elimiate the batch effectsber
batch <- as.factor(sample[,2])
data_ber <- ber(t(dat),batch)
data_ber <- t(data_ber)
write.table(data_ber,file = paste0(WhereAmI,"x_meta3.txt"), sep = "\t")

# Using ber_bg to elimiate the batch effects
batch <- as.factor(sample[,2])
data_ber <- ber_bg(t(dat),batch)
data_ber <- t(data_ber)
write.table(data_ber,file = paste0(WhereAmI,"x_meta4.txt"), sep = "\t")

