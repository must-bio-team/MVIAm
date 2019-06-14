setwd("/simulation")
library(sva)
library(ber)

data1<-read.table("Dataset1_matrix.txt",sep="\t",header=F)
data2<-read.table("Dataset2_matrix.txt",sep="\t",header=F)
data3<-read.table("Dataset3_matrix.txt",sep="\t",header=F)
sample<-read.table("Dataset_label.txt",sep = "\t",header = F)

dat<-t(rbind(data1,data2,data3))

# Combat_p #
data_combat<-ComBat(dat, batch=sample[,2],mod=NULL, par.prior=TRUE, prior.plots=TRUE)
write.table(data_combat,file="/Simu_data_meta/x_meta1.txt",sep="\t")

# Combat_np #
data_combat<-ComBat(dat, batch=sample[,2],mod=NULL, par.prior=FALSE, prior.plots=TRUE)
write.table(data_combat,file="Simu_data_meta/x_meta2.txt",sep="\t")

# ber #
batch=as.factor(sample[,2])
data_ber <- ber(t(dat),batch)
data_ber <- t(data_ber)
write.table(data_ber,file="/Simu_data_meta/x_meta3.txt",sep="\t")

# ber_bg #
batch=as.factor(sample[,2])
data_ber <- ber_bg(t(dat),batch)
data_ber <- t(data_ber)
write.table(data_ber,file="/Simu_data_meta/x_meta4.txt",sep="\t")

