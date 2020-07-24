## Power calculations 
## 23/04/2020

# Aim to compare power between bulk brain tissue and purified cell types.
# Steps:
# 1)     Calculate SD for each probe for total, double negative, neun+ve and sox10+ve.
# 2)     Initially calculate median SD for each sample type
# 3)     Calculate power for each sample type for different samples sizes (e.g. 100,200,…,5000) 
#        for fixed effect size with P = 9e-8 based on two sample t-test.  
# 4)     Plot sample size (x-axis) against power (y-axis) for each sample type.

## Load data
setwd("/gpfs/ts0/home/and202/Power_Calc")
load("/gpfs/mrc0/projects/Research_Project-MRC190311/DNAm/QCdMRC.rdata")
library(ggplot2)

# 1)     Calculate SD for each probe for total, double negative, neun+ve and sox10+ve.
# pull out required betas
betas_total <- as.data.frame(celltypenormbeta[,colnames(celltypenormbeta) %in% pheno[which(pheno$Cell.type == "Total"),"Basename"]])
betas_double_neg <- as.data.frame(celltypenormbeta[,colnames(celltypenormbeta) %in% pheno[which(pheno$Cell.type == "Double -ve"),"Basename"]])
betas_neun_pos <- as.data.frame(celltypenormbeta[,colnames(celltypenormbeta) %in% pheno[which(pheno$Cell.type == "NeuN +ve"),"Basename"]])
betas_sox10_pos <- as.data.frame(celltypenormbeta[,colnames(celltypenormbeta) %in% pheno[which(pheno$Cell.type == "Sox10 +ve"),"Basename"]])

#calculate SD for each dataframe
sd_betas_total <- apply(betas_total,1,sd)
sd_betas_double_neg <- apply(betas_double_neg,1,sd)
sd_betas_neun_pos <- apply(betas_neun_pos,1,sd)
sd_betas_sox10_pos <- apply(betas_sox10_pos,1,sd)

#combine data to visualise results 
#Boxplot
combined_sds <- cbind(sd_betas_total, sd_betas_double_neg, sd_betas_neun_pos, sd_betas_sox10_pos)
colnames(combined_sds) <- c(paste("Total", ncol(betas_total), sep = " - "), paste("Double -ve", ncol(betas_double_neg), sep = " - "),
                            paste("NeuN +ve", ncol(betas_neun_pos), sep = " - "), paste("Sox10 +ve", ncol(betas_sox10_pos), sep = " - "))

boxplot(combined_sds)

# Possible to do a line graph for each probe (4 lines)
library(reshape2)
colnames(combined_sds) <- c("Total", "DoubleNeg", "NeunPos", "Sox10Pos")
d <- melt(combined_sds, id.vars="row.names")
png("SD_CellTypes.png")
ggplot(d, aes(x = Var1, y = value, colour=Var2, group = 1)) +
  geom_line() 
dev.off()


# 2)     Initially calculate median SD for each sample type
md_sd_betas_total <- median(sd_betas_total)
md_sd_betas_double_neg <- median(sd_betas_double_neg)
md_sd_betas_neun_pos <- median(sd_betas_neun_pos)
md_sd_betas_sox10_pos <- median(sd_betas_sox10_pos)

combined_md_sd <- c(md_sd_betas_total, md_sd_betas_double_neg, md_sd_betas_neun_pos, md_sd_betas_sox10_pos)
combined_md_sd <- as.data.frame(combined_md_sd)
colnames(combined_md_sd) <- c("Total", "DoubleNeg", "NeunPos", "Sox10Pos")

barplot(combined_md_sd,  
        names.arg=c("Total", "DoubleNeg", "NeunPos", "Sox10Pos"), 
        cex.names=0.8, main = "Median Sd for each cell type")


# Save two plots to pdf's
pdf("CpG_Sd.pdf")
boxplot(combined_sds, main = "Sd of each CpG for each cell type")

barplot(combined_md_sd,  
        names.arg=c("Total", "DoubleNeg", "NeunPos", "Sox10Pos"), 
        cex.names=0.8, main = "Median Sd for each cell type")
dev.off()

# 3)     Calculate power for each sample type for different samples sizes (e.g. 100,200,…,5000) 
#        for fixed effect size with P = 9e-8 based on two sample t-test.  


# Difference is 5%
# P = 9e-8
# Run this for total cell types

library(pwr)
samples <- seq(10, 5000, 10)

pdf("Power_Calc_0.05.pdf", onefile = T)

# Total
power <- c()
for (i in samples){
    power <- c(power,pwr.t.test(n=i, d = 0.05/md_sd_betas_total, sig.level = 9e-8,type="two.sample",alternative="two.sided")$power)
}
df <- cbind(samples, power)
plot(df, main = "Total", xlab = "Sample Sizes", ylab = "Power")

# DoubleNeg
power <- c()
for (i in samples){
  power <- c(power,pwr.t.test(n=i, d = 0.05/md_sd_betas_double_neg, sig.level = 9e-8,type="two.sample",alternative="two.sided")$power)
}
df <- cbind(samples, power)
plot(df, main = "DoubleNeg", xlab = "Sample Sizes", ylab = "Power")

# NeunPos
power <- c()
for (i in samples){
  power <- c(power,pwr.t.test(n=i, d = 0.05/md_sd_betas_neun_pos, sig.level = 9e-8,type="two.sample",alternative="two.sided")$power)
}
df <- cbind(samples, power)
plot(df, main = "NeunPos", xlab = "Sample Sizes", ylab = "Power")

# Sox10Pos
power <- c()
for (i in samples){
  power <- c(power,pwr.t.test(n=i, d = 0.05/md_sd_betas_sox10_pos, sig.level = 9e-8,type="two.sample",alternative="two.sided")$power)
}
df <- cbind(samples, power)
plot(df,main = "Sox10Pos", xlab = "Sample Sizes", ylab = "Power", col = "red")

dev.off()


## Merge all sample sizes onto one plot

#make df with each column different cell type

sd_list <- c(md_sd_betas_total, md_sd_betas_double_neg, md_sd_betas_neun_pos, md_sd_betas_sox10_pos)
samples <- seq(100, 5000, 10)

df <- samples
for(sd in sd_list){ #4 different sd's
  power <- c()
    for (i in samples){
      power <- c(power,pwr.t.test(n=i, d = 0.05/sd, sig.level = 9e-8,type="two.sample",alternative="two.sided")$power)
    }
  df <- cbind(df, power)
}
colnames(df) <- c("Sample_Size", "Total", "DoubleNeg", "NeunPos", "Sox10Pos")

#clean up df
df <- as.data.frame(df)
head(df)
# rownames(df) <- df$Sample_Size
# df <- df[,-1]

# Plot
cellTypes<-c("Total", "Double -ve","NeuN +ve","Sox10 +ve")
col_pal<-c("darkgray", "darkgreen", "darkmagenta", "deeppink")

pdf("PowerCalc_0.05_100-5000.pdf")
plot(df$Sample_Size, df$Total, type="l",main = "5% mean difference", xlab = "Sample Sizes", ylab = "Power", col = "darkgray")
lines(df$Sample_Size,df$DoubleNeg,col="darkgreen")
lines(df$Sample_Size,df$NeunPos,col="darkmagenta")
lines(df$Sample_Size,df$Sox10Pos,col="deeppink")
abline(h=0.8, col="red", lty=2)
legend("topright", cellTypes, col = col_pal, pch = 15, bg='white')
dev.off()
