# setwd("/mnt/data1/Eilis/Projects/LRAP/")
setwd("/mnt/data1/Aisha/CellLine/")
load("/mnt/data1/Eilis/Projects/LRAP/Normalised.rda")

library(pwr)
library(ggplot2)
# 521 samples

table(sampleSheet$SentrixID)
table(sampleSheet$Cell)
table(sampleSheet$Dosage)
table(sampleSheet$Compound)
unique(sampleSheet$Compound)

table(sampleSheet$Cell, sampleSheet$Compound)


# Filter to water and dmso compounds
# "Water1" "Water2" "DMSO1" "DMSO2" "DMSO3" "DMSO4" 

pheno <- sampleSheet[which(sampleSheet$Compound %in% c("Water1","Water2","DMSO1","DMSO2","DMSO3","DMSO4")),] #105 samples
# make 3df for different cell lines
pheno_KELLY <- pheno[which(pheno$Cell %in% c("KELLY")),]
pheno_MCF7 <- pheno[which(pheno$Cell %in% c("MCF7")),]
pheno_SHSY5Y <- pheno[which(pheno$Cell %in% c( "SHSY5Y")),]

# betas is not a matrix but 3 lists - perhaps by cell lines?
identical(rownames(sampleSheet[which(sampleSheet$Cell %in% c("KELLY")),]), colnames(betas[[2]]))
identical(rownames(sampleSheet[which(sampleSheet$Cell %in% c("MCF7")),]), colnames(betas[[3]]))
identical(rownames(sampleSheet[which(sampleSheet$Cell %in% c("SHSY5Y")),]), colnames(betas[[1]]))
# convert lists to matrix to pull data out
betas_KELLY <- betas[[2]]
betas_MCF7 <- betas[[3]]
betas_SHSY5Y<- betas[[1]]
# filter the cell line matrices to water/dmso compounds
betas_KELLY <- as.data.frame(betas_KELLY[,which(colnames(betas_KELLY) %in% rownames(pheno_KELLY))])
betas_MCF7 <- as.data.frame(betas_MCF7[,which(colnames(betas_MCF7) %in% rownames(pheno_MCF7))])
betas_SHSY5Y <- as.data.frame(betas_SHSY5Y[,which(colnames(betas_SHSY5Y) %in% rownames(pheno_SHSY5Y))])

#calculate SD for each dataframe
sd_betas_KELLY <- apply(betas_KELLY,1,sd)
sd_betas_MCF7 <- apply(betas_MCF7,1,sd)
sd_betas_SHSY5Y <- apply(betas_SHSY5Y,1,sd)

combined_sds <- cbind(sd_betas_KELLY, sd_betas_MCF7, sd_betas_SHSY5Y)
colnames(combined_sds) <- c("KELLY", "MCF7", "SHSY5Y")

boxplot(combined_sds)

md_sd_betas_KELLY <- median(sd_betas_KELLY)
md_sd_betas_MCF7 <- median(sd_betas_MCF7)
md_sd_betas_SHSY5Y <- median(sd_betas_SHSY5Y)

combined_md_sd <- c(md_sd_betas_KELLY, md_sd_betas_MCF7, md_sd_betas_SHSY5Y)
combined_md_sd <- as.data.frame(combined_md_sd)
rownames(combined_md_sd) <- c("KELLY", "MCF7", "SHSY5Y")
print(combined_md_sd)

# Calculate power!

sd_list <- c(md_sd_betas_KELLY, md_sd_betas_MCF7, md_sd_betas_SHSY5Y)
samples <- seq(0, 100, 2)

df <- samples
for(sd in sd_list){ #4 different sd's
  power <- c()
  for (i in samples){
    power <- c(power,pwr.t.test(n=i, d = 0.05/sd, sig.level = 9e-8,type="two.sample",alternative="two.sided")$power)
  }
  df <- cbind(df, power)
}
colnames(df) <- c("Sample_Size", "KELLY", "MCF7", "SHSY5Y")

#clean up df
df <- as.data.frame(df)
head(df)

cellTypes<-c("KELLY", "MCF7", "SHSY5Y")
col_pal<-c("purple", "darkgreen", "blue")

pdf("PowerCalcCellLine_0.05.pdf")
plot(df$Sample_Size, df$KELLY, type="l",main = "5% mean difference", xlab = "Sample Sizes", ylab = "Power", col = "purple")
lines(df$Sample_Size,df$MCF7,col="darkgreen")
lines(df$Sample_Size,df$SHSY5Y, col="blue")
abline(h=0.8, col="red", lty=2)
legend("topright", cellTypes, col = col_pal, pch = 15, bg='white')
dev.off()



#### 2. Proportion of sites that have 80% power
library(doParallel)
cl<-makeCluster(16)
registerDoParallel(cl)
#define function to calculate power for each sample size
powerCpG <- function(sd_betas_total){
  library(pwr)
  power1 <- c()
  samples <- seq(1, 400, 10)
  
  for(sample in samples){
    power1 <- c(power1, pwr.t.test(n=sample, d = 0.05/(sd_betas_total), sig.level = 9e-8,type="two.sample",alternative="two.sided")$power)
  }
  power2 <- as.data.frame(power1)
  names(power2) <- c(labels(sd_betas_total))
  return(power2)
}
# use parSapply() to run function on each row (for each CpG)
# and for each cell type

#betas_KELLY
print("Running Parrallel Power Calculations for betas_KELLY")
start <- proc.time()
results_KELLY <- as.data.frame(parSapply(cl, sd_betas_KELLY, powerCpG))
samples <- seq(1, 400, 10)
results_KELLY <- cbind(samples, results_KELLY)
end <- proc.time()
print(end - start)
# results_KELLY2 <- results_KELLY
# colnames(results_KELLY2) <- c("samples", rownames(betas_KELLY))
# write.csv(t(results_KELLY2), "KELLY_cpgs_power.csv")

#betas_MCF7
print("Running Parrallel Power Calculations for betas_MCF7")
start <- proc.time()
results_MCF7 <- as.data.frame(parSapply(cl, sd_betas_MCF7, powerCpG))
samples <- seq(1, 400, 10)
results_MCF7 <- cbind(samples, results_MCF7)
end <- proc.time()
print(end - start)
# results_MCF72<- results_MCF7
# colnames(results_MCF72) <- c("samples", rownames(betas_MCF7))
# write.csv(t(results_MCF72), "MCF7_cpgs_power.csv")


#betas_SHSY5Y
print("Running Parrallel Power Calculations for betas_SHSY5Y")
start <- proc.time()
results_SHSY5Y <- as.data.frame(parSapply(cl, sd_betas_SHSY5Y, powerCpG))
samples <- seq(1, 400, 10)
results_SHSY5Y <- cbind(samples, results_SHSY5Y)
end <- proc.time()
print(end - start)
# results_SHSY5Y2 <- results_SHSY5Y
# colnames(results_SHSY5Y2) <- c("samples", rownames(betas_SHSY5Y))
# write.csv(t(results_SHSY5Y2), "SHSY5Y_cpgs_power.csv")
stopCluster(cl)


prop_fun <- function(row){
  prop <- sum(row > 0.8)/length(row)
}
# change df to matrix and apply prop_fun to each row
props_df_KELLY <- as.matrix(results_KELLY[,c(-1)])
props_KELLY <- apply(props_df_KELLY, 1, prop_fun)
# change df to matrix and apply prop_fun to each row
props_df_MCF7 <- as.matrix(results_MCF7[,c(-1)])
props_MCF7 <- apply(props_df_MCF7, 1, prop_fun)
# change df to matrix and apply prop_fun to each row
props_df_SHSY5Y <- as.matrix(results_SHSY5Y[,c(-1)])
props_SHSY5Y <- apply(props_df_SHSY5Y, 1, prop_fun)
sample_props <- cbind(samples, props_KELLY, props_MCF7, props_SHSY5Y)
write.csv(sample_props, "celllines_props_samples.csv")

library(ggplot2)
library(reshape2)
sample_props <- read.csv("celllines_props_samples.csv", row.names = 1, header = T, stringsAsFactors = F)
sample_props <- sample_props[-1,]
colnames(sample_props) <- c("samples","KELLY", "MCF7", "SHSY5Y")
plot_df <- melt(sample_props, id.vars = "samples")
pdf("Proportion of sites over 0.8.pdf")
ggplot(plot_df, aes(x = samples, y = value, colour = variable))+
  geom_line()+
  labs(x = "Sample Sizes", y = "Proportion of CpG's over 80% power", color = "Cell Lines")+
  ggtitle("5% mean difference")
dev.off()

