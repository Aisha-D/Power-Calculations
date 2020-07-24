## Cell line pltos etc
setwd("/gpfs/ts0/home/and202/Power_Calc/data")
## Fix the bug where csv file cannot be reloaded
load("/gpfs/mrc0/projects/Research_Project-MRC190311/DNAm/QCdMRC.rdata")
celltypenormbeta <- celltypenormbeta[-grep("rs", rownames(celltypenormbeta)),]


meandiff <- "5"
meandiffdec <- as.numeric(meandiff)/100

print("loading data")
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
#Setting up parallel processors
library(doParallel)
cl<-makeCluster(16)
registerDoParallel(cl)
#define function to calculate power for each sample size
powerCpG <- function(sd_betas_total){
  library(pwr)
  power1 <- c()
  samples <- seq(1, 320, 1)
  
  for(sample in samples){
    power1 <- c(power1, pwr.t.test(n=sample, d = meandiffdec/(sd_betas_total), sig.level = 9e-8,type="two.sample",alternative="two.sided")$power)
  }
  power2 <- as.data.frame(power1)
  names(power2) <- c(labels(sd_betas_total))
  return(power2)
}
# use parSapply() to run function on each row (for each CpG)
# and for each cell type

#betas_total
print("Running Parrallel Power Calculations for betas_total")
start <- proc.time()
results_total <- as.data.frame(parSapply(cl, sd_betas_total, powerCpG))
samples <- seq(1, 320, 1)
results_total <- cbind(samples, results_total)
end <- proc.time()
print(end - start)
results_total2 <- results_total
colnames(results_total2) <- c("samples", rownames(betas_total))
write.csv(t(results_total2), "Total_cpgs_power.csv")

#betas_double_neg
print("Running Parrallel Power Calculations for betas_double_neg")
start <- proc.time()
results_double_neg <- as.data.frame(parSapply(cl, sd_betas_double_neg, powerCpG))
samples <- seq(1, 320, 1)
results_double_neg <- cbind(samples, results_double_neg)
end <- proc.time()
print(end - start)
results_double_neg2 <- results_double_neg
colnames(results_double_neg2) <- c("samples", rownames(betas_double_neg))
write.csv(t(results_double_neg2), "double_neg_cpgs_power.csv")

#betas_neun_pos
print("Running Parrallel Power Calculations for betas_neun_pos")
start <- proc.time()
results_neun_pos <- as.data.frame(parSapply(cl, sd_betas_neun_pos, powerCpG))
samples <- seq(1, 320, 1)
results_neun_pos <- cbind(samples, results_neun_pos)
end <- proc.time()
print(end - start)
results_neun_pos2 <- results_neun_pos
colnames(results_neun_pos2) <- c("samples", rownames(betas_neun_pos))
write.csv(t(results_neun_pos2), "neun_pos_cpgs_power.csv")

#betas_sox10_pos
print("Running Parrallel Power Calculations for betas_sox10_pos")
start <- proc.time()
results_sox10_pos <- as.data.frame(parSapply(cl, sd_betas_sox10_pos, powerCpG))
samples <- seq(1, 320, 1)
results_sox10_pos <- cbind(samples, results_sox10_pos)
end <- proc.time()
print(end - start)
results_sox10_pos2 <- results_sox10_pos
colnames(results_sox10_pos2) <- c("samples", rownames(betas_sox10_pos))
write.csv(t(results_sox10_pos2), "sox10_pos_cpgs_power.csv")
stopCluster(cl)

###calculate proportions
# create function to calcualte proportion of of values (cpgs) in in row that are over 0.8
prop_fun <- function(row){
  prop <- sum(row > 0.8)/length(row)
}
# change df to matrix and apply prop_fun to each row
props_df_total <- as.matrix(results_total[,c(-1)])
props_total <- apply(props_df_total, 1, prop_fun)
# change df to matrix and apply prop_fun to each row
props_df_double_neg <- as.matrix(results_double_neg[,c(-1)])
props_double_neg <- apply(props_df_double_neg, 1, prop_fun)
# change df to matrix and apply prop_fun to each row
props_df_neun_pos <- as.matrix(results_neun_pos[,c(-1)])
props_neun_pos <- apply(props_df_neun_pos, 1, prop_fun)
# change df to matrix and apply prop_fun to each row
props_df_sox10_pos <- as.matrix(results_sox10_pos[,c(-1)])
props_sox10_pos <- apply(props_df_sox10_pos, 1, prop_fun)
sample_props <- cbind(samples, props_total, props_double_neg, props_neun_pos, props_sox10_pos)
write.csv(sample_props, "total_props_1_320samples.csv")


sample_props <- read.csv("total_props_1500.csv", row.names = 1, header = T, stringsAsFactors = F)

## Can you produce another plot which is the ratio of total:neun+ve and total:sox10 for the y axis values, 
## across the sample sizes? I'd like to be able to say that you need 10x as many samples with bulk tissue for the same power in purified cells.
## Also do you have a table so we can pull out some statistics.

library(ggplot2)
library(reshape2)
sample_props <- read.csv("total_props_1_320samples.csv", row.names = 1, header = T, stringsAsFactors = F)
sample_props$props_total <- format(sample_props$props_total, scientific=F)
sample_props$props_double_neg <- format(sample_props$props_double_neg , scientific=F)
sample_props <- sample_props[-1,]
colnames(sample_props) <- c("samples", "Total", "Double_Neg", "Neun_Pos", "Sox10_Pos")
plot_df <- melt(sample_props, id.vars = "samples")
plot_df$value <- as.numeric(plot_df$value)
ggplot(plot_df, aes(x = samples, y = value, colour = variable))+
  geom_line()+
  labs(x = "Sample Sizes", y = "% of CpG's over 80% power", color = "Sorted Samples")+
  ggtitle("5% mean difference")


# select at point which total is 80% power
powers <- seq(0, 0.9, 0.1) 
tots <- c()
sox <- c()
neun <- c()
double <- c()
for(i in 1:length(powers)){
  j = powers[i]
  t <- sample_props[which(sample_props$props_total >= j),]
  tots <- c(tots, t[1,1])
  s <- sample_props[which(sample_props$props_sox10_pos >= j),]
  sox <- c(sox, s[1,1])
  n <- sample_props[which(sample_props$props_neun_pos >= j),]
  neun <- c(neun, n[1,1])
  d <- sample_props[which(sample_props$props_double_neg >= j),]
  double <- c(double, d[1,1])
}
samplescell <- cbind(powers, tots, sox,neun, double)
samplescell <- as.data.frame(samplescell)

tots.sox <- samplescell$tots/samplescell$sox
tots.neun <- samplescell$tots/samplescell$neun
tots.double <- samplescell$tots/samplescell$double
samplescell <- cbind(samplescell, tots.sox, tots.neun, tots.double)



df <- samplescell[,c("powers", "tots.sox", "tots.neun", "tots.double")]
colnames(df) <- c("powers", 'Sox10+ve', 'NeuN+ve', 'Double-ve')
df_mlt <- melt(df, id.vars = "powers")
filename <- paste("SampleRatiotoTotal_0.0", meandiff, sep  = "")
pdf(paste(filename, "_.pdf", sep = ""))
ggplot(df_mlt, aes(x = powers, y = value, color = variable)) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.1)) +
  labs(x = "proportion of CpG's at 80% power", y = "ratio of samples to total") +
  theme_bw()
dev.off()



