## Emma Walker + Aisha Dahir
## This script can take in any percentage and output results 
## percentage of interest can be inputed as an argument in the command land. 
## args inputs: 1)mean difference 2)samples start 3)samples end 4)samples increments
args <- commandArgs(trailingOnly = TRUE)

library(ggplot2)
library(reshape2)
library(pwr)


# args <- c(4,2,12,1) #helps write skills
setwd("/gpfs/ts0/home/and202/Power_Calc/data")
load("/gpfs/mrc0/projects/Research_Project-MRC190311/DNAm/QCdMRC.rdata")
celltypenormbeta <- celltypenormbeta[-grep("rs", rownames(celltypenormbeta)),]



meandiff <- as.character(args[1])
meandiffdec <- as.numeric(meandiff)/100

# make folder for the results of this percentage difference
cmd <- paste("mkdir", meandiffdec, sep = " ")
system(cmd)
setwd(paste(getwd(), meandiffdec, sep = "/"))


################################### Sample Sizes per sample plots #############################################


# Calculate SD for each probe for total, double negative, neun+ve and sox10+ve.
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

md_sd_betas_total <- median(sd_betas_total)
md_sd_betas_double_neg <- median(sd_betas_double_neg)
md_sd_betas_neun_pos <- median(sd_betas_neun_pos)
md_sd_betas_sox10_pos <- median(sd_betas_sox10_pos)

# Calculate power for each sample type for different samples sizes (e.g. 100,200,â€¦,5000) 
# for fixed effect size with P = 9e-8 based on two sample t-test.  
start <- as.integer(args[2])
end <- as.integer(args[3])
increments <- as.integer(args[4])
samples <- seq(start, end, increments)

sd_list <- c(md_sd_betas_total, md_sd_betas_double_neg, md_sd_betas_neun_pos, md_sd_betas_sox10_pos)
df <- samples
for(sd in sd_list){ #4 different sd's
  power <- c()
  for (i in samples){
    power <- c(power,pwr.t.test(n=i, d = meandiffdec/sd, sig.level = 9e-8,type="two.sample",alternative="two.sided")$power)
  }
  df <- cbind(df, power)
}
colnames(df) <- c("Sample_Size", "Total", "DoubleNeg", "NeunPos", "Sox10Pos")
df <- as.data.frame(df)


## Plot figure
cellTypes<-c("Total", "Double -ve","NeuN +ve","Sox10 +ve")
col_pal<-c("darkgray", "darkgreen", "darkmagenta", "deeppink")


pdf("Power of Samples types.pdf", )

title <- paste(meandiff, "% mean difference", sep = "")
plot(df$Sample_Size, df$Total, type="l",main = title, xlab = "Sample Sizes", ylab = "Power", col = "darkgray")
lines(df$Sample_Size,df$DoubleNeg,col="darkgreen")
lines(df$Sample_Size,df$NeunPos,col="darkmagenta")
lines(df$Sample_Size,df$Sox10Pos,col="deeppink")
abline(h=0.8, col="red", lty=2)
legend("topright", cellTypes, col = col_pal, pch = 15, bg='white')

dev.off()

###################################### Proportion of CpG sites per sample types #################################################
library(doParallel)
cl<-makeCluster(16)
registerDoParallel(cl)
#define function to calculate power for each sample size
powerCpG <- function(var1, X, var2){
  library(pwr)
  power1 <- c()
  samples <- var2
  
  for(sample in samples){
    power1 <- c(power1, pwr.t.test(n=sample, d =var1/(X), sig.level = 9e-8,type="two.sample",alternative="two.sided")$power)
  }
  power2 <- as.data.frame(power1)
  names(power2) <- c(labels(X))
  return(power2)
}

# results_total <- as.data.frame(parSapply(cl, X=sd_betas, FUN=powerCpG, var1=meandiffdec, var2=seq(start,end,increments)))
# use parSapply() to run function on each row (for each CpG)
# and for each cell type

#betas_total
print("Running Parrallel Power Calculations for betas_total")
start.time <- proc.time()
results_total <- as.data.frame(parSapply(cl, X=sd_betas_total, FUN=powerCpG, var1=meandiffdec, var2=seq(start,end,increments)))
samples <- seq(start, end, increments)
results_total <- cbind(samples, results_total)
end.time <- proc.time()
print(end.time - start.time)
results_total2 <- results_total
colnames(results_total2) <- c("samples", rownames(betas_total))
# write.csv(t(results_total2), "Total_cpgs_power.csv")

#betas_double_neg
print("Running Parrallel Power Calculations for betas_double_neg")
start.time <- proc.time()
results_double_neg <- as.data.frame(parSapply(cl, X=sd_betas_double_neg, FUN=powerCpG, var1=meandiffdec, var2=seq(start,end,increments)))
samples <- seq(start, end, increments)
results_double_neg <- cbind(samples, results_double_neg)
end.time <- proc.time()
print(end.time - start.time)
results_double_neg2 <- results_double_neg
colnames(results_double_neg2) <- c("samples", rownames(betas_double_neg))
# write.csv(t(results_double_neg2), "double_neg_cpgs_power.csv")

#betas_neun_pos
print("Running Parrallel Power Calculations for betas_neun_pos")
start.time <- proc.time()
results_neun_pos <- as.data.frame(parSapply(cl, X=sd_betas_neun_pos, FUN=powerCpG, var1=meandiffdec, var2=seq(start,end,increments)))
samples <- seq(start, end, increments)
results_neun_pos <- cbind(samples, results_neun_pos)
end.time <- proc.time()
print(end.time - start.time)
results_neun_pos2 <- results_neun_pos
colnames(results_neun_pos2) <- c("samples", rownames(betas_neun_pos))
# write.csv(t(results_neun_pos2), "neun_pos_cpgs_power.csv")

#betas_sox10_pos
print("Running Parrallel Power Calculations for betas_sox10_pos")
start.time <- proc.time()
results_sox10_pos <- as.data.frame(parSapply(cl, X=sd_betas_sox10_pos, FUN=powerCpG, var1=meandiffdec, var2=seq(start,end,increments)))
samples <- seq(start, end, increments)
results_sox10_pos <- cbind(samples, results_sox10_pos)
end.time <- proc.time()
print(end.time - start.time)
results_sox10_pos2 <- results_sox10_pos
colnames(results_sox10_pos2) <- c("samples", rownames(betas_sox10_pos))
# write.csv(t(results_sox10_pos2), "sox10_pos_cpgs_power.csv")
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

filename <- paste("total_props_",paste(as.character(start), as.character(end), sep = "_"), sep = "")
write.csv(sample_props, paste(filename, ".csv", sep = ""))


sample_props <- read.csv(paste(filename, ".csv", sep = ""), row.names = 1, header = T, stringsAsFactors = F)
sample_props$props_total <- format(sample_props$props_total, scientific=F)
sample_props$props_double_neg <- format(sample_props$props_double_neg , scientific=F)
sample_props <- sample_props[-1,]
colnames(sample_props) <- c("samples", "Total", "Double_Neg", "Neun_Pos", "Sox10_Pos")
plot_df <- melt(sample_props, id.vars = "samples")
plot_df$value <- as.numeric(plot_df$value)

filename <- paste("Proportion of Cpg sites with over 0.8 power at meandiff 0.0", meandiff, sep  = "")
pdf(paste(filename, ".pdf", sep = ""))
ggplot(plot_df, aes(x = samples, y = value, colour = variable))+
  geom_line()+
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Sample Sizes", y = "% of CpG's over 80% power", color = "Sorted Samples")+
  ggtitle(paste(meandiff, "% mean difference", sep = ""))
dev.off()


################################### Ratio Proportion of CpG sites per sample types #############################################

# select at point which total is 80% power
powers <- seq(0, 0.9, 0.1) 
tots <- c()
sox <- c()
neun <- c()
double <- c()
for(i in 1:length(powers)){
  print(i)
  j = powers[i]
  print(j)
  t <- sample_props[which(sample_props$Total >= j),]
  tots <- c(tots, t[1,1])
  s <- sample_props[which(sample_props$Sox10_Pos >= j),]
  sox <- c(sox, s[1,1])
  n <- sample_props[which(sample_props$Neun_Pos >= j),]
  neun <- c(neun, n[1,1])
  d <- sample_props[which(sample_props$Double_Neg >= j),]
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
filename <- paste("Sample Ratios to Total 0.0", meandiff, sep  = "")
pdf(paste(filename, "_.pdf", sep = ""))
ggplot(df_mlt, aes(x = powers, y = value, color = variable)) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.1)) +
  labs(x = "proportion of CpG's at 80% power", y = "ratio of samples to total",
       main = paste("Sample Ratios to Total at ", paste(meandiff, "% mean difference", sep = ""), sep="") ) +
  theme_bw()
dev.off()