## Pull together all percentage differences for power calcs
library(reshape2)
library(ggplot2)

setwd("/gpfs/ts0/home/and202/Power_Calc/data")

sample_props_0.01 <- read.csv("0.01/total_props_1_8001.csv", row.names = 1, header = T, stringsAsFactors = F)
sample_props_0.02 <- read.csv("0.02/total_props_1_2000.csv", row.names = 1, header = T, stringsAsFactors = F)
sample_props_0.03 <- read.csv("0.03/total_props_1_8001.csv", row.names = 1, header = T, stringsAsFactors = F)
sample_props_0.04 <- read.csv("0.04/total_props_2_500.csv", row.names = 1, header = T, stringsAsFactors = F)
sample_props_0.05 <- read.csv("0.05/total_props_1_320.csv", row.names = 1, header = T, stringsAsFactors = F)


#some dfs have NAs
sample_props_0.01 <- sample_props_0.01[-1,]
sample_props_0.03 <- sample_props_0.03[-1,]
sample_props_0.05 <- sample_props_0.05[-1,]

#calculate ratio to total
ratiotototal <- function(sample_props){
  colnames(sample_props) <- c("samples", "Total", "Double_Neg", "Neun_Pos", "Sox10_Pos")
  
  powers <- seq(0, 0.9, 0.05) 
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
  
  return(samplescell)
}

sample_ratio_0.01 <- ratiotototal(sample_props_0.01)
sample_ratio_0.02 <- ratiotototal(sample_props_0.02)
sample_ratio_0.03 <- ratiotototal(sample_props_0.03)
sample_ratio_0.04 <- ratiotototal(sample_props_0.04)
sample_ratio_0.05 <- ratiotototal(sample_props_0.05)

# Plot each cell type seperately
pdf("AllCellTypesRatios.pdf")

TotaltoSoxPos <- as.data.frame(cbind(sample_ratio_0.01$powers, 
                                     sample_ratio_0.01$tots.sox, sample_ratio_0.02$tots.sox,
                                     sample_ratio_0.03$tots.sox, sample_ratio_0.04$tots.sox, sample_ratio_0.05$tots.sox))
colnames(TotaltoSoxPos) <- c("powers","0.01", "0.02", "0.03", "0.04", "0.05")
df_mlt <- melt(TotaltoSoxPos, id.vars = "powers")
ggplot(df_mlt, aes(x = powers, y = value, color = variable)) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(x = "proportion of CpG's at 80% power", y = "ratio of samples to total",
       title = "Ratio of Total to Sox10 Positive cell types") +
  guides(color=guide_legend(title="Mean Difference")) +
  theme_bw()



TotaltoNeun <- as.data.frame(cbind(sample_ratio_0.01$powers,
                                   sample_ratio_0.01$tots.neun, sample_ratio_0.02$tots.neun,
                                   sample_ratio_0.03$tots.neun, sample_ratio_0.04$tots.neun, sample_ratio_0.05$tots.neun))
colnames(TotaltoNeun) <- c("powers","0.01", "0.02","0.03", "0.04", "0.05")
df_mlt <- melt(TotaltoNeun, id.vars = "powers")
ggplot(df_mlt, aes(x = powers, y = value, color = variable)) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(x = "proportion of CpG's at 80% power", y = "ratio of samples to total",
       title = "Ratio of Total to Neun Positive cell types") +
  guides(color=guide_legend(title="Mean Difference")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()


TotaltoDouble <- as.data.frame(cbind(sample_ratio_0.01$powers, 
                                     sample_ratio_0.01$tots.double,  sample_ratio_0.02$tots.double,
                                     sample_ratio_0.03$tots.double, sample_ratio_0.04$tots.double, sample_ratio_0.05$tots.double))
colnames(TotaltoDouble) <- c("powers","0.01","0.02", "0.03", "0.04", "0.05")
df_mlt <- melt(TotaltoDouble, id.vars = "powers")
ggplot(df_mlt, aes(x = powers, y = value, color = variable)) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(x = "proportion of CpG's at 80% power", y = "ratio of samples to total",
      title = "Ratio of Total to Double Negative cell types") +
  guides(color=guide_legend(title="Mean Difference")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

dev.off()