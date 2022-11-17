library(ggplot2)
library(data.table)
library(dplyr)
library(pROC)

# load the subtype data
meth = data.table(fread("data/interim/analysis_dataset.tsv",sep="\t",header=T))
dat = meth %>% dplyr::select("Subtype_mRNA","ATF1_methylation","TFAP2B_methylation", "TFAP2B_expr")
dat$basal = ifelse(dat$Subtype_mRNA == "Basal",1,0)

dat %>% filter(basal == F) %>% summarize(median(TFAP2B_methylation))
dat %>% filter(basal == T) %>% summarize(median(TFAP2B_methylation))

dat %>% filter(basal == F) %>% summarize(median(TFAP2B_expr))
dat %>% filter(basal == T) %>% summarize(median(TFAP2B_expr))

# Make plot to answer Dawn's question
dat %>% filter(basal == T) %>% ggplot(aes(TFAP2B_expr)) + geom_histogram(bins=100)
dat %>% filter(basal == T) %>% summarize(nZeros = sum(TFAP2B_expr == 0))
dat %>% filter(basal == T) %>% summarize(mean(TFAP2B_expr))
# nZeros
# 1      5
# mean(TFAP2B_expr)
# 1          8.389515

mod1 = glm(factor(dat$basal) ~ dat$TFAP2B_methylation, family="binomial")
mod2 = glm(factor(dat$basal) ~ dat$TFAP2B_expr, family="binomial")
mod3 = glm(factor(dat$basal) ~ dat$TFAP2B_methylation*dat$TFAP2B_expr, family="binomial")

pointEstimates = exp(summary(mod3)$coef[,1])
lower = exp(summary(mod3)$coef[,1]-1.96*summary(mod3)$coef[,2])
upper = exp(summary(mod3)$coef[,1]+1.96*summary(mod3)$coef[,2])
mod3Res = data.frame("est" = pointEstimates,"lower"=lower,"upper"=upper)

prob1 = predict(mod1,type="response")
prob2 = predict(mod2,type="response")
prob3 = predict(mod3,type="response")

dat$prob1 = prob1
dat$prob2 = prob2
dat$prob3 = prob3

roc(basal ~ prob1,data = dat)
roc(basal ~ prob2,data = dat)
roc(basal ~ prob3,data = dat)

roc(basal ~ ATF1_methylation, data = dat)
roc(basal ~ TFAP2B_methylation, data = dat)
roc(basal ~ TFAP2B_expr, data = dat)

