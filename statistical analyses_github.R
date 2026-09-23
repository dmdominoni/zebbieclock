# load packages

library(tidyverse)
library(readxl)
library(ggplot2)
library(lme4)
library(janitor)
library(data.table)
library(lmerTest)
library(performance)
library(gtools)
library(furrr)
library(stargazer)
library(modelsummary)
library(glmnet)
library(caret)
library(ggeffects)



## load and prepare age predictions from clock models and metadata

meta_data<-read.delim(file.choose(), header=T,sep=",",stringsAsFactors = T) ## load files with predicted epigenetic ages. on github different example files are provided, for training and testing, as well as either single cpg clock or dbscan clocks. please email authors if different datasets of predicted ages are required.

sample_data<-read.delim(file.choose(), header=T,sep=",",stringsAsFactors = T) ## load "metadata_lot1_lot2_lot3.csv", it can then be subsetted depending on which lots need to be used (eg the lot3 can be excluded for these first steps)
names(sample_data)[1] <- "pred_age"
names(sample_data)[7] <- "Sample_Info" ## change the [] argument to 7 or 3 when using training or test files, respectively

full_data<-left_join(sample_data,meta_data,by="Sample_Info")


# AGE DISTRIBUTION PLOT

clock_data<-subset(meta_data,meta_data$Selection99=="clock set") ## this is for the distribution of ages in the training dataset

ggplot(clock_data, aes(x=Age.years, color=Sex, fill=Sex)) +
  geom_histogram(aes(y=..count..), position="dodge", alpha=1, bins=100)+
  #geom_density(alpha=0.2)+
  labs(
    x = "Age (years)",
    y = "Number of birds",
    color = "Sex") +
  scale_x_continuous(limits=c(-0.1,10),breaks=c(0,1,2,3,4,5,6,7,8,9,10),labels=c("0","1","2","3","4","5","6","7","8","9","10"))+
  scale_y_continuous(limits=c(0,9),breaks=c(1,2,3,4,5,6,7,8,9),labels=c("1","2","3","4","5","6","7","8","9"))+
  theme_classic()+
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.title=element_text(size=44),
        legend.text=element_text(size=44),
        axis.text=element_text(size=54),
        axis.title=element_text(size=54))


# STATISTICAL MODELS

# 1. MAE 

MAE<-round(MAE(full_data$pred_age, full_data$actual_age),0)

# 2. correlation test
cor<-cor.test(full_data$actual_age, full_data$pred_age, alternative = c("two.sided", "less", "greater"),
                    
                    method = c("pearson", "kendall", "spearman"),
                    
                    exact = NULL, conf.level = 0.95, continuity = FALSE)


cor$estimate
cor$p.value


# 3. models (linear mixed models for training, linear models for testing)

#LMMs for training datasets
full_data$ID <- substr(full_data$Sample_Info,1,4)
full_data$ID <- as.factor(full_data$ID)
model1<-lmer(pred_age~actual_age+Sex+Lot+(1|ID),data=full_data)
summary(model1)
r2train<-r2(model1)
r2train$R2_marginal
r2train$R2_conditional

modelsummary(model1, shape = term ~ model + statistic, fmt = fmt_decimal(2, 3), statistic = c("std.error", "statistic","p.value"),output = "modtest_linear_singlecpgs.docx") ## save summary table 

#LMs for testing datasets
model2<-lm(pred_age~actual_age+Sex,data=full_data)
summary(model2)
r_sq<-r2(model2)
r_sq_adj<-round(r_sq$R2_adjusted,2)

modelsummary(model2, shape = term ~ model + statistic, fmt = fmt_decimal(2, 3), statistic = c("std.error", "statistic","p.value"),output = "modtest_log_fixedcpg.docx") ## save summary table


# Plot results

ptrain<-print(ggplot(NULL, aes(x=full_data$actual_age, y=full_data$pred_age)) + ## train dataset
                geom_point(color="red",size=3) + 
                geom_smooth(se=FALSE)+
                xlab("Chronological age (days)") +
                ylab("Predicted age (days)") +
                geom_abline(intercept=0, slope=1) +
                xlim(NA,3500)+
                ylim(NA,3500) +
                annotate("text", x=800, y=3250, size=6, label=paste0(" MAE = ",MAE)) +
                annotate("text", x=800, y=2750, size=6, label=paste0(" R-squared = ",r_sq_adj)) +
                theme_classic()+
                theme(axis.text=element_text(size=16, font="arial"),
                      axis.title=element_text(size=20,face="bold",font="arial")))
ggsave(paste0("train_plot_ARcpgs.png"), height=7, width=10)



ptest<-print(ggplot(NULL, aes(x=sample_data_test$Age, y=predicted_ages_test)) + ## test dataset
               geom_point(color="red",size=3) + 
               geom_smooth(method=lm,se=FALSE)+
               xlab("Chronological age (days)") +
               ylab("Predicted age (days)") +
               geom_abline(intercept=0, slope=1) +
               xlim(NA,3500)+
               ylim(NA,3500) +
               annotate("text", x=800, y=3250, size=6, label=(" MAE = 197")) +
               annotate("text", x=800, y=2750, size=6, label=paste0(" R-squared = ",r_sq_test)) +
               theme_classic()+
               theme(axis.text=element_text(size=16),
                     axis.title=element_text(size=20,face="bold"))+
               ggtitle("Test data"))
ggsave(paste0("test_plot_ARcpgs.png"), height=7, width=10)



## Plot for sex differences

p1<-ggpredict(model1,c("Sex"))
ggplot(p1, aes(x = x, y = predicted)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, size = 3) +
  geom_point(size = 10) +
  labs(
    x = "Sex",
    y = "Predicted age (days)") +
 
  theme_classic() + 
  theme(
    legend.position = "none",
    axis.text = element_text(size = 60),
    axis.title = element_text(size = 50, face = "bold")
  )


## MODELS FOR STRESS EFFECT IN EARLY LIFE

stress_data<-read.delim(file.choose(), header=T,sep=",",stringsAsFactors = T) ## load predicted ages for lot3, the early life stress analyses, using file "predictedages_stress_models.csv"

stress_data$Age.diff<-stress_data$eAge_dbscan_730_4-stress_data$actual_age
  
model3<-lmer(Age.diff~treatment+sex+mat_age +(1|nest),data=stress_data)
summary(model3)

modelsummary(model3, shape = term ~ model + statistic, fmt = fmt_decimal(2, 3), statistic = c("std.error", "statistic","p.value"),output = "mod_stress.docx") ## save summary table

p2<-ggpredict(model3,c("Treatment","Sex"))


ggplot() +
  geom_jitter(
    data = stress_data,
    aes(x = Treatment, y = Age.diff, color = Sex),
    position = position_dodge(width = 0.3),
    size=2,alpha = 0.5
  ) +
  geom_point(
    data = p2,
    aes(x = x, y = predicted, color = group),
    size = 5,
    position = position_dodge(width = 0.3)
  ) +
  geom_errorbar(
    data = p2,
    aes(x = x, ymin = conf.low, ymax = conf.high, color = group),
    width = 0.2,size=2,
    position = position_dodge(width = 0.3)
  ) +
  labs(
    x = "Treatment",
    y = "Predicted age offset (days)",
    color = "Sex"
  ) +
  theme_classic()+
  theme(legend.title=element_text(size=30),legend.text=element_text(size=40,face="bold"),axis.text=element_text(size=40),axis.title=element_text(size=36,face="bold"))


## MODELS TO TEST THE DIFFERENCES BETWEEN DIFFERENT CLOCK OUTPUTS BASED ON THE DIFFERENT PARAMETERS

## linear clocks
all_clocks_linear<-read.delim(file.choose(), header=T,sep=",",stringsAsFactors = T) ## load "all_clocks_linear_results.xlsx"

model11<-lm(MAE_age_all_RRBS~partition+ar_filtering+scale(coverage)+scale(size),data=all_clocks_linear)
summary(model11)
anova(model11)
modelsummary(model11, shape = term ~ model + statistic, fmt = fmt_decimal(2, 3),output = "modtest_linear_allparameters.docx") ## save summary table 

## log-linear clocks
all_clocks_log<-read.delim(file.choose(), header=T,sep=",",stringsAsFactors = T) ## load "all_clocks_log_results.xlsx"

model12<-lm(MAE_age_all_RRBS~partition+ar_filtering+scale(coverage)+scale(size),data=all_clocks_log)
summary(model12)
anova(model12)
modelsummary(model12, shape = term ~ model + statistic, fmt = fmt_decimal(2, 3),output = "modtest_log_allparameters.docx") ## save summary table 



## FIND AGE-RELATED CPGS

# renaming age column
sample_data <- rename(sample_data, Age="Age days")

# exclude experimental samples, focus only on control sampels to build the clock
sample_test<-subset(sample_data,sample_data$Selection99=="test set")
sample_data<-subset(sample_data,sample_data$Selection99=="clock set")



##### IMPORT AND PREPARE METHYLATION DATA #####
meth_data <- read.delim("single cpgs/methylation files/Methylation_Report_min20inAll_CpGs_99_bTaeGut1.4.pri(27063).txt") ## change this to "Methylation_eps_5_min10pos_ALL(11972).txt" if you want to try with a regional partition matrix. These are just two example matrices provided on github, if you are interested in other clock matrices, please email corresponding authors

# Setting up the methylation datafile
colnames(meth_data)[13:111] <-  sub(".sorted_readname.dedup.bismark.cov.gz", "", colnames(meth_data)[13:111])
meth_data<-column_to_rownames(meth_data, var="Probe")
meth_data <- meth_data[,-c(1,2,3,4,5,6,6,7,8,9,10,11)]
#colnames(meth_data) <- sub("_[^_]+$", "", colnames(meth_data))

# Because the methods can't handle Nas, I've removed any regions containing any NAs for this study. There are alternative ways, for instance imputing or replacing NAs with 0s, but this is not my favourite way to go 
meth_data <- na.omit(meth_data)

# Transpose
meth_data <- t(meth_data)
meth_data <- meth_data[ order(row.names(meth_data)), ]

# We get errors if all the regions have zero methylation - so let's remove these
#all_zeros_index <- colSums(meth_data)==0 ## just store info on how many columns have zeros, not 100% necessary
meth_data <- meth_data[,colSums(meth_data)>0]



##### FIND CPG SITES RELATED TO AGE #####
meth_data_copy<-as.data.frame(meth_data)
meth_data_copy$Sample_Info<-rownames(meth_data_copy)
meta_age<-full_join(sample_data,meth_data_copy,by="Sample_Info")

# run loop of correlation tests to find cpg sites related to age
results_list <- lapply(20:ncol(meta_age), function(i) {
  # Get the current column
  current_var <- meta_age[[i]]
  
  # Run the correlation test
  cor_test <- cor.test(meta_age$Age, current_var)
  
  # Extract the correlation coefficient and p-value
  corr_coef <- cor_test$estimate
  p_value <- cor_test$p.value
  p_bonf<-round(p.adjust(p_value, "bonferroni"), 3)
  
  # Return a named list with results
  list(
    Variable = colnames(meta_age)[i],
    Correlation = corr_coef,
    P_Value = p_value,
    P_bonf= p_bonf
  )
})


# convert the list to a dataframe
AR_cpgs <- do.call(rbind, lapply(results_list, as.data.frame))
write.table(AR_cpgs,"AR_singlecpgs.csv") ## save full file if needed

# subset original methylation dataset based on cpg sites that have corrected pvalues < 0.05
cpg_bonf<-subset(AR_cpgs,AR_cpgs$P_bonf<0.001)
meth_bonf<- meth_data[, cpg_bonf$Variable]

# total proportion of age-related cpg sites over total number of cpg sites sequenced
nrow(cpg_bonf)/nrow(AR_cpgs) #0.034 single

# analyse and plot distribution of age-related cpg sites across chromosomes
AR_cpgs$chromosome<-substr(AR_cpgs$Variable, 1, 5)
AR_cpgs$chromosome<-gsub(":","",as.character(AR_cpgs$chromosome))
AR_cpgs$chromosome<-gsub("Chr","",as.character(AR_cpgs$chromosome))
AR_cpgs$Chromosome <- factor(AR_cpgs$chromosome, levels = mixedsort(unique(AR_cpgs$chromosome)))
results_summary <- as.data.frame(table(AR_cpgs$chromosome))
colnames(results_summary) <- c("Chromosome", "Count")

cpg_bonf$chromosome<-substr(cpg_bonf$Variable, 1, 5)
cpg_bonf$chromosome<-gsub(":","",as.character(cpg_bonf$chromosome))
cpg_bonf$chromosome<-gsub("Chr","",as.character(cpg_bonf$chromosome))
cpg_bonf<-subset(cpg_bonf,cpg_bonf$chromosome!="NW")
cpg_bonf$chromosome <- factor(cpg_bonf$chromosome, levels = mixedsort(unique(cpg_bonf$chromosome)))
cpg_bonf$chromosome <- factor(cpg_bonf$chromosome, 
                              levels = c("1", "1A", "2", "3", "4", "4A","5", "6", "7", "8", "9", "10","11", "12", "13", "14", "15", "16","17", "18", "19", "20", "21", "22","23", "24", "25", "26", "27", "28","29", "30", "31", "32", "33", "34","35", "36", "37", "W", "Z"))
cpg_bonf$sign<-ifelse(cpg_bonf$Correlation >0,"increasing","decreasing") 
cpg_summary <- as.data.frame(table(cpg_bonf$chromosome,cpg_bonf$sign))
colnames(cpg_summary) <- c("Chromosome", "Sign","Count")

merged_cpg<-merge(results_summary, cpg_summary, by = "Chromosome", suffixes = c("_total", "_age_related"))
merged_cpg$ratio<-merged_cpg$Count_age_related/merged_cpg$Count_total
merged_cpg$ratio<-round(merged_cpg$ratio*100, 1)
merged_cpg<-subset(merged_cpg,merged_cpg$Chromosome!="NW")
merged_cpg$Chromosome <- factor(merged_cpg$Chromosome, levels = mixedsort(unique(merged_cpg$Chromosome)))
merged_cpg$Chromosome <- factor(merged_cpg$Chromosome, 
                                levels = c("1", "1A", "2", "3", "4", "4A","5", "6", "7", "8", "9", "10","11", "12", "13", "14", "15", "16","17", "18", "19", "20", "21", "22","23", "24", "25", "26", "27", "28","29", "30", "31", "32", "33", "34","35", "36", "37", "W", "Z"))